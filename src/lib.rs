/// @brief Tools for working with RQTL2 format.
///
/// From
/// https://kbroman.org/qtl2/assets/vignettes/user_guide.html#Data_file_format:
///
/// The input data file formats for R/qtl cannot handle complex crosses, and so
/// for R/qtl2, we have defined a new format for the data files. We’ll describe
/// it here briefly; for details, see the separate vignette on the input file
/// format. QTL mapping data consists of a set of tables of data: marker
/// genotypes, phenotypes, marker maps, etc. In the new format, these different
/// tables are in separate comma-delimited (CSV) files. In each file, the first
/// column is a set of IDs for the rows, and the first row is a set of IDs for
/// the columns. For example, the phenotype data file will have individual IDs
/// in the first column and phenotype names in the first row.
pub mod util {
  use std::collections::HashMap;
  use std::fs::File;
  use std::io::BufRead;
  use std::io::BufReader;
  use std::io::Seek;
  use std::io::SeekFrom;
  /// @brief Batch size (number of lines to read).
  /// @brief R/QTL2 genotype data file parser.
  ///
  /// @note https://kbroman.org/qtl2/assets/vignettes/input_files.html
  pub struct GenoParser {
    file_reader: BufReader<File>,
    comments: Vec<String>,
    /// @note Markers names.
    markers: Vec<String>,
    /// @note Maps snps value to f64 values. E.g. A to 0.5, B to 1.0, etc.
    hab_mapper: HashMap<char, f64>,
    /// @note File cursor position where SNP records start.
    snp_pos_start: u64,
  }

  impl GenoParser {
    /// @brief Reads file at path.
    ///
    /// @param[in] path      path to R/QTl genotype data file.
    /// @param[in] strip_ids determines whether the first column (IDs) should be omitted.
    pub fn new(path: String, hab_mapper: HashMap<char, f64>) -> std::io::Result<Self> {
      let file = File::open(path)?;
      Self::new_with_file(file, hab_mapper)
    }

    pub fn new_with_file(file: File, hab_mapper: HashMap<char, f64>) -> std::io::Result<Self> {
      let mut file_reader = BufReader::new(file);
      let comments = Self::consume_comments(&mut file_reader)?;
      let markers = Self::consume_markers(&mut file_reader)?;
      Ok(GenoParser {
        snp_pos_start: file_reader.seek(SeekFrom::Current(0))?,
        file_reader: file_reader,
        comments: comments,
        markers: markers,
        hab_mapper: hab_mapper,
      })
    }

    pub fn iter(&mut self) -> std::io::Result<GenoParserIter> {
      self.file_reader.seek(SeekFrom::Start(self.snp_pos_start))?;
      GenoParserIter::new(&mut self.file_reader, &self.hab_mapper)
    }

    /// @brief Get comments from genotype file.
    pub fn get_comments(&self) -> &Vec<String> {
      &self.comments
    }

    /// @brief Returns vector of tuples (id, snps) parsed from file.
    ///
    /// @note Rewinds file cursor to the beginning of SNP lines after finishing
    /// reading.
    pub fn read_all(&mut self) -> std::io::Result<Vec<(String, Vec<f64>)>> {
      let snps_start_pos = self.file_reader.seek(SeekFrom::Current(0))?;
      let res = read_geno(&mut self.file_reader, &self.hab_mapper);
      self.file_reader.seek(SeekFrom::Start(snps_start_pos))?;
      res
    }

    fn parse_into(
      parsed_snp_buf: &mut [f64],
      snp_line: &String,
      hab_mapper: &HashMap<char, f64>,
    ) -> std::io::Result<()> {
      let io_err = |bad_str: String, msg: &str| {
        std::io::Error::new(
          std::io::ErrorKind::InvalidInput,
          format!("This line <{}> is an invalid SNP record: {}", &bad_str, msg),
        )
      };
      let snp = match snp_line.split('\t').skip(1).next() {
        Some(snp_str) => snp_str,
        None => {
          return Err(io_err(
            snp_line.clone(),
            "snp record and row id should be separated with tab.",
          ))
        }
      };
      if parsed_snp_buf.len() != snp.len() {
        return Err(io_err(
          snp_line.clone(),
          &format!(
            "Invalid record: there are {} markers, however {} SNPs were parsed.",
            parsed_snp_buf.len(),
            snp.len()
          ),
        ));
      }
      for (buf_slot, snp_char) in parsed_snp_buf.iter_mut().zip(snp.chars()) {
        *buf_slot = hab_mapper
          .get(&snp_char)
          .ok_or(io_err(
            String::from(snp),
            &format!(
              "failed to convert character <{}> to a float value.",
              snp_char
            ),
          ))?
          .clone();
      }
      Ok(())
    }

    fn fill_buffer(
      fill_buf: &mut Vec<f64>,
      lines_iter: &mut std::io::Lines<BufReader<&mut File>>,
      snp_line_size: usize,
      hab_mapper: HashMap<char, f64>,
    ) -> std::io::Result<usize> {
      let mut parsed_lines_counter: usize = 0;
      for (line_slice, snp_line) in fill_buf.chunks_mut(snp_line_size).zip(lines_iter) {
        Self::parse_into(line_slice, &snp_line?, &hab_mapper)?;
        parsed_lines_counter += 1;
      }
      Ok(parsed_lines_counter)
    }

    /// @brief Calculates kinship matrix for given geno data reading it in
    /// batches. The amount of buffer, so as memory consumption, depends on the
    /// amount of logical cores on the machine and amount of snps.
    ///
    /// @note The purpose of calculation of Kinship matrix in batches is to not
    /// load a complete SNPs dataset in memory.
    ///
    /// Given a genotype data file stripped to just containing SNPs called G, a
    /// Kinship matrix is a matrix product G.T (transposed) * G.
    ///
    /// The matrix times its transpose forms a symmetrical matrix, thus it is
    /// not needed to calculate the full matrix, just a triangular part of it.
    ///
    /// Matrix multiplication of transposed matrix by itself (or vice versa) can
    /// be performed via addition of non intersecting (1 column with 1 row, 2
    /// column with 2 row) products of columns of matrix's transpose by rows of
    /// matrix.
    ///
    /// [[1, 2],          [[1, 4],     [[1],  * [[1],[4]]     [[2],   * [[2],[5]]
    ///  [4, 5]]     *     [2, 5]]  =   [4]]               +   [5]]                =
    ///
    ///  = [[1,  4],        [[4,  10],      [[5,  14],
    ///     [4, 16]]    +    [10, 25]]  =    [14, 41]]
    ///
    /// When reading a row of genotype matrix, the column of transposed matrix
    /// can be obtained by transposing this row. Thus, it is possible to
    /// represent genotype matrix transpose times itself (G.T * G) product, as
    /// many products of it's columns by rows, and column can be obtained once
    /// the row is read.
    ///
    /// This is how Kinship matrix can be calculated: read batch from file, copy
    /// and transpose it, multiply transposed matrix by non transposed, add to
    /// result Kinship matrix. After each batch will be processed and added to
    /// the result Kinship matrix it will contain the full Kinship matrix.
    ///
    /// Actual algorithm does not involve matrix copying and transposing and
    /// instead just manipulates matrix indices calculation to achieve same
    /// result.
    ///
    /// Since processing of one batch does not depend on the others, the process
    /// of Kinship matrix calculation can be parallelized: each logical thread
    /// gets 2 buffer, first one contains read rows, and a second one stores the
    /// result of batch multiplication, it is done to not block a shared Kinship
    /// matrix buffer while the calculation is in process. When the thread is
    /// spawned, it locks the read and result buffer dispatched to him by a main
    /// thread and then starts the multiplication. Once the multiplication is
    /// finished and a result buffer contains the part of the resulting Kinship
    /// matrix, the thread locks shared Kinship matrix and merges the results
    /// simultaneously nullifying result buffer to not interfere with the
    /// results calculated by the next threads obtaining this buffer, then
    /// messaging the main thread that the buffer pair on this index is freed.
    ///
    /// Main thread works in a loop: loads data, parses it into a read buffer,
    /// dispatches read/result buffer pair to the thread. If all threads are
    /// busy performing calculations, it waits until one of them will put a
    /// freed buffer pair index to the concurrent queue.
    pub fn calc_kinship(&mut self, batch_size: usize) -> std::io::Result<Vec<f64>> {
      if batch_size < 1 {
        panic!("Batch size can't be less than 1.");
      }
      let ids_num = self.markers.len();
      // Kinship matrix is square.
      let mut common_kinship_matrix: Vec<f64> = vec![0.0; ids_num * ids_num];

      // This amount of snps will be parsed and processed on each iteration.
      let buf_size = ids_num * batch_size;

      let mut read_buf: Vec<f64> = vec![0.0; buf_size];

      let file_reader = self.file_reader.get_mut();
      let mut line_iter = BufReader::new(file_reader).lines();
      let mut total_snps_read: usize = 0;
      loop {
        let read_line_amount = match Self::fill_buffer(
          &mut read_buf,
          &mut line_iter,
          self.markers.len(),
          self.hab_mapper.clone(),
        )? {
          0 => break,
          n => {
            total_snps_read += n;
            n
          }
        };
        if read_line_amount < ids_num {
          let buf = &mut read_buf;
          buf.resize(ids_num * ids_num, 0.0);
        }
        calc_partial_kinship(&mut read_buf, &mut common_kinship_matrix, ids_num);
      }

      assert!(
        total_snps_read >= ids_num,
        format!(
          "Amount of SNPS (lines in file - (1+comments_lines_count)) should be \
           greater or equal to amount of ids \
           (amount of markers). SNP number: {}, IDS number: {}",
          total_snps_read, ids_num
        )
      );

      self.file_reader.seek(SeekFrom::Start(self.snp_pos_start))?;

      let mut res = common_kinship_matrix;

      // Mirror Kinship matrix, since only the upper part was calculated (the
      // Kinship matrix is symmetrical because it's formed from it's transpose times itself).
      for i in 0..ids_num {
        let row_length = ids_num;
        for j in 0..i + 1 {
          res[j * row_length + i] /= total_snps_read as f64;
          res[i * row_length + j] = res[j * row_length + i];
        }
      }
      Ok(res)
    }

    /// @brief Consumes comments lines from the stream. File cursor is left right
    /// after comments.
    pub fn consume_comments(file_reader: &mut BufReader<File>) -> std::io::Result<Vec<String>> {
      let mut buf_str = String::new();
      let mut res = Vec::<String>::new();
      let mut comments_bytes_count: u64 = 0;
      loop {
        let read_bytes_count: usize = file_reader.read_line(&mut buf_str)?;
        if buf_str.starts_with('#') {
          res.push(String::from(&buf_str[1..buf_str.len() - 1]));
        } else {
          // read_line returns Ok(0) when reached EOF.
          if read_bytes_count == 0 {
            return Err(std::io::Error::new(
              std::io::ErrorKind::InvalidInput,
              "File is empty.",
            ));
          }
          file_reader.seek(SeekFrom::Start(comments_bytes_count))?;
          return Ok(res);
        }
        buf_str.clear();
        comments_bytes_count += read_bytes_count as u64;
      }
    }

    /// @brief Consumes markers line from BufRead. File cursor is left right
    /// after comments.
    pub fn consume_markers(file_reader: &mut BufReader<File>) -> std::io::Result<Vec<String>> {
      let mut markers = String::new();
      let start_pos = file_reader.seek(SeekFrom::Current(0))?;
      let markers_len = file_reader.read_line(&mut markers)?;
      file_reader.seek(SeekFrom::Start(start_pos + markers_len as u64))?;
      markers.truncate(markers.len() - 1);
      Ok(
        markers
          .split('\t')
          .skip(1)
          .map(|marker| String::from(marker))
          .collect::<Vec<String>>(),
      )
    }
  }

  pub fn calc_partial_kinship(
    snps: &mut Vec<f64>,
    partial_matrix: &mut Vec<f64>,
    ids_num: usize,
  ) -> () {
    // Algorithm from BLAS dsyrk:
    // http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_gae0ba56279ae3fa27c75fefbc4cc73ddf.html#gae0ba56279ae3fa27c75fefbc4cc73ddf
    //
    // The BLAS Fortran stores array in a column-major format, but the R/qtl2
    // genotype data stored in a row-major format, so this algorithm corresponds
    // to the branch for non transposed, lower triangular part version, however
    // in fact it performs transposed, upper triangular part multiplication
    // (G.T*G).
    //
    // This algorithm branch (Lower, Non transposed) chosen based on CBLAS
    // http://www.netlib.org/blas/blast-forum/cblas.tgz code for dsyrk
    // (cblas_dsyrk.c), which transforms options (Upper, Transposed) to these
    // arguments when called for row-major matrixes.
    //
    // When the matrix stored in row-major way read in column major way,
    // obtained data is a transpose of this matrix:
    // https://en.wikipedia.org/wiki/Row-_and_column-major_order#Transposition
    extern crate cblas;
    extern crate openblas_src;

    let n = ids_num;
    let k = snps.len() / n;
    println!(
      "N = {:?}\nK = {:?}\nRES_MATR_SIZE = {}\nSNPS_SIZE = {}",
      n,
      k,
      partial_matrix.len(),
      snps.len()
    );
    // Only the upper triangular part of matrix will be calculated.
    unsafe {
      cblas::dsyrk(
        cblas::Layout::RowMajor,
        cblas::Part::Upper,
        cblas::Transpose::Ordinary,
        n as i32,
        k as i32,
        1.0,
        &snps[..],
        k as i32,
        1.0,
        &mut partial_matrix[..],
        n as i32,
      );
    }
  }

  /// @brief Reads snps from file.
  /// Returns vector of tuples (id, snps) parsed from file.
  pub fn parse_geno(
    file: &mut File,
    hab_mapper: &HashMap<char, f64>,
  ) -> std::io::Result<Vec<(String, Vec<f64>)>> {
    let mut file_reader = BufReader::new(file.try_clone()?);
    GenoParser::consume_comments(&mut file_reader)?;
    GenoParser::consume_markers(&mut file_reader)?;
    read_geno(&mut file_reader, hab_mapper)
  }

  /// @note Parse line with markers. File cursor is rewinded to the beginning of
  /// the file.
  /// Example: marker	10	12	38	39	42	54
  pub fn parse_markers(file: &mut File) -> std::io::Result<Vec<String>> {
    let mut buf_reader = BufReader::new(file.try_clone()?);
    GenoParser::consume_comments(&mut buf_reader)?;
    let res = GenoParser::consume_markers(&mut buf_reader);
    file.seek(SeekFrom::Start(0))?;
    res
  }

  /// @brief Parses comments from the beginning of the file. File cursor is
  /// rewinded to the beginning of the file.
  ///
  /// @note Comments lines example:
  /// # These are comments.
  /// # Only at the beginning of the R/QTL2 file geno file.
  pub fn parse_comments(file: &mut File) -> std::io::Result<Vec<String>> {
    let res = GenoParser::consume_comments(&mut BufReader::new(file.try_clone()?));
    file.seek(SeekFrom::Start(0))?;
    res
  }

  /// @brief Parses snps from BufRead.
  ///
  /// Returns vector of tuples (id, snps) parsed from file.
  pub fn read_geno(
    file_reader: &mut dyn BufRead,
    hab_mapper: &HashMap<char, f64>,
  ) -> std::io::Result<Vec<(String, Vec<f64>)>> {
    let mut contents = Vec::<(String, Vec<f64>)>::new();
    for line in file_reader.lines() {
      let id_snp_tuple: (String, Vec<f64>) = parse_snp_rec(line?, &hab_mapper)?;
      contents.push(id_snp_tuple);
    }
    Ok(contents)
  }

  /// @brief Parses snp geno record into tuple. Consumes line with record.
  /// <rs41245 AABH> to ("rs41245", Vec<f64>(0.0, 0.0, 1.0, 0.5))
  pub fn parse_snp_rec(
    line: String,
    hab_mapper: &HashMap<char, f64>,
  ) -> std::io::Result<(String, Vec<f64>)> {
    let line_str = line;
    let mut id_snp = line_str.split('\t');
    let id = id_snp.next().unwrap();
    let parse_snps = |snp_str: &str| {
      snp_str
        .chars()
        .map(|ch| {
          hab_mapper
            .get(&ch)
            .map(|v| *v)
            .clone()
            .expect(&format!("No key <{}> in SNP mapper.", ch)[..])
        })
        .collect::<Vec<f64>>()
    };
    let id_snp_tuple: (String, Vec<f64>) = (
      String::from(id),
      parse_snps(&id_snp.next().ok_or(std::io::Error::new(
        std::io::ErrorKind::InvalidInput,
        format!("This line <{}> is an invalid SNP record.", line_str),
      ))?),
    );
    Ok(id_snp_tuple)
  }

  /// @brief Parses lines from genotype file.
  pub struct GenoParserIter<'a> {
    lines_reader: std::io::Lines<&'a mut BufReader<File>>,
    hab_mapper: &'a HashMap<char, f64>,
  }

  impl<'a> GenoParserIter<'a> {
    /// @note File cursor must be located at the beginning of SNP records.
    fn new(
      file_reader: &'a mut BufReader<File>,
      hab_mapper: &'a HashMap<char, f64>,
    ) -> std::io::Result<Self> {
      Ok(Self {
        lines_reader: file_reader.lines(),
        hab_mapper: hab_mapper,
      })
    }
  }

  impl<'a> Iterator for GenoParserIter<'a> {
    type Item = (String, Vec<f64>);

    /// @brief Parse next line from genotype file. Returns tuple (row_id, snps).
    fn next(&mut self) -> Option<Self::Item> {
      // While EOF is not reached (and until buffer is filled).
      match self.lines_reader.next() {
        Some(line) => {
          return match parse_snp_rec(line.expect("Failed to read SNP line."), &self.hab_mapper) {
            Ok(val) => Some(val),
            Err(e) => {
              println!("Failed to parse the line. Error: {}", e);
              return self.next();
            }
          }
        }
        None => return None,
      }
    }
  }
}
