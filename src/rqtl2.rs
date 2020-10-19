use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Seek;
use std::io::SeekFrom;

pub fn io_err(bad_str: String, msg: &str) -> std::io::Error {
  std::io::Error::new(
    std::io::ErrorKind::InvalidInput,
    format!("This line <{}> is an invalid SNP record: {}", &bad_str, msg),
  )
}

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
    snp_line: &String,
    parsed_snp_buf: &mut [f64],
    hab_mapper: &HashMap<char, f64>,
  ) -> Result<(), crate::util::error::ParsingError> {
    let snp = match snp_line.split('\t').skip(1).next() {
      Some(snp_str) => snp_str,
      None => Err(io_err(
        snp_line.clone(),
        "snp record and row id must be separated with tab.",
      ))?,
    };
    if parsed_snp_buf.len() != snp.len() {
      Err(io_err(
        snp_line.clone(),
        &format!(
          "Invalid record: there are {} markers, however {} SNPs were parsed.",
          parsed_snp_buf.len(),
          snp.len()
        ),
      ))?;
    }
    for (buf_slot, snp_char) in parsed_snp_buf.iter_mut().zip(snp.chars()) {
      *buf_slot = hab_mapper
        .get(&snp_char)
        .ok_or(io_err(
          String::from(snp),
          &format!("failed to convert SNP <{}> to a float value.", snp_char),
        ))?
        .clone();
    }
    Ok(())
  }

  /// @brief Calculates kinship matrix for given geno data reading it in
  /// batches. The amount of buffer, so as memory consumption, depends on the
  /// amount of logical cores on the machine and amount of snps.
  pub fn calc_kinship(
    &mut self,
    batch_size: usize,
  ) -> Result<Vec<f64>, crate::util::error::ParsingError> {
    if batch_size < 1 {
      panic!("Batch size can't be less than 1.");
    }
    let ids_num = self.markers.len();
    // Kinship matrix is square.
    let mut common_kinship_matrix: Vec<f64> = vec![0.0; ids_num * ids_num];

    // This amount of snps will be parsed and processed on each iteration.
    let buf_size = ids_num * batch_size;

    let file_reader = self.file_reader.get_mut();
    let mut line_iter = BufReader::new(file_reader).lines();
    let mut total_snps_read: usize = 0;
    let markers_len = self.markers.len();
    let hab_mapper = self.hab_mapper.clone();

    use crate::util::kinship::WorkUnit;
    let kinship_merger_delegate =
      |work_unit: &mut WorkUnit| -> Result<bool, crate::util::error::ParsingError> {
        // Merge results and restore result buffer.
        for (buf_elem, common_matrix_elem) in work_unit
          .result_buf
          .iter_mut()
          .zip(common_kinship_matrix.iter_mut())
        {
          *common_matrix_elem += *buf_elem;
          *buf_elem = 0.0;
        }

        let parser =
          |source: &String, dest: &mut [f64]| Self::parse_into(source, dest, &hab_mapper);

        // Parse new batch (load work unit).
        let line_count = match crate::util::kinship::fill_buffer(
          &mut work_unit.input_buf.chunks_mut(markers_len),
          &mut line_iter,
          parser,
        )? {
          n => {
            total_snps_read += n;
            n
          }
        };
        // Reached EOF (there is not enough records to form a full batch.).
        // Resize buffer to discard data from previous iterations which was not
        // overwritten because there is not enough lines to fill the whole
        // buffer.
        if line_count < batch_size {
          work_unit.input_buf.resize(line_count * ids_num, 0.0);
          return Ok(true);
        }
        return Ok(false);
      };

    use crate::util::kinship::calc_kinship_parallel;
    calc_kinship_parallel(kinship_merger_delegate, buf_size, ids_num)
      .expect("Parallel processing failed.");

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

    crate::util::kinship::mirror_and_normalize_kinship(
      &mut common_kinship_matrix[..],
      ids_num,
      1.0 / (total_snps_read as f64),
    );

    Ok(common_kinship_matrix)
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
    // While EOF is not reached.
    let line = self.lines_reader.next()?;

    match parse_snp_rec(line.expect("Failed to read SNP line."), &self.hab_mapper) {
      Ok(val) => Some(val),
      Err(e) => {
        println!("Failed to parse the line. Error: {}", e);
        return self.next();
      }
    }
  }
}
