use crate::util::error::ParsingError;
use crate::util::tokenizers::tokenize_bimbam_or_rqtl2_line;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Seek;
use std::io::SeekFrom;

// The input files can be comma-delimited, space-delimited, tab-delimited,
// and semi-colon de-limited, or mixed use of those.  (i.e.  entries can be
// separated by commas, spaces, semi-colons or tabs).
pub const GENO_SEPARATORS: &[char] = &['\t', ' ', ',', ';'];

// Represents SNP id and alleles_types parsed from BIMBAM Mean Genotype File
// Format file.
// Example:
//  This is what this type alias object holds.
//  |   |  |
//  ∨   ∨  ∨
// rs1, A, T, 0.02, 0.80, 1.50 rs2, G, C, 0.98, 0.04, 1.00
type SNPAndAllelesTypes = (String, (String, String));

enum LineType {
  SNP,
  EMPTY,
  COMMENT,
}

pub struct GenoReader {
  file_reader: BufReader<File>,
  ids_num: usize,
}

impl GenoReader {
  pub fn new(path: &str) -> std::io::Result<Self> {
    let file = File::open(path)?;
    Self::with_file(file)
  }

  pub fn with_file(file: File) -> std::io::Result<Self> {
    let mut file_reader = BufReader::new(file);
    // Get number of individuals per line
    let mut buf_str = String::new();
    let read_bytes_count: usize = file_reader.read_line(&mut buf_str)?;

    if read_bytes_count == 0 {
      return Err(std::io::Error::new(
        std::io::ErrorKind::InvalidInput,
        "File is empty.",
      ));
    }

    let elems = buf_str
      .split(&GENO_SEPARATORS[..])
      .filter(|s| !s.is_empty());

    let len = elems.clone().count();
    if len < 4 {
      return Err(std::io::Error::new(
          std::io::ErrorKind::InvalidInput,
          format!(
            "File is incorrect.\n\
            This line is invalid: <{}>\n\
            (There have to be at least 4 entities at line (but only {} were read): SNP ID, 2 alleles types\
            and 1 mean genotype.\n\
            Example:\n\
            rs1, A, T, 0.02\n\
            rs2, G, C, 0.98\n",
            buf_str,
            len
          ),
        ));
    };
    // Skip SNP id and 2 alleles types
    let mean_genotypes = elems.skip(3);
    // There is mean genotype for each id
    let ids_num = mean_genotypes.count();

    file_reader.seek(std::io::SeekFrom::Start(0))?;
    Ok(GenoReader {
      file_reader,
      ids_num,
    })
  }

  /// @brief Parses mean geno line in form of bytestring into preallocated buffer.
  /// @note Only mean genos are filled into the buffer.
  fn parse_into(
    snp_name_from_snp_pos_iter: &str,
    source: &[u8],
    dest: &mut [f64],
  ) -> Result<(), ParsingError> {
    let mut tokens = tokenize_bimbam_or_rqtl2_line(source);

    let line_type = consume_and_check_id_consume_alleles(snp_name_from_snp_pos_iter, &mut tokens)?;

    match line_type {
      LineType::SNP => (),
      LineType::EMPTY | LineType::COMMENT => panic!(
        "Empty or comment lines should've been\
   covered inside fill_buff function."
      ),
    };

    // Buffer is splitted on chunks of size ids_num which are supplied to this
    // function
    let ids_num = dest.len();

    // https://thebird.nl/blog/work/rotate.html#org74d1758
    let mut geno_mean: f64 = 0.0;
    // let mut geno_var: f64 = 0.0;

    let mut parsed_tokens_num = 0;
    // Array containing indices of NA values in parsed SNP array.
    let mut na_genos = vec![false; ids_num];
    for ((i, buf_slot), token) in dest
      .iter_mut()
      .enumerate()
      .zip(tokens.by_ref().map(bstr::ByteSlice::to_str_lossy))
    {
      if token == "NA" {
        na_genos[i] = true;
        continue;
      }
      let mean_geno = token.parse::<f64>()?;
      geno_mean += mean_geno;
      // geno_var += mean_geno * mean_geno;
      *buf_slot = mean_geno;
      parsed_tokens_num += 1;
    }
    geno_mean /= ids_num as f64;
    // geno_var /= ids_num as f64;
    // geno_var -= geno_mean * geno_mean;

    // Replace NA values with SNP mean value.
    for (_, buf_slot) in na_genos
      .into_iter()
      .zip(dest.iter_mut())
      .filter(|(is_missing, _)| *is_missing)
    {
      *buf_slot = geno_mean;
    }
    // Transform row.
    for buf_slot in dest.iter_mut() {
      *buf_slot -= geno_mean;
      // *buf_slot /= geno_var.sqrt();
    }

    if tokens.next() != None || parsed_tokens_num != ids_num {
      use bstr::ByteSlice;
      return Err(ParsingError::from(std::io::Error::new(
        std::io::ErrorKind::InvalidInput,
        format!(
          "Error: buffer length is not equal to the amount of tokens.\n\
        This line is invalid: <{}>",
          source.to_str().unwrap()
        ),
      )));
    }

    Ok(())
  }

  /// @brief Calculates LOCO Kinship matrix for every chromosome.
  ///
  /// @note This function assumes that the SNP ids in both <Mean Genotype File
  /// Format> and <SNP Location File Format> are going in the same order and
  /// there is exactly one position entry for one geno entry. It also assumes
  /// that the SNP ids are sorted in the chromosome order.
  pub fn calc_kinship(
    &mut self,
    batch_size: usize,
    snp_pos_file: std::fs::File,
  ) -> std::io::Result<Vec<(String, Vec<f64>)>> {
    if batch_size < 1 {
      panic!("Batch size can't be less than 1.");
    }

    let snp_positions: Vec<SnpPosRecord> = SnpPosIter::with_file(snp_pos_file)?.collect();

    if snp_positions.is_empty() {
      panic!("File is empty.");
    }

    let chromosome_count = SnpPosIter::count_chromosomes(&snp_positions);

    // Kinship matrices for each chromosome. Kinship matrix is square.
    let mut chrs_kinship_matrices: Vec<Vec<f64>> =
      vec![vec![0.0; self.ids_num * self.ids_num]; chromosome_count];

    let ids_num = self.ids_num;
    // Maximum amount of SNPs that will be parsed and processed on each iteration.
    let buf_size = self.ids_num * batch_size;

    // let file_reader = self.file_reader.get_mut();
    // let mut line_iter = BufReader::new(file_reader).lines();

    // Amount of SNPs parsed for each chromosome. Since it is assumed that the
    // SNPs are ordered by chromosomes, it's simply a sequential array of values
    // (no need for String keys).
    let mut chrs_parsed_snps_nums = Vec::<usize>::new();
    // Maps chromosomes indices to their names.
    let mut chrs_names = Vec::<Option<String>>::new();

    use crate::util::kinship::WorkUnit;
    // Tracks the chromosome currently being processed.
    // All SNPs without specified chromosome are associated with None chromosome
    let mut cur_chromosome = snp_positions[0].2.clone();
    // Tracks amount of snps on current chromosome.
    let mut cur_chr_total_snps_parsed: usize = 0;
    let mut cur_snp_idx = 0;

    let kinship_merger_delegate = |work_unit: &mut WorkUnit| -> Result<bool, ParsingError> {
      // Merge results and restore result buffer.
      let merge_into = |chr_kinship_matrix: &mut Vec<f64>| {
        for (buf_elem, matrix_elem) in work_unit
          .result_buf
          .iter()
          .zip(chr_kinship_matrix.iter_mut())
        {
          *matrix_elem += *buf_elem;
        }
      };

      // LOCO
      for (i, res_matrix) in chrs_kinship_matrices.iter_mut().enumerate() {
        if i != work_unit.chr_num {
          merge_into(res_matrix);
        }
      }

      for buf_elem in work_unit.result_buf.iter_mut() {
        *buf_elem = 0.0;
      }

      let mut get_next_batch = || -> Option<(Vec<&String>, Option<Option<String>>)> {
        let mut snps_ids = Vec::new();
        loop {
          if cur_snp_idx == snp_positions.len() {
            // Parsed all records.
            return None;
          }
          let (snp_id, _, snp_chr) = &snp_positions[cur_snp_idx];
          if cur_chromosome != *snp_chr {
            // New chromosome encountered.
            return Some((snps_ids, Some(snp_chr.clone())));
          }
          if snps_ids.len() == batch_size {
            // Batch is full.
            break;
          }
          snps_ids.push(snp_id);
          cur_snp_idx += 1;

          if cur_snp_idx == snp_positions.len() {
            // No more records left.
            break;
          }
        }
        // Batch was filled. No new chromosome was encountered yet.
        Some((snps_ids, None))
      };

      // Batch of SNP ids on a single chromosome to be processed.
      let (snp_batch_on_chr, new_chr) = match get_next_batch() {
        Some(tup) => tup,
        // Reached EOF. Resize buffer to discard data from previous iterations
        // which was not overwritten because there is not enough lines to fill the
        // whole buffer.
        None => {
          if cur_chr_total_snps_parsed != 0 {
            chrs_parsed_snps_nums.push(cur_chr_total_snps_parsed);
            chrs_names.push(cur_chromosome.clone());
            cur_chr_total_snps_parsed = 0;
            cur_chromosome = None;
          }

          let mut str_buf = String::new();
          let err_not_enough_snps = std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "Error: BIMBAM Mean Genotype File Format file \
              has more records than SNP Location File Format file does.\n",
          );

          // Check if EOF is reached
          if self.file_reader.read_line(&mut str_buf)? != 0 {
            return Err(ParsingError::from(err_not_enough_snps));
          }
          return Ok(true);
        }
      };

      let mut snps_ids = snp_batch_on_chr.iter();
      let parser =
        |source: &[u8], dest: &mut [f64]| Self::parse_into(snps_ids.next().unwrap(), source, dest);

      // Resize the buffer so only the SNPs on this chromosome get parsed into it.
      work_unit
        .input_buf
        .resize(snp_batch_on_chr.len() * ids_num, 0.0);

      // Parse new batch (load work unit).
      let line_count = crate::util::kinship::fill_buffer_from_bytes(
        &mut work_unit.input_buf.chunks_mut(ids_num),
        &mut self.file_reader,
        parser,
      )?;

      cur_chr_total_snps_parsed += line_count;

      if line_count < snp_batch_on_chr.len() {
        return Err(ParsingError::from(std::io::Error::new(
          std::io::ErrorKind::InvalidInput,
          "Error: BIMBAM Mean Genotype File Format file\
            has fewer records than SNP Location File Format file does.\n",
        )));
      }

      // Set the current chromosome number. <chrs_parsed_snps_nums> contains
      // records for all chromosome parsed or in process of parsing up to this
      // moment. The id of the last chromosome in the array (or size of the
      // array) represents the current chromosome id, since the actual id,
      // which is represented as String does not matter (the records are
      // ordered by chromosome).
      work_unit.chr_num = chrs_parsed_snps_nums.len();

      // Resize buffer so old data gets discarded.
      if line_count < batch_size {
        work_unit.input_buf.resize(line_count * ids_num, 0.0);
      }

      if let Some(n_chr) = new_chr {
        chrs_parsed_snps_nums.push(cur_chr_total_snps_parsed);
        chrs_names.push(cur_chromosome.clone());
        cur_chromosome = n_chr;
        cur_chr_total_snps_parsed = 0;
      }

      Ok(false)
    };

    use crate::util::kinship::calc_kinship_parallel;
    calc_kinship_parallel(kinship_merger_delegate, buf_size, ids_num)
      .expect("Parallel processing failed.");

    let total: usize = chrs_parsed_snps_nums.iter().sum();

    for (chrs_kinship_matrix, chrs_parsed_snps_num) in chrs_kinship_matrices
      .iter_mut()
      .zip(chrs_parsed_snps_nums.into_iter())
    {
      crate::util::kinship::mirror_and_scale_kinship(
        &mut chrs_kinship_matrix[..],
        self.ids_num,
        1.0 / ((total - chrs_parsed_snps_num) as f64),
      );
    }

    self.file_reader.seek(SeekFrom::Start(0))?;

    Ok(
      chrs_names
        .into_iter()
        .map(|name| name.unwrap_or_else(|| "None".to_owned()))
        .zip(chrs_kinship_matrices.into_iter())
        .collect(),
    )
  }

  pub fn iter(&mut self) -> std::io::Result<GenoReaderIter> {
    self.file_reader.seek(SeekFrom::Start(0))?;
    Ok(GenoReaderIter {
      file_reader: &mut self.file_reader,
      ids_num: self.ids_num,
      read_buf: Vec::<u8>::new(),
    })
  }
}

/// @brief Parsed line of BIMBAM Mean Genotype File Format file and optionally
/// can contain data from SNP Location File Format file.
///
/// @note From BIMBAM manual
/// (https://www.haplotype.org/download/bimbam-manual.pdf):
///
/// The first column of the mean genotype files is the SNP ID, the second and
/// third columns are allele types with minor allele first.  The rest columns
/// are the mean genotypes of different individuals – numbers between 0 and 2
/// that represents the (posterior) mean genotype, or6 dosage of the minor
/// allele.  
///
/// An example of mean genotypes file of two SNPs and three individuals follows:
/// rs1, A, T, 0.02, 0.80, 1.50 rs2, G, C, 0.98, 0.04, 1.00
#[derive(Debug)]
pub struct MeanGenoLine {
  pub snp_id: String,
  pub alleles_types: (String, String),
  pub mean_genos: Vec<f64>,
  // Physical location of SNP
  pub pos: Option<usize>,
  // Chromosome ID
  pub chr: Option<usize>,
}

pub struct GenoReaderIter<'a> {
  file_reader: &'a mut BufReader<File>,
  ids_num: usize,
  read_buf: Vec<u8>,
}

impl<'a> Iterator for GenoReaderIter<'a> {
  type Item = MeanGenoLine;

  fn next(&mut self) -> Option<Self::Item> {
    // Dump old data
    self.read_buf.clear();

    let on_err = |it: &mut Self, e: ParsingError| {
      eprintln!("Failed to parse line. Parsing next line... Error: <{}>", e);
      it.next()
    };

    match self.file_reader.read_until(b'\n', &mut self.read_buf) {
      Ok(read_bytes_n) => {
        if read_bytes_n == 0 {
          // EOF is reached.
          return None;
        }
      }
      Err(e) => return on_err(self, ParsingError::Io(e)),
    };

    match parse_geno_line(&self.read_buf) {
      Ok(Some(mean_geno_line)) => {
        if mean_geno_line.mean_genos.len() != self.ids_num {
          eprintln!(
            "{}",
            format!(
              "Number of mean genotypes on this line is incorrect. \
          Expected {} genotypes, but {} were parsed. \
          Skipping this line and parsing next one...",
              self.ids_num,
              mean_geno_line.mean_genos.len()
            )
          );
          self.next()
        } else {
          Some(mean_geno_line)
        }
      }
      Ok(None) => {
        eprintln!("Skipped empty or comment line.");
        self.next()
      }
      Err(e) => on_err(self, e),
    }
  }
}

/// @brief Parses line from BIMBAM Mean Genotype File Format file.
/// If line is empty or contains a comment, None is returned.
///
/// @in line full line from file, e.g. "rs1, A, T, 0.02"
///
/// @note https://www.haplotype.org/download/bimbam-manual.pdf
pub fn parse_geno_line(line: &[u8]) -> Result<Option<MeanGenoLine>, ParsingError> {
  let mut tokens = tokenize_bimbam_or_rqtl2_line(line);

  let (snp_id, alleles_types) = match parse_id_and_alleles(&mut tokens)? {
    Some(id_and_alleles) => id_and_alleles,
    None => return Ok(None),
  };

  let mean_genos: Vec<f64> = tokens
    .map(bstr::ByteSlice::to_str_lossy)
    .map(|mean_geno| mean_geno.parse::<f64>())
    .collect::<Result<Vec<f64>, _>>()?;
  if mean_genos.is_empty() {
    use bstr::ByteSlice;
    return Err(ParsingError::from(std::io::Error::new(
      std::io::ErrorKind::InvalidInput,
      format!(
        "BIMBAM Mean Genotype File Format file\n\
        line must have at least 1 mean genotype.\n\
        This line is invalid: <{}>",
        line.to_str().unwrap()
      ),
    )));
  }
  Ok(Some(MeanGenoLine {
    snp_id,
    alleles_types,
    mean_genos,
    pos: None,
    chr: None,
  }))
}

/// @brief Attempts to extract snp id from token iterator. If the line
/// corresponding to the token iterator is empty or comment line, returns None.
fn extract_snp_if_not_empty_or_comment<'a>(
  tokens: &mut impl Iterator<Item = &'a [u8]>,
) -> Option<String> {
  use bstr::ByteSlice;
  let first_token = String::from(match tokens.next() {
    Some(token_str) => token_str.to_str_lossy(),
    // Empty line
    None => return None,
  });

  if first_token.starts_with('#') {
    // Comment line
    return None;
  }

  Some(first_token)
}

/// @brief Consumes token from token iterator over BIMBAM Mean Genotype File
/// Format file line. Determines whether the line contains SNP or is comment or
/// empty line.
pub fn check_is_empty_or_comment<'a>(tokens: &mut impl Iterator<Item = &'a [u8]>) -> bool {
  let first_token = match tokens.next() {
    Some(token_str) => token_str,
    // Empty line
    None => return true,
  };

  if *first_token.first().unwrap() == b'#' {
    // Comment line
    return true;
  }

  false
}

// @brief Consume snp id and alleles types from token iterator over BIMBAM Mean
// Genotype File format file line.
fn consume_and_check_id_consume_alleles<'a>(
  expected_snp_id: &str,
  tokens: &mut impl Iterator<Item = &'a [u8]>,
) -> Result<LineType, ParsingError> {
  let first_token = match tokens.next() {
    Some(token_str) => token_str,
    // Empty line
    None => return Ok(LineType::EMPTY),
  };

  if *first_token.first().unwrap() == b'#' {
    // Comment line
    return Ok(LineType::COMMENT);
  }

  let wrong_order_err = || {
    std::io::Error::new(
      std::io::ErrorKind::InvalidInput,
      "Error: BIMBAM Mean Genotype File Format file\
      has different SNP records order than SNP Location File Format file.\n",
    )
  };

  use bstr::ByteSlice;
  if first_token.to_str_lossy() != expected_snp_id {
    return Err(ParsingError::from(wrong_order_err()));
  }

  let no_alleles_types_err = || {
    std::io::Error::new(
      std::io::ErrorKind::InvalidInput,
      "BIMBAM Mean Genotype File Format file\n\
      line must have 2 alleles types.",
    )
  };

  // Alleles types.
  tokens.next().ok_or_else(no_alleles_types_err)?;
  tokens.next().ok_or_else(no_alleles_types_err)?;

  Ok(LineType::SNP)
}

/// @brief Parse snp id and alleles types from token iterator over BIMBAM Mean
/// Genotype File Format file line.
fn parse_id_and_alleles<'a>(
  tokens: &mut impl Iterator<Item = &'a [u8]>,
) -> Result<Option<SNPAndAllelesTypes>, ParsingError> {
  let snp_id = match extract_snp_if_not_empty_or_comment(tokens) {
    Some(snp_id_str) => snp_id_str,
    // Empty line
    None => return Ok(None),
  };

  let no_alleles_types_err = || {
    std::io::Error::new(
      std::io::ErrorKind::InvalidInput,
      "BIMBAM Mean Genotype File Format file\n\
      line must have 2 alleles types.",
    )
  };
  use bstr::ByteSlice;
  let alleles_types = (
    String::from(
      tokens
        .next()
        .ok_or_else(no_alleles_types_err)?
        .to_str_lossy(),
    ),
    String::from(
      tokens
        .next()
        .ok_or_else(no_alleles_types_err)?
        .to_str_lossy(),
    ),
  );

  Ok(Some((snp_id, alleles_types)))
}

/// @brief Parses mean genotypes values into new vector and returns it.
///
/// @in mean_geno_strs mean genotypes in a string form (result of splitting the
/// genotype line with separators).
/// @in ids_num amount of individuals per line.
pub fn parse_mean_geno(
  mean_geno_strs: &[&str],
  ids_num: usize,
) -> Result<Vec<f64>, <f64 as std::str::FromStr>::Err> {
  let mut mean_genos: Vec<f64> = vec![0.0; ids_num];
  parse_mean_geno_into(mean_geno_strs, &mut mean_genos[..])?;
  Ok(mean_genos)
}

/// @brief Parses mean genotypes values into buffer.
///
/// @in mean_geno_strs mean genotypes in a string form (result of splitting the
/// genotype line with separators).
/// @in buf buffer to fill. There must be e
pub fn parse_mean_geno_into(
  mean_geno_strs: &[&str],
  buf: &mut [f64],
) -> Result<(), <f64 as std::str::FromStr>::Err> {
  assert_eq!(mean_geno_strs.len(), buf.len());
  for (geno_str, buf_slot) in mean_geno_strs.iter().zip(buf.iter_mut()) {
    *buf_slot = geno_str.parse::<f64>()?;
  }
  Ok(())
}

// Represents SNP id, SNP physical location (position), and chromosome ID from  SNP Location File Format file.
// Example:
//  This is what this type alias object holds.
//  |     |    |
//  ∨     ∨    ∨
// rs1, 1200,  1
type SnpPosRecord = (String, usize, Option<String>);

/// @brief This iterator parses lines from SNP Location File Format file.
pub struct SnpPosIter {
  file_reader: BufReader<File>,
  read_buf: Vec<u8>,
}

impl SnpPosIter {
  pub fn new(path: String) -> std::io::Result<Self> {
    let file = File::open(path)?;
    Self::with_file(file)
  }

  pub fn with_file(file: File) -> std::io::Result<Self> {
    let file_reader = BufReader::new(file);
    Ok(SnpPosIter {
      file_reader,
      read_buf: Vec::<u8>::new(),
    })
  }

  pub fn count_chromosomes(snp_pos_list: &[SnpPosRecord]) -> usize {
    let mut chr_count = 0;
    let mut i = 0;
    loop {
      if i == snp_pos_list.len() {
        return chr_count;
      }
      chr_count += 1;
      let (_, _, cur_chr) = &snp_pos_list[i];
      loop {
        if i == snp_pos_list.len() {
          return chr_count;
        }
        let (_, _, next_chr) = &snp_pos_list[i];
        if next_chr == cur_chr {
          i += 1;
        } else {
          break;
        }
      }
    }
  }
}

impl Iterator for SnpPosIter {
  type Item = SnpPosRecord;

  fn next(&mut self) -> Option<Self::Item> {
    // Dump old data
    self.read_buf.clear();

    let on_err = |it: &mut Self, e: ParsingError| {
      eprintln!("Failed to parse line. Parsing next line... Error: <{}>", e);
      it.next()
    };

    match self.file_reader.read_until(b'\n', &mut self.read_buf) {
      Ok(read_bytes_n) => {
        if read_bytes_n == 0 {
          // EOF is reached.
          return None;
        }
      }
      Err(e) => return on_err(self, ParsingError::Io(e)),
    };

    let err_msg = |line: &str| {
      format!(
        "Failed to parse the line.\nThis line is invalid: <{}>.\
      \nParsing next line...",
        line
      )
    };

    let mut tokens = tokenize_bimbam_or_rqtl2_line(&self.read_buf);

    let snp_id = match extract_snp_if_not_empty_or_comment(&mut tokens) {
      Some(snp_id_str) => snp_id_str,
      // Empty line,
      None => return self.next(),
    };

    use bstr::ByteSlice;
    let pos = match tokens.next() {
      Some(token) => match token.to_str_lossy().parse::<usize>() {
        Ok(val) => val,
        Err(_) => {
          eprintln!("{}", err_msg(self.read_buf.to_str().unwrap()));
          return self.next();
        }
      },
      None => {
        eprintln!("{}", err_msg(self.read_buf.to_str().unwrap()));
        return self.next();
      }
    };

    // Chromosome field is optional.
    let chr = tokens
      .next()
      .map(bstr::ByteSlice::to_str_lossy)
      .map(String::from);

    Some((snp_id, pos, chr))
  }
}
