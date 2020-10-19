use crate::util::error;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Seek;
use std::io::SeekFrom;

// The input files can be comma-delimited, space-delimited, tab-delimited,
// and semi-colon de-limited, or mixed use of those.  (i.e.  entries can be
// separated by commas, spaces, semi-colons or tabs).
pub const GENO_SEPARATORS: &[char] = &['\t', ' ', ',', ';'];

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
      file_reader: file_reader,
      ids_num: ids_num,
    })
  }

  fn parse_into(source: &String, dest: &mut [f64]) -> Result<(), crate::util::error::ParsingError> {
    let mut tokens = tokenize_bimbam_line(source);
    match consume_id_and_alleles(&mut tokens)? {
      Some(s) => s,
      None => panic!(
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
    for ((i, buf_slot), token) in dest.iter_mut().enumerate().zip(tokens.by_ref()) {
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
      Err(std::io::Error::new(
        std::io::ErrorKind::InvalidInput,
        format!(
          "Error: buffer length is not equal to the amount of tokens.\n\
        This line is invalid: <{}>",
          source
        ),
      ))?
    }

    Ok(())
  }

  fn extract_geno_line(
    pos_snp_id: &String,
    line_iter: &mut dyn Iterator<Item = std::io::Result<String>>,
  ) -> std::io::Result<String> {
    let not_enough_snp_records_err = std::io::Error::new(
      std::io::ErrorKind::InvalidInput,
      "Error: BIMBAM Mean Genotype File Format file\
        has fewer records than SNP Location File Format file does.\n",
    );

    let line = match line_iter.next() {
      Some(Err(e)) => return Err(e),
      Some(Ok(s)) => s,
      None => return Err(not_enough_snp_records_err),
    };

    let mut tokens = tokenize_bimbam_line(&line);
    let no_snp_id_err = std::io::Error::new(
      std::io::ErrorKind::InvalidData,
      format!(
        "Error: Encountered empty or comment line in BIMBAM Mean Genotype File \
          Format file.\nLine: <{}>",
        line
      ),
    );

    let snp_id: String = match check_if_empty_or_comment(&mut tokens) {
      Some(s) => s,
      None => return Err(no_snp_id_err),
    };

    let wrong_order_err = std::io::Error::new(
      std::io::ErrorKind::InvalidInput,
      "Error: BIMBAM Mean Genotype File Format file\
        has different SNP records order than SNP Location File Format file.\n",
    );
    if snp_id != *pos_snp_id {
      return Err(wrong_order_err);
    }

    Ok(line)
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

    let snp_positions: Vec<(String, usize, String)> =
      SnpPosIter::with_file(snp_pos_file)?.collect();

    if snp_positions.len() == 0 {
      panic!("File is empty.");
    }

    let chromosome_count = SnpPosIter::count_chromosomes(&snp_positions);

    // Kinship matrices for each chromosome. Kinship matrix is square.
    let mut chrs_kinship_matrices: Vec<Vec<f64>> =
      vec![vec![0.0; self.ids_num * self.ids_num]; chromosome_count];

    let ids_num = self.ids_num;
    // Maximum amount of SNPs that will be parsed and processed on each iteration.
    let buf_size = self.ids_num * batch_size;

    let file_reader = self.file_reader.get_mut();
    let mut line_iter = BufReader::new(file_reader).lines();

    // Amount of SNPs parsed for each chromosome. Since it is assumed that the
    // SNPs are ordered by chromosomes, it's simply a sequential array of values
    // (no need for String keys).
    let mut chrs_parsed_snps_nums = Vec::<usize>::new();
    // Maps chromosomes indices to their names.
    let mut chrs_names = Vec::<String>::new();

    use crate::util::kinship::WorkUnit;
    // Tracks the chromosome currently being processed.
    let mut cur_chromosome = snp_positions[0].2.clone();
    // Tracks amount of snps on current chromosome.
    let mut cur_chr_total_snps_parsed: usize = 0;
    let mut cur_snp_idx = 0;

    let kinship_merger_delegate =
      |work_unit: &mut WorkUnit| -> Result<bool, crate::util::error::ParsingError> {
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

        let mut get_next_batch = || -> Option<(Vec<&String>, Option<String>)> {
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
          return Some((snps_ids, None));
        };

        // Batch of SNP ids on a single chromosome to be processed.
        let (snp_batch_on_chr, new_chr) = match get_next_batch() {
          Some(tup) => tup,
          // Reached EOF. Resize buffer to discard data from previous iterations
          // which was not overwritten because there is not enough lines to fill the
          // whole buffer.
          None => {
            if cur_chr_total_snps_parsed != 0 {
              chrs_parsed_snps_nums.push(cur_chr_total_snps_parsed.clone());
              chrs_names.push(cur_chromosome.clone());
              cur_chr_total_snps_parsed = 0;
              cur_chromosome.clear();
            }

            if line_iter.next().is_none() == false {
              Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "Error: BIMBAM Mean Genotype File Format file \
                has more records than SNP Location File Format file does.\n",
              ))?;
            }
            return Ok(true);
          }
        };

        // Parses line from Lines iterator over <Mean Genotype File Format> file,
        // checks it and returns String.
        let check_and_extract_line =
          |pos_snp_id: &&String| Self::extract_geno_line(*pos_snp_id, &mut line_iter);

        // This iter yields only the SNPs located on the same chromosome.
        let mut wrapped_line_iter = snp_batch_on_chr.iter().map(check_and_extract_line);

        let parser = |source: &String, dest: &mut [f64]| Self::parse_into(source, dest);

        // Parse new batch (load work unit).
        let line_count = match crate::util::kinship::fill_buffer(
          &mut work_unit.input_buf.chunks_mut(ids_num),
          &mut wrapped_line_iter,
          parser,
        )? {
          n => {
            cur_chr_total_snps_parsed += n;
            n
          }
        };

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

        if !new_chr.is_none() {
          chrs_parsed_snps_nums.push(cur_chr_total_snps_parsed.clone());
          chrs_names.push(cur_chromosome.clone());
          cur_chromosome = new_chr.unwrap();
          cur_chr_total_snps_parsed = 0;
        }

        return Ok(false);
      };

    use crate::util::kinship::calc_kinship_parallel;
    calc_kinship_parallel(kinship_merger_delegate, buf_size, self.ids_num)
      .expect("Parallel processing failed.");

    let total: usize = chrs_parsed_snps_nums.iter().sum();

    for (chrs_kinship_matrix, chrs_parsed_snps_num) in chrs_kinship_matrices
      .iter_mut()
      .zip(chrs_parsed_snps_nums.into_iter())
    {
      crate::util::kinship::mirror_and_normalize_kinship(
        &mut chrs_kinship_matrix[..],
        self.ids_num,
        1.0 / ((total - chrs_parsed_snps_num) as f64),
      );
    }

    self.file_reader.seek(SeekFrom::Start(0))?;

    Ok(
      chrs_names
        .into_iter()
        .zip(chrs_kinship_matrices.into_iter())
        .collect(),
    )
  }

  pub fn iter(&mut self) -> std::io::Result<GenoReaderIter> {
    self.file_reader.seek(SeekFrom::Start(0))?;
    Ok(GenoReaderIter {
      lines_reader: (&mut self.file_reader).lines(),
      ids_num: self.ids_num,
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
  lines_reader: std::io::Lines<&'a mut BufReader<File>>,
  ids_num: usize,
}

impl<'a> Iterator for GenoReaderIter<'a> {
  type Item = MeanGenoLine;

  fn next(&mut self) -> Option<Self::Item> {
    let on_err = |it: &mut Self, e: error::ParsingError| {
      eprintln!("Failed to parse line. Parsing next line... Error: <{}>", e);
      return it.next();
    };

    let line = match self.lines_reader.next()? {
      Ok(string) => string,
      Err(e) => return on_err(self, error::ParsingError::Io(e)),
    };
    match parse_geno_line(&line) {
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

/// @brief Tokenizes BIMBAM line (separates original line on tokens by the
/// separators from const SEPARATORS array)
pub fn tokenize_bimbam_line<'a>(
  line: &'a str,
) -> std::iter::Filter<std::str::Split<&'a [char]>, fn(&&str) -> bool> {
  // If a string contains multiple contiguous separators, you will end up with
  // empty strings in the output:
  // https://doc.rust-lang.org/std/primitive.str.html#method.split
  //
  // To specify the type of the iterator it should have a complete type. Since
  //  unboxed closures - a filter predicate in this case - are anonymous, this
  //  is not possible without explicitly defining function.
  //  https://stackoverflow.com/questions/30641167/figuring-out-return-type-of-closure
  let filter_fn = |token: &&str| !token.is_empty();

  line
    .split(&GENO_SEPARATORS[..])
    // Each named function has unique signature, converting to function pointer
    // via <as> keyword.
    .filter(filter_fn as fn(&&str) -> bool)
}

/// @brief Parses line from BIMBAM Mean Genotype File Format file.
/// If line is empty or contains a comment, None is returned.
///
/// @in line full line from file, e.g. "rs1, A, T, 0.02"
///
/// @note https://www.haplotype.org/download/bimbam-manual.pdf
pub fn parse_geno_line(line: &str) -> Result<Option<MeanGenoLine>, error::ParsingError> {
  let mut tokens = tokenize_bimbam_line(line);

  let (snp_id, alleles_types) = match consume_id_and_alleles(&mut tokens)? {
    Some(id_and_alleles) => id_and_alleles,
    None => return Ok(None),
  };

  let mean_genos: Vec<f64> = tokens
    .map(|mean_geno| mean_geno.parse::<f64>())
    .collect::<Result<Vec<f64>, _>>()?;
  if mean_genos.len() < 1 {
    return Err(std::io::Error::new(
      std::io::ErrorKind::InvalidInput,
      format!(
        "BIMBAM Mean Genotype File Format file\n\
        line must have at least 1 mean genotype.\n\
        This line is invalid: <{}>",
        line
      ),
    ))?;
  }
  Ok(Some(MeanGenoLine {
    snp_id: snp_id.to_owned(),
    alleles_types: alleles_types,
    mean_genos: mean_genos,
    pos: None,
    chr: None,
  }))
}

/// @brief Checks if line is comment and if it's empty. Returns Some(snp_id) if not.
fn check_if_empty_or_comment(
  tokens: &mut std::iter::Filter<std::str::Split<&[char]>, fn(&&str) -> bool>,
) -> Option<String> {
  let first_token = String::from(match tokens.next() {
    Some(token_str) => token_str,
    // Empty line
    None => return None,
  });

  if first_token.chars().nth(0).unwrap() == '#' {
    // Comment line
    return None;
  }

  Some(first_token)
}

/// @brief Consume snp id and alleles types from token iterator over BIMBAM Mean
/// Genotype File Format file line.
fn consume_id_and_alleles(
  tokens: &mut std::iter::Filter<std::str::Split<&[char]>, fn(&&str) -> bool>,
) -> Result<Option<(String, (String, String))>, error::ParsingError> {
  let snp_id = match check_if_empty_or_comment(tokens) {
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
  let alleles_types = (
    String::from(tokens.next().ok_or_else(no_alleles_types_err)?),
    String::from(tokens.next().ok_or_else(no_alleles_types_err)?),
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

/// @brief This iterator parses lines from SNP Location File Format file.
pub struct SnpPosIter {
  lines_reader: std::io::Lines<BufReader<File>>,
}

impl SnpPosIter {
  pub fn new(path: String) -> std::io::Result<Self> {
    let file = File::open(path)?;
    Self::with_file(file)
  }

  pub fn with_file(file: File) -> std::io::Result<Self> {
    let file_reader = BufReader::new(file);
    Ok(SnpPosIter {
      lines_reader: file_reader.lines(),
    })
  }

  pub fn count_chromosomes(snp_pos_list: &Vec<(String, usize, String)>) -> usize {
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
  type Item = (String, usize, String);

  fn next(&mut self) -> Option<Self::Item> {
    let on_err = |it: &mut Self, e: error::ParsingError| {
      eprintln!("Failed to parse line. Parsing next line... Error: <{}>", e);
      return it.next();
    };

    let line = match self.lines_reader.next()? {
      Ok(string) => string,
      Err(e) => return on_err(self, error::ParsingError::Io(e)),
    };

    let err_msg = format!(
      "Failed to parse the line.\nThis line is invalid: <{}>.\
      \nParsing next line...",
      line
    );

    let mut tokens = tokenize_bimbam_line(&line);

    let snp_id = match check_if_empty_or_comment(&mut tokens) {
      Some(snp_id_str) => snp_id_str,
      // Empty line,
      None => return self.next(),
    };

    let pos = match tokens.next() {
      Some(token) => match token.parse::<usize>() {
        Ok(val) => val,
        Err(_) => {
          eprintln!("{}", err_msg);
          return self.next();
        }
      },
      None => {
        eprintln!("{}", err_msg);
        return self.next();
      }
    };

    let chr = match tokens.next() {
      Some(token) => String::from(token),
      None => {
        eprintln!("{}", err_msg);
        return self.next();
      }
    };

    Some((snp_id, pos, chr))
  }
}
