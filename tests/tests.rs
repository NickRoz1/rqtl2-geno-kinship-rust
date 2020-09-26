#[cfg(test)]
mod tests {
  extern crate rqtl2;
  use std::collections::HashMap;
  use std::env;
  use std::fs;
  use std::io::prelude::*;
  use std::io::BufReader;
  use std::io::SeekFrom;

  fn create_test_file(name: &str, contents: &str) -> std::io::Result<fs::File> {
    let mut path = env::temp_dir();
    path.push(name);

    let mut f = std::fs::OpenOptions::new()
      .create(true)
      .write(true)
      .read(true)
      .open(&path)?;

    // @note Remove old test data (in case when new test data was provided).
    f.set_len(0)?;
    f.write(contents.as_bytes())?;
    f.seek(SeekFrom::Start(0))?;
    Ok(f)
  }
  #[test]
  fn parsers() {
    use rqtl2::util::GenoParser;
    let mut f = create_test_file(
      "test_geno_parsers_1.txt",
      "#test file\n#comment\nmarker	10	12	38	39\nrs31443144	ABAH",
    )
    .expect("Failed file creation");
    let comments = rqtl2::util::parse_comments(&mut f).expect("Parsing failed");
    let markers = rqtl2::util::parse_markers(&mut f).expect("Parsing failed");
    assert_eq!(["test file", "comment"], &comments[..]);
    assert_eq!(["10", "12", "38", "39"], &markers[..]);
    let consumed_comments = GenoParser::consume_comments(&mut BufReader::new(
      f.try_clone().expect("Failed cloning file handler"),
    ))
    .expect("Parsing failed");
    let consumed_markers = GenoParser::consume_markers(&mut BufReader::new(
      f.try_clone().expect("Failed cloning file handler"),
    ))
    .expect("Parsing failed");
    assert_eq!(["test file", "comment"], &consumed_comments[..]);
    assert_eq!(["10", "12", "38", "39"], &consumed_markers[..]);
    let mut snp_line = String::new();
    f.read_to_string(&mut snp_line).expect("Parsing failed");
    assert_eq!(snp_line.len(), b"rs31443144	ABAH".len());
  }

  #[test]
  fn read_snps() {
    let path = "test_geno_parsers_2.txt";
    let mut f = create_test_file(
      &path,
      "#test file\n#comment\nmarker	10	12	38	39\nrs31443144	ABAH\nrs31443154	ABHH\nrs31443157	BH--",
    )
    .expect("Failed to create test file.");

    let mut hab_mapper = HashMap::new();
    use std::f64::NAN;

    hab_mapper.insert('A', 0.0);
    hab_mapper.insert('H', 0.5);
    hab_mapper.insert('B', 1.0);
    hab_mapper.insert('-', NAN);

    let geno = rqtl2::util::parse_geno(&mut f, &hab_mapper).unwrap();
    let check_snps = |geno: &Vec<(String, Vec<f64>)>| {
      let ids = geno
        .iter()
        .map(|record| record.0.clone())
        .collect::<Vec<String>>();
      let snps = geno
        .iter()
        .map(|record| record.1.clone())
        .collect::<Vec<Vec<f64>>>();
      assert_eq!(["rs31443144", "rs31443154", "rs31443157"], &ids[..]);
      let check = |test_snps, snps: &Vec<f64>| assert_eq!(test_snps, &snps[..]);
      check([0.0, 1.0, 0.0, 0.5], &snps[0]);
      check([0.0, 1.0, 0.5, 0.5], &snps[1]);
      assert_eq!([1.0, 0.5], &snps[2][..2]);
      assert_eq!([true, true], [snps[2][2].is_nan(), snps[2][3].is_nan()]);
    };

    f.seek(SeekFrom::Start(0))
      .expect("Failed to rewind file cursor.");
    let geno_from_struct = rqtl2::util::GenoParser::new_with_file(f, hab_mapper)
      .expect("Failed geno parser creation")
      .read_all()
      .expect("Failed parsing snps");

    check_snps(&geno);
    check_snps(&geno_from_struct);
  }

  #[test]
  fn snp_iter() {
    let f = create_test_file(
      "test_geno_parsers_3.txt",
      "#test file\n#comment\nmarker	10	12	38	39\nrs31443144	ABAH\nrs31443154	ABHH",
    )
    .expect("Failed to create test file.");

    let mut hab_mapper = HashMap::new();
    use std::f64::NAN;

    hab_mapper.insert('A', 0.0);
    hab_mapper.insert('H', 0.5);
    hab_mapper.insert('B', 1.0);
    hab_mapper.insert('-', NAN);

    let mut geno_parser =
      rqtl2::util::GenoParser::new_with_file(f, hab_mapper).expect("Failed to create GenoParser");

    let test_recs: Vec<(String, Vec<f64>)> = vec![
      (String::from("rs31443144"), vec![0.0, 1.0, 0.0, 0.5]),
      (String::from("rs31443154"), vec![0.0, 1.0, 0.5, 0.5]),
    ];

    for (rec, test_rec) in geno_parser.iter().unwrap().zip(test_recs.iter()) {
      assert_eq!(&rec, test_rec);
    }
  }

  #[test]
  fn calc_kinsip() {
    let f = create_test_file(
      "test_geno_parsers_4.txt",
      "#test file\n#comment\nmarker	10	12	38\nrs31443144	ABH\nrs31443154	ABH\nrs31443144	BBA",
    )
    .expect("Failed to create test file.");

    let mut hab_mapper = HashMap::new();
    use std::f64::NAN;

    hab_mapper.insert('A', 0.0);
    hab_mapper.insert('H', 0.5);
    hab_mapper.insert('B', 1.0);
    hab_mapper.insert('-', NAN);

    let expected_kinship_matrix: Vec<f64> = vec![1.0, 1.0, 0.0, 1.0, 3.0, 1.0, 0.0, 1.0, 0.5];

    let mut geno_parser = rqtl2::util::GenoParser::new_with_file(f, hab_mapper.clone())
      .expect("Failed to create GenoParser");

    let mut matr = geno_parser.calc_kinship(10).unwrap();
    for e in matr.iter_mut() {
      *e *= 3.0;
    }
    assert_eq!(matr, expected_kinship_matrix);
  }
}
