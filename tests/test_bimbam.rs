mod common;

#[cfg(test)]
mod tests {
  extern crate rqtl2;
  use crate::common::create_test_file;

  fn compare_against_ref(
    mean_geno_line: &rqtl2::bimbam::MeanGenoLine,
    reference: &rqtl2::bimbam::MeanGenoLine,
  ) -> bool {
    mean_geno_line.snp_id == reference.snp_id
      && mean_geno_line.alleles_types == reference.alleles_types
      && mean_geno_line.mean_genos == reference.mean_genos
  }
  #[test]
  fn test_parse_mean_geno() {
    let line = "rs1, A, T, 0.02, 0.80, 1.50";

    let line_tokens: Vec<&str> = line
      .split(rqtl2::bimbam::GENO_SEPARATORS)
      .filter(|s| !s.is_empty())
      .collect();

    use rqtl2::bimbam::parse_mean_geno;

    // skip snp id and alleles types
    let geno_tokens = &line_tokens[3..];
    let parsed_geno = parse_mean_geno(&geno_tokens, 3).unwrap();

    use rqtl2::bimbam::parse_mean_geno_into;

    let mut buf: [f64; 3] = [0.0; 3];
    parse_mean_geno_into(&geno_tokens, &mut buf).unwrap();

    assert_eq!([0.02, 0.80, 1.50], &parsed_geno[..]);
    assert_eq!([0.02, 0.80, 1.50], &buf[..]);
  }

  #[test]
  fn test_parse_geno_line() {
    use rqtl2::bimbam::parse_geno_line;
    use rqtl2::bimbam::MeanGenoLine;

    let line_commas = b"rs1, A, T, 0.02, 0.80, 1.50";
    let line_spaces = b"rs1 A T 0.02 0.80 1.50";
    let line_tabs = b"rs1\tA\tT\t0.02\t0.80\t1.50";

    let comment_line_1 = b"#rs1, A, T, 0.02, 0.80, 1.50";
    let comment_line_2 = b"# rs1, A, T, 0.02, 0.80, 1.50";

    let empty_line = b"";

    // One allele type is missing
    let error_line_1 = b"rs1, A, 0.02, 0.80, 1.50";
    // Mean genos are missing
    let error_line_2 = b"rs1, A, T";

    let reference = MeanGenoLine {
      snp_id: "rs1".to_owned(),
      alleles_types: ("A".to_owned(), "T".to_owned()),
      mean_genos: vec![0.02, 0.80, 1.50],
      pos: None,
      chr: None,
    };

    let check = |line: &[u8]| {
      assert!(compare_against_ref(
        &parse_geno_line(&line).unwrap().unwrap(),
        &reference
      ))
    };

    check(line_commas);
    check(line_spaces);
    check(line_tabs);

    assert!(parse_geno_line(comment_line_1).unwrap().is_none());
    assert!(parse_geno_line(comment_line_2).unwrap().is_none());
    assert!(parse_geno_line(empty_line).unwrap().is_none());

    let expected_err = std::io::ErrorKind::InvalidInput;

    let check_err = |error_line: &[u8]| {
      if let Err(rqtl2::util::error::ProcessingError::Io(err)) = parse_geno_line(error_line) {
        assert_eq!(err.kind(), expected_err);
      }
    };

    check_err(error_line_1);
    check_err(error_line_2);
  }

  #[test]
  fn test_geno_reader() {
    use rqtl2::bimbam::MeanGenoLine;

    let file = create_test_file(
      "test_bimbam_tools_2.txt",
      "rs1, A, T, 0.02, 0.80, 1.50\n# comment line\n\
     \n\
     rs2, G, C, 0.98, 0.04, 1.00",
    )
    .unwrap();

    use rqtl2::bimbam::GenoReader;

    let mut geno_reader = GenoReader::with_file(file).unwrap();

    let expected_iter = vec![
      MeanGenoLine {
        snp_id: "rs1".to_owned(),
        alleles_types: ("A".to_owned(), "T".to_owned()),
        mean_genos: vec![0.02, 0.80, 1.50],
        pos: None,
        chr: None,
      },
      MeanGenoLine {
        snp_id: "rs2".to_owned(),
        alleles_types: ("G".to_owned(), "C".to_owned()),
        mean_genos: vec![0.98, 0.04, 1.00],
        pos: None,
        chr: None,
      },
    ];

    assert_eq!(geno_reader.iter().unwrap().count(), expected_iter.len());
    for (res, exp) in geno_reader.iter().unwrap().zip(expected_iter.iter()) {
      assert!(compare_against_ref(&res, &exp));
    }
  }

  #[test]
  fn test_snp_pos_iter() {
    let file = create_test_file(
      "test_bimbam_tools_3.txt",
      "rs31443144	3010274	1\n# comment line\n\
     \n\
     rs30991578	11793253	1",
    )
    .unwrap();

    let snp_pos_iter = rqtl2::bimbam::SnpPosIter::with_file(file).unwrap();

    let res: Vec<(String, usize, Option<String>)> = snp_pos_iter.collect();
    let expected: Vec<(String, usize, Option<String>)> = vec![
      ("rs31443144".to_owned(), 3010274, Some("1".to_owned())),
      ("rs30991578".to_owned(), 11793253, Some("1".to_owned())),
    ];

    assert_eq!(res.len(), expected.len());
    assert!(res.iter().eq(expected.iter()));
  }

  #[test]
  fn test_loco_kinship() {
    let input_file = create_test_file(
      "test_bimbam_tools_4_1.txt",
      "rs1, A, T, 0.02, 0.80, 1.50\n\
     rs2, G, C, 0.98, 0.04, 1.00\n\
     rs3    A K 1.2 4.5 1.2\n\
     rs4 G,   C, 0.12, 0.25, 1.45\n\
     rs5 G,   C, 0.12, 0.25, 1.45\n\
     rs6 G,   C, 0.12, 0.25, 1.45\n",
    )
    .unwrap();

    let snp_pos_file = create_test_file(
      "test_bimbam_tools_4_2.txt",
      "rs1, 1 0\n# comment line\n\
     \n\
     rs2, 2 0\n\
     rs3    3 0\n\
     rs4 4 9\n\
     rs5 1 9\n\
     rs6 2 9\n",
    )
    .unwrap();

    let mut geno_reader = rqtl2::bimbam::GenoReader::with_file(input_file).unwrap();
    let result = geno_reader.calc_kinship(1, snp_pos_file).unwrap();

    let expected = [
      [
        0.23684444444444439,
        0.1735777777777777,
        -0.4104222222222222,
        0.1735777777777777,
        0.12721111111111105,
        -0.30078888888888883,
        -0.4104222222222222,
        -0.30078888888888883,
        0.7112111111111111,
      ],
      [
        0.6238518518518521,
        -0.8781037037037037,
        0.25425185185185206,
        -0.8781037037037037,
        1.7472740740740738,
        -0.8691703703703705,
        0.25425185185185206,
        -0.8691703703703705,
        0.6149185185185188,
      ],
    ];
    assert_eq!(result.len(), expected.len());

    let epsilon = 0.000_000_000_1;
    for (res, exp) in result.iter().map(|(_, arr)| arr).zip(expected.iter()) {
      assert_eq!(res.len(), exp.len());
      for (e_1, e_2) in res.iter().zip(exp.iter()) {
        assert!((e_1 - e_2) * (e_1 - e_2) < epsilon * epsilon);
      }
    }
  }
}
