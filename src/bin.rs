fn main() {
  test_bimbam_kinship();
  // test_rqtl_kinship();
}

fn test_bimbam_kinship() -> () {
  let open_file_in_cur_dir = |name| {
    let mut path_input_file = std::env::current_dir().unwrap();
    path_input_file.push(name);
    std::fs::OpenOptions::new()
      .read(true)
      .open(path_input_file)
      .unwrap()
  };

  let mut geno_reader =
    rqtl2::bimbam::GenoReader::with_file(open_file_in_cur_dir("src/BXD_geno.txt")).unwrap();
  let test_results = geno_reader
    .calc_kinship(1000, open_file_in_cur_dir("src/BXD_snps.txt"))
    .unwrap();

  for (chr_name, res_str_repr) in test_results {
    let strings: Vec<String> = res_str_repr.iter().map(|n| n.to_string()).collect();
    let mut output_file =
      std::fs::File::create(format!("data/LOCO.Kinship.{}.txt", chr_name)).unwrap();
    let mut line = String::new();
    for val_group in strings.chunks(198) {
      let mut l = val_group.join("\t");
      l.push('\n');
      line.push_str(&l);
    }
    use std::io::prelude::*;
    writeln!(output_file, "{}", line).unwrap();
  }
}

// fn test_rqtl_kinship() -> () {
//   use std::collections::HashMap;
//   let mut hab_mapper = HashMap::new();
//   use std::f64::NAN;

//   hab_mapper.insert('A', 0.0);
//   hab_mapper.insert('H', 0.5);
//   hab_mapper.insert('B', 1.0);
//   hab_mapper.insert('-', NAN);
//   use std::env;
//   let name = "src/21487-pheno_geno.txt";
//   let mut path = env::current_dir().unwrap();
//   path.push(name);

//   let f2 = std::fs::OpenOptions::new()
//     .create(true)
//     .write(true)
//     .read(true)
//     .open(&path)
//     .unwrap();

//   let mut geno_parser_2 = rqtl2::rqtl2::GenoParser::new_with_file(f2, hab_mapper.clone())
//     .expect("Failed to create GenoParser");

//   println!("{:?}", geno_parser_2.calc_kinship(1000));
// }
