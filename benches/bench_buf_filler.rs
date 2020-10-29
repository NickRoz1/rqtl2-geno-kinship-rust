use bstr::ByteSlice;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

pub const GENO_SEPARATORS: &[char] = &['\t', ' ', ',', ';'];

fn bstr_tokenizer(source: &Vec<u8>, dest: &mut Vec<f64>) -> () {
  for (token, slot) in source.words().zip(dest.iter_mut()) {
    *slot = token.parse::<f64>().unwrap();
  }
}

fn std_tokenizer(source: &Vec<u8>, dest: &mut Vec<f64>) -> () {
  let parsed_line = source.to_str_lossy();
  let filter_fn = |token: &&str| !token.is_empty();

  let tokens = parsed_line
    .split(&GENO_SEPARATORS[..])
    // Each named function has unique signature, converting to function pointer
    // via <as> keyword.
    .filter(filter_fn as fn(&&str) -> bool);

  // eprintln!("\n{}\n", str_buf);
  for (token, slot) in tokens.zip(dest.iter_mut()) {
    // eprintln!("\n{}\n", token);
    *slot = token.parse::<f64>().unwrap();
  }
}

fn bstr_custom_tokenizer(source: &Vec<u8>, dest: &mut Vec<f64>) -> () {
  extern crate rqtl2;
  use rqtl2::util::tokenizers::tokenize_bimbam_or_rqtl2_line;

  let tokens = tokenize_bimbam_or_rqtl2_line(&source);
  for (slot, token) in dest
    .iter_mut()
    .zip(tokens.map(bstr::ByteSlice::to_str_lossy))
  {
    *slot = token.parse::<f64>().unwrap();
  }
}

fn bstr_fields_tokenizer(source: &Vec<u8>, dest: &mut Vec<f64>) -> () {
  let tokens = source.fields_with(|c| GENO_SEPARATORS.iter().any(|delim| c == *delim));
  for (slot, token) in dest
    .iter_mut()
    .zip(tokens.map(bstr::ByteSlice::to_str_lossy))
  {
    *slot = token.parse::<f64>().unwrap();
  }
}

fn criterion_benchmark(c: &mut Criterion) {
  let src_str: Vec<u8> = b"2352352353.1,\t ;;;,,77457455.1\t;           1.1".to_vec();
  let mut dest: Vec<f64> = vec![0.0; 3];
  let clean = |v: &mut Vec<f64>| {
    v.clear();
    v.resize(3, 0.0);
  };

  let mut group = c.benchmark_group("longer-warmup");
  group.warm_up_time(std::time::Duration::new(10, 0));

  group.bench_function("STD Tokenizer", |b| {
    b.iter(|| std_tokenizer(black_box(&src_str), &mut dest))
  });
  println!("{:?}", dest);
  clean(&mut dest);
  group.bench_function("BSTR custom Tokenizer", |b| {
    b.iter(|| bstr_custom_tokenizer(black_box(&src_str), &mut dest))
  });
  println!("{:?}", dest);
  clean(&mut dest);
  group.bench_function("BSTR Tokenizer", |b| {
    b.iter(|| bstr_tokenizer(black_box(&src_str), &mut dest))
  });
  println!("{:?}", dest);
  clean(&mut dest);
  group.bench_function("BSTR fields Tokenizer", |b| {
    b.iter(|| bstr_fields_tokenizer(black_box(&src_str), &mut dest))
  });
  println!("{:?}", dest);
  group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
