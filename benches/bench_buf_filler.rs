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
  // let mut str_buf: String = String::new();
  // View will shrink as we find new delimiters.
  let mut view = &source[..];

  let mut i = 0;
  loop {
    let end_idx = match view.find_byteset(b",; \t") {
      Some(idx) => idx,
      // Reached end of string
      None => {
        let view_len = view.len();
        if view_len == 0 {
          break;
        }
        view_len
      }
    };
    let parsed_num_str = view[0..end_idx].to_str_lossy();
    dest[i] = parsed_num_str.parse::<f64>().unwrap();
    // str_buf.clear();
    i += 1;
    view = &view[end_idx..];
    match view.find_not_byteset(b",; \t") {
      Some(start_idx) => view = &view[start_idx..],
      None => break,
    }
  }
}

fn criterion_benchmark(c: &mut Criterion) {
  let src_str: Vec<u8> = b"2352352353.1,\t ;;;,,77457455.1\t;           1.1".to_vec();
  // let src_std: String = "2352352353.1,\t ;;;,,77457455.1\t;           1.1".to_owned();
  // let mut str_buf = String::new();
  let mut dest: Vec<f64> = vec![0.0; 3];
  let clean = |v: &mut Vec<f64>| {
    v.clear();
    v.resize(3, 0.0);
  };
  c.bench_function("STD Tokenizer", |b| {
    b.iter(|| std_tokenizer(black_box(&src_str), &mut dest))
  });
  println!("{:?}", dest);
  clean(&mut dest);
  c.bench_function("BSTR custom Tokenizer", |b| {
    b.iter(|| bstr_custom_tokenizer(black_box(&src_str), &mut dest))
  });
  println!("{:?}", dest);
  clean(&mut dest);
  c.bench_function("BSTR Tokenizer", |b| {
    b.iter(|| bstr_tokenizer(black_box(&src_str), &mut dest))
  });
  println!("{:?}", dest);
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
