use std::env;
use std::fs;
use std::io::prelude::*;
use std::io::SeekFrom;

pub fn create_test_file(name: &str, contents: &str) -> std::io::Result<fs::File> {
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
