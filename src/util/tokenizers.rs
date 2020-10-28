/// @brief BimBam and R/atl2 line splitter. Tokenizes line by [';', ' ', ',', '\t'] delimiters.
/// @note The splitter does not convert byte tokens into strings.
pub struct BimbamAndRqtl2Split<'a> {
  byte_string: &'a [u8],
}

impl<'a> BimbamAndRqtl2Split<'a> {
  pub fn new(byte_str: &'a [u8]) -> Self {
    use bstr::ByteSlice;
    BimbamAndRqtl2Split {
      byte_string: byte_str.trim(),
    }
  }
}

impl<'a> Iterator for BimbamAndRqtl2Split<'a> {
  type Item = &'a [u8];

  fn next(&mut self) -> Option<Self::Item> {
    use bstr::ByteSlice;
    // View will shrink as we find new delimiters.
    let view = &mut self.byte_string;

    let end_idx = match view.find_byteset(b",; \t") {
      Some(idx) => idx,
      // Reached end of string
      None => {
        let view_len = view.len();
        if view_len == 0 {
          return None;
        }
        view_len
      }
    };

    let token = &view[0..end_idx];

    *view = &view[end_idx..];
    if let Some(start_idx) = view.find_not_byteset(b",; \t") {
      *view = &view[start_idx..];
    }
    Some(token)
  }
}

/// @brief Tokenizes BIMBAM or R/QTL2 line (bytestring).
pub fn tokenize_bimbam_or_rqtl2_line(line: &[u8]) -> BimbamAndRqtl2Split {
  BimbamAndRqtl2Split::new(line)
}
