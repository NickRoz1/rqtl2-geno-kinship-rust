use std::error;
use std::fmt;
use std::io;
use std::num;

pub type Result<T> = std::result::Result<T, ParsingError>;

// from https://blog.burntsushi.net/rust-error-handling/#composing-option-and-result
//
// We derive `Debug` because all types should probably derive `Debug`.
// This gives us a reasonable human readable description of `CliError` values.
#[derive(Debug)]
pub enum ParsingError {
  Io(io::Error),
  ParseFloat(num::ParseFloatError),
}

impl From<io::Error> for ParsingError {
  fn from(err: io::Error) -> ParsingError {
    ParsingError::Io(err)
  }
}

impl From<num::ParseFloatError> for ParsingError {
  fn from(err: num::ParseFloatError) -> ParsingError {
    ParsingError::ParseFloat(err)
  }
}

impl fmt::Display for ParsingError {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    match *self {
      // Both underlying errors already impl `Display`, so we defer to
      // their implementations.
      ParsingError::Io(ref err) => write!(f, "IO error: {}", err),
      ParsingError::ParseFloat(ref err) => write!(f, "Parse error: {}", err),
    }
  }
}

impl error::Error for ParsingError {
  fn source(&self) -> Option<&(dyn error::Error + 'static)> {
    match *self {
      // N.B. Both of these implicitly cast `err` from their concrete
      // types (either `&io::Error` or `&num::ParseIntError`)
      // to a trait object `&Error`. This works because both error types
      // implement `Error`.
      ParsingError::Io(ref err) => Some(err),
      ParsingError::ParseFloat(ref err) => Some(err),
    }
  }
}
