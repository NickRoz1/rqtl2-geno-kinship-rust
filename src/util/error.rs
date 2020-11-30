use std::error;
use std::fmt;
use std::io;
use std::num;

pub type Result<T> = std::result::Result<T, ProcessingError>;

// from https://blog.burntsushi.net/rust-error-handling/#composing-option-and-result
//
// We derive `Debug` because all types should probably derive `Debug`.
// This gives us a reasonable human readable description of `CliError` values.
#[derive(Debug)]
pub enum ProcessingError {
  Io(io::Error),
  ParseFloat(num::ParseFloatError),
  GPUerror(GPUerror),
}

impl From<io::Error> for ProcessingError {
  fn from(err: io::Error) -> ProcessingError {
    ProcessingError::Io(err)
  }
}

impl From<num::ParseFloatError> for ProcessingError {
  fn from(err: num::ParseFloatError) -> ProcessingError {
    ProcessingError::ParseFloat(err)
  }
}

impl From<GPUerror> for ProcessingError {
  fn from(err: GPUerror) -> ProcessingError {
    ProcessingError::GPUerror(err)
  }
}

impl fmt::Display for ProcessingError {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    match *self {
      // Both underlying errors already impl `Display`, so we defer to
      // their implementations.
      ProcessingError::Io(ref err) => write!(f, "IO error: {}", err),
      ProcessingError::ParseFloat(ref err) => write!(f, "Parse error: {}", err),
      ProcessingError::GPUerror(ref err) => write!(f, "GPU error: {}", err),
    }
  }
}

impl error::Error for ProcessingError {
  fn source(&self) -> Option<&(dyn error::Error + 'static)> {
    match *self {
      // N.B. Both of these implicitly cast `err` from their concrete
      // types (either `&io::Error` or `&num::ParseIntError`)
      // to a trait object `&Error`. This works because both error types
      // implement `Error`.
      ProcessingError::Io(ref err) => Some(err),
      ProcessingError::ParseFloat(ref err) => Some(err),
      ProcessingError::GPUerror(ref err) => Some(err),
    }
  }
}

#[derive(Debug)]
pub enum GPUerror {
  NoDynamicLibrary(libloading::Error),
  NoDevice,
}

impl From<libloading::Error> for GPUerror {
  fn from(err: libloading::Error) -> GPUerror {
    GPUerror::NoDynamicLibrary(err)
  }
}

impl error::Error for GPUerror {
  fn source(&self) -> Option<&(dyn error::Error + 'static)> {
    match *self {
      GPUerror::NoDynamicLibrary(ref err) => Some(err),
      GPUerror::NoDevice => None,
    }
  }
}

impl fmt::Display for GPUerror {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    match *self {
      GPUerror::NoDynamicLibrary(ref err) => write!(
        f,
        "Unable to load CUDA runtime library libcudart.so: {}",
        err
      ),
      GPUerror::NoDevice => write!(f, "No CUDA compatible device was found on this machine."),
    }
  }
}
