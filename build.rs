extern crate cc;
fn main() {
  #[cfg(feature = "elapsed")]
  let conditional_flag = "-Dbench";
  #[cfg(not(feature = "elapsed"))]
  let conditional_flag = "-DDummyVariableFromRustWorkAround";
  cc::Build::new()
    .cuda(true)
    .flag("-cudart=shared")
    .flag("-lcublas")
    .flag("-lineinfo")
    .flag(conditional_flag)
    .flag("-Wno-deprecated-gpu-targets")
    .files(&["./src/util/cublas_dsyrk.cu"])
    .compile("cublas_mult.a");
  println!("cargo:rustc-link-search=native=/CUDA_Toolkit_8.0/lib64/");
  println!("cargo:rustc-link-search=/CUDA_Toolkit_8.0/lib64/");
  println!("cargo:rustc-env=LD_LIBRARY_PATH=/CUDA_Toolkit_8.0/lib64/");
  println!("cargo:rustc-link-lib=dylib=cublas");
  println!("cargo:rustc-link-lib=dylib=cudart");
}
