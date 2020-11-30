# Installation

1. Install `CUDA Toolkit 8.0` (version deployed on the server). Instructions: https://rodrigodzf.github.io/setup/cuda/2019/04/15/cuda-setup.html.
2. Install gcc `sudo apt install gcc-4.8` (required by `nvcc` to work).
3. Create symlink `gcc` to `gcc-4.8` (`sudo ln -s /usr/bin/gcc-4.8 ./gcc`).
4. In the `build.rs` file set correct paths to `CUDA_Toolkit_8.0`.
5. Substitute the right paths and run `export PATH=path/to/CUDA_Toolkit_8.0/bin:$PATH` and `export CXX=path/to/gcc-4.8-for-cuda/`.
6. Run `cargo build`.

You can enable elapsed time trackers via `cargo build --features "elapsed"`.
