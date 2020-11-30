#include "cublas_v2.h"
#include <cuda_runtime.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h> // for clock_gettime()

void traceMemory() {
  size_t free, total;
  cudaMemGetInfo(&free, &total);
  printf("FREE MEMORY: %lu\nTOTAL MEMORY: %lu\n", free, total);
}

// from
// https://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api
#define gpuErrchk(ans)                                                         \
  { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line,
                      bool abort = true) {
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file,
            line);
    if (abort)
      exit(code);
  }
}

void print_elapsed_time_and_reset(timeval *start, timeval *end,
                                  const char *entity_name) {
  long secs_used, micros_used;

  gettimeofday(end, NULL);
  secs_used = (end->tv_sec - start->tv_sec); // avoid overflow by subtracting
  micros_used = ((secs_used * 1000000) + end->tv_usec) - (start->tv_usec);

  printf("micros_used by %s: %ld\n", entity_name, micros_used);

  *start = *end;
}

void cublas_dsyrk(double *dev_input, double *dev_res_buf, const size_t row_num,
                  const size_t ids_num) {

#ifdef bench
  struct timeval start, end;
  gettimeofday(&start, NULL);
#endif

  cublasStatus_t stat;
  cublasHandle_t handle;
  stat = cublasCreate(&handle);

  if (stat != CUBLAS_STATUS_SUCCESS) {
    printf("CUBLAS initialization failed\n");
    cudaFree(dev_input);
    cudaFree(dev_res_buf);
    exit(1);
  }

  double alpha = 1.0;
  double beta = 0.0;

  stat =
      cublasDsyrk(handle, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, ids_num, row_num,
                  &alpha, dev_input, ids_num, &beta, dev_res_buf, ids_num);

  if (stat != CUBLAS_STATUS_SUCCESS) {
    printf("CUBLAS DSYRK failed\n");
    cudaFree(dev_input);
    cudaFree(dev_res_buf);
    exit(1);
  }

  stat = cublasDestroy(handle);

  if (stat != CUBLAS_STATUS_SUCCESS) {
    printf("CUBLAS deconstruction failed\n");
    cudaFree(dev_input);
    cudaFree(dev_res_buf);
    exit(1);
  }

#ifdef bench
  print_elapsed_time_and_reset(&start, &end, "CUBLAS DSYRK");
#endif
}

extern "C" {
void call_cublas_dsyrk(const double *input, double *res_buf,
                       const long long ids_num, const long long row_count) {

#ifdef bench
  struct timeval start, end;
  gettimeofday(&start, NULL);
#endif
  // CUDA context creation takes a lot of time - several seconds. The first CUDA
  // function called in the application will trigger the context creation, hence
  // measuring elapsed time of this function call is inaccurate. To avoid that,
  // trigger CUDA context creation by this dummy function.
  // Source: forums.developer.nvidia.com/t/cudamalloc-slow/40238/2
  gpuErrchk(cudaFree(0));

#ifdef bench
  print_elapsed_time_and_reset(&start, &end, "CUDA init context.");
#endif
  double *dev_input, *dev_res;

  // Res matrix is square
  const size_t N = ids_num * ids_num;
  const size_t INPUT_ARR_SIZE = ids_num * row_count;

  gpuErrchk(cudaMalloc((void **)&dev_input, INPUT_ARR_SIZE * sizeof(double)));
  gpuErrchk(cudaMalloc((void **)&dev_res, N * sizeof(double)));

#ifdef bench
  print_elapsed_time_and_reset(&start, &end, "CUDA Malloc.");
#endif

  gpuErrchk(cudaMemcpy(dev_input, input, INPUT_ARR_SIZE * sizeof(double),
                       cudaMemcpyHostToDevice));

#ifdef bench
  print_elapsed_time_and_reset(&start, &end, "CUDA Memcpy to DEVICE.");
#endif

  cublas_dsyrk(dev_input, dev_res, row_count, ids_num);

  gpuErrchk(
      cudaMemcpy(res_buf, dev_res, N * sizeof(double), cudaMemcpyDeviceToHost));

#ifdef bench
  print_elapsed_time_and_reset(&start, &end, "CUDA Memcpy to HOST.");
#endif

  gpuErrchk(cudaFree(dev_input));
  gpuErrchk(cudaFree(dev_res));

#ifdef bench
  print_elapsed_time_and_reset(&start, &end, "CUDA Free.");
#endif
}
}

extern "C" {
bool check_gpu_device_availability() {
  int devices = 0;

  cudaError_t err = cudaGetDeviceCount(&devices);
  if (devices > 0 && err == cudaSuccess) {
    return true;
  } else {
    return false;
  }
}
}