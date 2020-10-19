// This module hosts Kinship calculator implementation.
// ----------------------------------------------------------------------------

// @brief Unit of work when processing geno file in parallel.
// @note Gets dispatched to a thread.
pub struct WorkUnit {
  pub sender: std::sync::mpsc::Sender<WorkUnit>,
  pub input_buf: Vec<f64>,
  pub result_buf: Vec<f64>,
  pub chr_num: usize,
}

/// @brief Calculates Kinship matrix in parallel. Uses <processor> delegate to
/// customize data parsing and merging behavior.
/// @note The purpose of calculation of Kinship matrix in batches is to not
/// load a complete SNPs dataset in memory.
///
/// Given a genotype data file stripped to just containing SNPs called G, a
/// Kinship matrix is a matrix product G.T (transposed) * G.
///
/// The matrix times its transpose forms a symmetrical matrix, thus it is
/// not needed to calculate the full matrix, just a triangular part of it.
///
/// Matrix multiplication of transposed matrix by itself (or vice versa) can
/// be performed via addition of non intersecting (1 column with 1 row, 2
/// column with 2 row) products of columns of matrix's transpose by rows of
/// matrix.
///
/// [[1, 2],          [[1, 4],     [[1],  * [[1],[4]]     [[2],   * [[2],[5]]
///  [4, 5]]     *     [2, 5]]  =   [4]]               +   [5]]                =
///
///  = [[1,  4],        [[4,  10],      [[5,  14],
///     [4, 16]]    +    [10, 25]]  =    [14, 41]]
///
/// When reading a row of genotype matrix, the column of transposed matrix
/// can be obtained by transposing this row. Thus, it is possible to
/// represent genotype matrix transpose times itself (G.T * G) product, as
/// many products of it's columns by rows, and column can be obtained once
/// the row is read.
///
/// This is how Kinship matrix can be calculated: read batch from file, copy
/// and transpose it, multiply transposed matrix by non transposed, add to
/// result Kinship matrix. After each batch will be processed and added to
/// the result Kinship matrix it will contain the full Kinship matrix.
///
/// Actual algorithm does not involve matrix copying and transposing and
/// instead just manipulates matrix indices calculation to achieve same
/// result.
///
/// Since processing of one batch does not depend on the others, the process
/// of Kinship matrix calculation can be parallelized: each logical thread
/// gets 2 buffers, first one contains read rows, and a second one stores the
/// result of batch multiplication, it is done to not block a shared Kinship
/// matrix buffer while the calculation is in process.
//
/// Main thread has several channels connecting him to each work thread,
/// this is done to emulate threadpool. Threads connected to the main thread
/// via mpsc queue (which is also implemented via channels). Notice, that
/// each thread connected to the main thread with two channels.
///
/// The MAIN -> WORKER channel messages a work unit to the worker thread
/// from the main thread (dispatches work, just like a threadpool
/// implementation would). WORKER -> MAIN channel messages work unit
/// (processed) back to the main thread, where it gets merged into the
/// common kinship matrix. The WorkUnit contains the sender field, which
/// contains a channel sender. This way the main threads knows which thread
/// completed its work and can be loaded with the new batch (dispatched
/// through this sender).
///
/// When the main thread hits EOF in the geno file, it starts dropping work
/// units. When the work unit is dropped, the unique sender inside it
/// associated with particular thread is also dropped. When the unique
/// sender is dropped, the consumer calling recv errors. In a thread, it
/// causes it to stop waiting for new tasks and terminate. Once all worker
/// threads finished executing, there is no more senders for the main thread
/// consumer (which receives, merges and dispatches the work units). This
/// causes recv to error, terminating the dispatching loop, finishing the
/// geno file processing.
pub fn calc_kinship_parallel<
  F: FnMut(&mut WorkUnit) -> Result<bool, crate::util::error::ParsingError>,
>(
  mut processor: F,
  read_buf_size: usize,
  ids_num: usize,
) -> Result<(), crate::util::error::ParsingError> {
  // For each physical thread a buffers will be created.
  let buf_num = num_cpus::get();
  use std::sync::mpsc::channel;
  use std::thread;

  // Channel which returns results of calculation to the main thread
  // to merge them with the end Kinship matrix.
  let (worker_thread_sender, main_thread_consumer) = channel::<WorkUnit>();

  let mut threads = Vec::<thread::JoinHandle<()>>::new();

  for _ in 0..buf_num {
    // Channel which sends work unit to the worker thread for processing.
    let (main_thread_sender, worker_thread_consumer) = channel::<WorkUnit>();
    let worker_thread_sender_clone = worker_thread_sender.clone();
    // Prefill the queue.
    let work_unit = WorkUnit {
      sender: main_thread_sender,
      input_buf: vec![0.0; read_buf_size],
      result_buf: vec![0.0; ids_num * ids_num],
      chr_num: 0,
    };
    worker_thread_sender.send(work_unit).unwrap();

    threads.push(std::thread::spawn(move || {
      while let Ok(mut work_unit) = worker_thread_consumer.recv() {
        calc_partial_kinship(&work_unit.input_buf, &mut work_unit.result_buf, ids_num);
        worker_thread_sender_clone.send(work_unit).unwrap();
      }
      // Worker thread terminates.
    }));
  }

  // Drop worker thread sender, so the main thread consumer automatically
  // terminates when the threads are destroyed (when threads are destroyed,
  // there no more senders, so the consumer fails on recv call.)
  drop(worker_thread_sender);

  while let Ok(mut work_unit) = main_thread_consumer.recv() {
    let is_eof_reached = processor(&mut work_unit)?;
    // Continue consuming work unit queue until every thread is finished.
    if is_eof_reached {
      continue;
    }
    // Work unit contains sender tied with receiver in thread to allow for
    // dispatching work to worker threads.
    work_unit.sender.clone().send(work_unit).unwrap();
  }
  threads.into_iter().for_each(|thread| {
    thread
      .join()
      .expect("The thread creating or execution failed!")
  });
  Ok(())
}

fn calc_partial_kinship(snps: &Vec<f64>, partial_matrix: &mut Vec<f64>, ids_num: usize) -> () {
  let n = ids_num;
  let k = snps.len() / n;
  // Algorithm from BLAS dsyrk:
  // http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_gae0ba56279ae3fa27c75fefbc4cc73ddf.html#gae0ba56279ae3fa27c75fefbc4cc73ddf
  //
  // The BLAS Fortran stores array in a column-major format, but the R/qtl2
  // genotype data stored in a row-major format, so this algorithm corresponds
  // to the branch for non transposed, lower triangular part version, however
  // in fact it performs transposed, upper triangular part multiplication
  // (G.T*G).
  //
  // This algorithm branch (Lower, Non transposed) chosen based on CBLAS
  // http://www.netlib.org/blas/blast-forum/cblas.tgz code for dsyrk
  // (cblas_dsyrk.c), which transforms options (Upper, Transposed) to these
  // arguments when called for row-major matrixes.
  //
  // When the matrix stored in row-major way read in column major way,
  // obtained data is a transpose of this matrix:
  // https://en.wikipedia.org/wiki/Row-_and_column-major_order#Transposition
  //
  // Since this is an exact copy of Fortran code, and Fortran utilizes double
  // index (i,j) to operate over single dimension array (which represents 2D
  // array), the code below performs index flattening for column-major
  // storages exactly how Fortran does. Normally, to flatten index in
  // row-major languages we will multiply row index i by row width and add
  // column index j, here, since this is a direct copy of Fortran code which
  // is a colum-major language, we flatten it as column index j *
  // column height + row index i.
  for j in 0..n {
    for l in 0..k {
      for i in j..n {
        partial_matrix[j * ids_num + i] += snps[l * ids_num + j] * snps[l * ids_num + i];
      }
    }
  }
}
// std::io::Lines<std::io::BufReader<&mut std::fs::File>>
/// @brief Fills preallocated buffer with parsed values.
pub fn fill_buffer<
  F: FnMut(&String, &mut [f64]) -> Result<(), crate::util::error::ParsingError>,
>(
  fill_buf: &mut std::slice::ChunksMut<f64>,
  lines_iter: &mut dyn Iterator<Item = std::io::Result<String>>,
  mut parser: F,
) -> Result<usize, crate::util::error::ParsingError> {
  let mut parsed_lines_counter: usize = 0;
  for (line_slice, snp_line) in fill_buf.zip(lines_iter) {
    parser(&snp_line?, line_slice)?;
    parsed_lines_counter += 1;
  }
  Ok(parsed_lines_counter)
}

/// @brief Mirrors and scales Kinship matrix, since only the upper part was
/// calculated (the Kinship matrix is symmetrical because it's formed from it's
/// transpose times itself).
pub fn mirror_and_normalize_kinship(
  common_kinship_matrix: &mut [f64],
  ids_num: usize,
  scale: f64,
) -> () {
  for i in 0..ids_num {
    let row_length = ids_num;
    for j in 0..i + 1 {
      common_kinship_matrix[j * row_length + i] *= scale as f64;
      common_kinship_matrix[i * row_length + j] = common_kinship_matrix[j * row_length + i];
    }
  }
}
