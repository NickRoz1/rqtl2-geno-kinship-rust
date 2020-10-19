/// @brief Tools for working with RQTL2 and BIMBAM format.
///
/// From
/// https://kbroman.org/qtl2/assets/vignettes/user_guide.html#Data_file_format:
///
/// The input data file formats for R/qtl cannot handle complex crosses, and so
/// for R/qtl2, we have defined a new format for the data files. We’ll describe
/// it here briefly; for details, see the separate vignette on the input file
/// format. QTL mapping data consists of a set of tables of data: marker
/// genotypes, phenotypes, marker maps, etc. In the new format, these different
/// tables are in separate comma-delimited (CSV) files. In each file, the first
/// column is a set of IDs for the rows, and the first row is a set of IDs for
/// the columns. For example, the phenotype data file will have individual IDs
/// in the first column and phenotype names in the first row.
///
/// For BIMBAM:
/// https://www.haplotype.org/download/bimbam-manual.pdf
pub mod bimbam;
pub mod rqtl2;
pub mod util;
