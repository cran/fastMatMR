% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fastMatMR-package.R
\name{write_fmm}
\alias{write_fmm}
\title{Convert Various Numeric Types to Matrix Market Format}
\usage{
write_fmm(input, filename = "out.mtx")
}
\arguments{
\item{input}{A numeric object to be converted. This can be a numeric vector,
a matrix, or a sparse matrix.}

\item{filename}{The name of the output file
where the Matrix Market formatted data will be saved.  It is recommended to
use a filename ending with ".mtx" for clarity.}
}
\value{
A boolean indicating success or failure. Writes a MTX file to disk.
}
\description{
This function takes different types of numeric inputs—vectors, matrices, and
sparse matrices— and converts them into Matrix Market files. The output file
is written to disk.
}
\examples{
vec <- c(1, 2, 3)
mat <- matrix(c(1, 2, 3, 4), nrow = 2)
sparse_mat_diag <- Matrix::Matrix(c(1, 0, 0, 2), nrow = 2, sparse = TRUE)
## Diagonal ^-
sparse_mat <- Matrix::Matrix(c(1, 1, 0, 2), nrow = 2, sparse = TRUE)
## And not diagonal -^
write_fmm(vec, tempfile(fileext = ".mtx"))
write_fmm(mat, tempfile(fileext = ".mtx"))
write_fmm(sparse_mat_diag, tempfile(fileext = ".mtx"))
write_fmm(sparse_mat, tempfile(fileext = ".mtx"))

}
