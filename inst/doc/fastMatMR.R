## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("ropensci/fastMatMR")

## -----------------------------------------------------------------------------
library(fastMatMR)

## -----------------------------------------------------------------------------
vec <- c(1, 2, 3)
temp_file_vec <- tempfile(fileext = ".mtx")
write_fmm(vec, temp_file_vec)

## -----------------------------------------------------------------------------
mat <- matrix(c(1, 2, 3, 4), nrow = 2)
temp_file_mat <- tempfile(fileext = ".mtx")
write_fmm(mat, temp_file_mat)

## -----------------------------------------------------------------------------
sp_mat <- Matrix::sparseMatrix(i = c(1, 3), j = c(2, 4), x = 7:8)
temp_file_sp_mat <- tempfile(fileext = ".mtx")
write_fmm(sp_mat, temp_file_sp_mat)

## -----------------------------------------------------------------------------
vec <- c(1, 2, 3.32, 225.61)
temp_file_vec_r <- tempfile(fileext = ".mtx")
vec_to_fmm(vec, temp_file_vec_r)
fmm_to_vec(temp_file_vec_r)

## ----eval=FALSE---------------------------------------------------------------
#  spmat <- Matrix::Matrix(c(1, 0, 3, NA), nrow = 2, sparse = TRUE)
#  temp_file_sp_na <- tempfile(fileext = ".mtx")
#  Matrix::writeMM(spmat, temp_file_sp_na)
#  Matrix::readMM(temp_file_sp_na)
#  ## NULL
#  ## 2 x 2 sparse Matrix of class "dgTMatrix"
#  ## [1,] 1  3e+00
#  ## [2,] . 1e+308

