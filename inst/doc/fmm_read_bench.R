## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "",
  fig.retina = 2,
  fig.width = 800 / 72,
  fig.height = 400 / 72,
  dev.args = list(type = "cairo-png")
)


## -----------------------------------------------------------------------------
library(Matrix)
library(fastMatMR)
library(microbenchmark)
library(ggplot2)


## -----------------------------------------------------------------------------
# Function to create a sparse matrix of given size
create_sparse_matrix <- function(n, sparsity = 0.7) {
      mat <- matrix(0, nrow = n, ncol = n)
      for (i in 1:n) {
        for (j in 1:n) {
          if (runif(1) > sparsity) {
            mat[i, j] <- rnorm(1)
          }
        }
      }
      return(Matrix(mat, sparse = TRUE))
    }
 
    # Define a range of matrix sizes
    sizes <- c(10, 100, 500, 1000, 2000, 3000)
 
    # Prepare data frame to store results
    results_fixed_sparsity <- data.frame()
 
    # Benchmarking
    for (n in sizes) {
      message("Benchmarking for matrix size: ", n, "x", n)
 
      # Generate a sparse matrix of size n x n
      testmat <- create_sparse_matrix(n)
      write_fmm(testmat, "sparse.mtx")
 
      # Run the benchmarks, we coerce to a sparse matrix for readMM for fairness
      bm <- microbenchmark(
        Matrix_readMM = as(readMM("sparse.mtx"), "CsparseMatrix"),
        fastMatMR_read_fmm = fmm_to_sparse_Matrix("sparse.mtx"),
        times = 10
      )
 
      bm$size <- n
      results_fixed_sparsity <- rbind(results_fixed_sparsity, bm)
    }


## ----fixed-sparse-read--------------------------------------------------------
# Plotting
suppressWarnings(print(
      ggplot(results_fixed_sparsity, aes(x = size, y = time, color = expr)) +
        geom_point() +
        geom_smooth(method = "loess") +
        ggtitle("Benchmarking reads with fixed sparsity for 70% sparsity") +
        xlab("Matrix Size") +
        ylab("Time (ns)")
))


## -----------------------------------------------------------------------------
# Sparsity levels to test
sparsity_levels <- seq(0.45, 0.95, by = 0.1)

# Prepare data frame to store results
results_varying_sparsity <- data.frame()

# Benchmarking
for (sparsity in sparsity_levels) {
  message("Benchmarking for sparsity level: ", sparsity)

  # Generate a sparse matrix of size 2000 x 2000 with varying sparsity
  testmat <- create_sparse_matrix(2000, sparsity)
  write_fmm(testmat, "sparse.mtx")

  # Run the benchmarks
  bm <- microbenchmark(
    Matrix_readMM = as(readMM("sparse.mtx"), "CsparseMatrix"),
    fastMatMR_read_fmm = fmm_to_sparse_Matrix("sparse.mtx"),
    times = 10
  )

  bm$sparsity <- sparsity
  results_varying_sparsity <- rbind(results_varying_sparsity, bm)
}


## ----varying-sparse-read------------------------------------------------------
ggplot(results_varying_sparsity, aes(x = sparsity, y = time, color = expr)) +
  geom_point() +
  geom_smooth(method = "loess") +
  scale_x_log10() +
  scale_y_log10() +
  ggtitle("Benchmarking reads with varying sparsity for 2000 entries") +
  xlab("Sparsity Level (log10)") +
  ylab("Time (ns, log10)")


## -----------------------------------------------------------------------------
sessionInfo()

