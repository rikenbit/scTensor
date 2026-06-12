library(scTensor)

context("BuildSignedCCI")

test_that("BuildSignedCCI basic construction", {
  # 3 cell types, 2 LR pairs
  lig <- matrix(c(1, 0, 0, 0, 1, 0), nrow = 3)
  rec <- matrix(c(0, 1, 0, 0, 0, 1), nrow = 3)
  signs <- c(1, -1)

  result <- BuildSignedCCI(lig, rec, signs)

  expect_identical(dim(result$A_pos), c(3L, 3L))
  expect_identical(dim(result$A_neg), c(3L, 3L))
  expect_identical(dim(result$X_pos), c(3L, 3L, 2L))
  expect_identical(dim(result$X_neg), c(3L, 3L, 2L))

  # LR pair 1 (+): cell 1 ligand, cell 2 receptor -> A_pos[1,2] = 1
  expect_identical(result$A_pos[1, 2], 1)
  # LR pair 2 (-): cell 2 ligand, cell 3 receptor -> A_neg[2,3] = 1
  expect_identical(result$A_neg[2, 3], 1)
})

test_that("BuildSignedCCI edge_table is correct", {
  lig <- matrix(c(1, 0, 0, 1), nrow = 2)
  rec <- matrix(c(0, 1, 1, 0), nrow = 2)
  signs <- c(1, -1)

  result <- BuildSignedCCI(lig, rec, signs)

  expect_true(is.data.frame(result$edge_table))
  expect_true(all(c("sender", "receiver", "lr_pair", "score", "sign") %in%
    names(result$edge_table)))
  expect_gt(nrow(result$edge_table), 0)
})

test_that("BuildSignedCCI all positive signs yields zero A_neg", {
  lig <- matrix(runif(12), nrow = 3)
  rec <- matrix(runif(12), nrow = 3)
  signs <- rep(1, 4)

  result <- BuildSignedCCI(lig, rec, signs)

  expect_identical(sum(result$A_neg), 0)
  expect_gt(sum(result$A_pos), 0)
})

test_that("BuildSignedCCI validates inputs", {
  expect_error(BuildSignedCCI("not a matrix", matrix(1), 1))
  expect_error(BuildSignedCCI(matrix(1, 2, 3), matrix(1, 3, 2), c(1, 1, 1)))
  expect_error(BuildSignedCCI(matrix(1, 2, 2), matrix(1, 2, 2), c(1, 0)))
})
