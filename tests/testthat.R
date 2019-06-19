library("testthat")
library("scTensor")

options(testthat.use_colours = FALSE)

test_file("testthat/test_GermMale.R")
test_file("testthat/test_tsneGermMale.R")
test_file("testthat/test_labelGermMale.R")
test_file("testthat/test_cellCellFunctions.R")
test_file("testthat/test_CCSParamsFunctions.R")
test_file("testthat/test_convertToNCBIGeneID.R")
