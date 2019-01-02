context("cellCellFunctions")

library(SingleCellExperiment)
library(LRBase.Hsa.eg.db)

data(GermMale)
data(labelGermMale)
data(tsneGermMale)

# SingleCellExperiment-class
sce <- SingleCellExperiment(assays = list(counts = GermMale))
reducedDims(sce) <- SimpleList(TSNE=tsneGermMale$Y)

expect_true(is.null(metadata(sce)$lrbase))
expect_true(is.null(metadata(sce)$color))
expect_true(is.null(metadata(sce)$label))

# Setting
cellCellSetting(sce, LRBase.Hsa.eg.db, labelGermMale, names(labelGermMale))

expect_false(is.null(metadata(sce)$lrbase))
expect_false(is.null(metadata(sce)$color))
expect_false(is.null(metadata(sce)$label))
