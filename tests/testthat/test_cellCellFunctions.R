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

# ranks
rks <- cellCellRanks(sce)
expect_true(length(rks$selected) == 3)
expect_true(is.vector(rks$mode1))
expect_true(is.vector(rks$mode2))
expect_true(is.vector(rks$mode3))

# Decomposition
cellCellDecomp(sce, ranks=rks$selected)
expect_equivalent(names(metadata(sce)),
    c("lrbase", "color", "label", "algorithm", "sctensor"))
expect_equivalent(metadata(sce)$algorithm, "ntd")

# # Report
# tmp <- tempdir()
# # cl <- snow::makeCluster(2, "SOCK")
# print(
# 	system.time(
# 		cellCellReport(sce, reducedDimNames="TSNE", out.dir=tmp,
# 	    title="Cell-cell interaction within Germline, Male, GSE86146",
# 	    author="Koki Tsuyuzaki", thr=10, top=10, html.open=TRUE)
# 	)
# )
# expect_true(file.exists(paste0(tmp, "/index.html")))
# # stopCluster(cl)
# # gc();gc()
