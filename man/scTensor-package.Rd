\name{scTensor-package}
\alias{scTensor-package}
\alias{scTensor}
\docType{package}
\title{
\packageTitle{scTensor}
}
\description{
\packageDescription{scTensor}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{scTensor}
\packageIndices{scTensor}
}
\author{
\packageAuthor{scTensor}

Maintainer: \packageMaintainer{scTensor}
}

\keyword{ package }
\seealso{\code{\link{GermMale}},\code{\link{labelGermMale}},
\code{\link{tsneGermMale}},\code{\link{cellCellSetting}},
\code{\link{cellCellDecomp}},\code{\link{cellCellReport}}}
\examples{
if(interactive()){
    # Package Loading
    library(SingleCellExperiment)
    library(LRBase.Hsa.eg.db)
    # Data Loading
    data(GermMale)
    data(labelGermMale)
    data(tsneGermMale)
    # SingleCellExperiment-class
    sce <- SingleCellExperiment(assays = list(counts = GermMale))
    reducedDims(sce) <- SimpleList(TSNE=tsneGermMale$Y)
    # Setting
    cellCellSetting(sce, LRBase.Hsa.eg.db, labelGermMale, names(labelGermMale))
    # Decomposition
    cellCellDecomp(sce, ranks=c(4,4,5))
    # Report
    tmp <- tempdir()
    cellCellReport(sce, reducedDimNames="TSNE", out.dir=tmp,
        html.open = TRUE,
        title="Cell-cell interaction within Germline, Male, GSE86146",
        author="Koki Tsuyuzaki", thr=5)
}else{
    ls("package:scTensor")
}
}