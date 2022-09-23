[![DOI](https://zenodo.org/badge/135140761.svg)](https://zenodo.org/badge/latestdoi/135140761)
[![](https://img.shields.io/badge/release%20version-2.6.0-green.svg)](https://www.bioconductor.org/packages/release/bioc/html/scTensor.html)

# scTensor
 R package for detection of cell-cell interaction using Non-negative Tensor Decomposition

Installation of Dependent Packages
======
```r
# CRAN
install.packages(c("RSQLite", "igraph", "plotly", "nnTensor",
    "rTensor", "abind", "plotrix", "heatmaply", "tagcloud",
    "rmarkdown", "knitr", "outliers", "crayon", "checkmate",
    "testthat", "Seurat", "BiocManager"),
    repos="http://cran.r-project.org")

# Bioconductor
library("BiocManager")
BiocManager::install(c("S4Vectors", "reactome.db", "AnnotationDbi",
    "SummarizedExperiment", "SingleCellExperiment", "BiocStyle",
    "biomaRt", "MeSHDbi", "Category", "meshr", "GOstats", "ReactomePA",
    "DOSE", "LRBase.Hsa.eg.db", "MeSH.Hsa.eg.db", "LRBase.Mmu.eg.db",
    "MeSH.Mmu.eg.db", "LRBaseDbi", "Homo.sapiens"),
    suppressUpdates=TRUE)
```

Installation
======
```r
git clone https://github.com/rikenbit/scTensor/
R CMD INSTALL scTensor
```
or type the code below in the R console window
```r
install.packages("devtools", repos="http://cran.r-project.org")
library(devtools)
devtools::install_github("rikenbit/scTensor")
```

## License
Copyright (c) 2018 Koki Tsuyuzaki and Laboratory for Bioinformatics Research, RIKEN Center for Biosystems Dynamics Reseach
Released under the [Artistic License 2.0](http://www.perlfoundation.org/artistic_license_2_0).

## Authors
- Koki Tsuyuzaki
- Manabu Ishii
- Itoshi Nikaido