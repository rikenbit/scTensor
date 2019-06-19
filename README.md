# scTensor
 R package for detection of cell-cell interaction using Non-negative Tensor Decomposition


Installation of Dependent Packages
======
```r
# CRAN
install.packages("RSQLite", repos="http://cran.r-project.org")
install.packages("igraph", repos="http://cran.r-project.org")
install.packages("plotly", repos="http://cran.r-project.org")
install.packages("nnTensor", repos="http://cran.r-project.org")
install.packages("rTensor", repos="http://cran.r-project.org")
install.packages("abind", repos="http://cran.r-project.org")
install.packages("plotrix", repos="http://cran.r-project.org")
install.packages("heatmaply", repos="http://cran.r-project.org")
install.packages("tagcloud", repos="http://cran.r-project.org")
install.packages("rmarkdown", repos="http://cran.r-project.org")
install.packages("knitr", repos="http://cran.r-project.org")
install.packages("outliers", repos="http://cran.r-project.org")
install.packages("crayon", repos="http://cran.r-project.org")
install.packages("checkmate", repos="http://cran.r-project.org")
install.packages("testthat", repos="http://cran.r-project.org")
install.packages("Seurat", repos="http://cran.r-project.org")
install.packages("BiocManager", repos="http://cran.r-project.org")

# Bioconductor
BiocManager::install("S4Vectors", suppressUpdates=TRUE)
BiocManager::install("reactome.db", suppressUpdates=TRUE)
BiocManager::install("AnnotationDbi", suppressUpdates=TRUE)
BiocManager::install("SummarizedExperiment", suppressUpdates=TRUE)
BiocManager::install("SingleCellExperiment", suppressUpdates=TRUE)
BiocManager::install("BiocStyle", suppressUpdates=TRUE)
BiocManager::install("biomaRt", suppressUpdates=TRUE)
BiocManager::install("MeSHDbi", suppressUpdates=TRUE)
BiocManager::install("Category", suppressUpdates=TRUE)
BiocManager::install("meshr", suppressUpdates=TRUE)
BiocManager::install("GOstats", suppressUpdates=TRUE)
BiocManager::install("ReactomePA", suppressUpdates=TRUE)
BiocManager::install("DOSE", suppressUpdates=TRUE)
BiocManager::install("LRBase.Hsa.eg.db", suppressUpdates=TRUE)
BiocManager::install("MeSH.Hsa.eg.db", suppressUpdates=TRUE)
BiocManager::install("LRBase.Mmu.eg.db", suppressUpdates=TRUE)
BiocManager::install("MeSH.Mmu.eg.db", suppressUpdates=TRUE)
BiocManager::install("LRBaseDbi", suppressUpdates=TRUE)
BiocManager::install("Homo.sapiens", suppressUpdates=TRUE)
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
