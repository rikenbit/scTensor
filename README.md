# scTensor
 R package for detection of cell-cell interaction using Non-negative Tensor Decomposition


Installation Dependent Packages
======
~~~~
install.packages("igraph", repos="http://cran.r-project.org")
install.packages("plotly", repos="http://cran.r-project.org")
install.packages("nnTensor", repos="http://cran.r-project.org")
install.packages("rTensor", repos="http://cran.r-project.org")
install.packages("abind", repos="http://cran.r-project.org")
install.packages("plotrix", repos="http://cran.r-project.org")
install.packages("heatmaply", repos="http://cran.r-project.org")
install.packages("tagcloud", repos="http://cran.r-project.org")
install.packages("RColorBrewer", repos="http://cran.r-project.org")
install.packages("rmarkdown", repos="http://cran.r-project.org")

source("http://bioconductor.org/biocLite.R")
biocLite("BiocStyle", suppressUpdates=TRUE)
biocLite("biomaRt", suppressUpdates=TRUE)
biocLite("reactome.db", suppressUpdates=TRUE)
biocLite("MeSHDbi", suppressUpdates=TRUE)
biocLite("MeSH.Hsa.eg.db", suppressUpdates=TRUE)
biocLite("MeSH.Mmu.eg.db", suppressUpdates=TRUE)
biocLite("MeSH.Ath.eg.db", suppressUpdates=TRUE)
biocLite("MeSH.Rno.eg.db", suppressUpdates=TRUE)
biocLite("MeSH.Bta.eg.db", suppressUpdates=TRUE)
biocLite("MeSH.Cel.eg.db", suppressUpdates=TRUE)
biocLite("MeSH.Dme.eg.db", suppressUpdates=TRUE)
biocLite("MeSH.Dre.eg.db", suppressUpdates=TRUE)
biocLite("MeSH.Gga.eg.db", suppressUpdates=TRUE)
biocLite("MeSH.Pab.eg.db", suppressUpdates=TRUE)
biocLite("MeSH.Xtr.eg.db", suppressUpdates=TRUE)
biocLite("MeSH.Ssc.eg.db", suppressUpdates=TRUE)
~~~~

Installation
======
~~~~
git clone https://github.com/rikenbit/LRBaseDbi/
git clone https://github.com/rikenbit/scTensor/
R CMD INSTALL LRBaseDbi
R CMD INSTALL scTensor
~~~~
or type the code below in the R console window
~~~~
library(devtools)
devtools::install_github("rikenbit/LRBaseDbi")
devtools::install_github("rikenbit/scTensor")
~~~~

## License
Copyright (c) 2018 Koki Tsuyuzaki and Laboratory for Bioinformatics Research, RIKEN Center for Biosystems Dynamics Reseach
Released under the [Artistic License 2.0](http://www.perlfoundation.org/artistic_license_2_0).

## Authors
- Koki Tsuyuzaki
- Manabu Ishii
- Itoshi Nikaido
