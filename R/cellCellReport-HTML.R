# Header Row
.XYZ_HEADER1 <- function(index, i){
    TITLE <- paste0("Details of (",
    paste(index[i, seq_len(3)], collapse=","),
    ") Pattern")
    HEADER <- paste0("---\n",
    "title: <font color='#87b13f'>XXXXX</font>\n",
    "output:\n",
    "    html_document:\n",
    "        toc: true\n",
    "        toc_float: true\n",
    "        toc_depth: 2\n",
    "---\n",
    "# <font color='#1881c2'>(", paste(index[i, seq_len(2)], collapse=","),
    ",:) Pattern</font>\n\n",
    "![](figures/CCIHypergraph_",
    paste(index[i, seq_len(2)], collapse="_"), ".png)\n\n")
    sub("XXXXX", TITLE, HEADER)
}

.XYZ_HEADER1_2 <- function(index, i){
    TITLE <- paste0("Details of (",
    paste(c(index[i, seq_len(2)], ""), collapse=","),
    ") -related L-R pairs")
    HEADER <- paste0("---\n",
    "title: <font color='#87b13f'>XXXXX</font>\n",
    "output:\n",
    "    html_document:\n",
    "        toc: true\n",
    "        toc_float: true\n",
    "        toc_depth: 2\n",
    "---\n",
    "# <font color='#1881c2'>(", paste(index[i, seq_len(2)], collapse=","),
    ",:) Pattern</font>\n\n",
    "![](figures/CCIHypergraph_",
    paste(index[i, seq_len(2)], collapse="_"), ".png)\n\n")
    sub("XXXXX", TITLE, HEADER)
}

.XYZ_HEADER2 <- function(i, x, top){
    paste0("# <font color='#1881c2'>(\\*,\\*,", i, ") Pattern",
    " (Top", top, " LR-pairs)</font>\n\n",
    "```{r}\n", # Top
    "library(\"scTensor\")\n",
    "load(\"reanalysis.RData\")\n",
    "col <- .setColor(\"many\")[", x, "]\n",
    "scTensor:::.myvisNetwork(g, col)\n",
    "```\n\n", # Bottom
    "|Rank|Ligand Gene|Receptor Gene|",
    "Ligand Expression (Log10(exp + 1))|",
    "Receptor Expression (Log10(exp + 1))|",
    "LR-pair factor value (Percentage)|",
    "P-value (Grubbs test)|",
    "Q-value (BH method)|",
    "L-R Evidence cf. [rikenbit/lrbase-workflow](https://github.com/rikenbit/lrbase-workflow)|",
    "PubMed|\n",
    "|----|----|----|----|----|----|----|----|----|----|\n")
}

.XYZ_HEADER2_2 <- function(index, i, top){
    paste0("# <font color='#1881c2'>(",
    paste(c(index[i, seq_len(2)], ""), collapse=","),
    ") -related L-R pairs (Top", top, " pairs)</font>\n\n",
    "```{r}\n", # Top
    "library(\"scTensor\")\n",
    "load(\"reanalysis.RData\")\n",
    "col <- .setColor(\"many\")[", i, "]\n",
    "scTensor:::.myvisNetwork(g, col)\n",
    "```\n\n", # Bottom
    "|Rank|Ligand Gene|Receptor Gene|",
    "Ligand Expression (Log10(exp + 1))|",
    "Receptor Expression (Log10(exp + 1))|",
    "LR-pair factor value (Percentage)|",
    "P-value (Grubbs test)|",
    "Q-value (BH method)|",
    "L-R Evidence cf. [rikenbit/lrbase-workflow](https://github.com/rikenbit/lrbase-workflow)|",
    "PubMed|\n",
    "|----|----|----|----|----|----|----|----|----|----|\n")
}

.MAINHEADER <- function(author, title){
    HEADER <- paste0("---\ntitle: XXXXX\n",
        "author: YYYYY\ndate:",
        " \"`r Sys.time()`\"\n",
        "output:\n",
        " BiocStyle::html_document:\n",
        "  toc_float: true\n",
        "vignette: >\n ",
        "%\\VignetteIndexEntry{Vignette Title}\n ",
        "%\\VignetteEngine{knitr::rmarkdown}\n ",
        "%\\VignetteEncoding{UTF-8}\n---\n")
    sub("YYYYY", author, sub("XXXXX", title, HEADER))
}

# 1. About Dataset and scTensor Algorithm
.BODY1 <- paste0("\n\n# About scTensor Algorithm\n\n",
    "![](Workflow.png)\n",
    "[scTensor](https://bioconductor.org/packages/release/",
    "bioc/html/scTensor.html) is the R/Bioconductor package",
    " for visualization of cell-cell interaction within ",
    "single-cell RNA-Seq data. ",
    "The calculation consists of several steps.\n\n",
    "Firstly, [LRBase.XXX.eg.db](https://bioconductor.org/",
    "packages/release/bioc/html/LRBase.Hsa.eg.db.html)-type",
    " package is loaded for retriving ligand-receptor gene ",
    "relationship (XXX is the abbreviation of some organisms ",
    "like \"Hsa\" as Homo sapiens). ",
    "scTensor searches the corresponding pair of genes in the ",
    "rownames of input data matrix and extracted as vector. ",
    "In this step, the gene identifier is limited as [NCBI ",
    "Gene ID](https://www.ncbi.nlm.nih.gov/gene) for now.\n\n",
    "Next, the outer product of two vectors are calculated ",
    "and summarized as a matrix. ",
    "Here, the multiple matrices can be represented as a three-order",
    " \"tensor\" (Ligand-Cell * Receptor-Cell * LR-Pair). ",
    "scTensor decomposes the tensor into a small tensor (core tensor)",
    " and three factor matrices. ",
    "Tensor decomposition is very similar to the matrix decomposition",
    " like PCA (principal component analysis). ",
    "The core tensor is similar to eigenvalue of PCA and means how",
    " much the pattern is outstanding. ",
    "Likewise, three matrices is similar to the PC scores and ",
    "loadings of PCA and represents which ligand-cell/receptor-cell",
    "/LR-pair are informative. ",
    "When the matrices have negative values, distinguishing that ",
    "which direction (+/-) is important and which is not, is ",
    "difficult to interpret and laboring task. ",
    "That's why, scTensor performs non-negative Tucker ",
    "decomposition (NTD), which is non-negative version ",
    "of Tucker decomposition (c.f. [nnTensor]",
    "(https://cran.r-project.org/package=nnTensor)). \n\n",
    "Finaly, the result of NTD is summarized as this report. ",
    "Most of plots belows are visualized by ",
    "[plotly](https://plot.ly/r/) package and interactively ",
    "search the presise information of the plot. ",
    "The three factor matrices can be interactively viewed ",
    "and which celltypes are responds to the cell-cell interaction.",
    "The mode-3 (LR-pair) direction of sum of the core tensor ",
    "is calculated and visualized as Ligand-Receptor Patterns. ",
    "Detail of (Ligand-Cell, Receptor-Cell, LR-pair) Patterns ",
    "are also visualized as below.\n\n",
    "For more detail, visit the [vignette](https://bioconductor.org",
    "/packages/devel/bioc/vignettes/scTensor/inst/doc/scTensor.html)",
    " of [scTensor](https://bioconductor.org/",
    "packages/release/bioc/html/scTensor.html)"
    )

.BODY1_2 <- paste0("\n\n# About scTensor Algorithm\n\n",
    "![](Workflow_2.png)\n",
    "[scTensor](https://bioconductor.org/packages/release/",
    "bioc/html/scTensor.html) is the R/Bioconductor package",
    " for visualization of cell-cell interaction within ",
    "single-cell RNA-Seq data. ",
    "The calculation consists of several steps.\n\n",
    "Firstly, [LRBase.XXX.eg.db](https://bioconductor.org/",
    "packages/release/bioc/html/LRBase.Hsa.eg.db.html)-type",
    " package is loaded for retriving ligand-receptor gene ",
    "relationship (XXX is the abbreviation of some organisms ",
    "like \"Hsa\" as Homo sapiens). ",
    "scTensor searches the corresponding pair of genes in the ",
    "rownames of input data matrix and extracted as vector. ",
    "In this step, the gene identifier is limited as [NCBI ",
    "Gene ID](https://www.ncbi.nlm.nih.gov/gene) for now.\n\n",
    "Next, the all elements of extracted two vectors are summed ",
    "up with all possible combination (Kronecker sum) and summarized",
    " as a matrix. ",
    "Here, the multiple matrices can be represented as a three-order",
    " \"tensor\" (Ligand-Cell * Receptor-Cell * LR-Pair). ",
    "scTensor decomposes the tensor into a small tensor (core tensor)",
    " and three factor matrices. ",
    "Tensor decomposition is very similar to the matrix decomposition",
    " like PCA (principal component analysis). ",
    "The core tensor is similar to eigenvalue of PCA and means how",
    " much the pattern is outstanding. ",
    "Likewise, three matrices is similar to the PC scores and ",
    "loadings of PCA and represents which ligand-cell/receptor-cell",
    "/LR-pair are informative. ",
    "When the matrices have negative values, distinguishing that ",
    "which direction (+/-) is important and which is not, is ",
    "difficult to interpret and laboring task. ",
    "That's why, scTensor performs non-negative Tucker ",
    "decomposition (NTD), which is non-negative version ",
    "of Tucker decomposition (c.f. [nnTensor]",
    "(https://cran.r-project.org/package=nnTensor)). \n\n",
    "Finaly, the result of NTD is summarized as this report. ",
    "Most of plots belows are visualized by ",
    "[plotly](https://plot.ly/r/) package and interactively ",
    "search the presise information of the plot. ",
    "The three factor matrices can be interactively viewed ",
    "and which celltypes are responds to the cell-cell interaction.",
    "The mode-3 (LR-pair) direction of sum of the core tensor ",
    "is calculated and visualized as Ligand-Receptor Patterns. ",
    "Detail of (Ligand-Cell, Receptor-Cell, LR-pair) Patterns ",
    "are also visualized as below.\n\n",
    "For more detail, visit the [vignette](https://bioconductor.org",
    "/packages/devel/bioc/vignettes/scTensor/inst/doc/scTensor.html)",
    " of [scTensor](https://bioconductor.org/",
    "packages/release/bioc/html/scTensor.html)"
    )

# 2. Global statistics and plots
.BODY2 <- paste0("\n\n# Global statistics and plots\n\n",
    "The result of scTensor is saved as a R binary file",
    " (reanalysis.RData).\n",
    "```{r}\n", # Top
    "load(\"reanalysis.RData\")\n\n",
    "# SingleCellExperiment object\n",
    "sce\n",
    "# Data size\n",
    "metadata(sce)$datasize\n",
    "# Reduced data size\n",
    "metadata(sce)$ranks\n",
    "# Reconstruction Error of NTD\n",
    "head(metadata(sce)$recerror)\n",
    "tail(metadata(sce)$recerror)\n",
    "# Relative Change of NTD\n",
    "head(metadata(sce)$relchange)\n",
    "tail(metadata(sce)$relchange)\n",
    "# Gene expression matrix\n",
    "is(input)\n",
    "dim(input)\n",
    "input[seq_len(2), seq_len(2)]\n",
    "# The result of 2D dimensional reduction (e.g. t-SNE)\n",
    "is(twoD)\n",
    "dim(twoD)\n",
    "head(twoD)\n",
    "# Ligand-Receptor corresponding table",
    " extracted from LRBase.XXX.eg.db\n",
    "is(LR)\n",
    "dim(LR)\n",
    "head(LR, 2)\n",
    "# Celltype label and color scheme\n",
    "is(celltypes)\n",
    "length(celltypes)\n",
    "head(celltypes)\n",
    "# Core tensor values\n",
    "is(index)\n",
    "dim(index)\n",
    "head(index)\n",
    "is(corevalue)\n",
    "length(corevalue)\n",
    "head(corevalue)\n",
    "# Selected corevalue position with thr threshold \"thr\"\n",
    "is(selected)\n",
    "length(selected)\n",
    "head(selected)\n",
    "# The result of 2-class clustering\n",
    "# The result of analysis in each L vector\n",
    "is(ClusterL)\n",
    "dim(ClusterL)\n",
    "# The result of analysis in each R vector\n",
    "is(ClusterR)\n",
    "dim(ClusterR)\n",
    "# The result of analysis in each LR vector\n",
    "is(out.vecLR)\n",
    "dim(out.vecLR)\n",
    "head(out.vecLR)\n",
    "head(out.vecLR[1,1][[1]])\n",
    "```\n\n", # Bottom
    "\n\n## Number of cells in each celltype\n\n",
    "```{r}\n", # Top
    "color <- names(celltypes)\n",
    "colors <- celltypes[vapply(unique(color), ",
    "function(x){which(color == x)[1]}, 0L)]\n",
    "numcells <- vapply(names(colors), function(x){",
    "length(which(color == x))",
    "}, 0L)\n",
    "numcells <- data.frame(numcells=numcells, ",
    "celltypes=names(numcells))\n\n",
    "suppressPackageStartupMessages(library(plotly))\n",
    "plot_ly(numcells, x=~celltypes, y=~numcells, type=\"bar\", ",
    "marker = list(color = colors))\n",
    "```\n\n", # Bottom
    "\n\n## Number of expressed genes in ",
    "each celltype (Non-zero genes)\n\n",
    "```{r}\n", # Top
    "expgenes <- apply(input, 2, ",
    "function(x){length(which(x != 0))})\n",
    "expgenes <- data.frame(expgenes=expgenes, ",
    "celltypes=names(celltypes))\n",
    "plot_ly(expgenes, y=~expgenes, color=color, ",
    "colors=colors, type=\"box\")\n",
    "```\n\n", # Bottom
    "\n\n## Two dimensional plot of all cells\n\n",
    "```{r}\n", # Top
    "plot_ly(x=twoD[,1], y=twoD[,2], ",
    "color = color, ",
    "colors = colors, ",
    "type = \"scatter\", ",
    "text = rownames(twoD), ",
    "mode = \"markers\")\n",
    "```\n\n", # Bottom
    "\n\n## Distribution of core tensor values\n\n",
    "```{r}\n", # Top
    "corenames <- vapply(seq_len(nrow(index)), ",
    "function(x){paste(index[x,seq_len(3)], collapse=\",\")}, \"\")\n",
    "plot_ly(x=seq_along(corevalue), y=corevalue, ",
    "type=\"bar\", color=names(corevalue), text=corenames, ",
    "colors = c(\"#999999\", \"#E41A1C\"))\n",
    "```\n" # Bottom
    )

.BODY2_2 <- paste0("\n\n# Global statistics and plots\n\n",
    "The result of scTensor is saved as a R binary file",
    " (reanalysis.RData).\n",
    "```{r}\n", # Top
    "load(\"reanalysis.RData\")\n\n",
    "# SingleCellExperiment object\n",
    "sce\n",
    "# Data size\n",
    "metadata(sce)$datasize\n",
    "# Reduced data size\n",
    "metadata(sce)$ranks\n",
    "# Reconstruction Error of NTD\n",
    "head(metadata(sce)$recerror)\n",
    "tail(metadata(sce)$recerror)\n",
    "# Relative Change of NTD\n",
    "head(metadata(sce)$relchange)\n",
    "tail(metadata(sce)$relchange)\n",
    "# Gene expression matrix\n",
    "is(input)\n",
    "dim(input)\n",
    "input[seq_len(2), seq_len(2)]\n",
    "# The result of 2D dimensional reduction (e.g. t-SNE)\n",
    "is(twoD)\n",
    "dim(twoD)\n",
    "head(twoD)\n",
    "# Ligand-Receptor corresponding table",
    " extracted from LRBase.XXX.eg.db\n",
    "is(LR)\n",
    "dim(LR)\n",
    "head(LR, 2)\n",
    "# Celltype label and color scheme\n",
    "is(celltypes)\n",
    "length(celltypes)\n",
    "head(celltypes)\n",
    "# Core tensor values\n",
    "is(index)\n",
    "dim(index)\n",
    "head(index)\n",
    "is(corevalue)\n",
    "length(corevalue)\n",
    "head(corevalue)\n",
    "# Selected corevalue position with thr threshold \"thr\"\n",
    "is(selected)\n",
    "length(selected)\n",
    "head(selected)\n",
    "# The result of 2-class clustering\n",
    "# The result of analysis in each L vector\n",
    "is(ClusterL)\n",
    "dim(ClusterL)\n",
    "# The result of analysis in each R vector\n",
    "is(ClusterR)\n",
    "dim(ClusterR)\n",
    "# The result of analysis in each LR vector\n",
    "is(out.vecLR)\n",
    "```\n\n", # Bottom
    "\n\n## Number of cells in each celltype\n\n",
    "```{r}\n", # Top
    "color <- names(celltypes)\n",
    "colors <- celltypes[vapply(unique(color), ",
    "function(x){which(color == x)[1]}, 0L)]\n",
    "numcells <- vapply(names(colors), function(x){",
    "length(which(color == x))",
    "}, 0L)\n",
    "numcells <- data.frame(numcells=numcells, ",
    "celltypes=names(numcells))\n\n",
    "suppressPackageStartupMessages(library(plotly))\n",
    "plot_ly(numcells, x=~celltypes, y=~numcells, type=\"bar\", ",
    "marker = list(color = colors))\n",
    "```\n\n", # Bottom
    "\n\n## Number of expressed genes in ",
    "each celltype (Non-zero genes)\n\n",
    "```{r}\n", # Top
    "expgenes <- apply(input, 2, ",
    "function(x){length(which(x != 0))})\n",
    "expgenes <- data.frame(expgenes=expgenes, ",
    "celltypes=names(celltypes))\n",
    "plot_ly(expgenes, y=~expgenes, color=color, ",
    "colors=colors, type=\"box\")\n",
    "```\n\n", # Bottom
    "\n\n## Two dimensional plot of all cells\n\n",
    "```{r}\n", # Top
    "plot_ly(x=twoD[,1], y=twoD[,2], ",
    "color = color, ",
    "colors = colors, ",
    "type = \"scatter\", ",
    "text = rownames(twoD), ",
    "mode = \"markers\")\n",
    "```\n\n", # Bottom
    "\n\n## Distribution of core tensor values\n\n",
    "```{r}\n", # Top
    "corenames <- vapply(seq_len(nrow(index)), ",
    "function(x){paste(index[x,seq_len(3)], collapse=\",\")}, \"\")\n",
    "plot_ly(x=seq_along(corevalue), y=corevalue, ",
    "type=\"bar\", color=names(corevalue), text=corenames, ",
    "colors = c(\"#999999\", \"#E41A1C\"))\n",
    "```\n" # Bottom
    )

# 3. L-Pattern
.BODY3 <- function(numLPattern, ClusterL){
    BODY3 <- paste0("\n\n# Ligand-Cell Patterns\n\n",
    "```{r}\n", # Top
    "suppressPackageStartupMessages(library(\"heatmaply\"))\n",
    "l <- metadata(sce)$sctensor$ligand\n",
    "if(nrow(l) >= 2){\n",
    "rownames(l) <- paste0(\"(\",seq_len(nrow(l)), \",*,*)\")\n",
    "heatmaply(l,",
    "xlab=\"Celltype\",",
    "ylab=\"Pattern\",",
    "fontsize_col=20,",
    "fontsize_row=20,",
    "subplot_widths=c(0.7, 0.1),",
    "subplot_heights=c(0.2, 0.7),",
    "labRow = rownames(l),",
    "labCol = colnames(l)",
    ")\n",
    "}\n",
    "```\n" # Bottom
    )
    BODY3 <- paste0(BODY3,
        paste(vapply(seq_len(numLPattern), function(i){
            ClusterNameL <- paste(names(which(ClusterL[i,] == "selected")),
                collapse=" & ")
            titleL <- paste0("(", i, ",*,*)-Pattern", " = ", ClusterNameL)
            paste0("\n\n## ", titleL, "\n![](figures/Pattern_", i, "__.png)")
        }, ""), collapse=""))
    BODY3 <- paste0(BODY3, "\n")

}

# 4. R-Pattern
.BODY4 <- function(numRPattern, ClusterR){
    BODY4 <- paste0("\n\n# Receptor-Cell Patterns\n\n",
        "```{r}\n", # Top
        "r <- metadata(sce)$sctensor$receptor\n",
        "if(nrow(r) >= 2){\n",
        "rownames(r) <- paste0(\"(*,\",seq_len(nrow(r)), \",*)\")\n",
        "heatmaply(r,",
        "xlab=\"Celltype\",",
        "ylab=\"Pattern\",",
        "fontsize_col=20,",
        "fontsize_row=20,",
        "subplot_widths=c(0.7, 0.1),",
        "subplot_heights=c(0.2, 0.7),",
        "labRow = rownames(r),",
        "labCol = colnames(r)",
        ")\n",
        "}\n",
        "```\n" # Bottom
        )
    BODY4 <- paste0(BODY4,
        paste(vapply(seq_len(numRPattern), function(i){
            ClusterNameR <- paste(names(which(ClusterR[i,] == "selected")),
                collapse=" & ")
            titleR <- paste0("(*,", i, ",*)-Pattern", " = ", ClusterNameR)
            paste0("\n\n## ", titleR, "\n![](figures/Pattern__", i, "_.png)")
        }, ""), collapse=""))
    BODY4 <- paste0(BODY4, "\n")
}

# 5. LR-Pattern
.BODY5 <- paste0("\n\n# LR-pair Patterns\n\n",
"```{r}\n", # Top
"lr <- metadata(sce)$sctensor$lrpair\n",
"rownames(lr) <- paste0(\"(*,*,\",seq_len(nrow(lr)), \")\")\n",
"if(!is.null(GeneInfo$GeneName)){\n",
"GeneName <- GeneInfo$GeneName\n",
"Ligand <- vapply(colnames(lr), function(x){\n",
"    vapply(strsplit(x, \"_\")[[1]][1], function(xx){\n",
"    GeneName[which(GeneName[,2] == xx), 1][1]\n",
"    }, \"\")\n",
"}, \"\")\n",
"Receptor <- vapply(colnames(lr), function(x){\n",
"    vapply(strsplit(x, \"_\")[[1]][2], function(xx){\n",
"    GeneName[which(GeneName[,2] == xx), 1][1]\n",
"    }, \"\")\n",
"}, \"\")\n",
"}else{\n",
"Ligand <- vapply(colnames(lr), function(x){strsplit(x, \"_\")[[1]][1]}, \"\")\n",
"Receptor <- vapply(colnames(lr), function(x){strsplit(x, \"_\")[[1]][2]}, \"\")\n",
"}\n",
"colnames(lr) <- vapply(seq_along(Ligand), function(x){\n",
"    paste(c(Ligand[x], Receptor[x]), collapse=\" - \")\n",
"}, \"\")\n",
"target <- which(rank(colMaxs(lr)) <= 100)\n",
"heatmaply(lr[, target],",
"xlab=\"LR-Pair\",",
"ylab=\"Pattern\",",
"fontsize_col=20,",
"fontsize_row=20,",
"subplot_widths=c(0.7, 0.1),",
"subplot_heights=c(0.2, 0.7),",
"labRow = rownames(lr[, target]),",
"labCol = colnames(lr[, target])",
")\n",
"```\n" # Bottom
)

# 6. CCI-wise Hypergraph
.BODY6 <- paste0("\n\n# CCI-wise Hypergraph\n\n",
    "![](figures/CCIHypergraph.png){ width=100% }\n")

# 7. Gene-wise Hypergraph
.BODY7 <- paste0(
    "\n\n# Gene-wise Hypergraph\n\n",
    "```{r}\n", # Top
    "library(\"scTensor\")\n",
    "load(\"reanalysis.RData\")\n",
    "scTensor:::.myvisNetwork(g)\n",
    "```\n\n", # Bottom
    "[Details of Ligand Gene-centric Overview (selected)](ligand.html)\n\n",
    "[Details of Ligand Gene-centric Overview (all)](ligand_all.html)\n\n",
    "[Details of Receptor Gene-centric Overview (selected)](receptor.html)\n\n",
    "[Details of Receptor Gene-centric Overview (all)](receptor_all.html)\n\n"
    )

# 8. (Ligand, Receptor, LR-pair)-Pattern
.BODY8 <- function(selected, rmdfiles, index, corevalue){
    if(length(selected) != 0){
        htmlfiles <- gsub("Rmd", "html", rmdfiles)
        BODY8 <- vapply(selected, function(x){
            i <- selected[x]
            paste0("\n\n## (", paste(index[i, seq_len(3)],
                collapse=","),
                ") Pattern : (", round(corevalue[i], 2), "%)\n",
                "[Details of (", paste(index[i, seq_len(3)],
                collapse=","),
                ") Pattern", "](", htmlfiles[i], ")\n")
            }, "")
        BODY8 <- paste(BODY8, collapse = "\n")
        BODY8 <- paste0("# (Ligand-Cell, Receptor-Cell, LR-pair)",
            " Patterns\n\n", BODY8)
    }else{
        BODY8 <- "# (Ligand-Cell, Receptor-Cell, LR-pair) Patterns\n\n"
    }
}

.BODY8_2 <- function(selected, rmdfiles, index, corevalue){
    if(length(selected) != 0){
        htmlfiles <- gsub("Rmd", "html", rmdfiles)
        BODY8 <- vapply(selected, function(i){
            paste0("\n\n## (",
                paste(c(index[i, seq_len(2)], ""), collapse=","),
                ") -related L-R Pairs : (", round(corevalue[i], 2), "%)\n",
                "[Details of (",
                paste(c(index[i, seq_len(2)], ""), collapse=","),
                ") -related L-R Pairs", "](", htmlfiles[i], ")\n")
            }, "")
        BODY8 <- paste(BODY8, collapse = "\n")
        BODY8 <- paste0("# (Ligand-Cell, Receptor-Cell, )",
            " -related L-R Pairs\n\n", BODY8)
    }else{
        BODY8 <- "# (Ligand-Cell, Receptor-Cell, ) -related L-R Pairs\n\n"
    }
}

# 9. Session Information
.BODY9 <- paste0("\n\n# Session Information\n\n",
    "\n```{r}\n",
    "sessionInfo()\n",
    "```\n")

# 10. License
.BODY10 <- paste0("\n\n# License\n\n",
    "Copyright (c) 2018 Koki Tsuyuzaki and Laboratory for ",
    "Bioinformatics Research, RIKEN Center for Biosystems Dynamics",
    " Reseach Released under the ",
    "[Artistic License 2.0](",
    "http://www.perlfoundation.org/artistic_license_2_0)\n")

.LIGAND_HEADER <- paste0(
    "---\n",
    "title: <font color='#1881c2'>Details of Ligand Gene-centric Overview (selected)",
    "</font>\n",
    "---\n\n",
    "```{r}\n", # Top
    "library(\"scTensor\")\n",
    "load(\"reanalysis.RData\")\n",
    "scTensor:::.myvisNetwork(g)\n",
    "```\n\n", # Bottom
    "<style type='text/css'>\n",
    "table,th,td {\n",
    "width: 20%;\n",
    "border: 1px solid #f0f0f0;\n",
    "}\n",
    "</style>\n\n",
    "|Rank|No. of related pairs|Ligand Gene|Receptor Genes|",
    "Related CCIs|\n",
    "|---------------|---------------|---------------|---------------|---------------|"
)

.RECEPTOR_HEADER <- paste0(
    "---\n",
    "title: <font color='#1881c2'>Details of Receptor Gene-centric Overview (selected)",
    "</font>\n",
    "---\n\n",
    "```{r}\n", # Top
    "library(\"scTensor\")\n",
    "load(\"reanalysis.RData\")\n",
    "scTensor:::.myvisNetwork(g)\n",
    "```\n\n", # Bottom
    "<style type='text/css'>\n",
    "table,th,td {\n",
    "width: 20%;\n",
    "border: 1px solid #f0f0f0;\n",
    "}\n",
    "</style>\n\n",
    "|Rank|No. of related pairs|Receptor Gene|Ligand Genes|",
    "Related CCIs|\n",
    "|---------------|---------------|---------------|---------------|---------------|"
)

.LIGAND_BODY_2 <- function(out.vecLR, GeneInfo, index, selected){
    # Setting
    convertGeneName <- function(geneid, GeneInfo){
        if(!is.null(GeneInfo$GeneName)){
            genename <- GeneInfo$GeneName[
                which(GeneInfo$GeneName$ENTREZID == geneid), "SYMBOL"]
            if(length(genename) == 0 || genename %in% c("", NA)){
                genename = geneid
            }
            if(length(genename) != 1){
                # Black list
                genename = setdiff(genename, "cxcl11.6")[1]
            }
            genename
        }else{
            geneid
        }
    }

    # Node
    nodes <- lapply(seq_len(length(out.vecLR)), function(x){
        names(out.vecLR[[x]]$TARGET)
    })
    LnodesGeneID <- lapply(nodes, function(x){
        vapply(x, function(xx){
            strsplit(xx, "_")[[1]][1]
        }, "")
    })
    RnodesGeneID <- lapply(nodes, function(x){
        vapply(x, function(xx){
            strsplit(xx, "_")[[1]][2]
        }, "")
    })
    LnodesGeneName <- lapply(LnodesGeneID, function(x){
        vapply(x, function(xx){
            convertGeneName(xx, GeneInfo)
        }, "")
    })
    RnodesGeneName <- lapply(RnodesGeneID, function(x){
        vapply(x, function(xx){
            convertGeneName(xx, GeneInfo)
        }, "")
    })
    uniqueLnodesGeneName <- unique(unlist(LnodesGeneName))
    Freq <- vapply(uniqueLnodesGeneName, function(x){
            length(which(unlist(LnodesGeneName) == x))
        }, 0L)
    Freq <- sort(Freq, decreasing=TRUE)
    Rank <- rank(-Freq)
    Ligand <- names(Freq)
    Receptor <- lapply(Ligand, function(x){
        ReceptorGeneName <- unique(
            unlist(
                lapply(seq_len(length(LnodesGeneName)), function(xx){
                    RnodesGeneName[[xx]][which(x == LnodesGeneName[[xx]])]
                })
            )
        )
        ReceptorGeneID <- unlist(lapply(ReceptorGeneName, function(y){
            target <- which(GeneInfo$GeneName[, "SYMBOL"] == y)
            if(length(target) == 0){
                y
            }else{
                GeneInfo$GeneName[target[1], "ENTREZID"]
            }
        }))
        paste(
            paste0("[", ReceptorGeneName,
                "](https://www.ncbi.nlm.nih.gov/gene/",
            ReceptorGeneID, ") ![](figures/Receptor/",
            ReceptorGeneID, ".png)"), collapse=" ")
    })

    CCI <- vapply(Ligand, function(x){
        hit <- vapply(seq_len(length(out.vecLR)), function(xx){
            length(which(LnodesGeneName[[xx]] == x))
        }, 0L)
        vec <- seq_len(length(out.vecLR))[which(hit != 0)]
        # Out
        cciout <- vapply(vec, function(xx){
            paste0("[See the details of (",
                paste(c(index[xx, paste0("Mode", 1:2)], ""), collapse=","),
                ")-related L-R patterns](",
                "pattern_",
                paste(index[xx, paste0("Mode", 1:2)], collapse="_"),
                ".html)")
        }, "")
        paste(unique(cciout), collapse=" ")
    }, "")
    Ligand <- vapply(Ligand, function(x){
        target <- which(GeneInfo$GeneName$SYMBOL == x)
        if(length(target) == 0){
            LigandGeneID <- x
        }else{
            LigandGeneID <- GeneInfo$GeneName[target[1], "ENTREZID"]
        }
        paste0("[", x, "](https://www.ncbi.nlm.nih.gov/gene/",
            LigandGeneID, ") ![](figures/Ligand/",
            LigandGeneID, ".png)")
    }, "")    
    paste0(
    "|", Rank,
    "|", Freq,
    "|", Ligand,
    "|", Receptor,
    "|", CCI, "|\n", collapse="")
}

.LIGAND_BODY <- function(out.vecLR, GeneInfo, index, selected){
    # Setting
    convertGeneName <- function(geneid, GeneInfo){
        if(!is.null(GeneInfo$GeneName)){
            genename <- GeneInfo$GeneName[
                which(GeneInfo$GeneName$ENTREZID == geneid), "SYMBOL"]
            if(length(genename) == 0 || genename %in% c("", NA)){
                genename = geneid
            }
            if(length(genename) != 1){
                # Black list
                genename = setdiff(genename, "cxcl11.6")[1]
            }
            genename
        }else{
            geneid
        }
    }

    # Node
    nodes <- lapply(seq_len(ncol(out.vecLR)), function(x){
        names(out.vecLR["TARGET", x][[1]])
    })
    LnodesGeneID <- lapply(nodes, function(x){
        vapply(x, function(xx){
            strsplit(xx, "_")[[1]][1]
        }, "")
    })
    RnodesGeneID <- lapply(nodes, function(x){
        vapply(x, function(xx){
            strsplit(xx, "_")[[1]][2]
        }, "")
    })
    LnodesGeneName <- lapply(LnodesGeneID, function(x){
        vapply(x, function(xx){
            convertGeneName(xx, GeneInfo)
        }, "")
    })
    RnodesGeneName <- lapply(RnodesGeneID, function(x){
        vapply(x, function(xx){
            convertGeneName(xx, GeneInfo)
        }, "")
    })
    uniqueLnodesGeneName <- unique(unlist(LnodesGeneName))
    Freq <- vapply(uniqueLnodesGeneName, function(x){
            length(which(unlist(LnodesGeneName) == x))
        }, 0L)
    Freq <- sort(Freq, decreasing=TRUE)
    Rank <- rank(-Freq)
    Ligand <- names(Freq)
    Receptor <- lapply(Ligand, function(x){
        ReceptorGeneName <- unique(
            unlist(
                lapply(seq_len(length(LnodesGeneName)), function(xx){
                    RnodesGeneName[[xx]][which(x == LnodesGeneName[[xx]])]
                })
            )
        )
        ReceptorGeneID <- unlist(lapply(ReceptorGeneName, function(y){
            target <- which(GeneInfo$GeneName[, "SYMBOL"] == y)
            if(length(target) == 0){
                y
            }else{
                GeneInfo$GeneName[target[1], "ENTREZID"]
            }
        }))
        paste(
            paste0("[", ReceptorGeneName,
                "](https://www.ncbi.nlm.nih.gov/gene/",
            ReceptorGeneID, ") ![](figures/Receptor/",
            ReceptorGeneID, ".png)"), collapse=" ")
    })

    CCI <- vapply(Ligand, function(x){
        hit <- vapply(seq_len(ncol(out.vecLR)), function(xx){
            length(which(LnodesGeneName[[xx]] == x))
        }, 0L)
        vec <- seq_len(ncol(out.vecLR))[which(hit != 0)]
        target <- unlist(lapply(vec, function(xx){
                p <- gsub("pattern", "", colnames(out.vecLR)[xx])
                which(index[selected, "Mode3"] == p)
            }))
        out <- index[target, paste0("Mode", seq_len(3))]
        # Out
        if(is.vector(out)){
            paste0("[See the details of (",
                paste(out, collapse=","),
                ")-Pattern](",
                "pattern_",
                paste(out, collapse="_"),
                ".html)")
        }else{
            paste(vapply(seq_len(nrow(out)), function(xx){
                paste0("[See the details of (",
                    paste(out[xx,], collapse=","),
                    ")-Pattern](",
                    "pattern_",
                    paste(out[xx,], collapse="_"),
                    ".html)")
            }, ""), collapse=" ")
        }
    }, "")
    Ligand <- vapply(Ligand, function(x){
        target <- which(GeneInfo$GeneName$SYMBOL == x)
        if(length(target) == 0){
            LigandGeneID <- x
        }else{
            LigandGeneID <- GeneInfo$GeneName[target[1], "ENTREZID"]
        }
        paste0("[", x, "](https://www.ncbi.nlm.nih.gov/gene/",
            LigandGeneID, ") ![](figures/Ligand/",
            LigandGeneID, ".png)")
    }, "")
    paste0(
    "|", Rank,
    "|", Freq,
    "|", Ligand,
    "|", Receptor,
    "|", CCI, "|\n", collapse="")
}

.RECEPTOR_BODY_2 <- function(out.vecLR, GeneInfo, index, selected){
    # Setting
    convertGeneName <- function(geneid, GeneInfo){
        if(!is.null(GeneInfo$GeneName)){
            genename <- GeneInfo$GeneName[
                which(GeneInfo$GeneName$ENTREZID == geneid),
                "SYMBOL"]
            if(length(genename) == 0 || genename %in% c("", NA)){
                genename = geneid
            }
            if(length(genename) != 1){
                # Black list
                genename = setdiff(genename, "cxcl11.6")[1]
            }
            genename
        }else{
            geneid
        }
    }

    # Node
    nodes <- lapply(seq_len(length(out.vecLR)), function(x){
        names(out.vecLR[[x]]$TARGET)
    })
    LnodesGeneID <- lapply(nodes, function(x){
        vapply(x, function(xx){
            strsplit(xx, "_")[[1]][1]
        }, "")
    })
    RnodesGeneID <- lapply(nodes, function(x){
        vapply(x, function(xx){
            strsplit(xx, "_")[[1]][2]
        }, "")
    })
    LnodesGeneName <- lapply(LnodesGeneID, function(x){
        vapply(x, function(xx){
            convertGeneName(xx, GeneInfo)
        }, "")
    })
    RnodesGeneName <- lapply(RnodesGeneID, function(x){
        vapply(x, function(xx){
            convertGeneName(xx, GeneInfo)
        }, "")
    })
    uniqueRnodesGeneName <- unique(unlist(RnodesGeneName))
    Freq <- vapply(uniqueRnodesGeneName, function(x){
            length(which(unlist(RnodesGeneName) == x))
        }, 0L)
    Freq <- sort(Freq, decreasing=TRUE)
    Rank <- rank(-Freq)
    Receptor <- names(Freq)
    Ligand <- lapply(Receptor, function(x){
        LigandGeneName <- unique(
            unlist(
                lapply(seq_len(length(RnodesGeneName)), function(xx){
                    LnodesGeneName[[xx]][which(x == RnodesGeneName[[xx]])]
                })
            )
        )
        LigandGeneID <- unlist(lapply(LigandGeneName, function(y){
            target <- which(GeneInfo$GeneName$SYMBOL == y)
            if(length(target) == 0){
                y
            }else{
                GeneInfo$GeneName[target[1], "ENTREZID"]
            }
        }))
        paste(
            paste0("[", LigandGeneName,
                "](https://www.ncbi.nlm.nih.gov/gene/",
            LigandGeneID, ") ![](figures/Ligand/",
            LigandGeneID, ".png)"), collapse=" ")
    })

    CCI <- vapply(Receptor, function(x){
        hit <- vapply(seq_len(length(out.vecLR)), function(xx){
            length(which(RnodesGeneName[[xx]] == x))
        }, 0L)
        vec <- seq_len(length(out.vecLR))[which(hit != 0)]
        # Out
        cciout <- vapply(vec, function(xx){
            paste0("[See the details of (",
                paste(c(index[xx, paste0("Mode", 1:2)], ""), collapse=","),
                ")-related L-R patterns](",
                "pattern_",
                paste(index[xx, paste0("Mode", 1:2)], collapse="_"),
                ".html)")
        }, "")
        paste(unique(cciout), collapse=" ")
    }, "")
    Receptor <- vapply(Receptor, function(x){
        target <- which(GeneInfo$GeneName$SYMBOL == x)
        if(length(target) == 0){
            ReceptorGeneID <- x
        }else{
            ReceptorGeneID <- GeneInfo$GeneName[target[1], "ENTREZID"]
        }
        paste0("[", x, "](https://www.ncbi.nlm.nih.gov/gene/",
            ReceptorGeneID, ") ![](figures/Receptor/",
            ReceptorGeneID, ".png)")
    }, "")
    paste0(
    "|", Rank,
    "|", Freq,
    "|", Receptor,
    "|", Ligand,
    "|", CCI, "|\n", collapse="")
}

.RECEPTOR_BODY <- function(out.vecLR, GeneInfo, index, selected){
    # Setting
    convertGeneName <- function(geneid, GeneInfo){
        if(!is.null(GeneInfo$GeneName)){
            genename <- GeneInfo$GeneName[
                which(GeneInfo$GeneName$ENTREZID == geneid),
                "SYMBOL"]
            if(length(genename) == 0 || genename %in% c("", NA)){
                genename = geneid
            }
            if(length(genename) != 1){
                # Black list
                genename = setdiff(genename, "cxcl11.6")[1]
            }
            genename
        }else{
            geneid
        }
    }

    # Node
    nodes <- lapply(seq_len(ncol(out.vecLR)), function(x){
        names(out.vecLR["TARGET", x][[1]])
    })
    LnodesGeneID <- lapply(nodes, function(x){
        vapply(x, function(xx){
            strsplit(xx, "_")[[1]][1]
        }, "")
    })
    RnodesGeneID <- lapply(nodes, function(x){
        vapply(x, function(xx){
            strsplit(xx, "_")[[1]][2]
        }, "")
    })
    LnodesGeneName <- lapply(LnodesGeneID, function(x){
        vapply(x, function(xx){
            convertGeneName(xx, GeneInfo)
        }, "")
    })
    RnodesGeneName <- lapply(RnodesGeneID, function(x){
        vapply(x, function(xx){
            convertGeneName(xx, GeneInfo)
        }, "")
    })
    uniqueRnodesGeneName <- unique(unlist(RnodesGeneName))
    Freq <- vapply(uniqueRnodesGeneName, function(x){
            length(which(unlist(RnodesGeneName) == x))
        }, 0L)
    Freq <- sort(Freq, decreasing=TRUE)
    Rank <- rank(-Freq)
    Receptor <- names(Freq)
    Ligand <- lapply(Receptor, function(x){
        LigandGeneName <- unique(
            unlist(
                lapply(seq_len(length(RnodesGeneName)), function(xx){
                    LnodesGeneName[[xx]][which(x == RnodesGeneName[[xx]])]
                })
            )
        )
        LigandGeneID <- unlist(lapply(LigandGeneName, function(y){
            target <- which(GeneInfo$GeneName$SYMBOL == y)
            if(length(target) == 0){
                y
            }else{
                GeneInfo$GeneName[target[1], "ENTREZID"]
            }
        }))
        paste(
            paste0("[", LigandGeneName,
                "](https://www.ncbi.nlm.nih.gov/gene/",
            LigandGeneID, ") ![](figures/Ligand/",
            LigandGeneID, ".png)"), collapse=" ")
    })

    CCI <- vapply(Receptor, function(x){
        hit <- vapply(seq_len(ncol(out.vecLR)), function(xx){
            length(which(RnodesGeneName[[xx]] == x))
        }, 0L)
        vec <- seq_len(ncol(out.vecLR))[which(hit != 0)]
        target <- unlist(lapply(vec, function(xx){
                p <- gsub("pattern", "", colnames(out.vecLR)[xx])
                which(index[selected, "Mode3"] == p)
            }))
        out <- index[target, c("Mode1", "Mode2", "Mode3")]
        # Out
        if(is.vector(out)){
            paste0("[See the details of (",
                paste(out, collapse=","),
                ")-Pattern](",
                "pattern_",
                paste(out, collapse="_"),
                ".html)")
        }else{
            paste(vapply(seq_len(nrow(out)), function(xx){
                paste0("[See the details of (",
                    paste(out[xx,], collapse=","),
                    ")-Pattern](",
                    "pattern_",
                    paste(out[xx,], collapse="_"),
                    ".html)")
            }, ""), collapse=" ")
        }
    }, "")
    Receptor <- vapply(Receptor, function(x){
        target <- which(GeneInfo$GeneName$SYMBOL == x)
        if(length(target) == 0){
            ReceptorGeneID <- x
        }else{
            ReceptorGeneID <- GeneInfo$GeneName[target[1], "ENTREZID"]
        }
        paste0("[", x, "](https://www.ncbi.nlm.nih.gov/gene/",
            ReceptorGeneID, ") ![](figures/Receptor/",
            ReceptorGeneID, ".png)")
    }, "")
    paste0(
    "|", Rank,
    "|", Freq,
    "|", Receptor,
    "|", Ligand,
    "|", CCI, "|\n", collapse="")
}

.LIGANDALL_HEADER <- paste0(
    "---\n",
    "title: <font color='#1881c2'>Details of Ligand Gene-centric Overview (all)",
    "</font>\n",
    "---\n\n",
    "<style type='text/css'>\n",
    "table,th,td {\n",
    "width: 50%;\n",
    "border: 1px solid #f0f0f0;\n",
    "}\n",
    "</style>\n\n",
    "|Ligand Gene|Receptor Gene|\n",
    "|---------------------------------|---------------------------------|"
)

.RECEPTORALL_HEADER <- paste0(
    "---\n",
    "title: <font color='#1881c2'>Details of Receptor Gene-centric Overview (all)",
    "</font>\n",
    "---\n\n",
    "<style type='text/css'>\n",
    "table,th,td {\n",
    "width: 50%;\n",
    "border: 1px solid #f0f0f0;\n",
    "}\n",
    "</style>\n\n",
    "|Receptor Gene|Ligand Gene|\n",
    "|---------------------------------|---------------------------------|"
)

.LIGANDALL_BODY <- function(GeneInfo, LR, input){
    if(!is.null(GeneInfo$GeneName)){
        GeneName <- GeneInfo$GeneName
        LigandGeneID <- unique(LR$GENEID_L)
        LigandGeneName <- vapply(LigandGeneID, function(x){
            GeneName[which(GeneName$ENTREZID == x)[1], "SYMBOL"]}, "")
        naposition <- which(vapply(LigandGeneName, is.na, TRUE))
        LigandGeneName[naposition] <- LigandGeneID[naposition]
        # Sort by Alphabet of the ligand genes
        orderLigand <- order(LigandGeneName)
        LigandGeneID <- LigandGeneID[orderLigand]
        LigandGeneName <- LigandGeneName[orderLigand]
    }else{
        LigandGeneName <- LigandGeneID  
    }

    # Ligand Link
    Ligand <- vapply(seq_along(LigandGeneID), function(x){
        target <- which(rownames(input) == LigandGeneID[x])
        if(length(target) != 0){
            paste0("[", LigandGeneName[x],
                "](https://www.ncbi.nlm.nih.gov/gene/",
                LigandGeneID[x], ") ",
                "![](figures/Ligand/", LigandGeneID[x],
                ".png)")
        }else{
            paste0("[", LigandGeneName[x],
                "](https://www.ncbi.nlm.nih.gov/gene/",
                LigandGeneID[x], ")")
        }
    }, "")

    # Receptor Link
    Receptor <- vapply(seq_along(LigandGeneID), function(x){
        target <- which(LR$GENEID_L == LigandGeneID[x])
        ReceptorGeneID <- unique(LR$GENEID_R[target])
        if(!is.null(GeneInfo$GeneName)){
            ReceptorGeneName <- vapply(ReceptorGeneID, function(xx){
                GeneName[which(GeneName$ENTREZID == xx)[1], "SYMBOL"]
            }, "")
            naposition <- which(vapply(ReceptorGeneName, is.na, TRUE))
            ReceptorGeneName[naposition] <- ReceptorGeneID[naposition]
        }else{
            ReceptorGeneName <- ReceptorGeneID
        }
        paste(paste0("[", ReceptorGeneName,
            "](https://www.ncbi.nlm.nih.gov/gene/",
            ReceptorGeneID, ")"),
            collapse=" ")
    }, "")
    # Output
    paste0(
    "|", Ligand,
    "|", Receptor, "|\n", collapse="")
}

.RECEPTORALL_BODY <- function(GeneInfo, LR, input){
    if(!is.null(GeneInfo$GeneName)){
        GeneName <- GeneInfo$GeneName
        ReceptorGeneID <- unique(LR$GENEID_R)
        ReceptorGeneName <- vapply(ReceptorGeneID, function(x){
            GeneName[which(GeneName$ENTREZID == x)[1], "SYMBOL"]}, "")
        naposition <- which(vapply(ReceptorGeneName, is.na, TRUE))
        ReceptorGeneName[naposition] <- ReceptorGeneID[naposition]
        # Sort by Alphabet of the Receptor genes
        orderReceptor <- order(ReceptorGeneName)
        ReceptorGeneID <- ReceptorGeneID[orderReceptor]
        ReceptorGeneName <- ReceptorGeneName[orderReceptor]
    }else{
        ReceptorGeneName <- ReceptorGeneID
    }

    # Receptor Link
    Receptor <- vapply(seq_along(ReceptorGeneID), function(x){
        target <- which(rownames(input) == ReceptorGeneID[x])
        if(length(target) != 0){
            paste0("[", ReceptorGeneName[x],
                "](https://www.ncbi.nlm.nih.gov/gene/",
                ReceptorGeneID[x], ") ",
                "![](figures/Receptor/", ReceptorGeneID[x],
                ".png)")
        }else{
            paste0("[", ReceptorGeneName[x],
                "](https://www.ncbi.nlm.nih.gov/gene/",
                ReceptorGeneID[x], ")")
        }
    }, "")

    # Ligand Link
    Ligand <- vapply(seq_along(ReceptorGeneID), function(x){
        target <- which(LR$GENEID_R == ReceptorGeneID[x])
        LigandGeneID <- unique(LR$GENEID_L[target])
        if(!is.null(GeneInfo$GeneName)){
            LigandGeneName <- vapply(LigandGeneID, function(xx){
                GeneName[which(GeneName$ENTREZID == xx)[1], "SYMBOL"]
            }, "")
            naposition <- which(vapply(LigandGeneName, is.na, TRUE))
            LigandGeneName[naposition] <- LigandGeneID[naposition]
        }else{
            LigandGeneName <- LigandGeneID
        }
        paste(paste0("[", LigandGeneName,
            "](https://www.ncbi.nlm.nih.gov/gene/",
            LigandGeneID, ")"),
            collapse=" ")
    }, "")

    # Output
    paste0(
    "|", Receptor,
    "|", Ligand, "|\n", collapse="")
}

.XYZ_ENRICH_2 <- function(out.vecLR, i){
    go_bp <- .check_NULL_GO_BP_2(out.vecLR, i)
    go_mf <- .check_NULL_GO_MF_2(out.vecLR, i)
    go_cc <- .check_NULL_GO_CC_2(out.vecLR, i)
    mesh_a <- .check_NULL_MeSH_A_2(out.vecLR, i)
    mesh_b <- .check_NULL_MeSH_B_2(out.vecLR, i)
    mesh_c <- .check_NULL_MeSH_C_2(out.vecLR, i)
    mesh_d <- .check_NULL_MeSH_D_2(out.vecLR, i)
    mesh_e <- .check_NULL_MeSH_E_2(out.vecLR, i)
    mesh_f <- .check_NULL_MeSH_F_2(out.vecLR, i)
    mesh_g <- .check_NULL_MeSH_G_2(out.vecLR, i)
    mesh_h <- .check_NULL_MeSH_H_2(out.vecLR, i)
    mesh_i <- .check_NULL_MeSH_I_2(out.vecLR, i)
    mesh_j <- .check_NULL_MeSH_J_2(out.vecLR, i)
    mesh_k <- .check_NULL_MeSH_K_2(out.vecLR, i)
    mesh_l <- .check_NULL_MeSH_L_2(out.vecLR, i)
    mesh_m <- .check_NULL_MeSH_M_2(out.vecLR, i)
    mesh_n <- .check_NULL_MeSH_N_2(out.vecLR, i)
    mesh_v <- .check_NULL_MeSH_V_2(out.vecLR, i)
    mesh_z <- .check_NULL_MeSH_Z_2(out.vecLR, i)
    reactome <- .check_NULL_Reactome_2(out.vecLR, i)
    do <- .check_NULL_DO_2(out.vecLR, i)
    ncg <- .check_NULL_NCG_2(out.vecLR, i)
    dgn <- .check_NULL_DGN_2(out.vecLR, i)
    # Output
    paste(c(go_bp, go_mf, go_cc,
    mesh_a, mesh_b, mesh_c, mesh_d, mesh_e, mesh_f,
    mesh_g, mesh_h, mesh_i, mesh_j, mesh_k, mesh_l,
    mesh_m, mesh_n, mesh_v, mesh_z,
    reactome,
    do, ncg, dgn), collapse="")
}

################ GO_BP ################
.check_NULL_GO_BP_2 <- function(out.vecLR, i){
    if(!is.null(out.vecLR[[i]]$Enrich$GO_BP$Pvalue)){
        paste0("## <font color='#1881c2'>GO-Enrichment Analysis (BP : Biological Process)</font>\n\n",
            "```{r}\n", # Top
            "negLogPval <- -log10(out.vecLR[['",
            i,
            "']]$Enrich$GO_BP$Pvalue + 1E-10)\n",
            "term <- out.vecLR[['",
            i,
            "']]$Enrich$GO_BP$Term\n",
            "target <- seq_len(min(100, length(negLogPval)))\n",
            "negLogPval <- negLogPval[target]\n",
            "term <- term[target]\n",
            "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
            "type=\"bar\", color=~negLogPval, text=term,\n",
            "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
            "layout(p, title=\"GO-Enrichment Analysis (BP)\",\n",
            "   xaxis=list(title=\"Term\"),\n",
            "   yaxis=list(title=\"-Log(P-value)\"),\n",
            "   legend=list(name=\"\"))\n",
            "```\n\n", # Bottom
            "![](figures/Tagcloud/GO_BP_",
            i,
            ".png)\n\n")
    }else{
        ""
    }
}

################ GO_MF ################
.check_NULL_GO_MF_2 <- function(out.vecLR, i){
    if(!is.null(out.vecLR[[i]]$Enrich$GO_MF$Pvalue)){
        paste0("## <font color='#1881c2'>GO-Enrichment Analysis (MF : Molecular Function)</font>\n\n",
            "```{r}\n", # Top
            "negLogPval <- -log10(out.vecLR[['",
            i,
            "']]$Enrich$GO_MF$Pvalue + 1E-10)\n",
            "term <- out.vecLR[['",
            i,
            "']]$Enrich$GO_MF$Term\n",
            "target <- seq_len(min(100, length(negLogPval)))\n",
            "negLogPval <- negLogPval[target]\n",
            "term <- term[target]\n",
            "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
            "type=\"bar\", color=~negLogPval, text=term,\n",
            "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
            "layout(p, title=\"GO-Enrichment Analysis (MF)\",\n",
            "   xaxis=list(title=\"Term\"),\n",
            "   yaxis=list(title=\"-Log(P-value)\"),\n",
            "   legend=list(name=\"\"))\n",
            "```\n\n", # Bottom
            "![](figures/Tagcloud/GO_MF_",
            i,
            ".png)\n\n")
    }else{
        ""
    }
}

################ GO_CC ################
.check_NULL_GO_CC_2 <- function(out.vecLR, i){
    if(!is.null(out.vecLR[[i]]$Enrich$GO_CC$Pvalue)){
        paste0("## <font color='#1881c2'>GO-Enrichment Analysis (CC : Cellular Component)</font>\n\n",
            "```{r}\n", # Top
            "negLogPval <- -log10(out.vecLR[['",
            i,
            "']]$Enrich$GO_CC$Pvalue + 1E-10)\n",
            "term <- out.vecLR[['",
            i,
            "']]$Enrich$GO_CC$Term\n",
            "target <- seq_len(min(100, length(negLogPval)))\n",
            "negLogPval <- negLogPval[target]\n",
            "term <- term[target]\n",
            "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
            "type=\"bar\", color=~negLogPval, text=term,\n",
            "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
            "layout(p, title=\"GO-Enrichment Analysis (CC)\",\n",
            "   xaxis=list(title=\"Term\"),\n",
            "   yaxis=list(title=\"-Log(P-value)\"),\n",
            "   legend=list(name=\"\"))\n",
            "```\n\n", # Bottom
            "![](figures/Tagcloud/GO_CC_",
            i,
            ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_A ################
.check_NULL_MeSH_A_2 <- function(out.vecLR, i){
    if(!is.null(out.vecLR[[i]]$Enrich$MeSH_A$Pvalue)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (A : Anatomy)</font>\n\n",
            "```{r}\n", # Top
            "negLogPval <- -log10(out.vecLR[['",
            i,
            "']]$Enrich$MeSH_A$Pvalue + 1E-10)\n",
            "term <- out.vecLR[['",
            i,
            "']]$Enrich$MeSH_A$Term\n",
            "target <- seq_len(min(100, length(negLogPval)))\n",
            "negLogPval <- negLogPval[target]\n",
            "term <- term[target]\n",
            "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
            "type=\"bar\", color=~negLogPval, text=term,\n",
            "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
            "layout(p, title=\"MeSH-Enrichment Analysis (A)\",\n",
            "   xaxis=list(title=\"Term\"),\n",
            "   yaxis=list(title=\"-Log(P-value)\"),\n",
            "   legend=list(name=\"\"))\n",
            "```\n\n", # Bottom
            "![](figures/Tagcloud/MeSH_A_",
            i,
            ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_B ################
.check_NULL_MeSH_B_2 <- function(out.vecLR, i){
    if(!is.null(out.vecLR[[i]]$Enrich$MeSH_B$Pvalue)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (B : Organisms)</font>\n\n",
            "```{r}\n", # Top
            "negLogPval <- -log10(out.vecLR[['",
            i,
            "']]$Enrich$MeSH_B$Pvalue + 1E-10)\n",
            "term <- out.vecLR[['",
            i,
            "']]$Enrich$MeSH_B$Term\n",
            "target <- seq_len(min(100, length(negLogPval)))\n",
            "negLogPval <- negLogPval[target]\n",
            "term <- term[target]\n",
            "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
            "type=\"bar\", color=~negLogPval, text=term,\n",
            "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
            "layout(p, title=\"MeSH-Enrichment Analysis (B)\",\n",
            "   xaxis=list(title=\"Term\"),\n",
            "   yaxis=list(title=\"-Log(P-value)\"),\n",
            "   legend=list(name=\"\"))\n",
            "```\n\n", # Bottom
            "![](figures/Tagcloud/MeSH_B_",
            i,
            ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_C ################
.check_NULL_MeSH_C_2 <- function(out.vecLR, i){
    if(!is.null(out.vecLR[[i]]$Enrich$MeSH_C$Pvalue)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (C : Diseases)</font>\n\n",
            "```{r}\n", # Top
            "negLogPval <- -log10(out.vecLR[['",
            i,
            "']]$Enrich$MeSH_C$Pvalue + 1E-10)\n",
            "term <- out.vecLR[['",
            i,
            "']]$Enrich$MeSH_C$Term\n",
            "target <- seq_len(min(100, length(negLogPval)))\n",
            "negLogPval <- negLogPval[target]\n",
            "term <- term[target]\n",
            "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
            "type=\"bar\", color=~negLogPval, text=term,\n",
            "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
            "layout(p, title=\"MeSH-Enrichment Analysis (C)\",\n",
            "   xaxis=list(title=\"Term\"),\n",
            "   yaxis=list(title=\"-Log(P-value)\"),\n",
            "   legend=list(name=\"\"))\n",
            "```\n\n", # Bottom
            "![](figures/Tagcloud/MeSH_C_",
            i,
            ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_D ################
.check_NULL_MeSH_D_2 <- function(out.vecLR, i){
    if(!is.null(out.vecLR[[i]]$Enrich$MeSH_D$Pvalue)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (D : Drugs)</font>\n\n",
            "```{r}\n", # Top
            "negLogPval <- -log10(out.vecLR[['",
            i,
            "']]$Enrich$MeSH_D$Pvalue + 1E-10)\n",
            "term <- out.vecLR[['",
            i,
            "']]$Enrich$MeSH_D$Term\n",
            "target <- seq_len(min(100, length(negLogPval)))\n",
            "negLogPval <- negLogPval[target]\n",
            "term <- term[target]\n",
            "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
            "type=\"bar\", color=~negLogPval, text=term,\n",
            "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
            "layout(p, title=\"MeSH-Enrichment Analysis (D)\",\n",
            "   xaxis=list(title=\"Term\"),\n",
            "   yaxis=list(title=\"-Log(P-value)\"),\n",
            "   legend=list(name=\"\"))\n",
            "```\n\n", # Bottom
            "![](figures/Tagcloud/MeSH_D_",
            i,
            ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_E ################
.check_NULL_MeSH_E_2 <- function(out.vecLR, i){
    if(!is.null(out.vecLR[[i]]$Enrich$MeSH_E$Pvalue)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (E : Analytical, Diagnostic and Therapeutic Techniques and Equipment)</font>\n\n",
            "```{r}\n", # Top
            "negLogPval <- -log10(out.vecLR[['",
            i,
            "']]$Enrich$MeSH_E$Pvalue + 1E-10)\n",
            "term <- out.vecLR[['",
            i,
            "']]$Enrich$MeSH_E$Term\n",
            "target <- seq_len(min(100, length(negLogPval)))\n",
            "negLogPval <- negLogPval[target]\n",
            "term <- term[target]\n",
            "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
            "type=\"bar\", color=~negLogPval, text=term,\n",
            "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
            "layout(p, title=\"MeSH-Enrichment Analysis (E)\",\n",
            "   xaxis=list(title=\"Term\"),\n",
            "   yaxis=list(title=\"-Log(P-value)\"),\n",
            "   legend=list(name=\"\"))\n",
            "```\n\n", # Bottom
            "![](figures/Tagcloud/MeSH_E_",
            i,
            ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_F ################
.check_NULL_MeSH_F_2 <- function(out.vecLR, i){
    if(!is.null(out.vecLR[[i]]$Enrich$MeSH_F$Pvalue)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (F : Psychiatry and Psychology)</font>\n\n",
            "```{r}\n", # Top
            "negLogPval <- -log10(out.vecLR[['",
            i,
            "']]$Enrich$MeSH_F$Pvalue + 1E-10)\n",
            "term <- out.vecLR[['",
            i,
            "']]$Enrich$MeSH_F$Term\n",
            "target <- seq_len(min(100, length(negLogPval)))\n",
            "negLogPval <- negLogPval[target]\n",
            "term <- term[target]\n",
            "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
            "type=\"bar\", color=~negLogPval, text=term,\n",
            "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
            "layout(p, title=\"MeSH-Enrichment Analysis (F)\",\n",
            "   xaxis=list(title=\"Term\"),\n",
            "   yaxis=list(title=\"-Log(P-value)\"),\n",
            "   legend=list(name=\"\"))\n",
            "```\n\n", # Bottom
            "![](figures/Tagcloud/MeSH_F_",
            i,
            ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_G ################
.check_NULL_MeSH_G_2 <- function(out.vecLR, i){
    if(!is.null(out.vecLR[[i]]$Enrich$MeSH_G$Pvalue)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (G : Phenomena and Processes)</font>\n\n",
            "```{r}\n", # Top
            "negLogPval <- -log10(out.vecLR[['",
            i,
            "']]$Enrich$MeSH_G$Pvalue + 1E-10)\n",
            "term <- out.vecLR[['",
            i,
            "']]$Enrich$MeSH_G$Term\n",
            "target <- seq_len(min(100, length(negLogPval)))\n",
            "negLogPval <- negLogPval[target]\n",
            "term <- term[target]\n",
            "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
            "type=\"bar\", color=~negLogPval, text=term,\n",
            "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
            "layout(p, title=\"MeSH-Enrichment Analysis (G)\",\n",
            "   xaxis=list(title=\"Term\"),\n",
            "   yaxis=list(title=\"-Log(P-value)\"),\n",
            "   legend=list(name=\"\"))\n",
            "```\n\n", # Bottom
            "![](figures/Tagcloud/MeSH_G_",
            i,
            ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_H ################
.check_NULL_MeSH_H_2 <- function(out.vecLR, i){
    if(!is.null(out.vecLR[[i]]$Enrich$MeSH_H$Pvalue)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (H : Disciplines and Occupations)</font>\n\n",
            "```{r}\n", # Top
            "negLogPval <- -log10(out.vecLR[['",
            i,
            "']]$Enrich$MeSH_H$Pvalue + 1E-10)\n",
            "term <- out.vecLR[['",
            i,
            "']]$Enrich$MeSH_H$Term\n",
            "target <- seq_len(min(100, length(negLogPval)))\n",
            "negLogPval <- negLogPval[target]\n",
            "term <- term[target]\n",
            "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
            "type=\"bar\", color=~negLogPval, text=term,\n",
            "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
            "layout(p, title=\"MeSH-Enrichment Analysis (H)\",\n",
            "   xaxis=list(title=\"Term\"),\n",
            "   yaxis=list(title=\"-Log(P-value)\"),\n",
            "   legend=list(name=\"\"))\n",
            "```\n\n", # Bottom
            "![](figures/Tagcloud/MeSH_H_",
            i,
            ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_I ################
.check_NULL_MeSH_I_2 <- function(out.vecLR, i){
    if(!is.null(out.vecLR[[i]]$Enrich$MeSH_I$Pvalue)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (I : Anthropology, Education, Sociology and Social Phenomena)</font>\n\n",
            "```{r}\n", # Top
            "negLogPval <- -log10(out.vecLR[['",
            i,
            "']]$Enrich$MeSH_I$Pvalue + 1E-10)\n",
            "term <- out.vecLR[['",
            i,
            "']]$Enrich$MeSH_I$Term\n",
            "target <- seq_len(min(100, length(negLogPval)))\n",
            "negLogPval <- negLogPval[target]\n",
            "term <- term[target]\n",
            "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
            "type=\"bar\", color=~negLogPval, text=term,\n",
            "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
            "layout(p, title=\"MeSH-Enrichment Analysis (I)\",\n",
            "   xaxis=list(title=\"Term\"),\n",
            "   yaxis=list(title=\"-Log(P-value)\"),\n",
            "   legend=list(name=\"\"))\n",
            "```\n\n", # Bottom
            "![](figures/Tagcloud/MeSH_I_",
            i,
            ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_J ################
.check_NULL_MeSH_J_2 <- function(out.vecLR, i){
    if(!is.null(out.vecLR[[i]]$Enrich$MeSH_J$Pvalue)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (J : Technology and Food and Beverages)</font>\n\n",
            "```{r}\n", # Top
            "negLogPval <- -log10(out.vecLR[['",
            i,
            "']]$Enrich$MeSH_J$Pvalue + 1E-10)\n",
            "term <- out.vecLR[['",
            i,
            "']]$Enrich$MeSH_J$Term\n",
            "target <- seq_len(min(100, length(negLogPval)))\n",
            "negLogPval <- negLogPval[target]\n",
            "term <- term[target]\n",
            "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
            "type=\"bar\", color=~negLogPval, text=term,\n",
            "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
            "layout(p, title=\"MeSH-Enrichment Analysis (J)\",\n",
            "   xaxis=list(title=\"Term\"),\n",
            "   yaxis=list(title=\"-Log(P-value)\"),\n",
            "   legend=list(name=\"\"))\n",
            "```\n\n", # Bottom
            "![](figures/Tagcloud/MeSH_J_",
            i,
            ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_K ################
.check_NULL_MeSH_K_2 <- function(out.vecLR, i){
    if(!is.null(out.vecLR[[i]]$Enrich$MeSH_K$Pvalue)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (K : Humanities)</font>\n\n",
            "```{r}\n", # Top
            "negLogPval <- -log10(out.vecLR[['",
            i,
            "']]$Enrich$MeSH_K$Pvalue + 1E-10)\n",
            "term <- out.vecLR[['",
            i,
            "']]$Enrich$MeSH_K$Term\n",
            "target <- seq_len(min(100, length(negLogPval)))\n",
            "negLogPval <- negLogPval[target]\n",
            "term <- term[target]\n",
            "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
            "type=\"bar\", color=~negLogPval, text=term,\n",
            "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
            "layout(p, title=\"MeSH-Enrichment Analysis (K)\",\n",
            "   xaxis=list(title=\"Term\"),\n",
            "   yaxis=list(title=\"-Log(P-value)\"),\n",
            "   legend=list(name=\"\"))\n",
            "```\n\n", # Bottom
            "![](figures/Tagcloud/MeSH_K_",
            i,
            ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_L ################
.check_NULL_MeSH_L_2 <- function(out.vecLR, i){
    if(!is.null(out.vecLR[[i]]$Enrich$MeSH_L$Pvalue)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (L : Information Science)</font>\n\n",
            "```{r}\n", # Top
            "negLogPval <- -log10(out.vecLR[['",
            i,
            "']]$Enrich$MeSH_L$Pvalue + 1E-10)\n",
            "term <- out.vecLR[['",
            i,
            "']]$Enrich$MeSH_L$Term\n",
            "target <- seq_len(min(100, length(negLogPval)))\n",
            "negLogPval <- negLogPval[target]\n",
            "term <- term[target]\n",
            "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
            "type=\"bar\", color=~negLogPval, text=term,\n",
            "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
            "layout(p, title=\"MeSH-Enrichment Analysis (L)\",\n",
            "   xaxis=list(title=\"Term\"),\n",
            "   yaxis=list(title=\"-Log(P-value)\"),\n",
            "   legend=list(name=\"\"))\n",
            "```\n\n", # Bottom
            "![](figures/Tagcloud/MeSH_L_",
            i,
            ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_M ################
.check_NULL_MeSH_M_2 <- function(out.vecLR, i){
    if(!is.null(out.vecLR[[i]]$Enrich$MeSH_M$Pvalue)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (M : Persons)</font>\n\n",
            "```{r}\n", # Top
            "negLogPval <- -log10(out.vecLR[['",
            i,
            "']]$Enrich$MeSH_M$Pvalue + 1E-10)\n",
            "term <- out.vecLR[['",
            i,
            "']]$Enrich$MeSH_M$Term\n",
            "target <- seq_len(min(100, length(negLogPval)))\n",
            "negLogPval <- negLogPval[target]\n",
            "term <- term[target]\n",
            "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
            "type=\"bar\", color=~negLogPval, text=term,\n",
            "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
            "layout(p, title=\"MeSH-Enrichment Analysis (M)\",\n",
            "   xaxis=list(title=\"Term\"),\n",
            "   yaxis=list(title=\"-Log(P-value)\"),\n",
            "   legend=list(name=\"\"))\n",
            "```\n\n", # Bottom
            "![](figures/Tagcloud/MeSH_M_",
            i,
            ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_N ################
.check_NULL_MeSH_N_2 <- function(out.vecLR, i){
    if(!is.null(out.vecLR[[i]]$Enrich$MeSH_N$Pvalue)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (N : Health Care)</font>\n\n",
            "```{r}\n", # Top
            "negLogPval <- -log10(out.vecLR[['",
            i,
            "']]$Enrich$MeSH_N$Pvalue + 1E-10)\n",
            "term <- out.vecLR[['",
            i,
            "']]$Enrich$MeSH_N$Term\n",
            "target <- seq_len(min(100, length(negLogPval)))\n",
            "negLogPval <- negLogPval[target]\n",
            "term <- term[target]\n",
            "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
            "type=\"bar\", color=~negLogPval, text=term,\n",
            "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
            "layout(p, title=\"MeSH-Enrichment Analysis (N)\",\n",
            "   xaxis=list(title=\"Term\"),\n",
            "   yaxis=list(title=\"-Log(P-value)\"),\n",
            "   legend=list(name=\"\"))\n",
            "```\n\n", # Bottom
            "![](figures/Tagcloud/MeSH_N_",
            i,
            ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_V ################
.check_NULL_MeSH_V_2 <- function(out.vecLR, i){
    if(!is.null(out.vecLR[[i]]$Enrich$MeSH_V$Pvalue)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (V : Publication Type)</font>\n\n",
            "```{r}\n", # Top
            "negLogPval <- -log10(out.vecLR[['",
            i,
            "']]$Enrich$MeSH_V$Pvalue + 1E-10)\n",
            "term <- out.vecLR[['",
            i,
            "']]$Enrich$MeSH_V$Term\n",
            "target <- seq_len(min(100, length(negLogPval)))\n",
            "negLogPval <- negLogPval[target]\n",
            "term <- term[target]\n",
            "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
            "type=\"bar\", color=~negLogPval, text=term,\n",
            "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
            "layout(p, title=\"MeSH-Enrichment Analysis (V)\",\n",
            "   xaxis=list(title=\"Term\"),\n",
            "   yaxis=list(title=\"-Log(P-value)\"),\n",
            "   legend=list(name=\"\"))\n",
            "```\n\n", # Bottom
            "![](figures/Tagcloud/MeSH_V_",
            i,
            ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_Z ################
.check_NULL_MeSH_Z_2 <- function(out.vecLR, i){
    if(!is.null(out.vecLR[[i]]$Enrich$MeSH_Z$Pvalue)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (Z : Geographical Locations)</font>\n\n",
            "```{r}\n", # Top
            "negLogPval <- -log10(out.vecLR[['",
            i,
            "']]$Enrich$MeSH_Z$Pvalue + 1E-10)\n",
            "term <- out.vecLR[['",
            i,
            "']]$Enrich$MeSH_Z$Term\n",
            "target <- seq_len(min(100, length(negLogPval)))\n",
            "negLogPval <- negLogPval[target]\n",
            "term <- term[target]\n",
            "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
            "type=\"bar\", color=~negLogPval, text=term,\n",
            "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
            "layout(p, title=\"MeSH-Enrichment Analysis (Z)\",\n",
            "   xaxis=list(title=\"Term\"),\n",
            "   yaxis=list(title=\"-Log(P-value)\"),\n",
            "   legend=list(name=\"\"))\n",
            "```\n\n", # Bottom
            "![](figures/Tagcloud/MeSH_Z_",
            i,
            ".png)\n\n")
    }else{
        ""
    }
}

################ Reactome ################
.check_NULL_Reactome_2 <- function(out.vecLR, i){
    if(!is.null(out.vecLR[[i]]$Enrich$Reactome$Pvalue)){
        paste0("## <font color='#1881c2'>Reactome Pathway-Enrichment Analysis</font>\n\n",
            "```{r}\n", # Top
            "negLogPval <- -log10(out.vecLR[['",
            i,
            "']]$Enrich$Reactome$Pvalue + 1E-10)\n",
            "term <- out.vecLR[['",
            i,
            "']]$Enrich$Reactome$Term\n",
            "target <- seq_len(min(100, length(negLogPval)))\n",
            "negLogPval <- negLogPval[target]\n",
            "term <- term[target]\n",
            "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
            "type=\"bar\", color=~negLogPval, text=term,\n",
            "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
            "layout(p, title=\"Reactome-Enrichment Analysis\",\n",
            "   xaxis=list(title=\"Term\"),\n",
            "   yaxis=list(title=\"-Log(P-value)\"),\n",
            "   legend=list(name=\"\"))\n",
            "```\n\n", # Bottom
            "![](figures/Tagcloud/Reactome_",
            i,
            ".png)\n\n")
    }else{
        ""
    }
}

################ DO ################
.check_NULL_DO_2 <- function(out.vecLR, i){
    if(!is.null(out.vecLR[[i]]$Enrich$DO$Pvalue)){
        paste0("## <font color='#1881c2'>DO-Enrichment Analysis</font>\n\n",
            "```{r}\n", # Top
            "negLogPval <- -log10(out.vecLR[['",
            i,
            "']]$Enrich$DO$Pvalue + 1E-10)\n",
            "term <- out.vecLR[['",
            i,
            "']]$Enrich$DO$Term\n",
            "target <- seq_len(min(100, length(negLogPval)))\n",
            "negLogPval <- negLogPval[target]\n",
            "term <- term[target]\n",
            "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
            "type=\"bar\", color=~negLogPval, text=term,\n",
            "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
            "layout(p, title=\"DO-Enrichment Analysis\",\n",
            "   xaxis=list(title=\"Term\"),\n",
            "   yaxis=list(title=\"-Log(P-value)\"),\n",
            "   legend=list(name=\"\"))\n",
            "```\n\n", # Bottom
            "![](figures/Tagcloud/DO_",
            i,
            ".png)\n\n")
    }else{
        ""
    }
}

################ NCG ################
.check_NULL_NCG_2 <- function(out.vecLR, i){
    if(!is.null(out.vecLR[[i]]$Enrich$NCG$Pvalue)){
        paste0("## <font color='#1881c2'>NCG-Enrichment Analysis</font>\n\n",
            "```{r}\n", # Top
            "negLogPval <- -log10(out.vecLR[['",
            i,
            "']]$Enrich$NCG$Pvalue + 1E-10)\n",
            "term <- out.vecLR[['",
            i,
            "']]$Enrich$NCG$Term\n",
            "target <- seq_len(min(100, length(negLogPval)))\n",
            "negLogPval <- negLogPval[target]\n",
            "term <- term[target]\n",
            "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
            "type=\"bar\", color=~negLogPval, text=term,\n",
            "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
            "layout(p, title=\"NCG-Enrichment Analysis\",\n",
            "   xaxis=list(title=\"Term\"),\n",
            "   yaxis=list(title=\"-Log(P-value)\"),\n",
            "   legend=list(name=\"\"))\n",
            "```\n\n", # Bottom
            "![](figures/Tagcloud/NCG_",
            i,
            ".png)\n\n")
    }else{
        ""
    }
}

################ DGN ################
.check_NULL_DGN_2 <- function(out.vecLR, i){
    if(!is.null(out.vecLR[[i]]$Enrich$DGN$Pvalue)){
        paste0("## <font color='#1881c2'>DGN-Enrichment Analysis</font>\n\n",
            "```{r}\n", # Top
            "negLogPval <- -log10(out.vecLR[['",
            i,
            "']]$Enrich$DGN$Pvalue + 1E-10)\n",
            "term <- out.vecLR[['",
            i,
            "']]$Enrich$DGN$Term\n",
            "target <- seq_len(min(100, length(negLogPval)))\n",
            "negLogPval <- negLogPval[target]\n",
            "term <- term[target]\n",
            "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
            "type=\"bar\", color=~negLogPval, text=term,\n",
            "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
            "layout(p, title=\"DGN-Enrichment Analysis\",\n",
            "   xaxis=list(title=\"Term\"),\n",
            "   yaxis=list(title=\"-Log(P-value)\"),\n",
            "   legend=list(name=\"\"))\n",
            "```\n\n", # Bottom
            "![](figures/Tagcloud/DGN_",
            i,
            ".png)\n\n")
    }else{
        ""
    }
}

.XYZ_ENRICH <- function(out.vecLR, i){
    patternName <- paste0("pattern", i)
    go_bp <- .check_NULL_GO_BP(out.vecLR, patternName)
    go_mf <- .check_NULL_GO_MF(out.vecLR, patternName)
    go_cc <- .check_NULL_GO_CC(out.vecLR, patternName)
    mesh_a <- .check_NULL_MeSH_A(out.vecLR, patternName)
    mesh_b <- .check_NULL_MeSH_B(out.vecLR, patternName)
    mesh_c <- .check_NULL_MeSH_C(out.vecLR, patternName)
    mesh_d <- .check_NULL_MeSH_D(out.vecLR, patternName)
    mesh_e <- .check_NULL_MeSH_E(out.vecLR, patternName)
    mesh_f <- .check_NULL_MeSH_F(out.vecLR, patternName)
    mesh_g <- .check_NULL_MeSH_G(out.vecLR, patternName)
    mesh_h <- .check_NULL_MeSH_H(out.vecLR, patternName)
    mesh_i <- .check_NULL_MeSH_I(out.vecLR, patternName)
    mesh_j <- .check_NULL_MeSH_J(out.vecLR, patternName)
    mesh_k <- .check_NULL_MeSH_K(out.vecLR, patternName)
    mesh_l <- .check_NULL_MeSH_L(out.vecLR, patternName)
    mesh_m <- .check_NULL_MeSH_M(out.vecLR, patternName)
    mesh_n <- .check_NULL_MeSH_N(out.vecLR, patternName)
    mesh_v <- .check_NULL_MeSH_V(out.vecLR, patternName)
    mesh_z <- .check_NULL_MeSH_Z(out.vecLR, patternName)
    reactome <- .check_NULL_Reactome(out.vecLR, patternName)
    do <- .check_NULL_DO(out.vecLR, patternName)
    ncg <- .check_NULL_NCG(out.vecLR, patternName)
    dgn <- .check_NULL_DGN(out.vecLR, patternName)
    # Output
    paste(c(go_bp, go_mf, go_cc,
    mesh_a, mesh_b, mesh_c, mesh_d, mesh_e, mesh_f,
    mesh_g, mesh_h, mesh_i, mesh_j, mesh_k, mesh_l,
    mesh_m, mesh_n, mesh_v, mesh_z,
    reactome,
    do, ncg, dgn), collapse="")
}


################ GO_BP ################
.check_NULL_GO_BP <- function(out.vecLR, patternName){
    pval <- eval(parse(text=paste0(
        "out.vecLR[\"Enrich\", \"",
        patternName, "\"][[1]]$GO_BP$Pvalue")))
    if(!is.null(pval)){
        paste0("## <font color='#1881c2'>GO-Enrichment Analysis (BP : Biological Process)</font>\n\n",
        "```{r}\n", # Top
        "negLogPval <- -log10(out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$GO_BP$Pvalue + 1E-10)\n",
        "term <- out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$GO_BP$Term\n",
        "target <- seq_len(min(100, length(negLogPval)))\n",
        "negLogPval <- negLogPval[target]\n",
        "term <- term[target]\n",
        "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
        "type=\"bar\", color=~negLogPval, text=term,\n",
        "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
        "layout(p, title=\"GO-Enrichment Analysis (BP)\",\n",
        "   xaxis=list(title=\"Term\"),\n",
        "   yaxis=list(title=\"-Log(P-value)\"),\n",
        "   legend=list(name=\"\"))\n",
        "```\n\n", # Bottom
        "![](figures/Tagcloud/GO_BP_",
        patternName,
        ".png)\n\n")
    }else{
        ""
    }
}

################ GO_MF ################
.check_NULL_GO_MF <- function(out.vecLR, patternName){
    pval <- eval(parse(text=paste0(
        "out.vecLR[\"Enrich\", \"",
        patternName, "\"][[1]]$GO_MF$Pvalue")))
    if(!is.null(pval)){
        paste0("## <font color='#1881c2'>GO-Enrichment Analysis (MF : Molecular Function)</font>\n\n",
        "```{r}\n", # Top
        "negLogPval <- -log10(out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$GO_MF$Pvalue + 1E-10)\n",
        "term <- out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$GO_MF$Term\n",
        "target <- seq_len(min(100, length(negLogPval)))\n",
        "negLogPval <- negLogPval[target]\n",
        "term <- term[target]\n",
        "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
        "type=\"bar\", color=~negLogPval, text=term,\n",
        "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
        "layout(p, title=\"GO-Enrichment Analysis (MF)\",\n",
        "   xaxis=list(title=\"Term\"),\n",
        "   yaxis=list(title=\"-Log(P-value)\"),\n",
        "   legend=list(name=\"\"))\n",
        "```\n\n", # Bottom
        "![](figures/Tagcloud/GO_MF_",
        patternName,
        ".png)\n\n")
    }else{
        ""
    }
}

################ GO_CC ################
.check_NULL_GO_CC <- function(out.vecLR, patternName){
    pval <- eval(parse(text=paste0(
        "out.vecLR[\"Enrich\", \"",
        patternName, "\"][[1]]$GO_CC$Pvalue")))
    if(!is.null(pval)){
        paste0("## <font color='#1881c2'>GO-Enrichment Analysis (CC : Cellular Component)</font>\n\n",
        "```{r}\n", # Top
        "negLogPval <- -log10(out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$GO_CC$Pvalue + 1E-10)\n",
        "term <- out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$GO_CC$Term\n",
        "target <- seq_len(min(100, length(negLogPval)))\n",
        "negLogPval <- negLogPval[target]\n",
        "term <- term[target]\n",
        "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
        "type=\"bar\", color=~negLogPval, text=term,\n",
        "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
        "layout(p, title=\"GO-Enrichment Analysis (CC)\",\n",
        "   xaxis=list(title=\"Term\"),\n",
        "   yaxis=list(title=\"-Log(P-value)\"),\n",
        "   legend=list(name=\"\"))\n",
        "```\n\n", # Bottom
        "![](figures/Tagcloud/GO_CC_",
        patternName,
        ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_A ################
.check_NULL_MeSH_A <- function(out.vecLR, patternName){
    pval <- eval(parse(text=paste0(
        "out.vecLR[\"Enrich\", \"",
        patternName, "\"][[1]]$MeSH_A$Pvalue")))
    if(!is.null(pval)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (A : Anatomy)</font>\n\n",
        "```{r}\n", # Top
        "negLogPval <- -log10(out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_A$Pvalue + 1E-10)\n",
        "term <- out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_A$Term\n",
        "target <- seq_len(min(100, length(negLogPval)))\n",
        "negLogPval <- negLogPval[target]\n",
        "term <- term[target]\n",
        "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
        "type=\"bar\", color=~negLogPval, text=term,\n",
        "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
        "layout(p, title=\"MeSH-Enrichment Analysis (A)\",\n",
        "   xaxis=list(title=\"Term\"),\n",
        "   yaxis=list(title=\"-Log(P-value)\"),\n",
        "   legend=list(name=\"\"))\n",
        "```\n\n", # Bottom
        "![](figures/Tagcloud/MeSH_A_",
        patternName,
        ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_B ################
.check_NULL_MeSH_B <- function(out.vecLR, patternName){
    pval <- eval(parse(text=paste0(
        "out.vecLR[\"Enrich\", \"",
        patternName, "\"][[1]]$MeSH_B$Pvalue")))
    if(!is.null(pval)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (B : Organisms)</font>\n\n",
        "```{r}\n", # Top
        "negLogPval <- -log10(out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_B$Pvalue + 1E-10)\n",
        "term <- out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_B$Term\n",
        "target <- seq_len(min(100, length(negLogPval)))\n",
        "negLogPval <- negLogPval[target]\n",
        "term <- term[target]\n",
        "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
        "type=\"bar\", color=~negLogPval, text=term,\n",
        "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
        "layout(p, title=\"MeSH-Enrichment Analysis (B)\",\n",
        "   xaxis=list(title=\"Term\"),\n",
        "   yaxis=list(title=\"-Log(P-value)\"),\n",
        "   legend=list(name=\"\"))\n",
        "```\n\n", # Bottom
        "![](figures/Tagcloud/MeSH_B_",
        patternName,
        ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_C ################
.check_NULL_MeSH_C <- function(out.vecLR, patternName){
    pval <- eval(parse(text=paste0(
        "out.vecLR[\"Enrich\", \"",
        patternName, "\"][[1]]$MeSH_C$Pvalue")))
    if(!is.null(pval)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (C : Diseases)</font>\n\n",
        "```{r}\n", # Top
        "negLogPval <- -log10(out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_C$Pvalue + 1E-10)\n",
        "term <- out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_C$Term\n",
        "target <- seq_len(min(100, length(negLogPval)))\n",
        "negLogPval <- negLogPval[target]\n",
        "term <- term[target]\n",
        "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
        "type=\"bar\", color=~negLogPval, text=term,\n",
        "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
        "layout(p, title=\"MeSH-Enrichment Analysis (C)\",\n",
        "   xaxis=list(title=\"Term\"),\n",
        "   yaxis=list(title=\"-Log(P-value)\"),\n",
        "   legend=list(name=\"\"))\n",
        "```\n\n", # Bottom
        "![](figures/Tagcloud/MeSH_C_",
        patternName,
        ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_D ################
.check_NULL_MeSH_D <- function(out.vecLR, patternName){
    pval <- eval(parse(text=paste0(
        "out.vecLR[\"Enrich\", \"",
        patternName, "\"][[1]]$MeSH_D$Pvalue")))
    if(!is.null(pval)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (D : Drugs)</font>\n\n",
        "```{r}\n", # Top
        "negLogPval <- -log10(out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_D$Pvalue + 1E-10)\n",
        "term <- out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_D$Term\n",
        "target <- seq_len(min(100, length(negLogPval)))\n",
        "negLogPval <- negLogPval[target]\n",
        "term <- term[target]\n",
        "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
        "type=\"bar\", color=~negLogPval, text=term,\n",
        "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
        "layout(p, title=\"MeSH-Enrichment Analysis (D)\",\n",
        "   xaxis=list(title=\"Term\"),\n",
        "   yaxis=list(title=\"-Log(P-value)\"),\n",
        "   legend=list(name=\"\"))\n",
        "```\n\n", # Bottom
        "![](figures/Tagcloud/MeSH_D_",
        patternName,
        ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_E ################
.check_NULL_MeSH_E <- function(out.vecLR, patternName){
    pval <- eval(parse(text=paste0(
        "out.vecLR[\"Enrich\", \"",
        patternName, "\"][[1]]$MeSH_E$Pvalue")))
    if(!is.null(pval)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (E : Analytical, Diagnostic and Therapeutic Techniques and Equipment)</font>\n\n",
        "```{r}\n", # Top
        "negLogPval <- -log10(out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_E$Pvalue + 1E-10)\n",
        "term <- out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_E$Term\n",
        "target <- seq_len(min(100, length(negLogPval)))\n",
        "negLogPval <- negLogPval[target]\n",
        "term <- term[target]\n",
        "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
        "type=\"bar\", color=~negLogPval, text=term,\n",
        "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
        "layout(p, title=\"MeSH-Enrichment Analysis (E)\",\n",
        "   xaxis=list(title=\"Term\"),\n",
        "   yaxis=list(title=\"-Log(P-value)\"),\n",
        "   legend=list(name=\"\"))\n",
        "```\n\n", # Bottom
        "![](figures/Tagcloud/MeSH_E_",
        patternName,
        ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_F ################
.check_NULL_MeSH_F <- function(out.vecLR, patternName){
    pval <- eval(parse(text=paste0(
        "out.vecLR[\"Enrich\", \"",
        patternName, "\"][[1]]$MeSH_F$Pvalue")))
    if(!is.null(pval)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (F : Psychiatry and Psychology)</font>\n\n",
        "```{r}\n", # Top
        "negLogPval <- -log10(out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_F$Pvalue + 1E-10)\n",
        "term <- out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_F$Term\n",
        "target <- seq_len(min(100, length(negLogPval)))\n",
        "negLogPval <- negLogPval[target]\n",
        "term <- term[target]\n",
        "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
        "type=\"bar\", color=~negLogPval, text=term,\n",
        "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
        "layout(p, title=\"MeSH-Enrichment Analysis (F)\",\n",
        "   xaxis=list(title=\"Term\"),\n",
        "   yaxis=list(title=\"-Log(P-value)\"),\n",
        "   legend=list(name=\"\"))\n",
        "```\n\n", # Bottom
        "![](figures/Tagcloud/MeSH_F_",
        patternName,
        ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_G ################
.check_NULL_MeSH_G <- function(out.vecLR, patternName){
    pval <- eval(parse(text=paste0(
        "out.vecLR[\"Enrich\", \"",
        patternName, "\"][[1]]$MeSH_G$Pvalue")))
    if(!is.null(pval)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (G : Phenomena and Processes)</font>\n\n",
        "```{r}\n", # Top
        "negLogPval <- -log10(out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_G$Pvalue + 1E-10)\n",
        "term <- out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_G$Term\n",
        "target <- seq_len(min(100, length(negLogPval)))\n",
        "negLogPval <- negLogPval[target]\n",
        "term <- term[target]\n",
        "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
        "type=\"bar\", color=~negLogPval, text=term,\n",
        "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
        "layout(p, title=\"MeSH-Enrichment Analysis (G)\",\n",
        "   xaxis=list(title=\"Term\"),\n",
        "   yaxis=list(title=\"-Log(P-value)\"),\n",
        "   legend=list(name=\"\"))\n",
        "```\n\n", # Bottom
        "![](figures/Tagcloud/MeSH_G_",
        patternName,
        ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_H ################
.check_NULL_MeSH_H <- function(out.vecLR, patternName){
    pval <- eval(parse(text=paste0(
        "out.vecLR[\"Enrich\", \"",
        patternName, "\"][[1]]$MeSH_H$Pvalue")))
    if(!is.null(pval)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (H : Disciplines and Occupations)</font>\n\n",
        "```{r}\n", # Top
        "negLogPval <- -log10(out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_H$Pvalue + 1E-10)\n",
        "term <- out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_H$Term\n",
        "target <- seq_len(min(100, length(negLogPval)))\n",
        "negLogPval <- negLogPval[target]\n",
        "term <- term[target]\n",
        "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
        "type=\"bar\", color=~negLogPval, text=term,\n",
        "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
        "layout(p, title=\"MeSH-Enrichment Analysis (H)\",\n",
        "   xaxis=list(title=\"Term\"),\n",
        "   yaxis=list(title=\"-Log(P-value)\"),\n",
        "   legend=list(name=\"\"))\n",
        "```\n\n", # Bottom
        "![](figures/Tagcloud/MeSH_H_",
        patternName,
        ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_I ################
.check_NULL_MeSH_I <- function(out.vecLR, patternName){
    pval <- eval(parse(text=paste0(
        "out.vecLR[\"Enrich\", \"",
        patternName, "\"][[1]]$MeSH_I$Pvalue")))
    if(!is.null(pval)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (I : Anthropology, Education, Sociology and Social Phenomena)</font>\n\n",
        "```{r}\n", # Top
        "negLogPval <- -log10(out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_I$Pvalue + 1E-10)\n",
        "term <- out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_I$Term\n",
        "target <- seq_len(min(100, length(negLogPval)))\n",
        "negLogPval <- negLogPval[target]\n",
        "term <- term[target]\n",
        "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
        "type=\"bar\", color=~negLogPval, text=term,\n",
        "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
        "layout(p, title=\"MeSH-Enrichment Analysis (I)\",\n",
        "   xaxis=list(title=\"Term\"),\n",
        "   yaxis=list(title=\"-Log(P-value)\"),\n",
        "   legend=list(name=\"\"))\n",
        "```\n\n", # Bottom
        "![](figures/Tagcloud/MeSH_I_",
        patternName,
        ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_J ################
.check_NULL_MeSH_J <- function(out.vecLR, patternName){
    pval <- eval(parse(text=paste0(
        "out.vecLR[\"Enrich\", \"",
        patternName, "\"][[1]]$MeSH_J$Pvalue")))
    if(!is.null(pval)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (J : Technology and Food and Beverages)</font>\n\n",
        "```{r}\n", # Top
        "negLogPval <- -log10(out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_J$Pvalue + 1E-10)\n",
        "term <- out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_J$Term\n",
        "target <- seq_len(min(100, length(negLogPval)))\n",
        "negLogPval <- negLogPval[target]\n",
        "term <- term[target]\n",
        "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
        "type=\"bar\", color=~negLogPval, text=term,\n",
        "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
        "layout(p, title=\"MeSH-Enrichment Analysis (J)\",\n",
        "   xaxis=list(title=\"Term\"),\n",
        "   yaxis=list(title=\"-Log(P-value)\"),\n",
        "   legend=list(name=\"\"))\n",
        "```\n\n", # Bottom
        "![](figures/Tagcloud/MeSH_J_",
        patternName,
        ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_K ################
.check_NULL_MeSH_K <- function(out.vecLR, patternName){
    pval <- eval(parse(text=paste0(
        "out.vecLR[\"Enrich\", \"",
        patternName, "\"][[1]]$MeSH_K$Pvalue")))
    if(!is.null(pval)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (K : Humanities)</font>\n\n",
        "```{r}\n", # Top
        "negLogPval <- -log10(out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_K$Pvalue + 1E-10)\n",
        "term <- out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_K$Term\n",
        "target <- seq_len(min(100, length(negLogPval)))\n",
        "negLogPval <- negLogPval[target]\n",
        "term <- term[target]\n",
        "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
        "type=\"bar\", color=~negLogPval, text=term,\n",
        "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
        "layout(p, title=\"MeSH-Enrichment Analysis (K)\",\n",
        "   xaxis=list(title=\"Term\"),\n",
        "   yaxis=list(title=\"-Log(P-value)\"),\n",
        "   legend=list(name=\"\"))\n",
        "```\n\n", # Bottom
        "![](figures/Tagcloud/MeSH_K_",
        patternName,
        ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_L ################
.check_NULL_MeSH_L <- function(out.vecLR, patternName){
    pval <- eval(parse(text=paste0(
        "out.vecLR[\"Enrich\", \"",
        patternName, "\"][[1]]$MeSH_L$Pvalue")))
    if(!is.null(pval)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (L : Information Science)</font>\n\n",
        "```{r}\n", # Top
        "negLogPval <- -log10(out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_L$Pvalue + 1E-10)\n",
        "term <- out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_L$Term\n",
        "target <- seq_len(min(100, length(negLogPval)))\n",
        "negLogPval <- negLogPval[target]\n",
        "term <- term[target]\n",
        "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
        "type=\"bar\", color=~negLogPval, text=term,\n",
        "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
        "layout(p, title=\"MeSH-Enrichment Analysis (L)\",\n",
        "   xaxis=list(title=\"Term\"),\n",
        "   yaxis=list(title=\"-Log(P-value)\"),\n",
        "   legend=list(name=\"\"))\n",
        "```\n\n", # Bottom
        "![](figures/Tagcloud/MeSH_L_",
        patternName,
        ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_M ################
.check_NULL_MeSH_M <- function(out.vecLR, patternName){
    pval <- eval(parse(text=paste0(
        "out.vecLR[\"Enrich\", \"",
        patternName, "\"][[1]]$MeSH_M$Pvalue")))
    if(!is.null(pval)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (M : Persons)</font>\n\n",
        "```{r}\n", # Top
        "negLogPval <- -log10(out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_M$Pvalue + 1E-10)\n",
        "term <- out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_M$Term\n",
        "target <- seq_len(min(100, length(negLogPval)))\n",
        "negLogPval <- negLogPval[target]\n",
        "term <- term[target]\n",
        "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
        "type=\"bar\", color=~negLogPval, text=term,\n",
        "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
        "layout(p, title=\"MeSH-Enrichment Analysis (M)\",\n",
        "   xaxis=list(title=\"Term\"),\n",
        "   yaxis=list(title=\"-Log(P-value)\"),\n",
        "   legend=list(name=\"\"))\n",
        "```\n\n", # Bottom
        "![](figures/Tagcloud/MeSH_M_",
        patternName,
        ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_N ################
.check_NULL_MeSH_N <- function(out.vecLR, patternName){
    pval <- eval(parse(text=paste0(
        "out.vecLR[\"Enrich\", \"",
        patternName, "\"][[1]]$MeSH_N$Pvalue")))
    if(!is.null(pval)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (N : Health Care)</font>\n\n",
        "```{r}\n", # Top
        "negLogPval <- -log10(out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_N$Pvalue + 1E-10)\n",
        "term <- out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_N$Term\n",
        "target <- seq_len(min(100, length(negLogPval)))\n",
        "negLogPval <- negLogPval[target]\n",
        "term <- term[target]\n",
        "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
        "type=\"bar\", color=~negLogPval, text=term,\n",
        "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
        "layout(p, title=\"MeSH-Enrichment Analysis (N)\",\n",
        "   xaxis=list(title=\"Term\"),\n",
        "   yaxis=list(title=\"-Log(P-value)\"),\n",
        "   legend=list(name=\"\"))\n",
        "```\n\n", # Bottom
        "![](figures/Tagcloud/MeSH_N_",
        patternName,
        ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_V ################
.check_NULL_MeSH_V <- function(out.vecLR, patternName){
    pval <- eval(parse(text=paste0(
        "out.vecLR[\"Enrich\", \"",
        patternName, "\"][[1]]$MeSH_V$Pvalue")))
    if(!is.null(pval)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (V : Publication Type)</font>\n\n",
        "```{r}\n", # Top
        "negLogPval <- -log10(out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_V$Pvalue + 1E-10)\n",
        "term <- out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_V$Term\n",
        "target <- seq_len(min(100, length(negLogPval)))\n",
        "negLogPval <- negLogPval[target]\n",
        "term <- term[target]\n",
        "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
        "type=\"bar\", color=~negLogPval, text=term,\n",
        "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
        "layout(p, title=\"MeSH-Enrichment Analysis (V)\",\n",
        "   xaxis=list(title=\"Term\"),\n",
        "   yaxis=list(title=\"-Log(P-value)\"),\n",
        "   legend=list(name=\"\"))\n",
        "```\n\n", # Bottom
        "![](figures/Tagcloud/MeSH_V_",
        patternName,
        ".png)\n\n")
    }else{
        ""
    }
}

################ MeSH_Z ################
.check_NULL_MeSH_Z <- function(out.vecLR, patternName){
    pval <- eval(parse(text=paste0(
        "out.vecLR[\"Enrich\", \"",
        patternName, "\"][[1]]$MeSH_Z$Pvalue")))
    if(!is.null(pval)){
        paste0("## <font color='#1881c2'>MeSH-Enrichment Analysis (Z : Geographical Locations)</font>\n\n",
        "```{r}\n", # Top
        "negLogPval <- -log10(out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_Z$Pvalue + 1E-10)\n",
        "term <- out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$MeSH_Z$Term\n",
        "target <- seq_len(min(100, length(negLogPval)))\n",
        "negLogPval <- negLogPval[target]\n",
        "term <- term[target]\n",
        "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
        "type=\"bar\", color=~negLogPval, text=term,\n",
        "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
        "layout(p, title=\"MeSH-Enrichment Analysis (Z)\",\n",
        "   xaxis=list(title=\"Term\"),\n",
        "   yaxis=list(title=\"-Log(P-value)\"),\n",
        "   legend=list(name=\"\"))\n",
        "```\n\n", # Bottom
        "![](figures/Tagcloud/MeSH_Z_",
        patternName,
        ".png)\n\n")
    }else{
        ""
    }
}

################ Reactome ################
.check_NULL_Reactome <- function(out.vecLR, patternName){
    pval <- eval(parse(text=paste0(
        "out.vecLR[\"Enrich\", \"",
        patternName, "\"][[1]]$Reactome$Pvalue")))
    if(!is.null(pval)){
        paste0("## <font color='#1881c2'>Reactome Pathway-Enrichment Analysis</font>\n\n",
        "```{r}\n", # Top
        "negLogPval <- -log10(out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$Reactome$Pvalue + 1E-10)\n",
        "term <- out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$Reactome$Term\n",
        "target <- seq_len(min(100, length(negLogPval)))\n",
        "negLogPval <- negLogPval[target]\n",
        "term <- term[target]\n",
        "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
        "type=\"bar\", color=~negLogPval, text=term,\n",
        "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
        "layout(p, title=\"Reactome-Enrichment Analysis\",\n",
        "   xaxis=list(title=\"Term\"),\n",
        "   yaxis=list(title=\"-Log(P-value)\"),\n",
        "   legend=list(name=\"\"))\n",
        "```\n\n", # Bottom
        "![](figures/Tagcloud/Reactome_",
        patternName,
        ".png)\n\n")
    }else{
        ""
    }
}

################ DO ################
.check_NULL_DO <- function(out.vecLR, patternName){
    pval <- eval(parse(text=paste0(
        "out.vecLR[\"Enrich\", \"",
        patternName, "\"][[1]]$DO$Pvalue")))
    if(!is.null(pval)){
        paste0("## <font color='#1881c2'>DO-Enrichment Analysis</font>\n\n",
        "```{r}\n", # Top
        "negLogPval <- -log10(out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$DO$Pvalue + 1E-10)\n",
        "term <- out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$DO$Term\n",
        "target <- seq_len(min(100, length(negLogPval)))\n",
        "negLogPval <- negLogPval[target]\n",
        "term <- term[target]\n",
        "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
        "type=\"bar\", color=~negLogPval, text=term,\n",
        "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
        "layout(p, title=\"DO-Enrichment Analysis\",\n",
        "   xaxis=list(title=\"Term\"),\n",
        "   yaxis=list(title=\"-Log(P-value)\"),\n",
        "   legend=list(name=\"\"))\n",
        "```\n\n", # Bottom
        "![](figures/Tagcloud/DO_",
        patternName,
        ".png)\n\n")
    }else{
        ""
    }
}

################ NCG ################
.check_NULL_NCG <- function(out.vecLR, patternName){
    pval <- eval(parse(text=paste0(
        "out.vecLR[\"Enrich\", \"",
        patternName, "\"][[1]]$NCG$Pvalue")))
    if(!is.null(pval)){
        paste0("## <font color='#1881c2'>NCG-Enrichment Analysis</font>\n\n",
        "```{r}\n", # Top
        "negLogPval <- -log10(out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$NCG$Pvalue + 1E-10)\n",
        "term <- out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$NCG$Term\n",
        "target <- seq_len(min(100, length(negLogPval)))\n",
        "negLogPval <- negLogPval[target]\n",
        "term <- term[target]\n",
        "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
        "type=\"bar\", color=~negLogPval, text=term,\n",
        "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
        "layout(p, title=\"NCG-Enrichment Analysis\",\n",
        "   xaxis=list(title=\"Term\"),\n",
        "   yaxis=list(title=\"-Log(P-value)\"),\n",
        "   legend=list(name=\"\"))\n",
        "```\n\n", # Bottom
        "![](figures/Tagcloud/NCG_",
        patternName,
        ".png)\n\n")
    }else{
        ""
    }
}

################ DGN ################
.check_NULL_DGN <- function(out.vecLR, patternName){
    pval <- eval(parse(text=paste0(
        "out.vecLR[\"Enrich\", \"",
        patternName, "\"][[1]]$DGN$Pvalue")))
    if(!is.null(pval)){
        paste0("## <font color='#1881c2'>DGN-Enrichment Analysis</font>\n\n",
        "```{r}\n", # Top
        "negLogPval <- -log10(out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$DGN$Pvalue + 1E-10)\n",
        "term <- out.vecLR[\"Enrich\", \"",
        patternName,
        "\"][[1]]$DGN$Term\n",
        "target <- seq_len(min(100, length(negLogPval)))\n",
        "negLogPval <- negLogPval[target]\n",
        "term <- term[target]\n",
        "p <- plot_ly(x=seq_along(negLogPval), y=~negLogPval,\n",
        "type=\"bar\", color=~negLogPval, text=term,\n",
        "colors=c(\"#4b61ba\", \"gray\", \"#a87963\", \"red\"))\n",
        "layout(p, title=\"DGN-Enrichment Analysis\",\n",
        "   xaxis=list(title=\"Term\"),\n",
        "   yaxis=list(title=\"-Log(P-value)\"),\n",
        "   legend=list(name=\"\"))\n",
        "```\n\n", # Bottom
        "![](figures/Tagcloud/DGN_",
        patternName,
        ".png)\n\n")
    }else{
        ""
    }
}

.XYZ_HEADER3 <- function(i){
    paste0("# <font color='#1881c2'>(\\*,\\*,", i,
        ") Pattern-related Enrichment Analysis</font>\n\n",
        "```{r}\n", # Top
        "load(\"reanalysis.RData\")\n",
        "library(\"plotly\")\n",
        "```\n\n" # Bottom
        )
}

.XYZ_HEADER3_2 <- function(index, i){
        paste0("# <font color='#1881c2'>(",
        paste(c(index[i, seq_len(2)], ""), collapse=","),
        ") -related Enrichment Analysis</font>\n\n",
        "```{r}\n", # Top
        "load(\"reanalysis.RData\")\n",
        "library(\"plotly\")\n",
        "```\n\n" # Bottom
        )
}
