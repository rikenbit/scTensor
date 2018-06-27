#
# cellCellSetting
#
setGeneric("cellCellSetting", function(sce, lrbase, color, label){
    standardGeneric("cellCellSetting")})
setMethod("cellCellSetting", signature(sce="SingleCellExperiment"),
    function(sce, lrbase, color, label){
        userobjects <- deparse(substitute(sce))
        .cellCellSetting(userobjects, lrbase, color, label)})
.cellCellSetting <- function(userobjects, lrbase, color, label){
    sce <- eval(parse(text=userobjects))
    # class-check
    if(is(lrbase)[1] != "LRBaseDb"){
        stop("Please specify the lrbase as LRBaseDbi-class object
        such as LRBase.Hsa.eg.db")
    }
    # size-check
    if(length(color) != ncol(assay(sce))){
        stop("Please specify the color as the vector which has
            same length of nrow(assay(sce))")
    }
    # size-check
    if(length(label) != ncol(assay(sce))){
        stop("Please specify the label as the vector which has
            same length of nrow(assay(sce))")
    }
    # size-check
    if(length(unique(color)) != length(unique(label))){
        stop("Please specify the kind of elements containing
            in color and label as same")
    }
    # Overwrite
    metadata(sce) <- list(lrbase=lrbase, color=color, label=label)
    assign(userobjects, sce, envir=.GlobalEnv)
}

#
# cellCellDecomp
#
setGeneric("cellCellDecomp", function(sce, algorithm="ntd", ranks=c(3,3,3),
    rank=3, thr1=log2(5), thr2=25, centering=TRUE, mergeas="sum",
    outer="+", comb="random", num.sampling=100, decomp=TRUE){
    standardGeneric("cellCellDecomp")})
setMethod("cellCellDecomp", signature(sce="SingleCellExperiment"),
    function(sce, algorithm, ranks, rank, thr1, thr2, centering,
        mergeas, outer, comb, num.sampling, decomp){
        userobjects <- deparse(substitute(sce))
        .cellCellDecomp(userobjects, algorithm, ranks, rank, thr1, thr2,
            centering, mergeas, outer, comb, num.sampling, decomp)})

.cellCellDecomp <- function(userobjects, algorithm, ranks, rank, thr1, thr2,
    centering, mergeas, outer, comb, num.sampling, decomp){
    # Import from sce object
    sce <- eval(parse(text=userobjects))
    input <- assay(sce)
    LR <- LRBaseDbi::select(metadata(sce)$lrbase,
        columns=c("GENEID_L", "GENEID_R"),
        keytype="GENEID_L",
        keys=LRBaseDbi::keys(metadata(sce)$lrbase, keytype="GENEID_L"))
    celltypes <- metadata(sce)$color
    names(celltypes) <- metadata(sce)$label

    # Gene Symbol-check
    genesymbols <- grep("^[A-Z]", rownames(input))
    if(length(genesymbols) != 0){
        message(paste0("Input data matrix may contains ",
            length(genesymbols),
            " gene symbols because the name contains some alphabets.\n",
            "scTensor uses only NCBI Gene IDs for now.\n",
            "Here, the gene symbols are removed and remaining ",
            nrow(input) - length(genesymbols),
            " NCBI Gene IDs are used for scTensor next step."
            ))
        input <- input[setdiff(seq_len(nrow(input)), genesymbols), ]
    }

    # class-check
    if(is(sce)[1] != "SingleCellExperiment"){
        stop("input data must be defined as SingleCellExperiment-class.
            Please use SingleCellExperiment package.")
    }
    # algorithm-check
    if(!algorithm %in% c("ntd", "nmf", "pearson", "spearman", "distance",
        "pearson.lr", "spearman.lr", "distance.lr", "pcomb")){
        stop(
            paste0("algorithm must be defined as 'ntd', 'nmf', 'pearson',",
            " 'spearman', 'distance', 'pearson.lr', 'spearman.lr', ",
            "'distance.lr', or 'pcomb'"))
    }
    # ranks-check
    max.rank <- length(unique(celltypes))
    if(algorithm == "ntd" && (ranks[1] > max.rank || ranks[2] > max.rank)){
        stop("ranks must be defined less than number of celltypes")
    }
    if(algorithm == "nmf" && rank > max.rank){
        stop("ranks must be defined less than number of celltypes")
    }
    # thr-check
    if(algorithm == "pcomb" && (thr1 <= 0 || thr2 <= 0)){
        warning("None of cell-cell interaction will be detected.")
    }
    # mergeas-check
    if(!mergeas %in% c("sum", "mean")){
        stop("mergeas must be defined as 'sum', or 'mean'")
    }
    # outer-check
    if(!outer %in% c("+", "*")){
        stop("outer must be defined as '+', or '*'")
    }
    # comb-check
    if(!comb %in% c("random", "all")){
        stop("comb must be defined as 'random', or 'all'")
    }
    # thr-check
    if(comb == "random" && (num.sampling <= 0)){
        warning("None of cell-cell interaction will be detected.")
    }
    # 3-Order = NTD
    if(algorithm == "ntd"){
        res.sctensor <- .cellCellDecomp.Third(input, LR, celltypes, ranks,
            centering, mergeas, outer, comb, num.sampling, decomp)
    }
    # 2-Order = NMF
    else if(algorithm == "nmf"){
        res.sctensor <- .cellCellDecomp.Second(input, LR, celltypes, rank,
            centering, mergeas, outer, comb, num.sampling, decomp)
    }
    # Other methods
    else if(algorithm == "pearson"){
        res.sctensor <- .cellCellDecomp.Pearson(input, celltypes)
    }
    else if(algorithm == "spearman"){
        res.sctensor <- .cellCellDecomp.Spearman(input, celltypes)
    }
    else if(algorithm == "distance"){
        res.sctensor <- .cellCellDecomp.Distance(input, celltypes)
    }
    else if(algorithm == "pearson.lr"){
        res.sctensor <- .cellCellDecomp.Pearson.LR(input, LR, celltypes)
    }
    else if(algorithm == "spearman.lr"){
        res.sctensor <- .cellCellDecomp.Spearman.LR(input, LR, celltypes)
    }
    else if(algorithm == "distance.lr"){
        res.sctensor <- .cellCellDecomp.Distance.LR(input, LR, celltypes)
    }
    else if(algorithm == "pcomb"){
        res.sctensor <- .cellCellDecomp.PossibleCombination(input, LR,
            celltypes, thr1=thr1, thr2=thr2)
    }else{
        stop("Please specify the appropriate algorithm\n")
    }
    # Overwrite
    metadata(sce) <- list(lrbase=metadata(sce)$lrbase,
        color=metadata(sce)$color, label=metadata(sce)$label,
        algorithm=algorithm, sctensor=res.sctensor)
    assign(userobjects, sce, envir=.GlobalEnv)
    # Output
}

#
# cellCellReport
#
setGeneric("cellCellReport", function(sce, reducedDimNames,
    out.dir=NULL, html.open=FALSE,
    title="The result of scTensor",
    author="The person who runs this script", thr=1){
    standardGeneric("cellCellReport")})
setMethod("cellCellReport", signature(sce="SingleCellExperiment"),
    function(sce, reducedDimNames, out.dir, html.open, title, author, thr){
        userobjects <- deparse(substitute(sce))
        .cellCellReport(userobjects, reducedDimNames, out.dir,
            html.open, title, author, thr)})
.cellCellReport <- function(userobjects, reducedDimNames,
    out.dir=NULL, html.open=FALSE,
    title="The result of scTensor",
    author="The person who runs this script", thr=1){
    # Import from sce object
    sce <- eval(parse(text=userobjects))
    # class-check
    if(is(sce)[1] != "SingleCellExperiment"){
        stop(paste0("input data must defined as SingleCellExperiment-class.",
            " Please use SingleCellExperiment package."))
    }
    # algorithm-check
    if(metadata(sce)$algorithm != "ntd"){
        stop(paste0("cellCellReport can be performed by the result of",
            " cellCellDecomp in which the algorithm is ",
            "specified as 'ntd' for now."))
    }
    # out.dir-check
    if(is.null(out.dir)){
        stop(paste0("Please specify the output directory for ",
            "saving your analysis result"))
    }else{
        if(!file.exists(out.dir)){
            dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
        }
    }

    # Data matrix
    input <- assay(sce)
    # Low dimensional data
    twoD <- eval(parse(text=paste0("reducedDims(sce)$", reducedDimNames)))
    # Ligand-Receptor, PMID
    LR <- LRBaseDbi::select(metadata(sce)$lrbase,
        columns=c("GENEID_L", "GENEID_R", "SOURCEID"),
        keytype="GENEID_L",
        keys=LRBaseDbi::keys(metadata(sce)$lrbase, keytype="GENEID_L"))
    # Species
    lrname <- LRBaseDbi::packageName(metadata(sce)$lrbase)
    spc <- substr(lrname, nchar(lrname) - 8, nchar(lrname))
    spc <- gsub(".eg.db", "", spc)

    # biomaRt Setting
    ens <- .ensembl(spc)
    # GeneName, Description, GO, Reactome, MeSH
    GeneInfo <- .geneinformation(sce, ens, spc, LR)

    # Cell Label
    celltypes <- metadata(sce)$color
    names(celltypes) <- metadata(sce)$label
    # Core Tensor
    index <- metadata(sce)$sctensor$index
    corevalue <- index[, "Value"]
    corevalue <- corevalue / sum(corevalue) * 100
    # Thresholding of the elements of core tensor
    selected <- which(corevalue > thr)
    if(length(selected) == 0){
        message(paste0("None of core tensor element is selected.\n",
        "Please specify the larger thr or perform cellCellDecomp\n",
        "with smaller ranks such as c(3,3,3)."))
    }
    names(corevalue) <- c(rep("selected", length=length(selected)),
        rep("not selected",
        length=length(corevalue) - length(selected)))
    # Save the result of scTensor
    save(sce, input, twoD, LR, celltypes, index, corevalue, selected,
        file=paste0(out.dir, "/reanalysis.RData"))

    # Tempolary Directory for saving the analytical result
    temp <- tempdir()
    dir.create(paste0(temp, "/figures"),
        showWarnings = FALSE, recursive = TRUE)
    dir.create(paste0(out.dir, "/figures"),
        showWarnings = FALSE, recursive = TRUE)

    # Table
    if(length(selected) != 0){
        rmdfiles <- paste0("pattern",
            apply(index[selected, seq_len(3)], 1,
                function(x){paste(x, collapse="_")}), ".Rmd")
        for(i in seq_along(rmdfiles)){
            cat(paste0(rmdfiles[i],
                " is created...(", i, " / ", length(rmdfiles), ")\n"))
            # Save Rmd
            ranking <- ncol(metadata(sce)$sctensor$ppi) -
                base::rank(metadata(sce)$sctensor$ppi[index[i,3], ])
            target <- sapply(seq_len(100), function(x){
                which(ranking == x)[1]
                })
            Value <- metadata(sce)$sctensor$ppi[index[i,3], target]
            Percentage <- Value /
                sum(metadata(sce)$sctensor$ppi[index[i,3], ]) * 100
            # Header Row
            XYZ <- paste0("# Details of (",
                paste(index[i, seq_len(3)], collapse=","),
                ") Pattern (Top100 LR-pairs)\n\n",
                "## (", paste(index[i, seq_len(2)], collapse=","),
                ",\\*) Pattern\n\n",
                "![](figures/LRNetwork_",
                paste(index[i, seq_len(2)], collapse="_"), ".png)\n\n",
                "## (\\*,\\*,", index[i, 3], ") Pattern\n\n",
                "|Rank|Ligand Gene|Receptor Gene|",
                "Ligand Expression in each cell (Log10(Exp+1))|",
                "Receptor Expression in each cell (Log10(Exp+1))|",
                "LR-pair factor value (Percentage)|",
                "PubMed|\n",
                "|----|----|----|----|----|----|----|\n")

            for(j in seq_along(target)){
                Ranking <- j
                L_R <- strsplit(colnames(metadata(sce)$sctensor$ppi)[j], "_")
                LigandGeneID <- L_R[[1]][1]
                ReceptorGeneID <- L_R[[1]][2]
                # Embed Hyper Link
                XYZ <- .hyperLinks(Ranking, XYZ, LigandGeneID,
                    ReceptorGeneID, LR, Value, Percentage, j,
                    spc, GeneInfo)

                # Plot (Ligand/Receptor)
                Ligandfile <- paste0(temp, "/figures/Ligand_",
                    LigandGeneID, ".png")
                png(filename=Ligandfile, width=1000, height=1000)
                .smallTwoDplot(input, LigandGeneID, twoD, "Reds")
                dev.off()
                Receptorfile <- paste0(temp, "/figures/Receptor_",
                    ReceptorGeneID, ".png")
                png(filename=Receptorfile, width=1000, height=1000)
                .smallTwoDplot(input, ReceptorGeneID, twoD, "Blues")
                dev.off()
            }
            # Write to Rmd
            sink(file = paste0(temp, "/", rmdfiles[i]))
            cat(XYZ)
            sink()
        }
    }

    # Plot (Mode-sum)
    png(filename=paste0(temp, "/figures/LRNetwork.png"), width=2000, height=950)
    .lrPlot(metadata(sce)$sctensor, twoDplot=twoD, label=celltypes)
    dev.off()
    # Plot (Each <L,R,*>)
    if(length(selected) != 0){
        for(i in seq_along(selected)){
            filenames <- paste0(temp, "/figures/LRNetwork_", index[i, 1],
                "_", index[i, 2], ".png")
            png(filename=filenames, width=2000, height=950)
            .lrPlot(metadata(sce)$sctensor, twoDplot=twoD, label=celltypes,
                emph=index[i, seq_len(2)])
            dev.off()
        }
    }
    # Number of Patterns
    numLPattern <- nrow(metadata(sce)$sctensor$ligand)
    numRPattern <- nrow(metadata(sce)$sctensor$receptor)
    col.ligand <- brewer.pal(9, "Reds")
    col.receptor <- brewer.pal(9, "Blues")

    for(i in seq_len(numLPattern)){
        label.ligand <- unlist(sapply(names(celltypes),
            function(x){
                metadata(sce)$sctensor$ligand[paste0("Dim", i), x]}))
        label.ligand[] <- smoothPalette(label.ligand,
            palfunc=colorRampPalette(col.ligand, alpha=TRUE))
        LPatternfile = paste0(temp, "/figures/Pattern_", i, "__", ".png")
        png(filename=LPatternfile, width=1000, height=1000)
        par(ps=20)
        plot(twoD, col=label.ligand, pch=16, cex=2, bty="n",
            xaxt="n", yaxt="n", xlab="", ylab="",
            main=paste0("(", i, ",*,*)-Pattern"))
        dev.off()
    }

    for(i in seq_len(numRPattern)){
        label.receptor <- unlist(sapply(names(celltypes),
            function(x){metadata(sce)$sctensor$receptor[paste0("Dim", i), x]}))
        label.receptor[] <- smoothPalette(label.receptor,
            palfunc=colorRampPalette(col.receptor, alpha=TRUE))
        RPatternfile = paste0(temp, "/figures/Pattern__", i, "_", ".png")
        png(filename=RPatternfile, width=1000, height=1000)
        par(ps=20)
        plot(twoD, col=label.receptor, pch=16, cex=2, bty="n",
            xaxt="n", yaxt="n", xlab="", ylab="",
            main=paste0("(*,", i, ",*)-Pattern"))
        dev.off()
    }

    # Header
    HEADER <- paste0("---\ntitle: XXXXX\n",
        "author: YYYYY\ndate:",
        " \"`r Sys.time()`\"\n",
        "output: ",
        "BiocStyle::html_document\nvignette: >\n ",
        "%\\VignetteIndexEntry{Vignette Title}\n ",
        "%\\VignetteEngine{knitr::rmarkdown}\n ",
        "%\\VignetteEncoding{UTF-8}\n---\n")
    HEADER <- sub("YYYYY", author, sub("XXXXX", title, HEADER))

    # 1. About Dataset and scTensor Algorithm
    BODY1 <- paste0("\n\n# About scTensor Algorithm\n\n",
        "![](Algorithm.png)\n",
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
    BODY2 <- paste0("\n\n# Global statistics and plots\n\n",
        "The result of scTensor is saved as a R binary file",
        " (reanalysis.RData).\n",
        "```{r}\n", # Top
        "load(\"reanalysis.RData\")\n\n",
        "# SingleCellExperiment object\n",
        "sce\n",
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
        "head(LR)\n",
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
        "```\n\n", # Bottom
        "\n\n## Number of cells in each celltype\n\n",
        "```{r}\n", # Top
        "color <- names(celltypes)\n",
        "colors <- celltypes[sapply(unique(color), ",
        "function(x){which(color == x)[1]})]\n",
        "numcells <- sapply(names(colors), function(x){",
        "length(which(color == x))",
        "})\n",
        "numcells <- data.frame(numcells=numcells, ",
        "celltypes=names(numcells))\n\n",
        "library(plotly)\n",
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
        "type = \"scatter\",",
        "text = rownames(twoD),",
        "mode = \"markers\")\n",
        "```\n\n", # Bottom
        "\n\n## Distribution of core tensor values\n\n",
        "```{r}\n", # Top
        "corenames <- sapply(seq_len(nrow(index)), ",
        "function(x){paste(index[x,seq_len(3)], collapse=\",\")})\n",
        "plot_ly(x=seq_along(corevalue), y=corevalue, ",
        "type=\"bar\", color=names(corevalue), text=corenames, ",
        "colors = c(\"#999999\", \"#E41A1C\"))\n",
        "```\n" # Bottom
        )

    # 3. L-Pattern
    BODY3 <- paste0("\n\n# Ligand-Cell Patterns\n\n",
        "```{r}\n", # Top
        "library(\"heatmaply\")\n",
        "l <- metadata(sce)$sctensor$ligand\n",
        "rownames(l) <- paste0(\"(\",seq_len(nrow(l)), \",*,*)\")\n",
        "heatmaply(l,",
        "xlab=\"Celltype\",",
        "ylab=\"Pattern\",",
        "fontsize_col=20,",
        "fontsize_row=20,",
        "subplot_widths=c(0.7, 0.1),",
        "subplot_heights=c(0.2, 0.7),",
        "labRow = rownames(l),",
        "labCol = colnames(l),",
        ")\n",
        "```\n" # Bottom
        )
    for(i in seq_len(numLPattern)){
        BODY3 <- paste0(BODY3, "![](figures/Pattern_", i, "__.png)")
    }
    BODY3 <- paste0(BODY3, "\n")


    # 4. R-Pattern
    BODY4 <- paste0("\n\n# Receptor-Cell Patterns\n\n",
        "```{r}\n", # Top
        "r <- metadata(sce)$sctensor$receptor\n",
        "rownames(r) <- paste0(\"(*,\",seq_len(nrow(r)), \",*)\")\n",
        "heatmaply(r,",
        "xlab=\"Celltype\",",
        "ylab=\"Pattern\",",
        "fontsize_col=20,",
        "fontsize_row=20,",
        "subplot_widths=c(0.7, 0.1),",
        "subplot_heights=c(0.2, 0.7),",
        "labRow = rownames(r),",
        "labCol = colnames(r),",
        ")\n",
        "```\n" # Bottom
        )
    for(i in seq_len(numRPattern)){
        BODY4 <- paste0(BODY4, "![](figures/Pattern__", i, "_.png)")
    }
    BODY4 <- paste0(BODY4, "\n")

    # 5. LR-Pattern
    BODY5 <- paste0("\n\n# LR-pair Patterns\n\n",
        "```{r}\n", # Top
        "lr <- metadata(sce)$sctensor$ppi\n",
        "rownames(lr) <- paste0(\"(*,*,\",seq_len(nrow(lr)), \")\")\n",
        "target <- which(rank(colMaxs(lr)) <= 100)\n",
        "heatmaply(lr[, target],",
        "xlab=\"LP-Pair\",",
        "ylab=\"Pattern\",",
        "fontsize_col=20,",
        "fontsize_row=20,",
        "subplot_widths=c(0.7, 0.1),",
        "subplot_heights=c(0.2, 0.7),",
        "labRow = rownames(lr[, target]),",
        "labCol = colnames(lr[, target]),",
        ")\n",
        "```\n" # Bottom
        )

    # 6. Ligand-Receptor Pattern
    BODY6 <- paste0("\n\n# Ligand-Receptor Patterns ",
        "(Mode-3/LR-pair direction sum of core tensor)\n\n",
        "![](figures/LRNetwork.png){ width=100% }\n")

    # 7. (Ligand, Receptor, LR-pair)-Pattern
    if(length(selected) != 0){
        htmlfiles <- gsub("Rmd", "html", rmdfiles)
        BODY7 <- selected
        for (i in selected){
            BODY7[i] <- paste0("\n\n## (", paste(index[i,seq_len(3)],
                collapse=","),
                ") Pattern : (", round(corevalue[i], 2), " %)\n",
                "[Details of (", paste(index[i,seq_len(3)],
                collapse=","),
                ") Pattern", "](", htmlfiles[i], ")\n")
        }
        BODY7 <- paste(BODY7, collapse = "\n")
        BODY7 <- paste0("# (Ligand-Cell, Receptor-Cell, LR-pair)",
            " Patterns\n\n", BODY7)
    }else{
        BODY7 <- "# (Ligand-Cell, Receptor-Cell, LR-pair) Patterns\n\n"
    }

    # 8. Session Information
    BODY8 <- paste0("\n\n# Session Information\n\n",
        "\n```{r}\n",
        "sessionInfo()\n",
        "```\n")

    # 9. License
    BODY9 <- paste0("\n\n# License\n\n",
        "Copyright (c) 2018 Koki Tsuyuzaki and Laboratory for ",
        "Bioinformatics Research, RIKEN Center for Biosystems Dynamics",
        " Reseach Released under the ",
        "[Artistic License 2.0](",
        "http://www.perlfoundation.org/artistic_license_2_0)\n")

    # Output
    cat("index.Rmd is created...\n")
    sink(file = paste0(temp, "/index.Rmd"))
    cat(paste0(HEADER, "\n\n"))
    cat(paste0(BODY1, "\n"))
    cat(paste0(BODY2, "\n"))
    cat(paste0(BODY3, "\n"))
    cat(paste0(BODY4, "\n"))
    cat(paste0(BODY5, "\n"))
    cat(paste0(BODY6, "\n"))
    cat(paste0(BODY7, "\n"))
    cat(paste0(BODY8, "\n"))
    cat(paste0(BODY9, "\n"))
    sink()

    # File Copy from Tempolary Directory
    if(temp != out.dir){
        file.copy(from = paste0(temp, "/index.Rmd"),
            to = paste0(out.dir, "/index.Rmd"), overwrite = TRUE)
        file.copy(from = paste0(temp, "/figures"),
            to = out.dir, recursive=TRUE)
        if(length(selected) != 0){
            for(i in seq_along(rmdfiles)){
                file.copy(from = paste0(temp, "/", rmdfiles[i]),
                to = out.dir, recursive=TRUE)
            }
        }
    }else{
        out.dir = temp
    }
    file.copy(
        from = system.file("extdata", "Algorithm.png", package = "scTensor"),
        to = paste0(out.dir, "/Algorithm.png"),
        overwrite = TRUE)

    # Rendering
    if(length(selected) != 0){
        for(i in seq_along(rmdfiles)){
            cat(paste0(rmdfiles[i], " is compiled to ",
                gsub(".Rmd", ".html", rmdfiles[i]),
                "...(", i, " / ", length(rmdfiles), ")\n"))
            e <- try(render(paste0(out.dir, "/", rmdfiles[i]),
                quiet=TRUE))
            if(is(e)[1] == "try-error"){
                e <- try(render(paste0(out.dir, "/", rmdfiles[i]),
                    quiet=TRUE))
                if(is(e)[1] == "try-error"){
                    e <- try(render(paste0(out.dir, "/", rmdfiles[i]),
                        quiet=TRUE))
                    if(is(e)[1] == "try-error"){
                        e <- try(render(paste0(out.dir, "/", rmdfiles[i]),
                            quiet=TRUE))
                    }
                }
            }
        }
    }
    cat("index.Rmd is compiled to index.html...\n")
    e <- try(render(paste0(out.dir, "/index.Rmd"),
        quiet=TRUE))
    if(is(e)[1] == "try-error"){
        e <- try(render(paste0(out.dir, "/index.Rmd"),
            quiet=TRUE))
        if(is(e)[1] == "try-error"){
            e <- try(render(paste0(out.dir, "/index.Rmd"),
                quiet=TRUE))
            if(is(e)[1] == "try-error"){
                e <- try(render(paste0(out.dir, "/index.Rmd"),
                    quiet=TRUE))
            }
        }
    }
    cat("cellCellReport is finished!!!\n")
    # HTML Open
    if(html.open){
        browseURL(paste0(out.dir, "/index.html"))
    }
}
