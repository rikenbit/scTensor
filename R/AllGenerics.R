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
# cellCellRanks
#
setGeneric("cellCellRanks", function(sce, centering=TRUE, mergeas="mean", outer="*", comb="random", num.sampling=100, ftt=TRUE, thr1=0.8, thr2=0.8, thr3=0.8){
    standardGeneric("cellCellRanks")})
setMethod("cellCellRanks",
    signature(sce="SingleCellExperiment"),
    function(sce, centering, mergeas, outer, comb, num.sampling, ftt, thr1, thr2, thr3){
        userobjects <- deparse(substitute(sce))
        .cellCellRanks(userobjects, centering, mergeas, outer, comb, num.sampling, ftt, thr1, thr2, thr3)
    })

.cellCellRanks <- function(userobjects, centering, mergeas, outer, comb, num.sampling, ftt, thr1, thr2, thr3){
    # Import from sce object
    sce <- eval(parse(text=userobjects))
    if(ftt){
        input <- .FTT(assay(sce))
    }else{
        input <- assay(sce)
    }
    LR <- LRBaseDbi::select(metadata(sce)$lrbase,
        columns=c("GENEID_L", "GENEID_R"),
        keytype="GENEID_L",
        keys=LRBaseDbi::keys(metadata(sce)$lrbase, keytype="GENEID_L"))
    celltypes <- metadata(sce)$color
    names(celltypes) <- metadata(sce)$label


    # Tensor is generated, and then matricised
    tnsr <- .cellCellDecomp.Third(input, LR, celltypes, ranks=c(3,3,3), centering, mergeas, outer, comb, num.sampling, decomp=FALSE)$cellcelllrpairpattern
    d1 <- svd(rs_unfold(tnsr, m=1)@data)$d
    d2 <- svd(rs_unfold(tnsr, m=2)@data)$d
    d3 <- svd(rs_unfold(tnsr, m=3)@data)$d
    cumd1 <- cumsum(d1) / sum(d1)
    cumd2 <- cumsum(d2) / sum(d2)
    cumd3 <- cumsum(d3) / sum(d3)
    # Output
    selected = c(
        length(which(cumd1 <= thr1)),
        length(which(cumd2 <= thr2)),
        length(which(cumd3 <= thr3))
        )

    list(selected=selected,
        mode1=d1,
        mode2=d2,
        mode3=d3)
}

#
# cellCellDecomp
#
setGeneric("cellCellDecomp", function(sce, algorithm="ntd", ranks=c(3,3,3), rank=3, thr1=log2(5), thr2=25, centering=TRUE, mergeas="mean", outer="*", comb="random", num.sampling=100, decomp=TRUE, ftt=TRUE){
    standardGeneric("cellCellDecomp")})
setMethod("cellCellDecomp", signature(sce="SingleCellExperiment"),
    function(sce, algorithm, ranks, rank, thr1, thr2, centering,
        mergeas, outer, comb, num.sampling, decomp, ftt){
        userobjects <- deparse(substitute(sce))
        .cellCellDecomp(userobjects, algorithm, ranks, rank, thr1, thr2, centering, mergeas, outer, comb, num.sampling, decomp, ftt)})

.cellCellDecomp <- function(userobjects, algorithm, ranks, rank, thr1, thr2, centering, mergeas, outer, comb, num.sampling, decomp, ftt){
    # Import from sce object
    sce <- eval(parse(text=userobjects))
    if(ftt){
        input <- .FTT(assay(sce))
    }else{
        input <- assay(sce)
    }
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
        res.sctensor <- .cellCellDecomp.Third(input, LR, celltypes, ranks, centering, mergeas, outer, comb, num.sampling, decomp)
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
    author="The person who runs this script", thr=80, top="full", cl=NULL){
    standardGeneric("cellCellReport")})
setMethod("cellCellReport", signature(sce="SingleCellExperiment"),
    function(sce, reducedDimNames, out.dir, html.open, title, author, thr, top, cl){
        userobjects <- deparse(substitute(sce))
        .cellCellReport(userobjects, reducedDimNames, out.dir,
            html.open, title, author, thr, top, cl)})
.cellCellReport <- function(userobjects, reducedDimNames,
    out.dir=NULL, html.open=FALSE,
    title="The result of scTensor",
    author="The person who runs this script", thr=80, top="full", cl=NULL){
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
    lrname <- LRBaseDbi::lrPackageName(metadata(sce)$lrbase)
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
    selected <- which(cumsum(corevalue) <= thr)
    if(length(selected) == 0){
        message(paste0("None of core tensor element is selected.\n",
        "Please specify the larger thr or perform cellCellDecomp\n",
        "with smaller ranks such as c(3,3,3)."))
    }
    names(corevalue) <- c(rep("selected", length=length(selected)),
        rep("not selected",
        length=length(corevalue) - length(selected)))

    # Tempolary Directory for saving the analytical result
    temp <- tempdir()
    dir.create(paste0(temp, "/figures"),
        showWarnings = FALSE, recursive = TRUE)

    # Table
    if(length(selected) != 0){
        # Plot (Each <L,R,*>)
        sapply(seq_along(selected), function(i){
            filenames <- paste0(temp,
                "/figures/CCIHypergraph_", index[i, 1],
                "_", index[i, 2], ".png")
            png(filename=filenames, width=2000, height=950)
            .CCIhyperGraphPlot(metadata(sce)$sctensor,
                twoDplot=twoD,
                label=celltypes,
                emph=index[i, seq_len(2)])
            dev.off()
        })
        SelectedLR <- sort(unique(index[selected, "Mode3"]))

        # Setting for Parallel Computing
        cat(paste0(length(SelectedLR),
            " LR vectors will be calculated :\n"))
        e <<- new.env()
        e$index <- index
        e$sce <- sce
        e$.HCLUST <- .HCLUST
        e$.OUTLIERS <- .OUTLIERS
        e$top <- top
        e$spc <- spc
        e$.sapply_pb <- .sapply_pb
        e$GeneInfo <- GeneInfo
        e$temp <- temp
        e$.smallTwoDplot <- .smallTwoDplot
        e$input <- input
        e$twoD <- twoD
        e$.hyperLinks <- .hyperLinks
        e$LR <- LR
        e$.eachVecLR <- .eachVecLR
        e$.eachRender <- .eachRender
        e$.ENRICHMENT <- .ENRICHMENT
        e$.XYZ_HEADER1 <- .XYZ_HEADER1
        e$.XYZ_HEADER2 <- .XYZ_HEADER2
        e$.XYZ_HEADER3 <- .XYZ_HEADER3
        e$.XYZ_ENRICH <- .XYZ_ENRICH

        if (!is.null(cl)) {
            ############ Parallel ############
            # Package Loading in each node
            invisible(clusterEvalQ(cl, {
                library("outliers")
                library("S4Vectors")
                library("tagcloud")
                library("RColorBrewer")
                library("plotrix")
                library("plotly")
                library("rmarkdown")
                library("meshr")
                library("GOstats")
                library("ReactomePA")
            }))
            clusterExport(cl, "e")
            out.vecLR <- parSapply(cl, SelectedLR,
                function(x, e){.eachVecLR(x, e)}, e=e)
            colnames(out.vecLR) <- paste0("pattern", SelectedLR)
            e$out.vecLR <- out.vecLR
            clusterExport(cl, "e")
            ############ Parallel ############
        }else{
            out.vecLR <- sapply(SelectedLR,
                function(x, e){.eachVecLR(x, e)}, e=e)
            colnames(out.vecLR) <- paste0("pattern", SelectedLR)
            e$out.vecLR <- out.vecLR
        }
    }

    # Tagcloud
    invisible(.tagCloud(out.vecLR, temp))

    # Plot（CCI Hypergraph）
    png(filename=paste0(temp, "/figures/CCIHypergraph.png"), width=2000, height=950)
    .CCIhyperGraphPlot(metadata(sce)$sctensor, twoDplot=twoD, label=celltypes)
    dev.off()

    # Plot（Gene-wise Hypergraph）
    .geneHyperGraphPlot(out.vecLR, GeneInfo, temp)

    # Rmd（ligand）
    cat("ligand.Rmd is created...\n")
    sink(file = paste0(temp, "/ligand.Rmd"))
    cat(.LIGAND_HEADER)
    cat(paste0(.LIGAND_BODY(out.vecLR, GeneInfo, index, selected)
        , "\n"))
    sink()

    # Rmd（receptor）
    cat("receptor.Rmd is created...\n")
    sink(file = paste0(temp, "/receptor.Rmd"))
    cat(.RECEPTOR_HEADER)
    cat(paste0(.RECEPTOR_BODY(out.vecLR, GeneInfo, index, selected)
        , "\n"))
    sink()

    # Number of Patterns
    vecL <- metadata(sce)$sctensor$ligand
    vecR <- metadata(sce)$sctensor$receptor
    numLPattern <- nrow(vecL)
    numRPattern <- nrow(vecR)
    col.ligand <- brewer.pal(9, "Reds")
    col.receptor <- brewer.pal(9, "Blues")
    # Clustering
    ClusterL <- t(apply(vecL, 1, .HCLUST))
    ClusterR <- t(apply(vecR, 1, .HCLUST))

    # Ligand Pattern
    sapply(seq_len(numLPattern), function(i){
        label.ligand <- unlist(sapply(names(celltypes),
            function(x){
                metadata(sce)$sctensor$ligand[paste0("Dim", i), x]}))
        label.ligand[] <- smoothPalette(label.ligand,
            palfunc=colorRampPalette(col.ligand, alpha=TRUE))
        ClusterNameL <- paste(names(which(ClusterL[i,] == "selected")), collapse=" & ")
        titleL <- paste0("(", i, ",*,*)-Pattern", " = ", ClusterNameL)
        titleL <- .shrink(titleL)
        LPatternfile <- paste0(temp, "/figures/Pattern_", i, "__", ".png")
        png(filename=LPatternfile, width=1000, height=1000)
        par(ps=20)
        plot(twoD, col=label.ligand, pch=16, cex=2, bty="n",
            xaxt="n", yaxt="n", xlab="", ylab="",
            main=titleL)
        dev.off()
    })

    # Receptor Pattern
    sapply(seq_len(numRPattern), function(i){
        label.receptor <- unlist(sapply(names(celltypes),
            function(x){metadata(sce)$sctensor$receptor[paste0("Dim", i), x]}))
        label.receptor[] <- smoothPalette(label.receptor,
            palfunc=colorRampPalette(col.receptor, alpha=TRUE))
        ClusterNameR <- paste(names(which(ClusterR[i,] == "selected")), collapse=" & ")
        titleR <- paste0("(*,", i, ",*)-Pattern", " = ", ClusterNameR)
        titleR <- .shrink(titleR)
        RPatternfile = paste0(temp, "/figures/Pattern__", i, "_", ".png")
        png(filename=RPatternfile, width=1000, height=1000)
        par(ps=20)
        plot(twoD, col=label.receptor, pch=16, cex=2, bty="n",
            xaxt="n", yaxt="n", xlab="", ylab="",
            main=titleR)
        dev.off()
    })

    # Save the result of scTensor
    save(sce, input, twoD, LR, celltypes, index, corevalue,
        selected, ClusterL, ClusterR, out.vecLR,
        file=paste0(temp, "/reanalysis.RData"))

    # Rendering
    cat("ligand.Rmd is compiled to index.html...\n")
    render(paste0(temp, "/ligand.Rmd"), quiet=TRUE)
    cat("receptor.Rmd is compiled to index.html...\n")
    render(paste0(temp, "/receptor.Rmd"), quiet=TRUE)
    if (!is.null(cl)) {
        cat(paste0(length(selected),
            " pattern_X_Y_Z.Rmd files will be created :\n"))
        parSapply(cl, selected,
            function(x, e, SelectedLR){.eachRender(x, e, SelectedLR)},
            e=e, SelectedLR=SelectedLR)
    }else{
        cat(paste0(length(selected),
            " pattern_X_Y_Z.Rmd files will be created :\n"))
        sapply(selected,
            function(x, e, SelectedLR){.eachRender(x, e, SelectedLR)}, e=e, SelectedLR=SelectedLR)
    }

    # File copy
    file.copy(
        from = system.file("extdata", "Workflow.jpeg", package = "scTensor"),
        to = paste0(temp, "/Workflow.jpeg"),
        overwrite = TRUE)

    # Output index.html
    if(length(selected) != 0){
        if(length(selected) == 1){
            RMDFILES <- apply(t(index[selected, seq_len(3)]), 1,
                function(x){
                paste0(c("pattern", x), collapse="_")
            })
        }else{
            RMDFILES <- apply(index[selected, seq_len(3)], 1,
                function(x){
                paste0(c("pattern", x), collapse="_")
            })
        }
        RMDFILES <- paste0(RMDFILES, ".Rmd")
    }
    cat("index.Rmd is created...\n")
    sink(file = paste0(temp, "/index.Rmd"))
    cat(paste0(.MAINHEADER(author, title), "\n\n"))
    cat(paste0(.BODY1, "\n"))
    cat(paste0(.BODY2, "\n"))
    cat(paste0(.BODY3(numLPattern), "\n"))
    cat(paste0(.BODY4(numRPattern), "\n"))
    cat(paste0(.BODY5, "\n"))
    cat(paste0(.BODY6, "\n"))
    cat(paste0(.BODY7, "\n"))
    if(length(selected) != 0){
        cat(paste0(.BODY8(selected, RMDFILES, index, corevalue), "\n"))
    }
    cat(paste0(.BODY9, "\n"))
    cat(paste0(.BODY10, "\n"))
    sink()

    # Rendering
    cat("index.Rmd is compiled to index.html...\n")
    render(paste0(temp, "/index.Rmd"), quiet=TRUE)

    # File Copy from Tempolary Directory
    if(temp != out.dir){
        file.copy(from = temp, to = out.dir,
            overwrite = TRUE, recursive = TRUE)
    }else{
        out.dir = temp
    }

    # HTML Open
    if(html.open){
        browseURL(paste0(out.dir, "/index.html"))
    }
}
