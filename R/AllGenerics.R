#
# cellCellSetting
#
setGeneric("cellCellSetting", function(sce, lrbase, color, label){
    standardGeneric("cellCellSetting")})
setMethod("cellCellSetting", signature(sce="SingleCellExperiment"),
    function(sce, lrbase, color, label){
        userobjects <- deparse(substitute(sce))
        .cellCellSetting(userobjects, lrbase, color, label, sce)})
.cellCellSetting <- function(userobjects, lrbase, color, label, ...){
    sce <- list(...)[[1]]
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
    # Rowname-check
    rn <- rownames(assay(sce))
    if(length(rn) != length(unique(rn))){
        stop("Please specify the row names of the input matrix is unique")
    }
    # Only matrix is accepted
    if(!is.matrix(assay(sce))){
        message("The input data is coverted to matrix format by as.matrix")
        assay(sce) <- as.matrix(assay(sce))
    }
    # NA-check
    NA1 <- length(which(is.na(colnames(assay(sce)))))
    NA2 <- length(which(is.na(label)))
    NA3 <- length(which(is.na(color)))
    if(NA1 != 0){
        stop("At least one NA element is in colnames(assay(sce))\nPlease remove it.")
    }
    if(NA2 != 0){
        stop("At least one NA element is in label\nPlease remove it.")
    }
    if(NA3 != 0){
        stop("At least one NA element is in color\nPlease remove it.")
    }
    # Overwrite
    metadata(sce) <- list(lrbase=lrbase$conn@dbname, color=color, label=label)
    assign(userobjects, sce, envir=.GlobalEnv)
}

#
# cellCellRanks
#
setGeneric("cellCellRanks", function(sce, centering=TRUE,
    mergeas=c("mean", "sum"), outer=c("*", "+"), comb=c("random", "all"),
    num.sampling=100, ftt=TRUE, thr1=0.9, thr2=0.9, thr3=0.9){
    standardGeneric("cellCellRanks")})
setMethod("cellCellRanks",
    signature(sce="SingleCellExperiment"),
    function(sce, centering=TRUE,
    mergeas=c("mean", "sum"), outer=c("*", "+"), comb=c("random", "all"),
    num.sampling=100, ftt=TRUE, thr1=0.9, thr2=0.9, thr3=0.9){
        # Argument Check
        mergeas <- match.arg(mergeas)
        outer <- match.arg(outer)
        comb <- match.arg(comb)
        .cellCellRanks(centering, mergeas, outer, comb, num.sampling,
            ftt, thr1, thr2, thr3, sce)
    })
.cellCellRanks <- function(centering, mergeas, outer, comb, num.sampling,
    ftt, thr1, thr2, thr3, ...){
    # Import from sce object
    sce <- list(...)[[1]]
    if(ftt){
        input <- .FTT(assay(sce))
    }else{
        input <- assay(sce)
    }
    con = dbConnect(SQLite(), metadata(sce)$lrbase)
    LR <- dbGetQuery(con, "SELECT * FROM DATA")[, c("GENEID_L", "GENEID_R")]
    LR <- unique(LR)
    dbDisconnect(con)
    celltypes <- metadata(sce)$color
    names(celltypes) <- metadata(sce)$label

    # Tensor is generated, and then matricised
    tnsr <- .cellCellDecomp.Third(input, LR, celltypes, ranks=c(1,1,1),
        rank=1, centering, mergeas, outer, comb, num.sampling, decomp=FALSE,
        thr1, thr2)$cellcelllrpairpattern
    d1 <- svd(rs_unfold(tnsr, m=1)@data)$d
    d2 <- svd(rs_unfold(tnsr, m=2)@data)$d
    d3 <- svd(rs_unfold(tnsr, m=3)@data)$d
    cumd1 <- cumsum(d1) / sum(d1)
    cumd2 <- cumsum(d2) / sum(d2)
    cumd3 <- cumsum(d3) / sum(d3)
    # Output
    rank1 <- length(which(cumd1 <= thr1))
    rank2 <- length(which(cumd2 <= thr2))
    rank3 <- length(which(cumd3 <= thr3))
    if(rank1 == 0){
        rank1 = 1
    }
    if(rank2 == 0){
        rank2 = 1
    }
    if(rank3 == 0){
        rank3 = 1
    }
    selected = c(rank1, rank2, rank3)
    list(selected=selected,
        mode1=d1,
        mode2=d2,
        mode3=d3)
}

#
# cellCellDecomp
#
setGeneric("cellCellDecomp", function(sce, algorithm=c("ntd", "nmf", "pearson",
    "spearman", "distance", "pearson.lr", "spearman.lr", "distance.lr",
    "pcomb", "label.permutation"), ranks=c(3,3,3), rank=3, thr1=log2(5), thr2=25,
    centering=TRUE, mergeas=c("mean", "sum"), outer=c("*", "+"),
    comb=c("random", "all"), num.sampling=100, decomp=TRUE, ftt=TRUE){
    standardGeneric("cellCellDecomp")})
setMethod("cellCellDecomp", signature(sce="SingleCellExperiment"),
    function(sce, algorithm=c("ntd", "nmf", "pearson", "spearman", "distance",
        "pearson.lr", "spearman.lr", "distance.lr", "pcomb", "label.permutation"), ranks=c(3,3,3),
        rank=3, thr1=log2(5), thr2=25, centering=TRUE,
        mergeas=c("mean", "sum"), outer=c("*", "+"), comb=c("random", "all"),
        num.sampling=100, decomp=TRUE, ftt=TRUE){
        # Argument Check
        algorithm = match.arg(algorithm)
        mergeas = match.arg(mergeas)
        outer = match.arg(outer)
        comb = match.arg(comb)

        userobjects <- deparse(substitute(sce))
        .cellCellDecomp(userobjects, algorithm, ranks, rank, thr1, thr2,
            centering, mergeas, outer, comb, num.sampling, decomp, ftt, sce)})

.cellCellDecomp <- function(userobjects, algorithm, ranks, rank, thr1,
    thr2, centering, mergeas, outer, comb, num.sampling, decomp, ftt, ...){
    # Import from sce object
    sce <- list(...)[[1]]
    if(ftt){
        input <- .FTT(assay(sce))
    }else{
        input <- assay(sce)
    }
    con = dbConnect(SQLite(), metadata(sce)$lrbase)
    LR <- dbGetQuery(con, "SELECT * FROM DATA")[, c("GENEID_L", "GENEID_R")]
    LR <- unique(LR)
    dbDisconnect(con)
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
    # thr-check
    if(comb == "random" && (num.sampling <= 0)){
        warning("None of cell-cell interaction will be detected.")
    }
    # Corresponding function of the algorithm
    f <- .flist[[algorithm]]
    if(is.null(f)){
        stop("Please specify the appropriate algorithm\n")
    }else{
        res.sctensor <- f(input, LR, celltypes, ranks, rank, centering,
            mergeas, outer, comb, num.sampling, decomp, thr1, thr2)
    }

    # Data size
    datasize <- c(ncol(res.sctensor[[2]]), ncol(res.sctensor[[3]]),
        ncol(res.sctensor[[4]]))
    recerror <- res.sctensor$recerror
    relchange <- res.sctensor$relchange

    # Overwrite
    metadata(sce) <- list(lrbase=metadata(sce)$lrbase,
        color=metadata(sce)$color, label=metadata(sce)$label,
        algorithm=algorithm, sctensor=res.sctensor, ranks=ranks,
        datasize=datasize, recerror=recerror, relchange=relchange)
    assign(userobjects, sce, envir=.GlobalEnv)
}

#
# cellCellReport
#
setGeneric("cellCellReport", function(sce, reducedDimNames,
    out.dir=tempdir(), html.open=FALSE,
    title="The result of scTensor",
    author="The person who runs this script", thr=40, top="full", p=0.05){
    standardGeneric("cellCellReport")})
setMethod("cellCellReport", signature(sce="SingleCellExperiment"),
    function(sce, reducedDimNames, out.dir, html.open, title, author,
        thr, top, p){
        .cellCellReport(reducedDimNames, out.dir,
            html.open, title, author, thr, top, p, sce)})
.cellCellReport <- function(reducedDimNames,
    out.dir=tempdir(), html.open=FALSE,
    title="The result of scTensor",
    author="The person who runs this script", thr=40, top="full", p=0.05, ...){
    # Import from sce object
    sce <- list(...)[[1]]
    # algorithm-check
    if(metadata(sce)$algorithm != "ntd"){
        stop(paste0("cellCellReport can be performed by the result of",
            " cellCellDecomp in which the algorithm is ",
            "specified as 'ntd' for now."))
    }
    # The Directory for saving the analytical result
    dir.create(paste0(out.dir, "/figures"),
        showWarnings = FALSE, recursive = TRUE)
    dir.create(paste0(out.dir, "/figures/Ligand"),
        showWarnings = FALSE, recursive = TRUE)
    dir.create(paste0(out.dir, "/figures/Receptor"),
        showWarnings = FALSE, recursive = TRUE)
    dir.create(paste0(out.dir, "/figures/Tagcloud"),
        showWarnings = FALSE, recursive = TRUE)
    # File copy
    file.copy(
        from = system.file("extdata", "Workflow.jpeg", package = "scTensor"),
        to = paste0(out.dir, "/Workflow.jpeg"),
        overwrite = TRUE)

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
    }else{
        names(corevalue) <- c(rep("selected", length=length(selected)),
            rep("not selected",
            length=length(corevalue) - length(selected)))
        # Data matrix
        input <- assay(sce)
        # Low dimensional data
        twoD <- eval(parse(text=paste0("reducedDims(sce)$", reducedDimNames)))
        # Ligand-Receptor, PMID
        con = dbConnect(SQLite(), metadata(sce)$lrbase)
        LR <- dbGetQuery(con, "SELECT * FROM DATA")[,
            c("GENEID_L", "GENEID_R", "SOURCEID")]
        LR <- unique(LR)
        dbDisconnect(con)
        # Species
        spc <- gsub(".eg.db.sqlite", "",
            strsplit(metadata(sce)$lrbase, "LRBase.")[[1]][3])
        # biomaRt Setting
        ens <- .ensembl[[spc]]()
        Sys.sleep(10)
        # GeneName, Description, GO, Reactome, MeSH
        GeneInfo <- .geneinformation(sce, ens, spc, LR)
        # Cell Label
        celltypes <- metadata(sce)$color
        names(celltypes) <- metadata(sce)$label

        # Plot Ligand/Receptor Genes
        invisible(.genePlot(input, twoD, out.dir, GeneInfo, LR))
        # Plot (Each <L,R,*>)
        vapply(seq_along(selected), function(i){
            filenames <- paste0(out.dir,
                "/figures/CCIHypergraph_", index[i, 1],
                "_", index[i, 2], ".png")
            png(filename=filenames, width=2000, height=950)
            invisible(.CCIhyperGraphPlot(metadata(sce)$sctensor,
                twoDplot=twoD,
                label=celltypes,
                emph=index[i, seq_len(2)]))
            dev.off()
        }, 0L)
        # <*,*,LR>
        SelectedLR <- sort(unique(index[selected, "Mode3"]))

        # Setting for Parallel Computing
        message(paste0(length(SelectedLR),
            " LR vectors will be calculated :"))
        e <<- new.env()
        e$p <- p
        e$index <- index
        e$sce <- sce
        e$.HCLUST <- .HCLUST
        e$.OUTLIERS <- .OUTLIERS
        e$top <- top
        e$spc <- spc
        e$GeneInfo <- GeneInfo
        e$out.dir <- out.dir
        e$.smallTwoDplot <- .smallTwoDplot
        e$input <- input
        e$twoD <- twoD
        e$.hyperLinks <- .hyperLinks
        e$LR <- LR
        e$.eachVecLR <- .eachVecLR
        e$.eachRender <- .eachRender
        e$.XYZ_HEADER1 <- .XYZ_HEADER1
        e$.XYZ_HEADER2 <- .XYZ_HEADER2
        e$.XYZ_HEADER3 <- .XYZ_HEADER3
        e$.XYZ_ENRICH <- .XYZ_ENRICH

        # EachVec（Heavy...）
        out.vecLR <- vapply(SelectedLR,
            function(x, e){.eachVecLR(x, e)},
            FUN.VALUE=rep(list(0L), 9), e=e)
        colnames(out.vecLR) <- paste0("pattern", SelectedLR)
        e$out.vecLR <- out.vecLR

        # Tagcloud
        invisible(.tagCloud(out.vecLR, out.dir))
        # Plot（CCI Hypergraph）
        png(filename=paste0(out.dir, "/figures/CCIHypergraph.png"),
            width=2000, height=950)
        invisible(.CCIhyperGraphPlot(metadata(sce)$sctensor,
            twoDplot=twoD, label=celltypes))
        dev.off()
        # Plot（Gene-wise Hypergraph）
        invisible(.geneHyperGraphPlot(out.vecLR, GeneInfo, out.dir))

        # Rmd（ligand, selected）
        message("ligand.Rmd is created...")
        outLg <- file(paste0(out.dir, "/ligand.Rmd"), "w")
        writeLines(.LIGAND_HEADER, outLg, sep="\n")
        writeLines(.LIGAND_BODY(out.vecLR, GeneInfo, index, selected),
            outLg, sep="\n")
        close(outLg)
        # Rmd（receptor, selected）
        message("receptor.Rmd is created...")
        outRp <- file(paste0(out.dir, "/receptor.Rmd"), "w")
        writeLines(.RECEPTOR_HEADER, outRp, sep="\n")
        writeLines(.RECEPTOR_BODY(out.vecLR, GeneInfo, index, selected),
            outRp, sep="\n")
        close(outRp)
        # Rmd（ligand, all）
        message("ligand_all.Rmd is created...")
        outLg_all <- file(paste0(out.dir, "/ligand_all.Rmd"), "w")
        writeLines(.LIGANDALL_HEADER, outLg_all, sep="\n")
        writeLines(.LIGANDALL_BODY(GeneInfo, LR, input),
            outLg_all, sep="\n")
        close(outLg_all)
        # Rmd（receptor, all）
        message("receptor_all.Rmd is created...")
        outRp_all <- file(paste0(out.dir, "/receptor_all.Rmd"), "w")
        writeLines(.RECEPTORALL_HEADER, outRp_all, sep="\n")
        writeLines(.RECEPTORALL_BODY(GeneInfo, LR, input),
            outRp_all, sep="\n")
        close(outRp_all)

        # Number of Patterns
        vecL <- metadata(sce)$sctensor$ligand
        vecR <- metadata(sce)$sctensor$receptor
        numLPattern <- nrow(vecL)
        numRPattern <- nrow(vecR)
        col.ligand <- .setColor("reds")
        col.receptor <- .setColor("blues")
        # Clustering
        ClusterL <- t(apply(vecL, 1, .HCLUST))
        ClusterR <- t(apply(vecR, 1, .HCLUST))
        # Ligand Pattern
        invisible(.ligandPatternPlot(numLPattern, celltypes, sce, col.ligand, ClusterL, out.dir, twoD))
        # Receptor Pattern
        invisible(.receptorPatternPlot(numRPattern, celltypes, sce,
            col.receptor, ClusterR, out.dir, twoD))
        # Save the result of scTensor
        save(sce, input, twoD, LR, celltypes, index, corevalue,
            selected, ClusterL, ClusterR, out.vecLR,
            file=paste0(out.dir, "/reanalysis.RData"))

        # Rendering
        message("ligand.Rmd is compiled to ligand.html...")
        render(paste0(out.dir, "/ligand.Rmd"), quiet=TRUE)
        message("ligand_all.Rmd is compiled to ligand_all.html...")
        render(paste0(out.dir, "/ligand_all.Rmd"), quiet=TRUE)
        message("receptor.Rmd is compiled to receptor.html...")
        render(paste0(out.dir, "/receptor.Rmd"), quiet=TRUE)
        message("receptor_all.Rmd is compiled to receptor_all.html...")
        render(paste0(out.dir, "/receptor_all.Rmd"), quiet=TRUE)
        message(paste0(length(selected),
            " pattern_X_Y_Z.Rmd files are compiled to pattern_X_Y_Z.html :"))
        vapply(selected,
            function(x, e, SelectedLR){
                .eachRender(x, e, SelectedLR)}, "", e=e, SelectedLR=SelectedLR)

        # Output index.html
        RMDFILES <- vapply(selected, function(x){
            paste0(paste(c("pattern", index[x, seq_len(3)]),
                collapse="_"), ".Rmd")
        }, "")
        message("index.Rmd is created...")
        outIdx <- file(paste0(out.dir, "/index.Rmd"), "w")
        writeLines(.MAINHEADER(author, title), outIdx, sep="\n")
        writeLines(.BODY1, outIdx, sep="\n")
        writeLines(.BODY2, outIdx, sep="\n")
        writeLines(.BODY3(numLPattern, ClusterL), outIdx, sep="\n")
        writeLines(.BODY4(numRPattern, ClusterR), outIdx, sep="\n")
        writeLines(.BODY5, outIdx, sep="\n")
        writeLines(.BODY6, outIdx, sep="\n")
        writeLines(.BODY7, outIdx, sep="\n")
        if(length(selected) != 0){
            writeLines(.BODY8(selected, RMDFILES, index, corevalue),
                outIdx, sep="\n")
        }
        writeLines(.BODY9, outIdx, sep="\n")
        writeLines(.BODY10, outIdx, sep="\n")
        close(outIdx)

        # Rendering
        message("index.Rmd is compiled to index.html...")
        render(paste0(out.dir, "/index.Rmd"), quiet=TRUE)
        if(out.dir == tempdir()){
            message(paste0("################################################\n",
                "Data files are saved in\n",
                out.dir, "\n################################################\n"))
        }
        # HTML Open
        if(html.open){
            browseURL(paste0(out.dir, "/index.html"))
        }
    }
}

#
# cellCellSimulate-related functions
#
setGeneric("getParam", function(object, name){
    standardGeneric("getParam")})
setGeneric("setParam<-", function(object, name, value){
    standardGeneric("setParam<-")})
setGeneric("show", function(object){
    standardGeneric("show")})

#
# cellCellSimulate
#
cellCellSimulate <- function(params = newCCSParams(), verbose = TRUE){
    # Class Check
    assertClass(params, classes = "CCSParams")
    # Get the parameters
    if(verbose){message("Getting the values of params...")}
    nGene <- getParam(params, "nGene")
    nCell <- getParam(params, "nCell")
    cciInfo <- getParam(params, "cciInfo")
    lambda <- getParam(params, "lambda")
    # Set random seed
    if(verbose){message("Setting random seed...")}
    seed <- getParam(params, "seed")
    # Simulation data
    if(verbose){message("Generating simulation data...")}
    # Generate Simulation data
    out <- .simulateDropoutCounts(nGene, nCell, cciInfo, lambda, seed)
    input <- out$simcount
    LR <- out$LR
    celltypes <- out$celltypes
    LR_CCI <- out$LR_CCI
    # Output
    if(verbose){message("Done!")}
    list(input=input, LR=LR, celltypes=celltypes, LR_CCI=LR_CCI)
}

#
# convertToNCBIGeneID
#
convertToNCBIGeneID <- function(input, rowID, LefttoRight){
    # Argument check
    if(dim(input)[1] != length(rowID)){
        stop("The number of rows of input and the length of rowID must be same.")
    }
    # NA check
    nr <- nrow(input)
    nc <- ncol(input)
    notNA <- which(!is.na(rowID))
    input <- input[notNA,]
    rowID <- rowID[notNA]
    LefttoRight <- LefttoRight[which(!is.na(LefttoRight[,1])), ]
    LefttoRight <- LefttoRight[which(!is.na(LefttoRight[,2])), ]
    # Statistics of input matrix
    original.var <- apply(input, 1, var)
    original.mean <- apply(input, 1, mean)
    score <- original.var * original.mean
    # Common ID
    names(original.var) <- rowID
    search.ID <- intersect(unique(rowID), unique(LefttoRight[,1]))
    if(length(search.ID) == 0){
        stop("There are none of common ID in rowID and LefttoRight[,1].")
    }
    # ID Conversion
    targetGeneID <- unlist(lapply(search.ID, function(x){
        target <- which(LefttoRight[,1] == x)[1]
        LefttoRight[target, 2]
    }))
    names(targetGeneID) <- search.ID
    position.input <- unlist(lapply(seq_along(targetGeneID), function(x){
        target <- which(rowID == names(targetGeneID)[x])
        if(length(target) != 1){
            targetID <- target[which(score[target] == max(score[target]))]
        }else{
            targetID <- target
        }
        targetID
    }))
    input <- input[position.input, ]
    rownames(input) <- targetGeneID
    input <- as.matrix(input)
    # Before <-> After
    dif <- nr - nrow(input)
    if(dif > 0){
        message(paste0(dif, " of genes are removed from input matrix (",
            nr, "*", nc, ")"))
    }
    input
}
