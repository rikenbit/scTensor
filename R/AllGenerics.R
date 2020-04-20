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
    for(i in names(assays(sce))){
        if(!is.matrix(assays(sce)[[i]])){
            message("The input data is coverted to matrix format by as.matrix")
            assays(sce)[[i]] <- as.matrix(assays(sce)[[i]])
        }
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
    mergeas=c("mean", "sum"), outerfunc=c("*", "+"), comb=c("random", "all"),
    num.sampling=100, num.perm=1000, assayNames="counts", verbose=FALSE,
    num.iter1=5, num.iter2=5, num.iter3=NULL){
    standardGeneric("cellCellRanks")})

setMethod("cellCellRanks",
    signature(sce="SingleCellExperiment"),
    function(sce, centering=TRUE,
    mergeas=c("mean", "sum"), outerfunc=c("*", "+"), comb=c("random", "all"),
    num.sampling=100, num.perm=1000, assayNames="counts", verbose=FALSE,
    num.iter1=5, num.iter2=5, num.iter3=NULL){
        # Argument Check
        mergeas <- match.arg(mergeas)
        outerfunc <- match.arg(outerfunc)
        comb <- match.arg(comb)
        .cellCellRanks(centering, mergeas, outerfunc, comb, num.sampling, num.perm,
            assayNames, verbose, num.iter1, num.iter2, num.iter3, sce)
    })

.cellCellRanks <- function(centering, mergeas, outerfunc, comb, num.sampling, num.perm,
    assayNames, verbose, num.iter1, num.iter2, num.iter3=NULL, ...){
    # value-check
    if(num.iter1 < 0){
        stop("Please specify the num.iter1 as positive integer")
    }
    if(num.iter2 < 0){
        stop("Please specify the num.iter2 as positive integer")
    }

    # Import from sce object
    sce <- list(...)[[1]]
    # Import expression matrix
    input <- .importAssays(sce, assayNames)
    con = dbConnect(SQLite(), metadata(sce)$lrbase)
    LR <- dbGetQuery(con, "SELECT * FROM DATA")[, c("GENEID_L", "GENEID_R")]
    LR <- .uniqueLR(LR)
    dbDisconnect(con)
    celltypes <- metadata(sce)$color
    names(celltypes) <- metadata(sce)$label
    l <- length(unique(celltypes))

    # Tensor is generated
    tnsr <- .cellCellDecomp.Third(input, LR, celltypes, ranks=c(1,1,1),
        rank=1, centering, mergeas, outerfunc, comb, num.sampling,
        num.perm, decomp=FALSE, thr1=log2(5), thr2=25, verbose)$cellcelllrpairpattern

    # Limit
    l1 <- min(dim(tnsr)[1], dim(tnsr)[2]*dim(tnsr)[3])
    l2 <- min(dim(tnsr)[2], dim(tnsr)[3]*dim(tnsr)[1])

    # NMF in matricised tensors in each mode
    out1 <- NMF(cs_unfold(tnsr, m=1)@data, runtime=num.iter1, rank.method="rss", J=1:l1)
    out2 <- NMF(cs_unfold(tnsr, m=2)@data, runtime=num.iter2, rank.method="rss", J=1:l2)

    # Reconsturction Error
    rss1 <- unlist(lapply(seq(l1), function(x, out1){
        eval(parse(text=paste0("mean(out1$Trial$Rank", x, "$original)")))
    }, out1=out1))
    rss2 <- unlist(lapply(seq(l2), function(x, out2){
        eval(parse(text=paste0("mean(out2$Trial$Rank", x, "$original)")))
    }, out2=out2))

    # Estimated rank
    rank1 <- min(which((max(rss1) - rss1) / (max(rss1) - min(rss1)) > 0.8))
    rank2 <- min(which((max(rss2) - rss2) / (max(rss2) - min(rss2)) > 0.8))

    if(!is.null(num.iter3)){
        # Limit
        l3 <- min(30, dim(tnsr)[3], dim(tnsr)[1]*dim(tnsr)[2])
        # NMF in matricised tensors in each mode
        out3 <- NMF(cs_unfold(tnsr, m=3)@data, runtime=num.iter3, rank.method="rss", J=1:l3)
        # Reconsturction Error
        rss3 <- unlist(lapply(seq(l3), function(x, out3){
            eval(parse(text=paste0("mean(out3$Trial$Rank", x, "$original)")))
        }, out3=out3))        
        # Estimated rank
        rank3 <- min(which((max(rss3) - rss3) / (max(rss3) - min(rss3)) > 0.8))
        list(selected=c(rank1, rank2, rank3))
    }else{
        list(selected=c(rank1, rank2))
    }
}

#
# cellCellDecomp
#
setGeneric("cellCellDecomp", function(sce,
    algorithm=c("ntd2", "ntd", "nmf", "pearson",
    "spearman", "distance", "pearson.lr", "spearman.lr", "distance.lr",
    "pcomb", "label.permutation", "cabello.aguilar", "halpern"), ranks=c(3,3),
    rank=3, thr1=log2(5), thr2=25, verbose=FALSE,
    centering=TRUE, mergeas=c("mean", "sum"), outerfunc=c("*", "+"),
    comb=c("random", "all"), num.sampling=100, num.perm=1000,
    assayNames="counts", decomp=TRUE){
    standardGeneric("cellCellDecomp")})

setMethod("cellCellDecomp", signature(sce="SingleCellExperiment"),
    function(sce,
        algorithm=c("ntd2", "ntd", "nmf", "pearson", "spearman", "distance",
        "pearson.lr", "spearman.lr", "distance.lr", "pcomb",
        "label.permutation", "cabello.aguilar", "halpern"), ranks=c(3,3),
        rank=3, thr1=log2(5), thr2=25, verbose=FALSE, centering=TRUE,
        mergeas=c("mean", "sum"), outerfunc=c("*", "+"), comb=c("random", "all"),
        num.sampling=100, num.perm=1000, assayNames="counts", decomp=TRUE){
        # Argument Check
        algorithm = match.arg(algorithm)
        mergeas = match.arg(mergeas)
        outerfunc = match.arg(outerfunc)
        comb = match.arg(comb)

        userobjects <- deparse(substitute(sce))
        .cellCellDecomp(userobjects, algorithm, ranks, rank, thr1, thr2,
            verbose, centering, mergeas, outerfunc, comb, num.sampling,
            num.perm, assayNames, decomp, sce)})

.cellCellDecomp <- function(userobjects, algorithm, ranks, rank, thr1,
    thr2, verbose, centering, mergeas, outerfunc, comb, num.sampling,
    num.perm, assayNames, decomp,  ...){
    # Import from sce object
    sce <- list(...)[[1]]
    # Import expression matrix
    input <- .importAssays(sce, assayNames)
    con = dbConnect(SQLite(), metadata(sce)$lrbase)
    LR <- dbGetQuery(con, "SELECT * FROM DATA")[, c("GENEID_L", "GENEID_R")]
    LR <- .uniqueLR(LR)
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
            mergeas, outerfunc, comb, num.sampling, num.perm, decomp,
            thr1, thr2, verbose)
    }

    # Data size
    if (algorithm == "ntd" || algorithm == "ntd2"){
        datasize <- c(ncol(res.sctensor[[2]]), ncol(res.sctensor[[3]]),
            ncol(res.sctensor[[4]]))
        recerror <- res.sctensor$recerror
        relchange <- res.sctensor$relchange
    }else{
        datasize <- NULL
        recerror <- NULL
        relchange <- NULL
    }

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
    author="The person who runs this script", assayNames="counts", thr=100,
    top="full", p=0.05, upper=20,
    goenrich=TRUE, meshenrich=TRUE, reactomeenrich=TRUE,
    doenrich=TRUE, ncgenrich=TRUE, dgnenrich=TRUE, nbins=40){
    standardGeneric("cellCellReport")})

setMethod("cellCellReport", signature(sce="SingleCellExperiment"),
    function(sce, reducedDimNames, out.dir, html.open, title, author,
        assayNames, thr, top, p, upper, goenrich, meshenrich,
        reactomeenrich, doenrich, ncgenrich, dgnenrich, nbins){
        .cellCellReport(reducedDimNames, out.dir,
            html.open, title, author, assayNames, thr, top, p, upper,
            goenrich, meshenrich, reactomeenrich,
            doenrich, ncgenrich, dgnenrich, nbins, sce)})

.cellCellReport <- function(reducedDimNames,
    out.dir=tempdir(), html.open=FALSE,
    title="The result of scTensor",
    author="The person who runs this script", assayNames="counts",
    thr=100, top="full", p=0.05, upper=20,
    goenrich=TRUE, meshenrich=TRUE, reactomeenrich=TRUE,
    doenrich=TRUE, ncgenrich=TRUE, dgnenrich=TRUE, nbins=40, ...){
    # Import from sce object
    sce <- list(...)[[1]]
    # algorithm-check
    if(metadata(sce)$algorithm != "ntd" && metadata(sce)$algorithm != "ntd2"){
        stop(paste0("cellCellReport can be performed by the result of",
            " cellCellDecomp in which the algorithm is ",
            "specified as 'ntd' or 'ntd2' for now."))
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
    if(metadata(sce)$algorithm == "ntd"){
        file.copy(
            from = system.file("extdata", "Workflow.png", package = "scTensor"),
            to = paste0(out.dir, "/Workflow.png"),
            overwrite = TRUE)
    }
    if(metadata(sce)$algorithm == "ntd2"){
        file.copy(
            from = system.file("extdata", "Workflow_2.png", package = "scTensor"),
            to = paste0(out.dir, "/Workflow_2.png"),
            overwrite = TRUE)
    }

    # HTML Report
    if(metadata(sce)$algorithm == "ntd"){
        .cellCellReport.Third(sce, thr, upper, assayNames, reducedDimNames, out.dir, author, title, p, top,
            goenrich, meshenrich, reactomeenrich,
            doenrich, ncgenrich, dgnenrich, nbins)
    }
    if(metadata(sce)$algorithm == "ntd2"){
        .cellCellReport.Third_2(sce, thr, upper, assayNames, reducedDimNames, out.dir, author, title, p, top,
            goenrich, meshenrich, reactomeenrich,
            doenrich, ncgenrich, dgnenrich, nbins)
    }
    # HTML Open
    message(paste0("################################################\n",
        "Data files are saved in\n",
        out.dir, "\n################################################\n"))
    if(html.open){
        browseURL(paste0(out.dir, "/index.html"))
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
            targetID <- target[which(score[target] == max(score[target]))[1]]
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
            nr, "*", nc, "),\n",
            "and only ", nrow(input), " of genes are kept."))
    }
    input
}
