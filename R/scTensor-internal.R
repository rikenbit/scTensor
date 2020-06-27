.frontal.normalization <- function(tnsr, total=1){
    original.dim <- dim(tnsr)
    denom <- apply(tnsr, 3, sum)
    dim(tnsr) <- c(prod(original.dim[1:2]), original.dim[3])
    tnsr2 <- vapply(seq(ncol(tnsr)), function(x, y){
        y$data[,x] / y$denom[x] * y$total
        }, y=list(data=tnsr, denom=denom, total=total),
        FUN.VALUE=tnsr[,1])
    tnsr2[which(is.nan(tnsr2))] <- 0
    dim(tnsr2) <- original.dim
    tnsr2
}

.uniqueLR <- function(LR){
  targetL <- grep("GENEID", LR$GENEID_L, invert = TRUE)
  targetR <- grep("GENEID", LR$GENEID_R, invert = TRUE)
  targetNotNAL <- which(!is.na(LR$GENEID_L))
  targetNotNAR <- which(!is.na(LR$GENEID_R))
  target <- intersect(targetL,
    intersect(targetR,
      intersect(targetNotNAL, targetNotNAR)))
  unique(LR[target, ])
}

.myvisNetwork <- function(g, col=NULL){
    # Edges
    edges = as.data.frame(get.edgelist(g), stringsAsFactors=FALSE)
    colnames(edges) <- c("from", "to")
    edges <- data.frame(edges,
        color=edge.attributes(g)$color,
        width=edge.attributes(g)$width)
    if(!is.null(col)){
        edges <- edges[which(edges$color == col), ]
    }
    # Nodes
    vattr <- get.vertex.attribute(g)
    uniqueid <- unique(c(edges$from, edges$to))
    color <- sapply(uniqueid, function(x){
        pos <- which(vattr$name == x)
        if(length(pos) == 1){
            vattr$color[pos]
        }else{
            "purple"
        }
    })
    size <- sapply(uniqueid, function(x){
        sum(vattr$size[which(vattr$name == x)])
    })
    nodes <- data.frame(
        id=uniqueid,
        type=TRUE,
        color=color,
        size=size*5,
        label=uniqueid
    )
    # Plot
    visNetwork(nodes, edges)
}

.cellCellReport.Third_2 <- function(sce, thr, upper, assayNames, reducedDimNames, out.dir, author, title, p, top,
    goenrich, meshenrich, reactomeenrich,
    doenrich, ncgenrich, dgnenrich, nbins){
    # Core Tensor
    index <- metadata(sce)$sctensor$index
    corevalue <- index[, "Value"]
    corevalue <- corevalue / sum(corevalue) * 100
    # Thresholding of the elements of core tensor
    selected <- which(cumsum(corevalue) <= thr)
    if(length(selected) > upper){
        selected <- seq_len(upper)
    }
    if(length(selected) == 0){
        stop(paste0("None of core tensor element is selected.\n",
        "Please specify the larger thr or perform cellCellDecomp\n",
        "with larger ranks such as c(# celltypes,# celltypes,1)."))
    }else{
        names(corevalue) <- c(rep("selected", length=length(selected)),
            rep("not selected",
            length=length(corevalue) - length(selected)))
        # Import expression matrix
        input <- .importAssays(sce, assayNames)
        # Low dimensional data
        twoD <- eval(parse(text=paste0("reducedDims(sce)$", reducedDimNames)))
        # Ligand-Receptor, PMID
        con = dbConnect(SQLite(), metadata(sce)$lrbase)
        LR <- dbGetQuery(con, "SELECT * FROM DATA")[,
            c("GENEID_L", "GENEID_R", "SOURCEID")]
        LR <- .uniqueLR(LR)
        dbDisconnect(con)
        # Species
        spc <- gsub(".eg.db.sqlite", "",
            strsplit(metadata(sce)$lrbase, "LRBase.")[[1]][3])
        # biomaRt Setting
        ah <- .annotationhub[[spc]]()
        # GeneName, Description, GO, Reactome, MeSH
        GeneInfo <- .geneinformation(sce, ah, spc, LR)
        # Cell Label
        celltypes <- metadata(sce)$color
        names(celltypes) <- metadata(sce)$label

        # Setting of schex
        sce <- make_hexbin(sce, nbins=nbins,
            dimension_reduction=reducedDimNames)
        # Plot Ligand/Receptor Genes
        suppressMessages(
            invisible(.genePlot(sce, assayNames, input, out.dir, GeneInfo, LR)))
        # Plot (Each <L,R,*>)
        out <- vapply(seq_along(selected), function(i){
            filenames <- paste0(out.dir,
                "/figures/CCIHypergraph_", index[i, 1],
                "_", index[i, 2], ".png")
            png(filename=filenames, width=2000, height=950)
            invisible(.CCIhyperGraphPlot_2(metadata(sce)$sctensor,
                twoDplot=twoD,
                label=celltypes,
                emph=index[i, seq_len(2)]))
            dev.off()
        }, 0L)
        # <L,R,>
        SelectedLR <- index[selected, seq_len(2)]

        # Setting for Parallel Computing
        message(paste0(length(selected),
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
        e$.eachRender_2 <- .eachRender_2
        e$.XYZ_HEADER1_2 <- .XYZ_HEADER1_2
        e$.XYZ_HEADER2_2 <- .XYZ_HEADER2_2
        e$.XYZ_HEADER3_2 <- .XYZ_HEADER3_2
        e$.XYZ_ENRICH_2 <- .XYZ_ENRICH_2
        e$algorithm <- metadata(sce)$algorithm
        e$goenrich <- goenrich
        e$meshenrich <- meshenrich
        e$reactomeenrich <- reactomeenrich
        e$doenrich <- doenrich
        e$ncgenrich <- ncgenrich
        e$dgnenrich <- dgnenrich

        # EachVec（Heavy...）
        if(is.vector(SelectedLR)){
            SelectedLR <- t(SelectedLR)
        }
        out.vecLR <- apply(SelectedLR, 1,
            function(x, e){.eachVecLR(x, e)}, e=e)
        names(out.vecLR) <- apply(SelectedLR, 1, function(x){
            paste(c("pattern", x[1], x[2]), collapse="_")
        })
        e$out.vecLR <- out.vecLR

        # Tagcloud
        invisible(.tagCloud_2(out.vecLR, out.dir))
        # Plot（CCI Hypergraph）
        par(ask=FALSE)
        png(filename=paste0(out.dir, "/figures/CCIHypergraph.png"),
            width=2000, height=950)
        invisible(.CCIhyperGraphPlot_2(metadata(sce)$sctensor,
            twoDplot=twoD, label=celltypes))
        dev.off()
        # Plot（Gene-wise Hypergraph）
        par(ask=FALSE)
        invisible(g <- .geneHyperGraphPlot_2(out.vecLR, GeneInfo, out.dir))

        # Rmd（ligand, selected）
        message("ligand.Rmd is created...")
        outLg <- file(paste0(out.dir, "/ligand.Rmd"), "w")
        writeLines(.LIGAND_HEADER, outLg, sep="\n")
        writeLines(.LIGAND_BODY_2(out.vecLR, GeneInfo, index, selected),
            outLg, sep="\n")
        close(outLg)
        # Rmd（receptor, selected）
        message("receptor.Rmd is created...")
        outRp <- file(paste0(out.dir, "/receptor.Rmd"), "w")
        writeLines(.RECEPTOR_HEADER, outRp, sep="\n")
        writeLines(.RECEPTOR_BODY_2(out.vecLR, GeneInfo, index, selected),
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
            selected, ClusterL, ClusterR, out.vecLR, g,
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
            " pattern_X_Y.Rmd files are compiled to pattern_X_Y.html :"))
        out <- vapply(selected,
            function(x, e, SelectedLR){
                .eachRender_2(x, e, SelectedLR)}, "", e=e, SelectedLR=SelectedLR)

        # Output index.html
        RMDFILES <- vapply(selected, function(x){
            paste0(paste(c("pattern", index[x, seq_len(2)]),
                collapse="_"), ".Rmd")
        }, "")
        message("index.Rmd is created...")
        outIdx <- file(paste0(out.dir, "/index.Rmd"), "w")
        writeLines(.MAINHEADER(author, title), outIdx, sep="\n")
        writeLines(.BODY1_2, outIdx, sep="\n")
        writeLines(.BODY2_2, outIdx, sep="\n")
        writeLines(.BODY3(numLPattern, ClusterL), outIdx, sep="\n")
        writeLines(.BODY4(numRPattern, ClusterR), outIdx, sep="\n")
        writeLines(.BODY6, outIdx, sep="\n")
        writeLines(.BODY7, outIdx, sep="\n")
        if(length(selected) != 0){
            writeLines(.BODY8_2(selected, RMDFILES, index, corevalue),
                outIdx, sep="\n")
        }
        writeLines(.BODY9, outIdx, sep="\n")
        writeLines(.BODY10, outIdx, sep="\n")
        close(outIdx)

        # Rendering
        message("index.Rmd is compiled to index.html...")
        render(paste0(out.dir, "/index.Rmd"), quiet=TRUE)
    }
}

.cellCellReport.Third <- function(sce, thr, upper, assayNames, reducedDimNames, out.dir, author, title, p, top,
    goenrich, meshenrich, reactomeenrich, doenrich,
    ncgenrich, dgnenrich, nbins){
    # Core Tensor
    index <- metadata(sce)$sctensor$index
    corevalue <- index[, "Value"]
    corevalue <- corevalue / sum(corevalue) * 100
    # Thresholding of the elements of core tensor
    selected <- which(cumsum(corevalue) <= thr)
    if(length(selected) > upper){
        selected <- seq_len(upper)
    }
    if(length(selected) == 0){
        stop(paste0("None of core tensor element is selected.\n",
        "Please specify the larger thr or perform cellCellDecomp\n",
        "with smaller ranks such as c(3,3,3)."))
    }else{
        names(corevalue) <- c(rep("selected", length=length(selected)),
            rep("not selected",
            length=length(corevalue) - length(selected)))
        # Import expression matrix
        input <- .importAssays(sce, assayNames)
        # Low dimensional data
        twoD <- eval(parse(text=paste0("reducedDims(sce)$", reducedDimNames)))
        # Ligand-Receptor, PMID
        con = dbConnect(SQLite(), metadata(sce)$lrbase)
        LR <- dbGetQuery(con, "SELECT * FROM DATA")[,
            c("GENEID_L", "GENEID_R", "SOURCEID")]
        LR <- .uniqueLR(LR)
        dbDisconnect(con)
        # Species
        spc <- gsub(".eg.db.sqlite", "",
            strsplit(metadata(sce)$lrbase, "LRBase.")[[1]][3])
        # biomaRt Setting
        ah <- .annotationhub[[spc]]()
        # GeneName, Description, GO, Reactome, MeSH
        GeneInfo <- .geneinformation(sce, ah, spc, LR)
        # Cell Label
        celltypes <- metadata(sce)$color
        names(celltypes) <- metadata(sce)$label

        # Setting of schex
        sce <- make_hexbin(sce, nbins=nbins,
            dimension_reduction=reducedDimNames)
        # Plot Ligand/Receptor Genes
        suppressMessages(
            invisible(.genePlot(sce, assayNames, input, out.dir, GeneInfo, LR)))
        # Plot (Each <L,R,*>)
        out <- vapply(seq_along(selected), function(i){
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
        e$algorithm <- metadata(sce)$algorithm
        e$goenrich <- goenrich
        e$meshenrich <- meshenrich
        e$reactomeenrich <- reactomeenrich
        e$doenrich <- doenrich
        e$ncgenrich <- ncgenrich
        e$dgnenrich <- dgnenrich

        # EachVec（Heavy...）
        out.vecLR <- vapply(SelectedLR,
            function(x, e){.eachVecLR(x, e)},
            FUN.VALUE=rep(list(0L), 9), e=e)
        colnames(out.vecLR) <- paste0("pattern", SelectedLR)
        e$out.vecLR <- out.vecLR

        # Tagcloud
        invisible(.tagCloud(out.vecLR, out.dir))
        # Plot（CCI Hypergraph）
        par(ask=FALSE)
        png(filename=paste0(out.dir, "/figures/CCIHypergraph.png"),
            width=2000, height=950)
        invisible(.CCIhyperGraphPlot(metadata(sce)$sctensor,
            twoDplot=twoD, label=celltypes))
        dev.off()
        # Plot（Gene-wise Hypergraph）
        par(ask=FALSE)
        invisible(g <- .geneHyperGraphPlot(out.vecLR, GeneInfo, out.dir))

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
            selected, ClusterL, ClusterR, out.vecLR, g,
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
        out <- vapply(selected,
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
    }
}

.importAssays <- function(sce, assayNames){
    if(assayNames %in% names(assays(sce))){
        eval(parse(text=paste0("assays(sce)$", assayNames)))
    }else{
        stop("Please specify the valid assayNames (cf. names(assays(sce)))")
    }
}

.CPMED <- function(input){
    libsize <- colSums(input)
    median(libsize) * t(t(input) / libsize)
}

.CPM <- function(input){
    libsize <- colSums(input)
    1e6 * t(t(input) / libsize)
}

.CPT <- function(input){
    libsize <- colSums(input)
    1e4 * t(t(input) / libsize)
}

.FTT <- function(input){
    sqrt(input) + sqrt(input + 1)
}

.cellType <- function(input, celltypes, mergeas="mean"){
    all.celltypes = unique(names(celltypes))
    if(mergeas == "mean"){
        ct <- vapply(all.celltypes, function(x){
            rowMeans(input[, which(names(celltypes) == x), drop=FALSE])
        }, rep(0.0, nrow(input)))
    }else if(mergeas == "sum"){
        ct <- vapply(all.celltypes, function(x){
            rowSums(input[, which(names(celltypes) == x), drop=FALSE])
        }, rep(0.0, nrow(input)))
    }else{
        stop("Please specify mergeas as mean or sum")
    }
    ct[which(is.nan(ct))] <- 0
    return(ct)
}

.divLR <- function(input, LR){
    input_L <- as.matrix(input[as.character(
        intersect(rownames(input), LR$GENEID_L)), ])
    input_R <- as.matrix(input[as.character(
        intersect(rownames(input), LR$GENEID_R)), ])
    return(list(input_L, input_R))
}

.cellCellReconst <- function(A1, A2, A3, index,
    type="matrix", names, collapse=FALSE){
    if(type == "matrix"){
        out <- lapply(seq_len(ncol(A1)), function(x){
            outer(A1[, x], A2[, x])})
    }else if(type == "tensor"){
        out <- list()
        counter <- 1
        for(i in seq_len(nrow(index))){
            element <- as.tensor(array(1, dim=c(1,1,1)))
            element@data[] <- index[i, "Value"]
            A <- as.matrix(A1[index[i, "Mode1"], ])
            B <- as.matrix(A2[index[i, "Mode2"], ])
            C <- as.matrix(A3[index[i, "Mode3"], ])
            rec_out <- recTensor(element, list(t(A), t(B), t(C)))
            if(collapse){
                out[[counter]] <- modeSum(rec_out, m=3, drop=TRUE)@data
                dimnames(out[[counter]]) <- list(names[[1]], names[[2]])
            }else{
                out[[counter]] <- rec_out
                dimnames(out[[counter]]@data) <- names
            }
            counter <- counter + 1
        }
    }
    out
}

.core2table <- function(core){
    d <- dim(core)
    out <- matrix(NA, nrow=prod(d), ncol=5)
    out[seq_len(nrow(out)), seq_len(3)] <- as.matrix(
        expand.grid(seq_len(d[1]), seq_len(d[2]), seq_len(d[3])))
    colnames(out) <- c("Mode1", "Mode2", "Mode3", "Value", "Rank")
    veccore <- c()
    counter <- 1
    vapply(seq_len(nrow(out)), function(x){
        i <- out[x, 1]
        j <- out[x, 2]
        k <- out[x, 3]
        out[x, 4] <<- core[i,j,k]
    }, 0.0)
    out[, 5] <- nrow(out) - base::rank(out[, 4]) + 1
    out[order(out[,4], decreasing=TRUE), ]
}

.core2table_2 <- function(core){
    d <- dim(core)
    out <- matrix(NA, nrow=d[1]*d[2], ncol=4)
    out[seq_len(nrow(out)), seq_len(2)] <- as.matrix(
        expand.grid(seq_len(d[1]), seq_len(d[2])))
    colnames(out) <- c("Mode1", "Mode2", "Value", "Rank")
    veccore <- c()
    counter <- 1
    vapply(seq_len(nrow(out)), function(x){
        i <- out[x, 1]
        j <- out[x, 2]
        out[x, 3] <<- core[i,j]
    }, 0.0)
    out[, 4] <- nrow(out) - base::rank(out[, 3]) + 1
    out[order(out[,3], decreasing=TRUE), ]
}

.fastPossibleCombination <- function(input, LR, celltypes, num.sampling, outerfunc){
    lr <- .divLR(input, LR)
    L <- as.matrix(lr[[1]])
    R <- as.matrix(lr[[2]])
    tnsr <- array(0, dim=c(length(unique(names(celltypes))),
        length(unique(names(celltypes))), 0))
    Pair.name <- c()
    comb <- as.matrix(expand.grid(unique(names(celltypes)),
        unique(names(celltypes))))
    for(i in seq_len(nrow(LR))){
        l <- which(LR$GENEID_L[i] == rownames(L))
        r <- which(LR$GENEID_R[i] == rownames(R))
        if(length(l)==1 && length(r) == 1){
            pre_tnsr <- array(0, dim=c(length(unique(names(celltypes))),
                length(unique(names(celltypes))), 1))
            for(j in seq_len(num.sampling)){
                target <- vapply(unique(names(celltypes)), function(x){
                    sample(which(names(celltypes) == x), 1)}, 0L)
                pre_tnsr[,,1] <- pre_tnsr[,,1] + outer(as.vector(L[l, target]),
                as.vector(R[r, target]), outerfunc)
            }
            tnsr <- abind(tnsr, pre_tnsr, along=3)
            Pair.name <- c(Pair.name, paste(LR$GENEID_L[i],
                LR$GENEID_R[i], sep="_"))
        }
    }
    tnsr <- .frontal.normalization(tnsr)
    dimnames(tnsr) <- list(unique(names(celltypes)),
        unique(names(celltypes)), Pair.name)
    list(tnsr=tnsr, pairname=Pair.name)
}

.slowPossibleCombination <- function(input, LR, celltypes, outerfunc){
    lr <- .divLR(input, LR)
    L <- as.matrix(lr[[1]])
    R <- as.matrix(lr[[2]])
    tnsr <- array(0, dim=c(length(unique(names(celltypes))),
        length(unique(names(celltypes))), 0))
    Pair.name <- c()
    comb <- as.matrix(expand.grid(unique(names(celltypes)),
        unique(names(celltypes))))
    for(i in seq_len(nrow(LR))){
        l <- which(LR$GENEID_L[i] == rownames(L))
        r <- which(LR$GENEID_R[i] == rownames(R))
        if(length(l)==1 && length(r) == 1){
            pre_tnsr <- outer(as.vector(L[l,]), as.vector(R[r,]))
            pre_tnsr <- apply(comb, 1, function(x){
                sum(pre_tnsr[which(names(celltypes) == x[1]),
                which(names(celltypes) == x[2])])})
            dim(pre_tnsr) <- c(length(unique(celltypes)),
            length(unique(celltypes)))
            tnsr <- abind(tnsr, pre_tnsr, along=3)
            Pair.name <- c(Pair.name, paste(LR$GENEID_L[i],
                LR$GENEID_R[i], sep="_"))
        }
    }
    tnsr <- .frontal.normalization(tnsr)
    dimnames(tnsr) <- list(unique(names(celltypes)),
        unique(names(celltypes)), Pair.name)
    list(tnsr=tnsr, pairname=Pair.name)
}

.cellCellDecomp.Third_2 <- function(input, LR, celltypes, ranks, rank, centering,
        mergeas, outerfunc, comb, num.sampling, num.perm, decomp, thr1, thr2, verbose){
    # ranks-check
    max.rank <- length(unique(celltypes))
    if(ranks[1] > max.rank || ranks[2] > max.rank){
        stop("ranks must be defined less than number of celltypes")
    }
    if(centering){
        fout <- .celltypemergedtensor(input, LR, celltypes, mergeas, outerfunc)
    }else{
        if(comb == "random"){
            fout <- .fastPossibleCombination(input, LR, celltypes,
                num.sampling, outerfunc)
        }else if(comb == "all"){
            fout <- .slowPossibleCombination(input, LR, celltypes, outerfunc)
        }
    }
    # rank check
    check1 <- dim(fout$tnsr)[2] * dim(fout$tnsr)[3] < ranks[1]
    check2 <- dim(fout$tnsr)[3] * dim(fout$tnsr)[1] < ranks[2]
    if(check1){
        stop("Please specify the ranks[1] as an smaller value")
    }
    if(check2){
        stop("Please specify the ranks[2] as an smaller value")
    }
    ranks <- ranks[1:2]
    if(dim(fout$tnsr)[1] >= 100 || dim(fout$tnsr)[2] >= 100){
        stop(paste0("Extreamly Large tensor will be generated!\n",
            "This problem will be solved by the next version of scTensor."))
    }
    tnsr <- as.tensor(fout$tnsr)
    Pair.name <- fout$pairname
    if(decomp){
        message(paste0(paste(dim(tnsr), collapse=" * "), " Tensor is created"))
        out <- try(NTD(X=tnsr, rank=ranks, modes=1:2, num.iter=30,
            verbose=verbose, algorithm="Frobenius"))
        if(is(out)[1] == "try-error"){
            out <- NTD(X=tnsr, rank=ranks, modes=1:2, num.iter=30, algorithm="Frobenius")
        }
        A1 <- out$A[[1]]
        A2 <- out$A[[2]]
        rownames(A1) <- paste0("Dim", seq_len(ranks[1]))
        colnames(A1) <- unique(names(celltypes))
        rownames(A2) <- paste0("Dim", seq_len(ranks[2]))
        colnames(A2) <- unique(names(celltypes))
        dimnames(out$S@data) <- list(
            Lpattern=paste0("Dim", seq_len(ranks[1])),
            Rpattern=paste0("Dim", seq_len(ranks[2])),
            LR=Pair.name
        )
        if (nrow(A1) * nrow(A2) == 1){
            core <- t(sum(out$S@data))
        }else{
            core <- modeSum(out$S, m=3, drop=TRUE)@data
            if(ranks[1] == 1 || ranks[2] == 1){
                dim(core) <- c(ranks[1], ranks[2])
            }
        }
        index <- .core2table_2(core)
        if(is.vector(index)){
            index <- t(index)
        }
        return(list(index=index, ligand=A1,
            receptor=A2, lrpair=out$S,
            recerror=out$RecError, relchange=out$RelChange))
    }else{
        tnsr_cc <- modeSum(tnsr, m=3, drop=TRUE)@data
        dimnames(tnsr_cc) <- list(unique(names(celltypes)),
            unique(names(celltypes)))
        return(list(cellcellpattern=tnsr_cc, cellcelllrpairpattern=tnsr))
    }
}

.cellCellDecomp.Third <- function(input, LR, celltypes, ranks, rank, centering,
        mergeas, outerfunc, comb, num.sampling, num.perm, decomp, thr1, thr2, verbose){
    # ranks-check
    max.rank <- length(unique(celltypes))
    if(ranks[1] > max.rank || ranks[2] > max.rank){
        stop("ranks must be defined less than number of celltypes")
    }
    if(centering){
        fout <- .celltypemergedtensor(input, LR, celltypes, mergeas, outerfunc)
    }else{
        if(comb == "random"){
            fout <- .fastPossibleCombination(input, LR, celltypes,
                num.sampling, outerfunc)
        }else if(comb == "all"){
            fout <- .slowPossibleCombination(input, LR, celltypes, outerfunc)
        }
    }
    # rank check
    check1 <- dim(fout$tnsr)[2] * dim(fout$tnsr)[3] < ranks[1]
    check2 <- dim(fout$tnsr)[3] * dim(fout$tnsr)[1] < ranks[2]
    check3 <- dim(fout$tnsr)[1] * dim(fout$tnsr)[2] < ranks[3]
    if(check1){
        stop("Please specify the ranks[1] as an smaller value")
    }
    if(check2){
        stop("Please specify the ranks[2] as an smaller value")
    }
    if(check3){
        stop("Please specify the ranks[3] as an smaller value")
    }
    if(dim(fout$tnsr)[1] >= 100 || dim(fout$tnsr)[2] >= 100){
        stop(paste0("Extreamly Large tensor will be generated!\n",
            "This problem will be solved by the next version of scTensor."))
    }
    tnsr <- as.tensor(fout$tnsr)
    Pair.name <- fout$pairname
    if(decomp){
        message(paste0(paste(dim(tnsr), collapse=" * "), " Tensor is created"))
        out <- try(NTD(X=tnsr, rank=ranks, num.iter=30,
            verbose=verbose, algorithm="Frobenius"))
        if(is(out)[1] == "try-error"){
            out <- try(NTD(X=tnsr, rank=ranks, num.iter=30,
                verbose=verbose, algorithm="Frobenius"))
        }
        A1 <- out$A[[1]]
        A2 <- out$A[[2]]
        A3 <- out$A[[3]]
        rownames(A1) <- paste0("Dim", seq_len(ranks[1]))
        colnames(A1) <- unique(names(celltypes))
        rownames(A2) <- paste0("Dim", seq_len(ranks[2]))
        colnames(A2) <- unique(names(celltypes))
        rownames(A3) <- paste0("Dim", seq_len(ranks[3]))
        colnames(A3) <- Pair.name
        core <- out$S@data
        index <- .core2table(core)
        if(is.vector(index)){
            index <- t(index)
        }
        lrpairout <- .cellCellReconst(A1, A2, A3, index, type="tensor",
            names=list(colnames(A1), colnames(A2), colnames(A3)),
            collapse=TRUE)
        return(list(index=index, ligand=A1,
            receptor=A2, lrpair=A3, cellcellpattern=lrpairout,
            recerror=out$RecError, relchange=out$RelChange))
    }else{
        tnsr_cc <- modeSum(tnsr, m=3, drop=TRUE)@data
        dimnames(tnsr_cc) <- list(unique(names(celltypes)),
            unique(names(celltypes)))
        return(list(cellcellpattern=tnsr_cc, cellcelllrpairpattern=tnsr))
    }
}

.cellCellDecomp.Second <- function(input, LR, celltypes, ranks, rank,
    centering, mergeas, outerfunc, comb, num.sampling, num.perm,
    decomp, thr1, thr2, verbose){
    # ranks-check
    max.rank <- length(unique(celltypes))
    if(rank > max.rank){
        stop("ranks must be defined less than number of celltypes")
    }
    if(centering){
        fout <- .celltypemergedtensor(input, LR, celltypes, mergeas, outerfunc)
    }else{
        if(comb == "random"){
            fout <- .fastPossibleCombination(input, LR, celltypes,
                num.sampling, outerfunc)
        }else if(comb == "all"){
            fout <- .slowPossibleCombination(input, LR, celltypes, outerfunc)
        }
    }
    tnsr <- as.tensor(fout$tnsr)
    Pair.name <- fout$pairname
    mat.tnsr <- modeSum(tnsr, m=3, drop=TRUE)@data
    dimnames(mat.tnsr) <- list(unique(names(celltypes)),
        unique(names(celltypes)))
    if(decomp){
        out <- try(NMF(mat.tnsr, J=rank, num.iter=30, algorithm="Frobenius"))
        colnames(out$U) <- paste0("Dim", seq_len(rank))
        colnames(out$V) <- paste0("Dim", seq_len(rank))
        sizeU <- apply(out$U, 2, function(x){sqrt(sum(x^2))})
        sizeV <- apply(out$V, 2, function(x){sqrt(sum(x^2))})
        value <- sizeU * sizeV
        index <- rank - base::rank(value) + 1
        return(list(ligand=out$U[,index], receptor=out$V[,index],
        cellcellpattern=.cellCellReconst(out$U[,index], out$V[,index])))
    }else{
        return(list(cellcellpattern=mat.tnsr))
    }
}

.cellCellDecomp.Pearson <- function(input, LR, celltypes, ranks, rank,
    centering, mergeas, outerfunc, comb, num.sampling, num.perm,
    decomp, thr1, thr2, verbose){
    out <- .cellType(input, celltypes)
    out <- cor(out, method="pearson")
    out <- as.matrix(out)
    return(list(cellcellpattern=out))
}

.cellCellDecomp.Spearman <- function(input, LR, celltypes, ranks, rank,
    centering, mergeas, outerfunc, comb, num.sampling, num.perm,
    decomp, thr1, thr2, verbose){
    out <- .cellType(input, celltypes)
    out <- cor(out, method="spearman")
    out <- as.matrix(out)
    return(list(cellcellpattern=out))
}

.cellCellDecomp.Distance <- function(input, LR, celltypes, ranks, rank,
    centering, mergeas, outerfunc, comb, num.sampling, num.perm,
    decomp, thr1, thr2, verbose){
    out <- .cellType(input, celltypes)
    out <- 1 / dist(t(out))
    out <- as.matrix(out)
    diag(out) <- 1
    return(list(cellcellpattern=out))
}

.cellCellDecomp.Pearson.LR <- function(input, LR, celltypes, ranks, rank,
    centering, mergeas, outerfunc, comb, num.sampling, num.perm,
    decomp, thr1, thr2, verbose){
    ct <- .cellType(input, celltypes)
    ct <- .divLR(ct, LR)
    L <- ct[[1]]
    R <- ct[[2]]
    Adj <- matrix(0, nrow=nrow(L), ncol=nrow(R))
    rownames(Adj) <- rownames(L)
    colnames(Adj) <- rownames(R)
    for(i in seq_len(nrow(L))){
        l <- as.character(rownames(L)[i])
        r <- as.character(LR[which(LR$GENEID_L == l), "GENEID_R"])
        Adj[l, intersect(colnames(Adj), r)] <- 1
    }
    out <- t(L) %*% Adj %*% R
    out <- cor(out, method="pearson")
    out <- as.matrix(out)
    return(list(cellcellpattern=out))
}

.cellCellDecomp.Spearman.LR <- function(input, LR, celltypes, ranks, rank,
    centering, mergeas, outerfunc, comb, num.sampling, num.perm,
    decomp, thr1, thr2, verbose){
    ct <- .cellType(input, celltypes)
    ct <- .divLR(ct, LR)
    L <- ct[[1]]
    R <- ct[[2]]
    Adj <- matrix(0, nrow=nrow(L), ncol=nrow(R))
    rownames(Adj) <- rownames(L)
    colnames(Adj) <- rownames(R)
    for(i in seq_len(nrow(L))){
        l <- as.character(rownames(L)[i])
        r <- as.character(LR[which(LR$GENEID_L == l), "GENEID_R"])
        Adj[l, intersect(colnames(Adj), r)] <- 1
    }
    out <- t(L) %*% Adj %*% R
    out <- cor(out, method="spearman")
    out <- as.matrix(out)
    return(list(cellcellpattern=out))
}

.cellCellDecomp.Distance.LR <- function(input, LR, celltypes, ranks, rank,
    centering, mergeas, outerfunc, comb, num.sampling, num.perm,
    decomp, thr1, thr2, verbose){
    ct <- .cellType(input, celltypes)
    ct <- .divLR(ct, LR)
    L <- ct[[1]]
    R <- ct[[2]]
    Adj <- matrix(0, nrow=nrow(L), ncol=nrow(R))
    rownames(Adj) <- rownames(L)
    colnames(Adj) <- rownames(R)
    for(i in seq_len(nrow(L))){
        l <- as.character(rownames(L)[i])
        r <- as.character(LR[which(LR$GENEID_L == l), "GENEID_R"])
        Adj[l, intersect(colnames(Adj), r)] <- 1
    }
    out <- t(L) %*% Adj %*% R
    out <- 1 / dist(t(out))
    out <- as.matrix(out)
    diag(out) <- 1
    return(list(cellcellpattern=out))
}

.cellCellDecomp.PossibleCombination <- function(input, LR, celltypes, ranks,
    rank, centering, mergeas, outerfunc, comb, num.sampling, num.perm,
    decomp, thr1, thr2, verbose){
    if(thr1 <= 0 || thr2 <= 0){
        warning("None of cell-cell interaction will be detected.")
    }
    lr <- .divLR(input, LR)
    L <- as.matrix(lr[[1]])
    R <- as.matrix(lr[[2]])

    L[which(L < thr1)] <- 0
    R[which(R < thr1)] <- 0
    L[which(L >= thr1)] <- 1
    R[which(R >= thr1)] <- 1

    tnsr_bin <- array(0, dim=c(length(unique(names(celltypes))),
        length(unique(names(celltypes))), 0))
    Pair.name <- c()
    comb <- as.matrix(expand.grid(unique(names(celltypes)),
        unique(names(celltypes))))
    for(i in seq_len(nrow(LR))){
        l <- which(LR$GENEID_L[i] == rownames(L))
        r <- which(LR$GENEID_R[i] == rownames(R))
        if(length(l)==1 && length(r) == 1){
            pre_tnsr_bin <- outer(as.vector(L[l,]), as.vector(R[r,]))
            pre_tnsr_bin <- apply(comb, 1, function(x){
                sum(pre_tnsr_bin[which(names(celltypes) == x[1]),
                    which(names(celltypes) == x[2])])
                })
            dim(pre_tnsr_bin) <- c(length(unique(celltypes)),
                length(unique(celltypes)))
            pre_tnsr_bin[which(pre_tnsr_bin < thr2)] <- 0
            pre_tnsr_bin[which(pre_tnsr_bin >= thr2)] <- 1
            tnsr_bin <- abind(tnsr_bin, pre_tnsr_bin, along=3)
            Pair.name <- c(Pair.name,
                paste(LR$GENEID_L[i], LR$GENEID_R[i], sep="_"))
        }
    }
    dimnames(tnsr_bin) <- list(unique(names(celltypes)),
        unique(names(celltypes)), Pair.name)
    tnsr_cc <- modeSum(as.tensor(tnsr_bin), m=3, drop=TRUE)@data
    dimnames(tnsr_cc) <- list(unique(names(celltypes)),
        unique(names(celltypes)))
    return(list(cellcellpattern=tnsr_cc, cellcelllrpairpattern=tnsr_bin))
}


.celltypemergedtensor <- function(input, LR, celltypes, mergeas, outerfunc){
    ct <- .cellType(input, celltypes, mergeas)
    ct <- .divLR(ct, LR)
    L <- ct[[1]]
    R <- ct[[2]]
    tnsr <- array(0, dim=c(ncol(L), ncol(R), 0))
    Pair.name <- c()
    for(i in seq_len(nrow(LR))){
        l <- which(LR$GENEID_L[i] == rownames(L))
        r <- which(LR$GENEID_R[i] == rownames(R))
        if(length(l)==1 && length(r) == 1){
            tnsr <- abind(tnsr,
                outer(as.vector(L[l,]), as.vector(R[r,]), outerfunc), along=3)
            Pair.name <- c(Pair.name,
                paste(LR$GENEID_L[i], LR$GENEID_R[i], sep="_"))
        }
    }
    tnsr <- .frontal.normalization(tnsr)
    list(tnsr=tnsr, pairname=Pair.name)
}

.celltypemergedtensor_2 <- function(input, LR, celltypes, mergeas, outerfunc){
    ct <- .cellType(input, celltypes, mergeas)
    ct <- .divLR(ct, LR)
    L <- ct[[1]]
    R <- ct[[2]]
    outer(as.vector(L), as.vector(R), outerfunc)
}

.labelpermutation <- function(observed, input, l, r, num.perm,
    celltypes, mergeas, outerfunc){
    perm <- array(0, dim=c(dim(observed), 0))
    for(j in seq_len(num.perm)){
        tmp <- .celltypemergedtensor_2(input,
            data.frame(GENEID_L=l, GENEID_R=r, stringsAsFactors=FALSE),
            sample(celltypes), mergeas, outerfunc)
        tmp <- observed - tmp
        tmp[which(tmp >= 0)] <- 0
        tmp[which(tmp < 0)] <- 1
        perm <- abind(perm, tmp)
    }
    modeSum(as.tensor(perm), m=3, drop=TRUE)@data / num.perm
}

.cellCellDecomp.LabelPerm.LR <- function(input, LR, celltypes, ranks, rank,
    centering, mergeas, outerfunc, comb, num.sampling, num.perm,
    decomp, thr1, thr2, verbose){
    fout <- .celltypemergedtensor(input, LR, celltypes,
        mergeas, outerfunc)
    tnsr <- as.tensor(fout$tnsr)
    Pair.name <- fout$pairname
    pval <- array(0, dim=c(dim(tnsr)[1:2], 0))
    for(i in seq_along(Pair.name)){
        cat(paste0(i, " / ", length(Pair.name), "\r"))
        l <- strsplit(Pair.name[i], "_")[[1]][1]
        r <- strsplit(Pair.name[i], "_")[[1]][2]
        tmp <- .labelpermutation(tnsr[,,i]@data, input,
            l, r, num.perm, celltypes, mergeas, outerfunc)
        pval <- abind(pval, tmp)
    }
    pval[which(is.nan(pval))] <- 1
    dimnames(pval) <- list(unique(names(celltypes)),
        unique(names(celltypes)), Pair.name)
    tnsr_cc <- modeSum(tnsr, m=3, drop=TRUE)@data
    dimnames(tnsr_cc) <- list(unique(names(celltypes)),
        unique(names(celltypes)))
    return(list(cellcellpattern=tnsr_cc, cellcelllrpairpattern=tnsr,
        pval=pval))
}

# SingleCellSignalR
.celltypemergedtensor.ca <- function(input, LR, celltypes, mergeas, outerfunc){
    ct <- .cellType(input, celltypes, mergeas)
    ct <- .divLR(ct, LR)
    L <- ct[[1]]
    R <- ct[[2]]
    tnsr <- array(0, dim=c(ncol(L), ncol(R), 0))
    Pair.name <- c()
    for(i in seq_len(nrow(LR))){
        l <- which(LR$GENEID_L[i] == rownames(L))
        r <- which(LR$GENEID_R[i] == rownames(R))
        if(length(l)==1 && length(r) == 1){
            tnsr <- abind(tnsr,
                outer(as.vector(L[l,]), as.vector(R[r,]), outerfunc), along=3)
            Pair.name <- c(Pair.name,
                paste(LR$GENEID_L[i], LR$GENEID_R[i], sep="_"))
        }
    }
    tnsr[which(is.nan(tnsr))] <- 0
    tnsr <- sqrt(tnsr) / (mean(rbind(L, R)) + sqrt(tnsr))
    list(tnsr=tnsr, pairname=Pair.name)
}

.celltypemergedtensor.ca_2 <- function(input, LR, meanC, celltypes,
    mergeas, outerfunc){
    ct <- .cellType(input, celltypes, mergeas)
    ct <- .divLR(ct, LR)
    L <- ct[[1]]
    R <- ct[[2]]
    tmp <- outer(as.vector(L), as.vector(R), outerfunc)
    sqrt(tmp) / (meanC + sqrt(tmp))
}

.meanC <- function(input, celltypes, mergeas, LR){
    ct <- .cellType(input, celltypes, mergeas)
    ct <- .divLR(ct, LR)
    L <- ct[[1]]
    R <- ct[[2]]
    mean(rbind(L, R))
}

.labelpermutation.ca <- function(observed, input, LR, l, r, num.perm,
    celltypes, mergeas, outerfunc){
    meanC <- .meanC(input, celltypes, mergeas, LR)
    perm <- array(0, dim=c(dim(observed), 0))
    for(j in seq_len(num.perm)){
        tmp <- .celltypemergedtensor.ca_2(input,
            data.frame(GENEID_L=l, GENEID_R=r, stringsAsFactors=FALSE),
            meanC, sample(celltypes), mergeas, outerfunc)
        tmp <- observed - tmp
        tmp[which(tmp >= 0)] <- 0
        tmp[which(tmp < 0)] <- 1
        perm <- abind(perm, tmp)
    }
    modeSum(as.tensor(perm), m=3, drop=TRUE)@data / num.perm
}

.cellCellDecomp.CabelloAguilar <- function(input, LR, celltypes, ranks,
    rank, centering, mergeas, outerfunc, comb, num.sampling, num.perm,
    decomp, thr1, thr2, verbose){
    fout <- .celltypemergedtensor.ca(input, LR, celltypes,
        mergeas, outerfunc)
    tnsr <- as.tensor(fout$tnsr)
    Pair.name <- fout$pairname
    pval <- array(0, dim=c(dim(tnsr)[1:2], 0))
    for(i in seq_along(Pair.name)){
        cat(paste0(i, " / ", length(Pair.name), "\r"))
        l <- strsplit(Pair.name[i], "_")[[1]][1]
        r <- strsplit(Pair.name[i], "_")[[1]][2]
        tmp <- .labelpermutation.ca(tnsr[,,i]@data, input,
            LR, l, r, num.perm, celltypes, mergeas, outerfunc)
        pval <- abind(pval, tmp)
    }
    pval[which(is.nan(pval))] <- 1
    dimnames(pval) <- list(unique(names(celltypes)),
        unique(names(celltypes)), Pair.name)
    tnsr_cc <- modeSum(tnsr, m=3, drop=TRUE)@data
    dimnames(tnsr_cc) <- list(unique(names(celltypes)),
        unique(names(celltypes)))
    return(list(cellcellpattern=tnsr_cc, cellcelllrpairpattern=tnsr,
        pval=pval))
}

# Halpern et al.（LECs）
.Zscaling <- function(x){
    (x - mean(x)) / sd(x)
}

.celltypemergedtensor.hl <- function(input, LR, celltypes, mergeas, outerfunc){
    ct <- .cellType(input, celltypes, mergeas)
    ct <- .divLR(ct, LR)
    L <- ct[[1]]
    R <- ct[[2]]
    # Z-scaling
    Z_L <- t(apply(L, 1, .Zscaling))
    Z_R <- t(apply(R, 1, .Zscaling))
    tnsr <- array(0, dim=c(ncol(L), ncol(R), 0))
    Pair.name <- c()
    for(i in seq_len(nrow(LR))){
        l <- which(LR$GENEID_L[i] == rownames(L))
        r <- which(LR$GENEID_R[i] == rownames(R))
        if(length(l)==1 && length(r) == 1){
            tnsr <- abind(tnsr,
                sqrt(outer(as.vector(Z_L[l,])^2, as.vector(Z_R[r,])^2, outerfunc)), along=3)
            Pair.name <- c(Pair.name,
                paste(LR$GENEID_L[i], LR$GENEID_R[i], sep="_"))
        }
    }
    tnsr[which(is.nan(tnsr))] <- 0
    list(tnsr=tnsr, pairname=Pair.name)
}

.celltypemergedtensor.hl_2 <- function(input, LR, celltypes, mergeas, outerfunc){
    ct <- .cellType(input, celltypes, mergeas)
    ct <- .divLR(ct, LR)
    Z_L <- .Zscaling(ct[[1]])
    Z_R <- .Zscaling(ct[[2]])
    sqrt(outer(as.vector(Z_L)^2, as.vector(Z_R)^2, outerfunc))
}

.labelpermutation.hl <- function(observed, input, l, r, num.perm,
    celltypes, mergeas, outerfunc){
    perm <- array(0, dim=c(dim(observed), 0))
    for(j in seq_len(num.perm)){
        tmp <- .celltypemergedtensor.hl_2(input,
            data.frame(GENEID_L=l, GENEID_R=r, stringsAsFactors=FALSE),
            sample(celltypes), mergeas, outerfunc)
        tmp <- observed - tmp
        tmp[which(tmp >= 0)] <- 0
        tmp[which(tmp < 0)] <- 1
        perm <- abind(perm, tmp)
    }
    modeSum(as.tensor(perm), m=3, drop=TRUE)@data / num.perm
}

.cellCellDecomp.Halpern <- function(input, LR, celltypes, ranks, rank,
    centering, mergeas, outerfunc, comb, num.sampling, num.perm,
    decomp, thr1, thr2, verbose){
    fout <- .celltypemergedtensor.hl(input, LR, celltypes,
        mergeas, outerfunc)
    fout$tnsr[which(is.nan(fout$tnsr))] <- 0
    tnsr <- as.tensor(fout$tnsr)
    Pair.name <- fout$pairname
    pval <- array(0, dim=c(dim(tnsr)[1:2], 0))
    for(i in seq_along(Pair.name)){
        cat(paste0(i, " / ", length(Pair.name), "\r"))
        l <- strsplit(Pair.name[i], "_")[[1]][1]
        r <- strsplit(Pair.name[i], "_")[[1]][2]
        tmp <- .labelpermutation.hl(tnsr[,,i]@data, input,
            l, r, num.perm, celltypes, mergeas, outerfunc)
        pval <- abind(pval, tmp)
    }
    pval[which(is.nan(pval))] <- 1
    dimnames(pval) <- list(unique(names(celltypes)),
        unique(names(celltypes)), Pair.name)
    tnsr_cc <- modeSum(tnsr, m=3, drop=TRUE)@data
    dimnames(tnsr_cc) <- list(unique(names(celltypes)),
        unique(names(celltypes)))
    return(list(cellcellpattern=tnsr_cc, cellcelllrpairpattern=tnsr,
        pval=pval))
}

.flist <- list(
    # 3-Order, 2-mode decomposition = NTD2
    "ntd2" = .cellCellDecomp.Third_2,
    # 3-Order = NTD
    "ntd" = .cellCellDecomp.Third,
    # 2-Order = NMF
    "nmf" = .cellCellDecomp.Second,
    # Other method (Possible Combination)
    "pcomb" = .cellCellDecomp.PossibleCombination,
    # Other method (Pearson's Correlation Coefficient without LR pair)
    "pearson" = .cellCellDecomp.Pearson,
    # Other method (Spearman's Correlation Coefficient without LR pair)
    "spearman" = .cellCellDecomp.Spearman,
    # Other method (Euclidian Distance without LR pair)
    "distance" = .cellCellDecomp.Distance,
    # Other method (Pearson's Correlation Coefficient with LR pair)
    "pearson.lr" = .cellCellDecomp.Pearson.LR,
    # Other method (Spearman's Correlation Coefficient with LR pair)
    "spearman.lr" = .cellCellDecomp.Spearman.LR,
    # Other method (Euclidian Distance with LR pair)
    "distance.lr" = .cellCellDecomp.Distance.LR,
    # Other method (Celltype Label Permutation with LR pair)
    "label.permutation" = .cellCellDecomp.LabelPerm.LR,
    # Other method (SingleCellSignalR)
    "cabello.aguilar" = .cellCellDecomp.CabelloAguilar,
    # Other method（LECs）
    "halpern" = .cellCellDecomp.Halpern
)

.CCIhyperGraphPlot <- function(outobj, twoDplot=NULL, vertex.size=18,
    xleft=1.75, ybottom=-0.5, xright=1.85, ytop=0.5, label="", emph=NULL, algorithm=""){
    # Number of Patterns
    numLPattern <- nrow(outobj$ligand)
    numRPattern <- nrow(outobj$receptor)

    #
    # Step.1 : Background Network
    #
    edgewd_L <- as.vector(vapply(seq_len(numLPattern), function(x){
        rep(x, numRPattern)
        }, rep(0L, numRPattern)))
    edgewd_R <- rep(seq_len(numRPattern), numLPattern)
    edgewd_Strength <- vapply(
        seq_len(numLPattern*numRPattern), function(x){
            targetL <- which(
                outobj$index[, "Mode1"] == edgewd_L[x])
            targetR <- which(
                outobj$index[, "Mode2"] == edgewd_R[x])
            sum(outobj$index[intersect(targetL, targetR), 4])
        }, 0.0)
    edgewd <- cbind(edgewd_L, edgewd_R, edgewd_Strength)
    colnames(edgewd) <- c("L", "R", "Strength")

    # Node name (Top<Ligand> and Bottom<Receptor>)
    nodesSetTop <- paste0("L", seq_len(numLPattern))
    nodesSetBottom <- paste0("R", seq_len(numRPattern))

    # Empty Graph
    g <- graph.empty()

    # Add nodes
    g <- add.vertices(g, nv=length(nodesSetTop),
        attr=list(name=nodesSetTop,
            type=rep(TRUE, length(nodesSetTop))))
    g <- add.vertices(g, nv=length(nodesSetBottom),
        attr=list(name=nodesSetBottom,
            type=rep(TRUE, length(nodesSetBottom))))

    # Add edges
    edgeListVec <- as.vector(t(as.matrix(
        data.frame(
            S1=paste0('L', edgewd[,1]),
            S2=paste0('R', edgewd[,2])
    ))))
    g <- add.edges(g, edgeListVec)

    # Edge weghts
    E(g)$weight <- edgewd[,3]

    # Edge color
    weight <- E(g)$weight
    E(g)$weight <- weight / max(weight) * 20
    mycolor <- smoothPalette(E(g)$weight,
        palfunc=colorRampPalette(.setColor("greens"), alpha=TRUE))
    if(!is.null(emph)){
        target <- intersect(
            which(get.edgelist(g)[, 1] == paste0("L", emph[1])),
            which(get.edgelist(g)[, 2] == paste0("R", emph[2])))
        mycolor[target] <- rgb(1,0,0,0.5)
    }

    # Layout
    x <- c(seq_along(nodesSetTop), seq_along(nodesSetBottom))
    y <- c(rep(1, length=length(nodesSetTop)),
        rep(0, length=length(nodesSetBottom)))
    mylayout <- cbind(x, y)

    # Network Plot
    par(ask=FALSE)
    par(oma=c(2,2,2,2))
    plot.igraph(g,
        layout=mylayout,
        vertex.size=18,
        vertex.label="",
        vertex.color="white",
        vertex.shape="square",
        edge.color=mycolor,
        vertex.frame.color="gray",
        edge.width=E(g)$weight)

    # Gradient
    par(ask=FALSE)
    gradient.rect(xleft, ybottom, xright, ytop,
        col=smoothPalette(sort(weight),
        palfunc=colorRampPalette(.setColor("greens"), alpha=TRUE)),
        gradient="y")
    text(2.2, ybottom+(ytop-ybottom)*0/4, round(quantile(weight)[1]))
    text(2.2, ybottom+(ytop-ybottom)*1/4, round(quantile(weight)[2]))
    text(2.2, ybottom+(ytop-ybottom)*2/4, round(quantile(weight)[3]))
    text(2.2, ybottom+(ytop-ybottom)*3/4, round(quantile(weight)[4]))
    text(2.2, ybottom+(ytop-ybottom)*4/4, round(quantile(weight)[5]))
    text(1.8, ybottom+(ytop-ybottom)*4.5/4, "CCI-Strength", cex=2)
    text(-1.5, -1, "Receptor Patterns", cex=2)
    text(-1.5, 1, "Ligand Patterns", cex=2)

    if(!is.null(twoDplot)){
        # Setting
        maLR <- max(numLPattern, numRPattern)
        if(1 <= maLR && maLR <= 16){
            omasize <- .omasize(numLPattern, numRPattern)
            oma4 <- .oma4(numLPattern, numRPattern)
            #
            # Step.2 : Ligand Plot
            #
            # Color
            col.ligand <- .setColor("reds")
            # Constant
            LOMA_1 = 48.85
            LOMA_2 = 42
            LOMA_3 = 4
            out <- vapply(seq_len(numLPattern), function(i){
                label.ligand <- unlist(vapply(names(label),
                    function(x){
                        outobj$ligand[paste0("Dim", i), x]
                    }, 0.0))
                label.ligand[] <- smoothPalette(label.ligand,
                    palfunc=colorRampPalette(col.ligand, alpha=TRUE))
                par(new=TRUE)
                par(oma = c(LOMA_1, LOMA_2+(i-1)*omasize,
                    LOMA_3, oma4-omasize*i))
                plot(twoDplot, col=label.ligand, pch=16, cex=0.5, bty="n",
                    xaxt="n", yaxt="n", xlab="", ylab="",
                    main=paste0("(", i, ",*,*)"))
                0L
                }, 0L)

            #
            # Step.3 : Receptor Plot
            #
            # Color
            col.receptor <- .setColor("blues")
            # Constant
            ROMA_1 = 4
            ROMA_2 = 42
            ROMA_3 = 48.85
            out <- vapply(seq_len(numRPattern), function(i){
                label.receptor <- unlist(vapply(names(label),
                    function(x){
                        outobj$receptor[paste0("Dim", i), x]
                    }, 0.0))
                label.receptor[] <- smoothPalette(label.receptor,
                    palfunc=colorRampPalette(col.receptor, alpha=TRUE))
                par(ask=FALSE)
                par(new=TRUE)
                par(oma = c(ROMA_1, ROMA_2+(i-1)*omasize,
                    ROMA_3, oma4-omasize*i))
                plot(twoDplot, col=label.receptor, pch=16, cex=0.5,
                    bty="n", xaxt="n", yaxt="n", xlab="", ylab="",
                    main=paste0("(*,", i, ",*)"))
                0L
                }, 0L)
        }else{
            warning(paste0("LR plot can be performed when \n",
                "the maximum number of Ligand/Receptor patterns are \n",
                "higher than 1 and smaller than 12 for now."))
        }
    }
}

.CCIhyperGraphPlot_2 <- function(outobj, twoDplot=NULL, vertex.size=18,
    xleft=1.75, ybottom=-0.5, xright=1.85, ytop=0.5, label="", emph=NULL, algorithm=""){
    # Number of Patterns
    numLPattern <- nrow(outobj$ligand)
    numRPattern <- nrow(outobj$receptor)

    #
    # Step.1 : Background Network
    #
    edgewd_L <- as.vector(vapply(seq_len(numLPattern), function(x){
        rep(x, numRPattern)
        }, rep(0L, numRPattern)))
    edgewd_R <- rep(seq_len(numRPattern), numLPattern)
    edgewd_Strength <- vapply(
        seq_len(numLPattern*numRPattern), function(x){
            targetL <- which(
                outobj$index[, "Mode1"] == edgewd_L[x])
            targetR <- which(
                outobj$index[, "Mode2"] == edgewd_R[x])
            sum(outobj$index[intersect(targetL, targetR), 3])
        }, 0.0)
    edgewd <- cbind(edgewd_L, edgewd_R, edgewd_Strength)
    colnames(edgewd) <- c("L", "R", "Strength")

    # Node name (Top<Ligand> and Bottom<Receptor>)
    nodesSetTop <- paste0("L", seq_len(numLPattern))
    nodesSetBottom <- paste0("R", seq_len(numRPattern))

    # Empty Graph
    g <- graph.empty()

    # Add nodes
    g <- add.vertices(g, nv=length(nodesSetTop),
        attr=list(name=nodesSetTop,
            type=rep(TRUE, length(nodesSetTop))))
    g <- add.vertices(g, nv=length(nodesSetBottom),
        attr=list(name=nodesSetBottom,
            type=rep(TRUE, length(nodesSetBottom))))

    # Add edges
    edgeListVec <- as.vector(t(as.matrix(
        data.frame(
            S1=paste0('L', edgewd[,1]),
            S2=paste0('R', edgewd[,2])
    ))))
    g <- add.edges(g, edgeListVec)

    # Edge weghts
    E(g)$weight <- edgewd[,3]

    # Edge color
    weight <- E(g)$weight
    E(g)$weight <- weight / max(weight) * 20
    mycolor <- smoothPalette(E(g)$weight,
        palfunc=colorRampPalette(.setColor("greens"), alpha=TRUE))
    if(!is.null(emph)){
        target <- intersect(
            which(get.edgelist(g)[, 1] == paste0("L", emph[1])),
            which(get.edgelist(g)[, 2] == paste0("R", emph[2])))
        mycolor[target] <- rgb(1,0,0,0.5)
    }

    # Layout
    x <- c(seq_along(nodesSetTop), seq_along(nodesSetBottom))
    y <- c(rep(1, length=length(nodesSetTop)),
        rep(0, length=length(nodesSetBottom)))
    mylayout <- cbind(x, y)

    # Network Plot
    par(ask=FALSE)
    par(oma=c(2,2,2,2))
    plot.igraph(g,
        layout=mylayout,
        vertex.size=18,
        vertex.label="",
        vertex.color="white",
        vertex.shape="square",
        edge.color=mycolor,
        vertex.frame.color="gray",
        edge.width=E(g)$weight)

    # Gradient
    par(ask=FALSE)
    gradient.rect(xleft, ybottom, xright, ytop,
        col=smoothPalette(sort(weight),
        palfunc=colorRampPalette(.setColor("greens"), alpha=TRUE)),
        gradient="y")
    text(2.2, ybottom+(ytop-ybottom)*0/4, round(quantile(weight)[1]))
    text(2.2, ybottom+(ytop-ybottom)*1/4, round(quantile(weight)[2]))
    text(2.2, ybottom+(ytop-ybottom)*2/4, round(quantile(weight)[3]))
    text(2.2, ybottom+(ytop-ybottom)*3/4, round(quantile(weight)[4]))
    text(2.2, ybottom+(ytop-ybottom)*4/4, round(quantile(weight)[5]))
    text(1.8, ybottom+(ytop-ybottom)*4.5/4, "CCI-Strength", cex=2)
    text(-1.5, -1, "Receptor Patterns", cex=2)
    text(-1.5, 1, "Ligand Patterns", cex=2)

    if(!is.null(twoDplot)){
        # Setting
        maLR <- max(numLPattern, numRPattern)
        if(1 <= maLR && maLR <= 16){
            omasize <- .omasize(numLPattern, numRPattern)
            oma4 <- .oma4(numLPattern, numRPattern)
            #
            # Step.2 : Ligand Plot
            #
            # Color
            col.ligand <- .setColor("reds")
            # Constant
            LOMA_1 = 48.85
            LOMA_2 = 42
            LOMA_3 = 4
            out <- vapply(seq_len(numLPattern), function(i){
                label.ligand <- unlist(vapply(names(label),
                    function(x){
                        outobj$ligand[paste0("Dim", i), x]
                    }, 0.0))
                label.ligand[] <- smoothPalette(label.ligand,
                    palfunc=colorRampPalette(col.ligand, alpha=TRUE))
                par(new=TRUE)
                par(oma = c(LOMA_1, LOMA_2+(i-1)*omasize,
                    LOMA_3, oma4-omasize*i))
                plot(twoDplot, col=label.ligand, pch=16, cex=0.5, bty="n",
                    xaxt="n", yaxt="n", xlab="", ylab="",
                    main=paste0("(", i, ",*,*)"))
                0L
                }, 0L)

            #
            # Step.3 : Receptor Plot
            #
            # Color
            col.receptor <- .setColor("blues")
            # Constant
            ROMA_1 = 4
            ROMA_2 = 42
            ROMA_3 = 48.85
            out <- vapply(seq_len(numRPattern), function(i){
                label.receptor <- unlist(vapply(names(label),
                    function(x){
                        outobj$receptor[paste0("Dim", i), x]
                    }, 0.0))
                label.receptor[] <- smoothPalette(label.receptor,
                    palfunc=colorRampPalette(col.receptor, alpha=TRUE))
                par(ask=FALSE)
                par(new=TRUE)
                par(oma = c(ROMA_1, ROMA_2+(i-1)*omasize,
                    ROMA_3, oma4-omasize*i))
                plot(twoDplot, col=label.receptor, pch=16, cex=0.5,
                    bty="n", xaxt="n", yaxt="n", xlab="", ylab="",
                    main=paste0("(*,", i, ",*)"))
                0L
                }, 0L)
        }else{
            warning(paste0("LR plot can be performed when \n",
                "the maximum number of Ligand/Receptor patterns are \n",
                "higher than 1 and smaller than 12 for now."))
        }
    }
}

.omasize <- function(l, r){
    # 1 - 16
    c(44.7, 44.7, 22.4, 14.9, 11.15,
    8.95, 7.45, 6.4, 5.6, 4.98,
    4.47, 4.06, 3.75, 3.52, 3.25,
    3.05)[max(l, r)]
}

.oma4 <- function(l, r){
    # 1 - 16
    c(131.4, 131.4, 109.2, 101.5, 97.75,
    95.5, 94.1, 93.1, 92.3, 91.7,
    91.2, 90.7, 90.7, 91.70, 91.000,
    91.000)[max(l, r)]
}

.smallTwoDplot <- function(sce, assayNames, geneid, genename, color){
    g <- schex::plot_hexbin_gene(sce, type=assayNames, gene=geneid,
        action="mean", xlab="Dim1", ylab="Dim2",
        title=paste0("Mean of ", genename))
    if(color == "reds"){
        g <- g + scale_fill_gradient(low = 'gray', high = 'red')
    }
    if(color == "blues"){
        g <- g + scale_fill_gradient(low = 'gray', high = 'blue')
    }
    g
}

.annotationhub <- list(
    "Hsa" = function(){
        query(AnnotationHub(), c("OrgDb", "Homo sapiens", "org.Hs.eg.db.sqlite"))[[1]]
    },
    "Mmu" = function(){
        query(AnnotationHub(), c("OrgDb", "Mus musculus", "org.Mm.eg.db.sqlite"))[[1]]
    },
    "Ath" = function(){
        query(AnnotationHub(), c("OrgDb", "Arabidopsis thaliana", "org.At.tair.db.sqlite"))[[1]]
    },
    "Rno" = function(){
        query(AnnotationHub(), c("OrgDb", "Rattus norvegicus", "org.Rn.eg.db.sqlite"))[[1]]
    },
    "Bta" = function(){
        query(AnnotationHub(), c("OrgDb", "Bos taurus", "org.Bt.eg.db.sqlite"))[[1]]
    },
    "Cel" = function(){
        query(AnnotationHub(), c("OrgDb", "Caenorhabditis elegans", "org.Ce.eg.db.sqlite"))[[1]]
    },
    "Dme" = function(){
        query(AnnotationHub(), c("OrgDb", "Drosophila melanogaster", "org.Dm.eg.db.sqlite"))[[1]]
    },
    "Dre" = function(){
        query(AnnotationHub(), c("OrgDb", "Danio rerio", "org.Dr.eg.db.sqlite"))[[1]]
    },
    "Gga" = function(){
        query(AnnotationHub(), c("OrgDb", "Gallus gallus", "org.Gg.eg.db.sqlite"))[[1]]
    },
    "Pab" = function(){
        query(AnnotationHub(), c("OrgDb", "Pongo abelii", "org.Pongo_abelii.eg.sqlite"))[[1]]
    },
    "Xtr" = function(){
        query(AnnotationHub(), c("OrgDb", "Xenopus", "Silurana", "org.Xenopus_\\(Silurana\\)_tropicalis.eg.sqlite"))[[1]]
    },
    "Ssc" = function(){
        query(AnnotationHub(), c("OrgDb", "Sus scrofa", "org.Ss.eg.db.sqlite"))[[1]]
    }
)

'%ni%' <- Negate('%in%')

.geneinformation <- function(sce, ah, spc, LR){
    targetGeneID <- as.character(unique(c(LR$GENEID_L, LR$GENEID_R)))
    if(!is.null(ah)){
        # Gene Name
        message("Related gene names are retrieved from AnnotationHub...")
        GeneName <- AnnotationDbi::select(ah, columns=c("SYMBOL", "ENTREZID"),
            keytype="ENTREZID", keys=targetGeneID)
        GeneName <- unique(GeneName)
        nonna1 <- which(!is.na(GeneName[,1]))
        nonna2 <- which(!is.na(GeneName[,2]))
        GeneName <- GeneName[intersect(nonna1, nonna2), ]
        # Bipartite Matching
        g <- graph.data.frame(as.data.frame(GeneName), directed=FALSE)
        V(g)$type <- bipartite_mapping(g)$type
        g <- max_bipartite_match(g)
        target <- as.character(unique(GeneName[,2]))
        GeneName <- data.frame(
            ENTREZID=as.character(g$matching[target]),
            SYMBOL=target,
            stringsAsFactors = FALSE)
        GeneName <- GeneName[, c("ENTREZID", "SYMBOL")]

        # Description
        message("Related gene descriptions are retrieved from AnnotationHub...")
        Description <- AnnotationDbi::select(ah, columns=c("GENENAME", "ENTREZID"),
            keytype="ENTREZID", keys=targetGeneID)
        Description <- Description[, c("ENTREZID", "GENENAME")]

        # GO
        message("Related GO IDs are retrieved from AnnotationHub...")
        GO <- AnnotationDbi::select(ah, columns=c("GO", "ENTREZID"),
            keytype="ENTREZID", keys=targetGeneID)
        GO <- GO[, c("ENTREZID", "GO")]

        # ENSG
        if(spc %in% c("Hsa", "Mmu", "Rno", "Cel", "Dre")){
            message("Related Ensembl Gene IDs are retrieved from AnnotationHub...")
            ENSG <- AnnotationDbi::select(ah, columns=c("ENSEMBL", "ENTREZID"),
                keytype="ENTREZID", keys=targetGeneID)
            ENSG <- ENSG[, c("ENTREZID", "ENSEMBL")]
        }else{
            ENSG <- NULL
        }

        # ENSP
        if(spc %ni% c("Ath", "Pab", "Xtr", "Ssc")){
            message("Related Ensembl Protein IDs are retrieved from AnnotationHub...")
            ENSP <- AnnotationDbi::select(ah, columns=c("ENSEMBLPROT", "ENTREZID"),
                keytype="ENTREZID", keys=targetGeneID)
            ENSP <- ENSP[, c("ENTREZID", "ENSEMBLPROT")]
        }else{
            ENSP <- NULL
        }

        # UniProtKB
        if(spc %ni% c("Ath", "Xtr")){
            message("Related UniProtKB IDs are retrieved from AnnotationHub...")
            UniProtKB <- AnnotationDbi::select(ah, columns=c("UNIPROT", "ENTREZID"),
                keytype="ENTREZID", keys=targetGeneID)
            UniProtKB <- UniProtKB[, c("ENTREZID", "UNIPROT")]
        }else{
            UniProtKB <- NULL
        }
    }else{
        GeneName <- NULL
        Description <- NULL
        GO <- NULL
        ENSG <- NULL
        ENSP <- NULL
        UniProtKB <- NULL
    }
    # MeSH
    message(paste0("Related MeSH IDs are retrieved from ",
        "MeSH.XXX.eg.db-type package..."))
    MeSHname <- paste0("MeSH.", gsub(".eg.db.sqlite", "",
        strsplit(metadata(sce)$lrbase, "LRBase.")[[1]][3]), ".eg.db")
    MeSHobj <- eval(parse(text=MeSHname))
    MeSH <- MeSHDbi::select(MeSHobj, columns=c("MESHID", "GENEID"),
        keytype="GENEID",
        keys=targetGeneID)
    if(spc != "Pab"){
        # Reactome
        message(paste0("Related Reactome IDs are retrieved from ",
            "reactome.db package..."))
        Reactome <- toTable(reactomeEXTID2PATHID)
        targetReactome <- unlist(lapply(targetGeneID,
            function(x){which(Reactome$gene_id == x)}))
        Reactome <- Reactome[targetReactome, ]
    }else{
        Reactome <- NULL
    }
    # Output
    list(GeneName=GeneName, Description=Description, GO=GO,
        ENSG=ENSG, ENSP=ENSP, UniProtKB=UniProtKB,
        Reactome=Reactome, MeSH=MeSH)
}

.TAXID <- c(
    "Hsa" = 9606,
    "Mmu" = 10090,
    "Ath" = 3702,
    "Rno" = 10116,
    "Bta" = 9913,
    "Cel" = 6239,
    "Dme" = 7227,
    "Dre" = 7955,
    "Gga" = 9031,
    "Pab" = 9601,
    "Xtr" = 8364,
    "Ssc" = 9823
)

.hyperLinks <- function(ranking, ligandGeneID, receptorGeneID,
    lr, value, percentage, spc, geneInfo, pvalue, qvalue){
    ## helper for vector dividing
    div <- function(x, d=1) {""
        delta <- ceiling(length(x) / d)
        y <- lapply(seq_len(d), function(i){
            as.vector(na.omit(x[((i-1)*delta+1):(i*delta)]))
        })
        return(y)
    }

    embedLink <- function(spc, genename1, geneid1, description1,
                go1, reactome1, mesh1,
                uni1, string1, refex1,
                ea1, sea1, scdb1, panglao1, cmap1,
                genename2, geneid2, description2,
                go2, reactome2, mesh2,
                uni2, string2, refex2,
                ea2, sea2, scdb2, panglao2, cmap2){
            paste0(
                "[", genename1,
                "](https://www.ncbi.nlm.nih.gov/gene/",
                geneid1, ")<br>",
                "Description: ", description1, "<br>",
                "GO: ", go1, "<br>",
                "Reactome: ", reactome1, "<br>",
                "UniProtKB: ", uni1, "<br>",
                "STRING: ", string1, "<br>",
                "RefEx: [", genename1, "](", refex1, ")<br>",
                "Expression Atlas: [", genename1, "](", ea1, ")<br>",
                "Single Cell Expression Atlas: [", genename1, "](", sea1, ")<br>",
                "scRNASeqDB: [", genename1, "](", scdb1, ")<br>",
                "PanglaoDB: [", genename1, "](", panglao1, ")<br>",
                "CMap: [", genename1, "](", cmap1, ")<br>",
                "MeSH: ", mesh1,
                "|",
                "[", genename2,
                "](https://www.ncbi.nlm.nih.gov/gene/",
                geneid2, ")<br>",
                "Description: ", description2, "<br>",
                "GO: ", go2, "<br>",
                "Reactome: ", reactome2, "<br>",
                "UniProtKB: ", uni2, "<br>",
                "STRING: ", string2, "<br>",
                "RefEx: [", genename2, "](", refex2, ")<br>",
                "Expression Atlas: [", genename2, "](", ea2, ")<br>",
                "Single Cell Expression Atlas: [", genename2, "](", sea2, ")<br>",
                "scRNASeqDB: [", genename2, "](", scdb2, ")<br>",
                "PanglaoDB: [", genename2, "](", panglao2, ")<br>",
                "CMap: [", genename2, "](", cmap2, ")<br>",
                "MeSH: ", mesh2,
                "|"
            )
    }

    convertGeneName <- function(geneid, geneInfo){
        genename <- geneInfo$GeneName[
            which(geneInfo$GeneName$ENTREZID == geneid),
            "SYMBOL"][1]
        if(is.null(genename)){
            genename = geneid
        }else{
            if(is.na(genename)){
                genename = geneid
            }
            if(length(genename) == 0 || genename == ""){
                genename = geneid
            }
        }
        genename
    }

    convertGeneDescription <- function(geneid, geneInfo){
        description <- geneInfo$Description[
            which(geneInfo$Description$ENTREZID ==
                geneid), "GENENAME"][1]
        if(length(description) == 0 || is.na(description)){
            description = ""
        }
        description
    }

    convertGeneOntology <- function(geneid, geneInfo){
        GO <- unique(unlist(geneInfo$GO[
            which(geneInfo$GO$ENTREZID == geneid),
            "GO"]))
        GO <- GO[which(GO != "")]
        GO <- gsub(":", "%3A", GO)
        if(length(GO) != 0){
            GO_loc <- div(seq_along(GO), ceiling(length(GO) / 100))
            GO <- lapply(GO_loc, function(x){
                mi <- min(x)
                ma <- max(x)
                paste0("[", mi, "-", ma, "](",
                    "http://amigo.geneontology.org/",
                    "goose?query=SELECT+*+FROM+term+WHERE+acc%3D%27",
                    paste0(GO[mi:ma], collapse="%27+OR+acc%3D%27"),
                    "%27%3B&mirror=bbop)")
                })
            GO <- paste(unlist(GO), collapse=" ")
        }else{
            GO <- ""
        }
        GO
    }

    convertReactome <- function(geneid, geneInfo){
        Reactome <- unique(unlist(geneInfo$Reactome[
            which(geneInfo$Reactome$gene_id ==
            geneid), "DB_ID"]))
        Reactome <- Reactome[which(Reactome != "")]
        if(length(Reactome) != 0){
            Reactome_loc <- div(seq_along(Reactome),
                ceiling(length(Reactome) / 100))
            Reactome <- lapply(Reactome_loc, function(x){
                mi <- min(x)
                ma <- max(x)
                paste0("[", mi, "-", ma, "](",
                    "https://reactome.org/content/query?q=",
                    paste0(Reactome[mi:ma], collapse="+"),
                    "&types=Pathway&cluster=true)")
            })
            Reactome = paste(unlist(Reactome), collapse=" ")
        }else{
            Reactome = ""
        }
        Reactome
    }

    convertMeSH <- function(geneid, geneInfo){
        MeSH <- geneInfo$MeSH[which(geneInfo$MeSH$GENEID == geneid),
            "MESHID"]
        MeSH <- MeSH[which(MeSH != "")]
        if(length(MeSH) != 0){
            MeSH_loc <- div(seq_along(MeSH), ceiling(length(MeSH) / 100))
            MeSH <- lapply(MeSH_loc, function(x){
                mi <- min(x)
                ma <- max(x)
                paste0("[", mi, "-", ma, "](",
                    "https://www.ncbi.nlm.nih.gov/mesh?term=",
                    paste0(MeSH[mi:ma], collapse="%20OR%20"), ")")
                })
            MeSH = paste(unlist(MeSH), collapse=" ")
        }else{
            MeSH = ""
        }
        MeSH
    }

    convertUniProtKB <- function(geneid, geneInfo){
        UniProtKB <- unique(unlist(geneInfo$UniProtKB[
            which(geneInfo$UniProtKB$ENTREZID == geneid),
            "UNIPROT"]))
        UniProtKB <- UniProtKB[which(UniProtKB != "")]
        if(length(UniProtKB) != 0){
            UniProtKB_loc <- div(seq_along(UniProtKB),
                ceiling(length(UniProtKB) / 100))
            UniProtKB <- lapply(UniProtKB_loc, function(x){
                mi <- min(x)
                ma <- max(x)
                paste0("[", mi, "-", ma, "](",
                    "https://www.uniprot.org/uniprot/?query=",
                    paste0(UniProtKB[mi:ma], collapse="+OR+"),
                    "&sort=score)")
                })
            UniProtKB <- paste(unlist(UniProtKB), collapse=" ")
        }else{
            UniProtKB <- ""
        }
        UniProtKB
    }

    convertSTRING <- function(geneid, geneInfo, spc){
        ENSP <- unique(unlist(geneInfo$ENSP[
            which(geneInfo$ENSP$ENTREZID == geneid),
            "ENSEMBLPROT"]))
        ENSP <- ENSP[which(ENSP != "")]
        TAXID <- .TAXID[spc]
        if(length(ENSP) != 0 && !is.na(TAXID)){
            STRING <- paste0("https://string-db.org/network/",
                TAXID, ".", ENSP)
            STRING <- paste(
                paste0("[", seq_along(STRING), "](", STRING, ")"),
                collapse=" ")
        }else{
            STRING <- ""
        }
        STRING
    }

    convertRefEx <- function(geneid, geneInfo, spc){
        genename <- geneInfo$GeneName[
            which(geneInfo$GeneName$ENTREZID == geneid),
            "SYMBOL"][1]
        if(length(genename) != 0){
            ref <- "http://refex.dbcls.jp/genelist.php?gene_name%5B%5D="
            if(spc == "Hsa"){
                paste0(ref, genename, "&lang=en&db=human")
            }else if(spc == "Mmu"){
                paste0(ref, genename, "&lang=en&db=mouse")
            }else if(spc == "Rno"){
                paste0(ref, genename, "&lang=en&db=rat")
            }else{
                ""
            }
        }else{
            ""
        }
    }

    convertEA <- function(geneid, geneInfo){
        ENSG <- geneInfo$ENSG[
            which(geneInfo$GeneName$ENTREZID == geneid),
            "ENSEMBL"][1]
        if(length(ENSG) != 0){
            paste0("https://www.ebi.ac.uk/gxa/genes/", tolower(ENSG))
        }else{
            ""
        }
    }

    convertSEA <- function(geneid, geneInfo){
        genename <- geneInfo$GeneName[
            which(geneInfo$GeneName$ENTREZID == geneid),
            "SYMBOL"][1]
        if(length(genename) != 0){
            paste0("https://www.ebi.ac.uk/gxa/sc/search?species=&q=",
                genename)
        }else{
            ""
        }
    }

    convertSCDB <- function(geneid, geneInfo, spc){
        genename <- geneInfo$GeneName[
            which(geneInfo$GeneName$ENTREZID == geneid),
            "SYMBOL"][1]
        if(length(genename) != 0 && spc == "Hsa"){
            paste0("https://bioinfo.uth.edu/scrnaseqdb/",
                "index.php?r=site/rankGene&gene=",
                genename,
                "&check=0")
        }else{
            ""
        }
    }

    convertPANGLAO <- function(geneid, geneInfo, spc){
        genename <- geneInfo$GeneName[
            which(geneInfo$GeneName$ENTREZID == geneid),
            "SYMBOL"][1]
        if(length(genename) != 0){
            ref <- "https://panglaodb.se/search.html?query="
            if(spc == "Hsa"){
                paste0(ref, genename, "&species=3")
            }else if(spc == "Mmu"){
                paste0(ref, genename, "&species=2")
            }else{
                ""
            }
        }else{
            ""
        }
    }

    convertCMAP <- function(geneid, geneInfo, spc){
        genename <- geneInfo$GeneName[
            which(geneInfo$GeneName$ENTREZID == geneid),
            "SYMBOL"][1]
        if(length(genename) != 0){
            ref <- "https://clue.io/command?q="
            if(spc == "Hsa"){
                paste0(ref, genename)
            }else{
                ""
            }
        }else{
            ""
        }
    }

    convertPubMed <- function(geneid1, geneid2, lr){
        target <- intersect(
            which(lr$GENEID_L == geneid1),
            which(lr$GENEID_R == geneid2))
        PubMed <- lr$SOURCEID[target]
        LEN <- length(target) != 0
        NUL <- !is.null(lr$SOURCEID[target])
        HYP <- length(grep("-", PubMed)) == 0
        if(LEN && NUL && HYP){
            PubMed <- unique(strsplit(lr$SOURCEID[target], "\\|")[[1]])
            PubMed_loc <- div(seq_along(PubMed),
                ceiling(length(PubMed) / 100))
            PubMed <- lapply(PubMed_loc, function(x){
                mi <- min(x)
                ma <- max(x)
                paste0("[", mi, "-", ma, "](",
                    "https://www.ncbi.nlm.nih.gov/pubmed/?term=",
                    paste0(PubMed[mi:ma], collapse="%20OR%20"), ")")
            })
            PubMed = paste(unlist(PubMed), collapse=" ")
        }else{
            PubMed = ""
        }
        PubMed
    }

    # ranking
    XYZ <- paste0("|", ranking, "|")
    # Gene Name (Ligand)
    GeneName_L <- convertGeneName(ligandGeneID, geneInfo)
    # Gene Name (Receptor)
    GeneName_R <- convertGeneName(receptorGeneID, geneInfo)
    # Description (Ligand)
    Description_L <- convertGeneDescription(ligandGeneID, geneInfo)
    # Description (Receptor)
    Description_R <- convertGeneDescription(receptorGeneID, geneInfo)
    # Gene Ontology (Ligand)
    GO_L <- convertGeneOntology(ligandGeneID, geneInfo)
    # Gene Ontology (Receptor)
    GO_R <- convertGeneOntology(receptorGeneID, geneInfo)
    # Reactome (Ligand)
    Reactome_L <- convertReactome(ligandGeneID, geneInfo)
    # Reactome (Receptor)
    Reactome_R <- convertReactome(receptorGeneID, geneInfo)
    # MeSH (Ligand)
    MeSH_L <- convertMeSH(ligandGeneID, geneInfo)
    # MeSH (Receptor)
    MeSH_R <- convertMeSH(receptorGeneID, geneInfo)
    # UniProtKB（Ligand）
    UniProtKB_L <- convertUniProtKB(ligandGeneID, geneInfo)
    # UniProtKB（Receptor）
    UniProtKB_R <- convertUniProtKB(receptorGeneID, geneInfo)
    # STRING（Ligand）
    STRING_L <- convertSTRING(ligandGeneID, geneInfo, spc)
    # STRING（Receptor）
    STRING_R <- convertSTRING(receptorGeneID, geneInfo, spc)
    # RefEx（Ligand）
    RefEx_L <- convertRefEx(ligandGeneID, geneInfo, spc)
    # RefEx（Receptor）
    RefEx_R <- convertRefEx(receptorGeneID, geneInfo, spc)
    # EA（Ligand）
    EA_L <- convertEA(ligandGeneID, geneInfo)
    # EA（Receptor）
    EA_R <- convertEA(receptorGeneID, geneInfo)
    # SEA（Ligand）
    SEA_L <- convertSEA(ligandGeneID, geneInfo)
    # SEA（Receptor）
    SEA_R <- convertSEA(receptorGeneID, geneInfo)
    # SCDB（Ligand）
    SCDB_L <- convertSCDB(ligandGeneID, geneInfo, spc)
    # SCDB（Receptor）
    SCDB_R <- convertSCDB(receptorGeneID, geneInfo, spc)
    # PANGLAO（Ligand）
    PANGLAO_L <- convertPANGLAO(ligandGeneID, geneInfo, spc)
    # PANGLAO（Receptor）
    PANGLAO_R <- convertPANGLAO(receptorGeneID, geneInfo, spc)
    # CMAP（Ligand）
    CMAP_L <- convertCMAP(ligandGeneID, geneInfo, spc)
    # CMAP（Receptor）
    CMAP_R <- convertCMAP(receptorGeneID, geneInfo, spc)
    # PubMed (L and R)
    PubMed <- convertPubMed(ligandGeneID, receptorGeneID, lr)

    # Embedding
    paste0(XYZ,
        embedLink(spc,
            GeneName_L, ligandGeneID, Description_L,
            GO_L, Reactome_L, MeSH_L,
            UniProtKB_L, STRING_L, RefEx_L,
            EA_L, SEA_L, SCDB_L, PANGLAO_L, CMAP_L,
            GeneName_R, receptorGeneID, Description_R,
            GO_R, Reactome_R, MeSH_R,
            UniProtKB_R, STRING_R, RefEx_R,
            EA_R, SEA_R, SCDB_R, PANGLAO_R, CMAP_R
            ),
        "![](figures/Ligand/", ligandGeneID, ".png)", "|",
        "![](figures/Receptor/", receptorGeneID, ".png)", "|",
        round(value, 3), " (", round(percentage, 3), "%)",
        "|", round(pvalue, 3),
        "|", round(qvalue, 3),
        "|", PubMed, "|\n")
}

.HCLUST <- function(x){
    out <- hclust(dist(x), method="ward.D")
    cluster <- cutree(out, 2)
    max1 <- max(x[which(cluster == 1)])
    max2 <- max(x[which(cluster == 2)])
    if(max1 > max2){
        cluster[which(cluster == 1)] <- "selected"
        cluster[which(cluster == 2)] <- "not selected"
    }else{
        cluster[which(cluster == 1)] <- "not selected"
        cluster[which(cluster == 2)] <- "selected"
    }
    cluster
}

.OUTLIERS <- function(x, method="grubbs"){
    ranking = length(x) - rank(x) + 1
    vapply(ranking, function(thr){
        remain = which(ranking > thr)
        if(length(remain) > 2){
            if(method == "grubbs"){
                grubbs.test(x[remain])$p
            }else if(method == "chisq"){
                chisq.out.test(x[remain])$p
            }
        }else{
            1
        }
    }, 0.0)
}

.shrink <- function(x){
    block <- strsplit(x , " & ")[[1]]
    l <- length(block)
    nc <- vapply(block, nchar, 0L)
    if(l >= 2){
        out <- block[1]
        for(i in 2:l){
            end <- length(out)
            tmp <- block[i]
            if(nchar(out[end])+nchar(tmp) <= 45){
                out[end] <- paste(c(out[end], tmp), collapse=" & ")
            }else{
                out[end+1] <- tmp
            }
        }
    paste(out, collapse="\n")
    }else{
        x
    }
}

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
    "load(\"reanalysis.RData\")\n\n",
    "col <- .setColor(\"many\")[", x, "]\n",
    "scTensor:::.myvisNetwork(g, col)\n",
    "```\n\n", # Bottom
    "|Rank|Ligand Gene|Receptor Gene|",
    "Ligand Expression (Log10(exp + 1))|",
    "Receptor Expression (Log10(exp + 1))|",
    "LR-pair factor value (Percentage)|",
    "P-value (Grubbs test)|",
    "Q-value (BH method)|",
    "PubMed|\n",
    "|----|----|----|----|----|----|----|----|----|\n")
}

.XYZ_HEADER2_2 <- function(index, i, top){
    paste0("# <font color='#1881c2'>(",
    paste(c(index[i, seq_len(2)], ""), collapse=","),
    ") -related L-R pairs (Top", top, " pairs)</font>\n\n",
    "```{r}\n", # Top
    "library(\"scTensor\")\n",
    "load(\"reanalysis.RData\")\n\n",
    "col <- .setColor(\"many\")[", i, "]\n",
    "scTensor:::.myvisNetwork(g, col)\n",
    "```\n\n", # Bottom
    "|Rank|Ligand Gene|Receptor Gene|",
    "Ligand Expression (Log10(exp + 1))|",
    "Receptor Expression (Log10(exp + 1))|",
    "LR-pair factor value (Percentage)|",
    "P-value (Grubbs test)|",
    "Q-value (BH method)|",
    "PubMed|\n",
    "|----|----|----|----|----|----|----|----|----|\n")
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
    "```\n", # Bottom
    "\n\n## Distribution of mode-1 matricised tensor (Ligand-Cell Direction)\n\n",
    "```{r}\n", # Top
    "rks <- cellCellRanks(sce)\n",
    "mode1value <- rks$mode1\n",
    "names(mode1value)[seq_len(length(mode1value))] <- \"not selected\"\n",
    "names(mode1value)[seq_len(max(index[, \"Mode1\"]))] <- \"selected\"\n",
    "plot_ly(x=seq_along(rks$mode1), y=rks$mode1, type=\"bar\",\n",
    "    color=names(mode1value), colors = c(\"#999999\", \"#E41A1C\"))\n",
    "```\n", # Bottom
    "\n\n## Distribution of mode-2 matricised tensor (Receptor-Cell Direction)\n\n",
    "```{r}\n", # Top
    "mode2value <- rks$mode2\n",
    "names(mode2value)[seq_len(length(mode2value))] <- \"not selected\"\n",
    "names(mode2value)[seq_len(max(index[, \"Mode2\"]))] <- \"selected\"\n",
    "plot_ly(x=seq_along(rks$mode2), y=rks$mode2, type=\"bar\",\n",
    "    color=names(mode2value), colors = c(\"#999999\", \"#E41A1C\"))\n",
    "```\n", # Bottom
    "\n\n## Distribution of mode-3 matricised tensor (LR-pair Direction)\n\n",
    "```{r}\n", # Top
    "mode3value <- rks$mode3\n",
    "names(mode3value)[seq_len(length(mode3value))] <- \"not selected\"\n",
    "names(mode3value)[seq_len(max(index[, \"Mode3\"]))] <- \"selected\"\n",
    "plot_ly(x=seq_along(rks$mode3), y=rks$mode3, type=\"bar\",\n",
    "    color=names(mode3value), colors = c(\"#999999\", \"#E41A1C\"))\n",
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
    "load(\"reanalysis.RData\")\n\n",
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

.GOANNOTATION <- list(
    "Hsa" = "org.Hs.eg.db",
    "Mmu" = "org.Mm.eg.db",
    "Ath" = "org.At.tair.db",
    "Rno" = "org.Rn.eg.db",
    "Bta" = "org.Bt.eg.db",
    "Cel" = "org.Ce.eg.db",
    "Dme" = "org.Dm.eg.db",
    "Dre" = "org.Dr.eg.db",
    "Gga" = "org.Gg.eg.db",
    "Ssc" = "org.Sc.sgd.db"
)

.MESHANNOTATION <- list(
    "Hsa" = "MeSH.Hsa.eg.db",
    "Mmu" = "MeSH.Mmu.eg.db",
    "Ath" = "MeSH.Ath.eg.db",
    "Rno" = "MeSH.Rno.eg.db",
    "Bta" = "MeSH.Bta.eg.db",
    "Cel" = "MeSH.Cel.eg.db",
    "Dme" = "MeSH.Dme.eg.db",
    "Dre" = "MeSH.Dre.eg.db",
    "Gga" = "MeSH.Gga.eg.db",
    "Pab" = "MeSH.Pab.eg.db",
    "Xtr" = "MeSH.Xtr.eg.db",
    "Ssc" = "MeSH.Ssc.eg.db"
)

.REACTOMESPC <- list(
    "Aga" = "anopheles",
    "Ath" = "arabidopsis",
    "Bta" = "bovine",
    "Cfa" = "canine",
    "Cel" = "celegans",
    "Gga" = "chicken",
    "Ptr" = "chimp",
    "Dme" = "fly",
    "Tgo" = "gondii",
    "Hsa" = "human",
    "Pfa" = "malaria",
    "Mmu" = "mouse",
    "Ssc" = "pig",
    "Rno" = "rat",
    "Xla" = "xenopus",
    "Dre" = "zebrafish"
)

.eachVecLR <- function(x, e){
    p <- e$p
    index <- e$index
    sce <- e$sce
    .HCLUST <- e$.HCLUST
    .OUTLIERS <- e$.OUTLIERS
    top <- e$top
    spc <- e$spc
    GeneInfo <- e$GeneInfo
    out.dir <- e$out.dir
    .smallTwoDplot <- e$.smallTwoDplot
    input <- e$input
    twoD <- e$twoD
    .hyperLinks <- e$.hyperLinks
    LR <- e$LR
    .eachRender <- e$.eachRender
    .XYZ_HEADER1 <- e$.XYZ_HEADER1
    .XYZ_HEADER2 <- e$.XYZ_HEADER2
    .XYZ_HEADER3 <- e$.XYZ_HEADER3
    .XYZ_ENRICH <- e$.XYZ_ENRICH
    out.vecLR <- e$out.vecLR
    algorithm <- e$algorithm
    goenrich <- e$goenrich
    meshenrich <- e$meshenrich
    reactomeenrich <- e$reactomeenrich
    doenrich <- e$doenrich
    ncgenrich <- e$ncgenrich
    dgnenrich <- e$dgnenrich

    # Each LR-Pattern Vector
    if(algorithm == "ntd2"){
        vecLR <- metadata(sce)$sctensor$lrpair@data[x[1], x[2],]
    }
    if(algorithm == "ntd"){
        vecLR <- metadata(sce)$sctensor$lrpair[x, ]
    }
    # Clustering
    ClusterLR <- .HCLUST(vecLR)
    # P-value of Grubbs test
    PvalueLR <- suppressWarnings(.OUTLIERS(vecLR))
    # FDR control by BH method
    QvalueLR <- p.adjust(PvalueLR, "BH")
    # TARGET elements by Clustering
    TARGET <- which(ClusterLR == "selected")
    TARGET <- TARGET[order(vecLR[TARGET], decreasing=TRUE)]
    # TOP elements
    if(top != "full"){
        TOP <- min(top, length(TARGET))
        TARGET <- TARGET[seq_len(TOP)]
    }else{
        TOP <- "full"
    }
    # Enrichment (Too Heavy)
    all <- unique(unlist(strsplit(names(vecLR), "_")))
    sig <- unique(unlist(strsplit(names(TARGET), "_")))
    goannotation <- .GOANNOTATION[[spc]]
    meshannotation <- .MESHANNOTATION[[spc]]
    reactomespc <- .REACTOMESPC[[spc]]
    if(!requireNamespace(goannotation, quietly=TRUE)){
        install(goannotation)
    }
    if(!requireNamespace(meshannotation, quietly=TRUE)){
        install(meshannotation)
    }
    dospc <- 0L
    ncgspc <- 0L
    dgnspc <- 0L
    names(dospc) <- "Hsa"
    names(ncgspc) <- "Hsa"
    names(dgnspc) <- "Hsa"
    dospc <- dospc[spc]
    ncgspc <- ncgspc[spc]
    dgnspc <- dgnspc[spc]
    Enrich <- suppressWarnings(.ENRICHMENT(all, sig, goannotation,
        meshannotation, reactomespc, dospc, ncgspc, dgnspc,
        goenrich, meshenrich, reactomeenrich, doenrich, ncgenrich, dgnenrich,
        p))
    # Eigen Value
    # Each LR-Pattern Vector
    if(algorithm == "ntd2"){
        Value <- metadata(sce)$sctensor$lrpair@data[x[1], x[2], TARGET]
        Percentage <- Value / sum(metadata(sce)$sctensor$lrpair@data[x[1], x[2], ]) * 100
    }
    if(algorithm == "ntd"){
        Value <- metadata(sce)$sctensor$lrpair[x, TARGET]
        Percentage <- Value / sum(metadata(sce)$sctensor$lrpair[x, ]) * 100
    }
    # Hyper Link (Too Heavy)
    cat("Hyper-links are embedded...\n")
    GeneName <- GeneInfo$GeneName[, c("SYMBOL", "ENTREZID")]
    LINKS <- vapply(seq_along(TARGET), function(xx){
            # IDs
            Ranking <- xx
            L_R <- strsplit(names(TARGET[xx]), "_")
            LigandGeneID <- L_R[[1]][1]
            ReceptorGeneID <- L_R[[1]][2]
            LigandGeneName <- GeneName[which(GeneName[,2] == LigandGeneID), 1]
            ReceptorGeneName <- GeneName[which(GeneName[,2] == ReceptorGeneID), 1]
            if(length(LigandGeneName) == 0 || LigandGeneName == ""){
                LigandGeneName <- LigandGeneID
            }
            if(length(ReceptorGeneName) == 0 || ReceptorGeneName == ""){
                ReceptorGeneName <- ReceptorGeneID
            }
            # Return Hyper links
            .hyperLinks(Ranking, LigandGeneID,
            ReceptorGeneID, LR, Value[xx],
            Percentage[xx],
            spc, GeneInfo, PvalueLR[xx],
            QvalueLR[xx])
        }, "")
    LINKS <- paste(LINKS, collapse="")

    # Output object
    list(
        ClusterLR=ClusterLR,
        PvalueLR=PvalueLR,
        QvalueLR=QvalueLR,
        TARGET=TARGET,
        TOP=TOP,
        Enrich=Enrich,
        Value=Value,
        Percentage=Percentage,
        LINKS=LINKS
    )
}

.eachRender_2 <- function(x, e, SelectedLR){
    index <- e$index
    out.dir <- e$out.dir
    out.vecLR <- e$out.vecLR
    .XYZ_HEADER1_2 <- e$.XYZ_HEADER1_2
    .XYZ_HEADER2_2 <- e$.XYZ_HEADER2_2
    .XYZ_HEADER3_2 <- e$.XYZ_HEADER3_2
    .XYZ_ENRICH_2 <- e$.XYZ_ENRICH_2

    indexLR <- paste0(c("pattern", SelectedLR[x, ]), collapse="_")
    TARGET <- out.vecLR[[indexLR]]$TARGET
    LINKS <- out.vecLR[[indexLR]]$LINKS

    # Bottom part of Rmarkdown
    XYZ_BOTTOM <- paste(
        c(.XYZ_HEADER2_2(index, x, length(TARGET)),
        LINKS,
        .XYZ_HEADER3_2(index, x),
        .XYZ_ENRICH_2(out.vecLR, indexLR)),
        collapse="")

    # Each (x,y,z)-rmdfile
    RMDFILE <- paste0(c("pattern", index[x, seq_len(2)]), collapse="_")
    RMDFILE <- paste0(RMDFILE, ".Rmd")
    cat(paste0("\n", RMDFILE, " is created...\n"))
    sink(file = paste0(out.dir, "/", RMDFILE))
    cat(paste(
        c(.XYZ_HEADER1_2(index, x), XYZ_BOTTOM),
        collapse=""))
    sink()

    # Rendering
    message(paste0(RMDFILE, " is compiled to ",
        gsub(".Rmd", ".html", RMDFILE)))
    render(paste0(out.dir, "/", RMDFILE), quiet=TRUE)
}

.eachRender <- function(x, e, SelectedLR){
    index <- e$index
    out.dir <- e$out.dir
    out.vecLR <- e$out.vecLR
    .XYZ_HEADER1 <- e$.XYZ_HEADER1
    .XYZ_HEADER2 <- e$.XYZ_HEADER2
    .XYZ_HEADER3 <- e$.XYZ_HEADER3
    .XYZ_ENRICH <- e$.XYZ_ENRICH

    indexLR <- index[x, "Mode3"]
    TARGET <- out.vecLR[, paste0("pattern", indexLR)]$TARGET
    LINKS <- out.vecLR[, paste0("pattern", indexLR)]$LINKS

    # Bottom part of Rmarkdown
    XYZ_BOTTOM <- paste(
        c(.XYZ_HEADER2(indexLR, x, length(TARGET)),
        LINKS,
        .XYZ_HEADER3(indexLR),
        .XYZ_ENRICH(out.vecLR, indexLR)),
        collapse="")

    # Each (x,y,z)-rmdfile
    RMDFILE <- paste0(c("pattern", index[x, seq_len(3)]), collapse="_")
    RMDFILE <- paste0(RMDFILE, ".Rmd")
    cat(paste0("\n", RMDFILE, " is created...\n"))
    sink(file = paste0(out.dir, "/", RMDFILE))
    cat(paste(
        c(.XYZ_HEADER1(index, x), XYZ_BOTTOM),
        collapse=""))
    sink()

    # Rendering
    message(paste0(RMDFILE, " is compiled to ",
        gsub(".Rmd", ".html", RMDFILE)))
    render(paste0(out.dir, "/", RMDFILE), quiet=TRUE)
}

.palf <- colorRampPalette(c("#4b61ba", "gray", "#a87963", "red"))

.shrink2 <- function(x, thr=7){
    block <- strsplit(as.character(x) , " ")[[1]]
    l <- length(block)
    nc <- vapply(block, nchar, 0L)
    if(l >= 2){
        out <- block[1]
        for(i in 2:l){
            end <- length(out)
            tmp <- block[i]
            if(nchar(out[end])+nchar(tmp) <= thr){
                out[end] <- paste(c(out[end], tmp), collapse=" ")
            }else{
                out[end+1] <- tmp
            }
        }
    paste(out, collapse="\n")
    }else{
        x
    }
}

.setColor <- function(col){
    if(col == "reds"){
        c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C",
            "#CB181D", "#A50F15", "#67000D")
    }else if(col == "blues"){
        c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6",
            "#2171B5", "#08519C", "#08306B")
    }else if(col == "greens"){
        c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D",
            "#238B45", "#006D2C", "#00441B")
    }else if(col == "many"){
        c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
            "#A65628", "#F781BF", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
            "#66A61E", "#E6AB02", "#A6761D", "#66C2A5", "#FC8D62", "#8DA0CB",
            "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3E2CD", "#FDCDAC",
            "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#7FC97F",
            "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17",
            "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
            "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F",
            "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC",
            "#E5D8BD", "#FDDAEC")
    }else{
        stop("Wrong col is specified in .setColor")
    }
}

.geneHyperGraphPlot_2 <- function(out.vecLR, GeneInfo, out.dir){
    # Setting
    convertGeneName <- function(geneid, geneInfo){
        genename <- geneInfo$GeneName[
            which(geneInfo$GeneName$ENTREZID == geneid),
            "SYMBOL"][1]
        if(length(genename) == 0){
            genename = geneid
        }
        genename
    }

    # Node
    nodes <- lapply(seq_len(length(out.vecLR)), function(x){
        names(out.vecLR[[x]]$TARGET)})
    Lnodes <- lapply(nodes, function(x){
        vapply(x, function(xx){
            strsplit(xx, "_")[[1]][1]
        }, "")
    })
    Rnodes <-lapply(nodes, function(x){
        vapply(x, function(xx){
            strsplit(xx, "_")[[1]][2]
        }, "")
    })
    Lnodes <- lapply(Lnodes, function(x){
        vapply(x, function(xx){
            convertGeneName(xx, GeneInfo)
        }, "")
    })
    Rnodes <- lapply(Rnodes, function(x){
        vapply(x, function(xx){
            convertGeneName(xx, GeneInfo)
        }, "")
    })
    uniqueLnodes <- unique(unlist(Lnodes))
    uniqueRnodes <- unique(unlist(Rnodes))

    # Empty Graph
    g <- graph.empty(directed=FALSE)
    # Add nodes
    g <- add.vertices(g, nv=length(uniqueLnodes),
        attr=list(name=uniqueLnodes,
            type=rep(TRUE, length(uniqueLnodes)),
            color=rgb(1,0,0,0.5)))
    g <- add.vertices(g, nv=length(uniqueRnodes),
        attr=list(name=uniqueRnodes,
            type=rep(TRUE, length(uniqueRnodes)),
            color=rgb(0,0,1,0.5)))

    # Nodes Weight
    freqLnodes <- vapply(uniqueLnodes, function(x){
        length(which(unlist(Lnodes) == x))
        }, 0L)
    freqRnodes <- vapply(uniqueRnodes, function(x){
        length(which(unlist(Rnodes) == x))
        }, 0L)
    freq <- c(freqLnodes, freqRnodes)
    freq <- freq / max(freq) * 10

    # Add edges
    edgeListVec <- as.vector(t(as.matrix(
        data.frame(
            L=unlist(Lnodes),
            R=unlist(Rnodes)
    ))))
    g <- add.edges(g, edgeListVec)

    # Plot
    cols <- .setColor("many")
    edge.cols <- unlist(lapply(seq_len(length(out.vecLR)), function(x){
            rep(cols[x], length(out.vecLR[[x]]$TARGET))
        }))

    # Setting
    V(g)$size <- freq
    E(g)$color <- edge.cols
    E(g)$width <- 0.7
    l <- layout_with_dh(g)

    # All Pattern
    png(filename=paste0(out.dir, "/figures/GeneHypergraph.png"),
        width=2500, height=2500)
    par(ask=FALSE)
    plot.igraph(g, layout=l)
    par(ask=FALSE)
    legend("topleft",
        legend=c("ligand", "receptor",
            names(out.vecLR)),
        col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5),
            cols[seq_len(length(out.vecLR))]),
        pch=16, cex=2.2)
    dev.off()

    # Each Pattern
    out <- vapply(seq_len(length(out.vecLR)), function(x){
        tmp_edgecolor <- edge.cols
        tmp_edgecolor[which(tmp_edgecolor != cols[x])] <- rgb(0,0,0,0.1)
        tmp_nodecolor <- V(g)$color
        grayout <- setdiff(
            setdiff(
                names(V(g)),
                Lnodes[[x]]
                ), Rnodes[[x]]
            )
        target <- unlist(lapply(grayout, function(xx){
            which(names(V(g)) == xx)
        }))
        tmp_nodecolor[target] <- rgb(0,0,0,0.1)

        # Plot
        png(filename=paste0(
            out.dir, "/figures/GeneHypergraph",
            gsub("pattern", "", names(out.vecLR)[x]),
            ".png"),
            width=2500, height=2500)
        par(ask=FALSE)
        plot.igraph(g,
            vertex.color=tmp_nodecolor,
            edge.color=tmp_edgecolor, layout=l)
        par(ask=FALSE)
        legend("topleft",
            legend=c("ligand", "receptor",
                names(out.vecLR)[x]),
            col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5),
                cols[x]),
            pch=16, cex=2.2)
        dev.off()
    }, 0L)
    return(g)
}

.geneHyperGraphPlot <- function(out.vecLR, GeneInfo, out.dir){
    # Setting
    convertGeneName <- function(geneid, geneInfo){
        genename <- geneInfo$GeneName[
            which(geneInfo$GeneName$ENTREZID == geneid),
            "SYMBOL"][1]
        if(length(genename) == 0){
            genename = geneid
        }
        genename
    }

    # Node
    nodes <- lapply(seq_len(ncol(out.vecLR)), function(x){
        names(out.vecLR["TARGET", x][[1]])})
    Lnodes <- lapply(nodes, function(x){
        vapply(x, function(xx){
            strsplit(xx, "_")[[1]][1]
        }, "")
    })
    Rnodes <-lapply(nodes, function(x){
        vapply(x, function(xx){
            strsplit(xx, "_")[[1]][2]
        }, "")
    })
    Lnodes <- lapply(Lnodes, function(x){
        vapply(x, function(xx){
            convertGeneName(xx, GeneInfo)
        }, "")
    })
    Rnodes <- lapply(Rnodes, function(x){
        vapply(x, function(xx){
            convertGeneName(xx, GeneInfo)
        }, "")
    })
    uniqueLnodes <- unique(unlist(Lnodes))
    uniqueRnodes <- unique(unlist(Rnodes))

    # Empty Graph
    g <- graph.empty(directed=FALSE)
    # Add nodes
    g <- add.vertices(g, nv=length(uniqueLnodes),
        attr=list(name=uniqueLnodes,
            type=rep(TRUE, length(uniqueLnodes)),
            color=rgb(1,0,0,0.5)))
    g <- add.vertices(g, nv=length(uniqueRnodes),
        attr=list(name=uniqueRnodes,
            type=rep(TRUE, length(uniqueRnodes)),
            color=rgb(0,0,1,0.5)))

    # Nodes Weight
    freqLnodes <- vapply(uniqueLnodes, function(x){
        length(which(unlist(Lnodes) == x))
        }, 0L)
    freqRnodes <- vapply(uniqueRnodes, function(x){
        length(which(unlist(Rnodes) == x))
        }, 0L)
    freq <- c(freqLnodes, freqRnodes)
    freq <- freq / max(freq) * 10

    # Add edges
    edgeListVec <- as.vector(t(as.matrix(
        data.frame(
            L=unlist(Lnodes),
            R=unlist(Rnodes)
    ))))
    g <- add.edges(g, edgeListVec)

    # Plot
    cols <- .setColor("many")
    edge.cols <- unlist(lapply(seq_len(ncol(out.vecLR)), function(x){
            rep(cols[x], length(out.vecLR["TARGET", x][[1]]))
        }))

    # Setting
    V(g)$size <- freq
    E(g)$color <- edge.cols
    E(g)$width <- 0.7
    l <- layout_with_dh(g)

    # All Pattern
    png(filename=paste0(out.dir, "/figures/GeneHypergraph.png"),
        width=2500, height=2500)
    par(ask=FALSE)
    plot.igraph(g, layout=l)
    par(ask=FALSE)
    legend("topleft",
        legend=c("ligand", "receptor",
            colnames(out.vecLR)),
        col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5),
            cols[seq_len(ncol(out.vecLR))]),
        pch=16, cex=2.2)
    dev.off()

    # Each Pattern
    out <- vapply(seq_len(ncol(out.vecLR)), function(x){
        tmp_edgecolor <- edge.cols
        tmp_edgecolor[which(tmp_edgecolor != cols[x])] <- rgb(0,0,0,0.1)
        tmp_nodecolor <- V(g)$color
        grayout <- setdiff(
            setdiff(
                names(V(g)),
                Lnodes[[x]]
                ), Rnodes[[x]]
            )
        target <- unlist(lapply(grayout, function(xx){
            which(names(V(g)) == xx)
        }))
        tmp_nodecolor[target] <- rgb(0,0,0,0.1)

        # Plot
        png(filename=paste0(
            out.dir, "/figures/GeneHypergraph_",
            gsub("pattern", "", colnames(out.vecLR)[x]),
            ".png"),
            width=2500, height=2500)
        par(ask=FALSE)
        plot.igraph(g,
            vertex.color=tmp_nodecolor,
            edge.color=tmp_edgecolor, layout=l)
        par(ask=FALSE)
        legend("topleft",
            legend=c("ligand", "receptor",
                colnames(out.vecLR)[x]),
            col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5),
                cols[x]),
            pch=16, cex=2.2)
        dev.off()
    }, 0L)
    return(g)
}

.LIGAND_HEADER <- paste0(
    "# <font color='#1881c2'>Details of Ligand Gene-centric Overview (selected)",
    "</font>\n\n",
    "```{r}\n", # Top
    "library(\"scTensor\")\n",
    "load(\"reanalysis.RData\")\n\n",
    "scTensor:::.myvisNetwork(g)\n",
    "```\n\n", # Bottom
    "<style type='text/css'>\n",
    "table,th,td {\n",
    "width: 20%;\n",
    "border: 1px solid #f0f0f0;\n",
    "}\n",
    "</style>\n\n",
    "|Rank|Frequency|Ligand Gene|Receptor Genes|",
    "Related CCIs|\n",
    "|---------------|---------------|---------------|---------------|---------------|"
)

.RECEPTOR_HEADER <- paste0(
    "# <font color='#1881c2'>Details of Receptor Gene-centric Overview (selected)",
    "</font>\n\n",
    "```{r}\n", # Top
    "library(\"scTensor\")\n",
    "load(\"reanalysis.RData\")\n\n",
    "scTensor:::.myvisNetwork(g)\n",
    "```\n\n", # Bottom
    "<style type='text/css'>\n",
    "table,th,td {\n",
    "width: 20%;\n",
    "border: 1px solid #f0f0f0;\n",
    "}\n",
    "</style>\n\n",
    "|Rank|Frequency|Receptor Gene|Ligand Genes|",
    "Related CCIs|\n",
    "|---------------|---------------|---------------|---------------|---------------|"
)

.LIGAND_BODY_2 <- function(out.vecLR, GeneInfo, index, selected){
    # Setting
    convertGeneName <- function(geneid, geneInfo){
        genename <- geneInfo$GeneName[
            which(geneInfo$GeneName$ENTREZID == geneid), "SYMBOL"]
        if(length(genename) == 0){
            genename = geneid
        }
        if(length(genename) != 1){
            # Black list
            genename = setdiff(genename, "cxcl11.6")[1]
        }
        genename
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
    convertGeneName <- function(geneid, geneInfo){
        genename <- geneInfo$GeneName[
            which(geneInfo$GeneName$ENTREZID == geneid), "SYMBOL"]
        if(length(genename) == 0){
            genename = geneid
        }
        if(length(genename) != 1){
            # Black list
            genename = setdiff(genename, "cxcl11.6")[1]
        }
        genename
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
    convertGeneName <- function(geneid, geneInfo){
        genename <- geneInfo$GeneName[
            which(geneInfo$GeneName$ENTREZID == geneid),
            "SYMBOL"]
        if(length(genename) == 0){
            genename = geneid
        }
        if(length(genename) != 1){
            # Black list
            genename = setdiff(genename, "cxcl11.6")[1]
        }
        genename
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
    convertGeneName <- function(geneid, geneInfo){
        genename <- geneInfo$GeneName[
            which(geneInfo$GeneName$ENTREZID == geneid),
            "SYMBOL"]
        if(length(genename) == 0){
            genename = geneid
        }
        if(length(genename) != 1){
            # Black list
            genename = setdiff(genename, "cxcl11.6")[1]
        }
        genename
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
    "# <font color='#1881c2'>Details of Ligand Gene-centric Overview (all)",
    "</font>\n\n",
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
    "# <font color='#1881c2'>Details of Receptor Gene-centric Overview (all)",
    "</font>\n\n",
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
        ReceptorGeneName <- vapply(ReceptorGeneID, function(xx){
            GeneName[which(GeneName$ENTREZID == xx)[1], "SYMBOL"]
        }, "")
        naposition <- which(vapply(ReceptorGeneName, is.na, TRUE))
        ReceptorGeneName[naposition] <- ReceptorGeneID[naposition]
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
        LigandGeneName <- vapply(LigandGeneID, function(xx){
            GeneName[which(GeneName$ENTREZID == xx)[1], "SYMBOL"]
        }, "")
        naposition <- which(vapply(LigandGeneName, is.na, TRUE))
        LigandGeneName[naposition] <- LigandGeneID[naposition]
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

.GOENRICHMENT <- function(all, sig, goannotation, category, p){
    if(is.na(goannotation)){
        list(Term=NULL, Pvalue=NULL)
    }else{
        goParams <- new("GOHyperGParams",
            geneIds=sig,
            universeGeneIds=all,
            annotation=goannotation,
            ontology=category,
            pvalueCutoff=p,
            conditional=FALSE,
            testDirection="over")
        # Hyper geometric p-value
        out <- try(summary(hyperGTest(goParams)), silent=TRUE)
        if(is(out)[1] == "try-error"){
            list(Term=NULL, Pvalue=NULL)
        }else{
            list(Term=out$Term, Pvalue=out$Pvalue)
        }
    }
}

.MeSHENRICHMENT <- function(all, sig, meshannotation, category, p){
    if(is.na(meshannotation)){
        list(Term=NULL, Pvalue=NULL)
    }else{
        meshParams <- new("MeSHHyperGParams",
            geneIds = sig,
            universeGeneIds = all,
            annotation = meshannotation,
            category = category,
            database = "gene2pubmed",
            pvalueCutoff = p,
            pAdjust = "none")
        # Hyper geometric p-value
        out <- try(summary(meshHyperGTest(meshParams)), silent=TRUE)
        if(is(out)[1] == "try-error"){
            list(Term=NULL, Pvalue=NULL)
        }else{
            outTerm <- unique(out$MESHTERM)
            outPvalue <- sapply(outTerm, function(x){
                out$Pvalue[which(out$MESHTERM == x)[1]]
            })
            list(Term=outTerm, Pvalue=outPvalue)
        }
    }
}

.ReactomeENRICHMENT <- function(all, sig, reactomespc, p){
    if(is.na(reactomespc)){
        list(Term=NULL, Pvalue=NULL)
    }else{
        out <- try(enrichPathway(gene=sig,
          organism=reactomespc,
          pvalueCutoff=p, readable=TRUE), silent=TRUE)
        if(is(out)[1] == "try-error"){
            list(Term=NULL, Pvalue=NULL)
        }else if(is.null(out)){
            list(Term=NULL, Pvalue=NULL)
        }else{
            list(Term=out@result$Description,
                Pvalue=out@result$pvalue)
        }
    }
}

.DOENRICHMENT <- function(all, sig, dospc, p){
    if(is.na(dospc)){
        list(Term=NULL, Pvalue=NULL)
    }else{
        out <- try(enrichDO(gene=sig,
          pvalueCutoff=p, readable=TRUE), silent=TRUE)
        if(is(out)[1] == "try-error"){
            list(Term=NULL, Pvalue=NULL)
        }else if(is.null(out)){
            list(Term=NULL, Pvalue=NULL)
        }else{
            list(Term=out@result$Description,
                Pvalue=out@result$pvalue)
        }
    }
}

.NCGENRICHMENT <- function(all, sig, ncgspc, p){
    if(is.na(ncgspc)){
        list(Term=NULL, Pvalue=NULL)
    }else{
        out <- try(enrichNCG(gene=sig,
          pvalueCutoff=p, readable=TRUE), silent=TRUE)
        if(is(out)[1] == "try-error"){
            list(Term=NULL, Pvalue=NULL)
        }else if(is.null(out)){
            list(Term=NULL, Pvalue=NULL)
        }else{
            list(Term=out@result$Description,
                Pvalue=out@result$pvalue)
        }
    }
}

.DGNENRICHMENT <- function(all, sig, dgnspc, p){
    if(is.na(dgnspc)){
        list(Term=NULL, Pvalue=NULL)
    }else{
        out <- try(enrichDGN(gene=sig,
          pvalueCutoff=p, readable=TRUE), silent=TRUE)
        if(is(out)[1] == "try-error"){
            list(Term=NULL, Pvalue=NULL)
        }else if(is.null(out)){
            list(Term=NULL, Pvalue=NULL)
        }else{
            list(Term=out@result$Description,
                Pvalue=out@result$pvalue)
        }
    }
}

.NOSIG <- list(Term=NULL, PValue=NULL)

.ENRICHMENT <- function(all, sig, goannotation, meshannotation, reactomespc,
    dospc, ncgspc, dgnspc,
    goenrich, meshenrich, reactomeenrich, doenrich, ncgenrich, dgnenrich, p){
    # GO
    if(goenrich){
        cat("GO-Enrichment Analysis is running...(1/3)\n")
        BP <- .GOENRICHMENT(all, sig, goannotation, "BP", p)
        cat("GO-Enrichment Analysis is running...(2/3)\n")
        MF <- .GOENRICHMENT(all, sig, goannotation, "MF", p)
        cat("GO-Enrichment Analysis is running...(3/3)\n")
        CC <- .GOENRICHMENT(all, sig, goannotation, "CC", p)
    }else{
        BP <- .NOSIG
        MF <- .NOSIG
        CC <- .NOSIG
    }
    # MeSH
    if(meshenrich){
        cat("MeSH-Enrichment Analysis is running...(1/16)\n")
        A <- .MeSHENRICHMENT(all, sig, meshannotation, "A", p)
        cat("MeSH-Enrichment Analysis is running...(2/16)\n")
        B <- .MeSHENRICHMENT(all, sig, meshannotation, "B", p)
        cat("MeSH-Enrichment Analysis is running...(3/16)\n")
        C <- .MeSHENRICHMENT(all, sig, meshannotation, "C", p)
        cat("MeSH-Enrichment Analysis is running...(4/16)\n")
        D <- .MeSHENRICHMENT(all, sig, meshannotation, "D", p)
        cat("MeSH-Enrichment Analysis is running...(5/16)\n")
        E <- .MeSHENRICHMENT(all, sig, meshannotation, "E", p)
        cat("MeSH-Enrichment Analysis is running...(6/16)\n")
        F <- .MeSHENRICHMENT(all, sig, meshannotation, "F", p)
        cat("MeSH-Enrichment Analysis is running...(7/16)\n")
        G <- .MeSHENRICHMENT(all, sig, meshannotation, "G", p)
        cat("MeSH-Enrichment Analysis is running...(8/16)\n")
        H <- .MeSHENRICHMENT(all, sig, meshannotation, "H", p)
        cat("MeSH-Enrichment Analysis is running...(9/16)\n")
        I <- .MeSHENRICHMENT(all, sig, meshannotation, "I", p)
        cat("MeSH-Enrichment Analysis is running...(10/16)\n")
        J <- .MeSHENRICHMENT(all, sig, meshannotation, "J", p)
        cat("MeSH-Enrichment Analysis is running...(11/16)\n")
        K <- .MeSHENRICHMENT(all, sig, meshannotation, "K", p)
        cat("MeSH-Enrichment Analysis is running...(12/16)\n")
        L <- .MeSHENRICHMENT(all, sig, meshannotation, "L", p)
        cat("MeSH-Enrichment Analysis is running...(13/16)\n")
        M <- .MeSHENRICHMENT(all, sig, meshannotation, "M", p)
        cat("MeSH-Enrichment Analysis is running...(14/16)\n")
        N <- .MeSHENRICHMENT(all, sig, meshannotation, "N", p)
        cat("MeSH-Enrichment Analysis is running...(15/16)\n")
        V <- .MeSHENRICHMENT(all, sig, meshannotation, "V", p)
        cat("MeSH-Enrichment Analysis is running...(16/16)\n")
        Z <- .MeSHENRICHMENT(all, sig, meshannotation, "Z", p)
    }else{
        A <- .NOSIG
        B <- .NOSIG
        C <- .NOSIG
        D <- .NOSIG
        E <- .NOSIG
        F <- .NOSIG
        G <- .NOSIG
        H <- .NOSIG
        I <- .NOSIG
        J <- .NOSIG
        K <- .NOSIG
        L <- .NOSIG
        M <- .NOSIG
        N <- .NOSIG
        V <- .NOSIG
        Z <- .NOSIG
    }
    # Reactome
    if(reactomeenrich){
        cat("Reactome-Enrichment Analysis is running...(1/1)\n")
        Reactome <- .ReactomeENRICHMENT(all, sig, reactomespc, p)
    }else{
        Reactome <- .NOSIG
    }
    # DO
    if(doenrich){
        cat("DO-Enrichment Analysis is running...(1/1)\n")
        DO <- .DOENRICHMENT(all, sig, dospc, p)
    }else{
        DO <- .NOSIG
    }
    # NCG
    if(ncgenrich){
        cat("NCG-Enrichment Analysis is running...(1/1)\n")
        NCG <- .NCGENRICHMENT(all, sig, ncgspc, p)
    }else{
        NCG <- .NOSIG
    }
    # DGN
    if(dgnenrich){
        cat("DGN-Enrichment Analysis is running...(1/1)\n")
        DGN <- .DGNENRICHMENT(all, sig, dgnspc, p)
    }else{
        DGN <- .NOSIG
    }

    # Output
    out <- list(BP, MF, CC,
            A, B, C, D, E, F, G, H, I, J, K, L, M, N, V, Z,
            Reactome, DO, NCG, DGN)
    # Exception
    out <- lapply(out, function(x){
            if(length(x$Term) == 0 || length(x$Pvalue) == 0){
                .NOSIG
            }else{
                x
            }
        })
    names(out) <- c(
        "GO_BP", "GO_MF", "GO_CC",
        "MeSH_A", "MeSH_B", "MeSH_C", "MeSH_D",
        "MeSH_E", "MeSH_F", "MeSH_G", "MeSH_H",
        "MeSH_I", "MeSH_J", "MeSH_K", "MeSH_L",
        "MeSH_M", "MeSH_N", "MeSH_V", "MeSH_Z",
        "Reactome", "DO", "NCG", "DGN"
        )
    out
}

.eachCircleColor <- c(
    rep(rgb(0, 1, 0, 5E-3), 3),
    rep(rgb(0.5, 0, 1, 5E-3), 16),
    rgb(1, 0.2, 0.4, 5E-3),
    rep(rgb(1, 1, 0, 1E-2), 3))

names(.eachCircleColor) <- c(
    "GO_BP", "GO_MF", "GO_CC",
    "MeSH_A", "MeSH_B", "MeSH_C", "MeSH_D", "MeSH_E", "MeSH_F",
    "MeSH_G", "MeSH_H", "MeSH_I", "MeSH_J", "MeSH_K", "MeSH_L",
    "MeSH_M", "MeSH_N", "MeSH_V", "MeSH_Z",
    "Reactome",
    "DO", "NCG", "DGN")

.tagCloud_2 <- function(out.vecLR, out.dir){
    sapply(seq_len(length(out.vecLR)), function(x){
        # Pvalue
        Pvalues <- list(
            GO_BP=out.vecLR[[x]]$Enrich$GO_BP$Pvalue,
            GO_MF=out.vecLR[[x]]$Enrich$GO_MF$Pvalue,
            GO_CC=out.vecLR[[x]]$Enrich$GO_CC$Pvalue,
            MeSH_A=out.vecLR[[x]]$Enrich$MeSH_A$Pvalue,
            MeSH_B=out.vecLR[[x]]$Enrich$MeSH_B$Pvalue,
            MeSH_C=out.vecLR[[x]]$Enrich$MeSH_C$Pvalue,
            MeSH_D=out.vecLR[[x]]$Enrich$MeSH_D$Pvalue,
            MeSH_E=out.vecLR[[x]]$Enrich$MeSH_E$Pvalue,
            MeSH_F=out.vecLR[[x]]$Enrich$MeSH_F$Pvalue,
            MeSH_G=out.vecLR[[x]]$Enrich$MeSH_G$Pvalue,
            MeSH_H=out.vecLR[[x]]$Enrich$MeSH_H$Pvalue,
            MeSH_I=out.vecLR[[x]]$Enrich$MeSH_I$Pvalue,
            MeSH_J=out.vecLR[[x]]$Enrich$MeSH_J$Pvalue,
            MeSH_K=out.vecLR[[x]]$Enrich$MeSH_K$Pvalue,
            MeSH_L=out.vecLR[[x]]$Enrich$MeSH_L$Pvalue,
            MeSH_M=out.vecLR[[x]]$Enrich$MeSH_M$Pvalue,
            MeSH_N=out.vecLR[[x]]$Enrich$MeSH_N$Pvalue,
            MeSH_V=out.vecLR[[x]]$Enrich$MeSH_V$Pvalue,
            MeSH_Z=out.vecLR[[x]]$Enrich$MeSH_Z$Pvalue,
            Reactome=out.vecLR[[x]]$Enrich$Reactome$Pvalue,
            DO=out.vecLR[[x]]$Enrich$DO$Pvalue,
            NCG=out.vecLR[[x]]$Enrich$NCG$Pvalue,
            DGN=out.vecLR[[x]]$Enrich$DGN$Pvalue
            )
        # Term
        Terms <- list(
            GO_BP=out.vecLR[[x]]$Enrich$GO_BP$Term,
            GO_MF=out.vecLR[[x]]$Enrich$GO_MF$Term,
            GO_CC=out.vecLR[[x]]$Enrich$GO_CC$Term,
            MeSH_A=out.vecLR[[x]]$Enrich$MeSH_A$Term,
            MeSH_B=out.vecLR[[x]]$Enrich$MeSH_B$Term,
            MeSH_C=out.vecLR[[x]]$Enrich$MeSH_C$Term,
            MeSH_D=out.vecLR[[x]]$Enrich$MeSH_D$Term,
            MeSH_E=out.vecLR[[x]]$Enrich$MeSH_E$Term,
            MeSH_F=out.vecLR[[x]]$Enrich$MeSH_F$Term,
            MeSH_G=out.vecLR[[x]]$Enrich$MeSH_G$Term,
            MeSH_H=out.vecLR[[x]]$Enrich$MeSH_H$Term,
            MeSH_I=out.vecLR[[x]]$Enrich$MeSH_I$Term,
            MeSH_J=out.vecLR[[x]]$Enrich$MeSH_J$Term,
            MeSH_K=out.vecLR[[x]]$Enrich$MeSH_K$Term,
            MeSH_L=out.vecLR[[x]]$Enrich$MeSH_L$Term,
            MeSH_M=out.vecLR[[x]]$Enrich$MeSH_M$Term,
            MeSH_N=out.vecLR[[x]]$Enrich$MeSH_N$Term,
            MeSH_V=out.vecLR[[x]]$Enrich$MeSH_V$Term,
            MeSH_Z=out.vecLR[[x]]$Enrich$MeSH_Z$Term,
            Reactome=out.vecLR[[x]]$Enrich$Reactome$Term,
            DO=out.vecLR[[x]]$Enrich$DO$Term,
            NCG=out.vecLR[[x]]$Enrich$NCG$Term,
            DGN=out.vecLR[[x]]$Enrich$DGN$Term
            )
        lapply(names(Pvalues), function(xx){
            # Pvalue
            pval <- eval(parse(text=paste0("Pvalues$", xx)))
            # Term
            t <- as.character(eval(parse(text=paste0("Terms$", xx))))
            # Plot
            png(filename=paste0(out.dir, "/figures/Tagcloud/", xx,
                "_", names(out.vecLR)[x],
                ".png"), width=1000, height=1000)
            if(is.null(pval)){
                .NULLPlot()
            }else if(length(pval) == 1){
                # background circle
                for(i in 1:120){
                    par(ask=FALSE)
                    plot(1,1, cex=(120:1)[i], pch=16, col=.eachCircleColor[xx], xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
                    par(new=TRUE)
                }
                t <- sapply(t, function(x){.shrink2(x, thr=7)})
                .SinglePlot(t)
            }else{
                # background circle
                for(i in 1:120){
                    par(ask=FALSE)
                    plot(1,1, cex=(120:1)[i], pch=16, col=.eachCircleColor[xx], xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
                    par(new=TRUE)
                }
                negLogPval <- -log10(pval+1E-10)
                target <- seq_len(min(100, length(pval)))
                if(length(pval) <= 5){
                    t <- sapply(t, function(x){.shrink2(x, thr=10)})
                }else{
                    t <- sapply(t, function(x){.shrink2(x, thr=25)})
                }
                par(ask=FALSE)
                tagcloud(t[target], weights = negLogPval[target],
                col = smoothPalette(negLogPval[target], palfunc = .palf),
                order = "size", algorithm = "oval",
                scale.multiplier=0.8)
            }
            dev.off()
        })
    })
}

.tagCloud <- function(out.vecLR, out.dir){
    sapply(seq_len(ncol(out.vecLR)), function(x){
        # Pvalue
        Pvalues <- list(
            GO_BP=out.vecLR["Enrich", x][[1]]$GO_BP$Pvalue,
            GO_MF=out.vecLR["Enrich", x][[1]]$GO_MF$Pvalue,
            GO_CC=out.vecLR["Enrich", x][[1]]$GO_CC$Pvalue,
            MeSH_A=out.vecLR["Enrich", x][[1]]$MeSH_A$Pvalue,
            MeSH_B=out.vecLR["Enrich", x][[1]]$MeSH_B$Pvalue,
            MeSH_C=out.vecLR["Enrich", x][[1]]$MeSH_C$Pvalue,
            MeSH_D=out.vecLR["Enrich", x][[1]]$MeSH_D$Pvalue,
            MeSH_E=out.vecLR["Enrich", x][[1]]$MeSH_E$Pvalue,
            MeSH_F=out.vecLR["Enrich", x][[1]]$MeSH_F$Pvalue,
            MeSH_G=out.vecLR["Enrich", x][[1]]$MeSH_G$Pvalue,
            MeSH_H=out.vecLR["Enrich", x][[1]]$MeSH_H$Pvalue,
            MeSH_I=out.vecLR["Enrich", x][[1]]$MeSH_I$Pvalue,
            MeSH_J=out.vecLR["Enrich", x][[1]]$MeSH_J$Pvalue,
            MeSH_K=out.vecLR["Enrich", x][[1]]$MeSH_K$Pvalue,
            MeSH_L=out.vecLR["Enrich", x][[1]]$MeSH_L$Pvalue,
            MeSH_M=out.vecLR["Enrich", x][[1]]$MeSH_M$Pvalue,
            MeSH_N=out.vecLR["Enrich", x][[1]]$MeSH_N$Pvalue,
            MeSH_V=out.vecLR["Enrich", x][[1]]$MeSH_V$Pvalue,
            MeSH_Z=out.vecLR["Enrich", x][[1]]$MeSH_Z$Pvalue,
            Reactome=out.vecLR["Enrich", x][[1]]$Reactome$Pvalue,
            DO=out.vecLR["Enrich", x][[1]]$DO$Pvalue,
            NCG=out.vecLR["Enrich", x][[1]]$NCG$Pvalue,
            DGN=out.vecLR["Enrich", x][[1]]$DGN$Pvalue
            )
        # Term
        Terms <- list(
            GO_BP=out.vecLR["Enrich", x][[1]]$GO_BP$Term,
            GO_MF=out.vecLR["Enrich", x][[1]]$GO_MF$Term,
            GO_CC=out.vecLR["Enrich", x][[1]]$GO_CC$Term,
            MeSH_A=out.vecLR["Enrich", x][[1]]$MeSH_A$Term,
            MeSH_B=out.vecLR["Enrich", x][[1]]$MeSH_B$Term,
            MeSH_C=out.vecLR["Enrich", x][[1]]$MeSH_C$Term,
            MeSH_D=out.vecLR["Enrich", x][[1]]$MeSH_D$Term,
            MeSH_E=out.vecLR["Enrich", x][[1]]$MeSH_E$Term,
            MeSH_F=out.vecLR["Enrich", x][[1]]$MeSH_F$Term,
            MeSH_G=out.vecLR["Enrich", x][[1]]$MeSH_G$Term,
            MeSH_H=out.vecLR["Enrich", x][[1]]$MeSH_H$Term,
            MeSH_I=out.vecLR["Enrich", x][[1]]$MeSH_I$Term,
            MeSH_J=out.vecLR["Enrich", x][[1]]$MeSH_J$Term,
            MeSH_K=out.vecLR["Enrich", x][[1]]$MeSH_K$Term,
            MeSH_L=out.vecLR["Enrich", x][[1]]$MeSH_L$Term,
            MeSH_M=out.vecLR["Enrich", x][[1]]$MeSH_M$Term,
            MeSH_N=out.vecLR["Enrich", x][[1]]$MeSH_N$Term,
            MeSH_V=out.vecLR["Enrich", x][[1]]$MeSH_V$Term,
            MeSH_Z=out.vecLR["Enrich", x][[1]]$MeSH_Z$Term,
            Reactome=out.vecLR["Enrich", x][[1]]$Reactome$Term,
            DO=out.vecLR["Enrich", x][[1]]$DO$Term,
            NCG=out.vecLR["Enrich", x][[1]]$NCG$Term,
            DGN=out.vecLR["Enrich", x][[1]]$DGN$Term
            )
        lapply(names(Pvalues), function(xx){
            # Pvalue
            pval <- eval(parse(text=paste0("Pvalues$", xx)))
            # Term
            t <- as.character(eval(parse(text=paste0("Terms$", xx))))
            # Plot
            png(filename=paste0(out.dir, "/figures/Tagcloud/", xx,
                "_", colnames(out.vecLR)[x],
                ".png"), width=1000, height=1000)
            if(is.null(pval)){
                .NULLPlot()
            }else if(length(pval) == 1){
                # background circle
                for(i in 1:120){
                    par(ask=FALSE)
                    plot(1,1, cex=(120:1)[i], pch=16, col=.eachCircleColor[xx], xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
                    par(new=TRUE)
                }
                t <- sapply(t, function(x){.shrink2(x, thr=7)})
                .SinglePlot(t)
            }else{
                # background circle
                for(i in 1:120){
                    par(ask=FALSE)
                    plot(1,1, cex=(120:1)[i], pch=16, col=.eachCircleColor[xx], xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
                    par(new=TRUE)
                }
                negLogPval <- -log10(pval+1E-10)
                target <- seq_len(min(100, length(pval)))
                if(length(pval) <= 5){
                    t <- sapply(t, function(x){.shrink2(x, thr=10)})
                }else{
                    t <- sapply(t, function(x){.shrink2(x, thr=25)})
                }
                par(ask=FALSE)
                tagcloud(t[target], weights = negLogPval[target],
                col = smoothPalette(negLogPval[target], palfunc = .palf),
                order = "size", algorithm = "oval",
                scale.multiplier=0.8)
            }
            dev.off()
        })
    })
}

.NULLPlot <- function(){
    par(ask=FALSE)
    plot(1, 1, col="white", ann=FALSE, xaxt="n", yaxt="n", axes=FALSE)
    par(ps=100)
    text(1, 1, "None", col=rgb(0,0,0,0.5))
}

.SinglePlot <- function(x){
    par(ask=FALSE)
    plot(1, 1, col="white", ann=FALSE, xaxt="n", yaxt="n", axes=FALSE)
    par(ps=100)
    text(1, 1, x, col="red")
}

.XYZ_ENRICH_2 <- function(out.vecLR, i){
    paste0(
    ################ GO_BP ################
    "## <font color='#1881c2'>GO-Enrichment Analysis (BP : Biological Process)</font>\n\n",
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
    ".png)\n\n",
    ################ GO_MF ################
    "## <font color='#1881c2'>GO-Enrichment Analysis (MF : Molecular Function)</font>\n\n",
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
    ".png)\n\n",
    ################ GO_CC ################
    "## <font color='#1881c2'>GO-Enrichment Analysis (CC : Cellular Component)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_A ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (A : Anatomy)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_B ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (B : Organisms)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_C ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (C : Diseases)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_D ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (D : Drugs)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_E ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (E : Analytical, Diagnostic and Therapeutic Techniques and Equipment)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_F ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (F : Psychiatry and Psychology)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_G ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (G : Phenomena and Processes)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_H ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (H : Disciplines and Occupations)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_I ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (I : Anthropology, Education, Sociology and Social Phenomena)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_J ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (J : Technology and Food and Beverages)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_K ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (K : Humanities)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_L ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (L : Information Science)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_M ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (M : Persons)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_N ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (N : Health Care)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_V ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (V : Publication Type)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_Z ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (Z : Geographical Locations)</font>\n\n",
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
    ".png)\n\n",
    ################ Reactome ################
    "## <font color='#1881c2'>Reactome Pathway-Enrichment Analysis</font>\n\n",
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
    ".png)\n\n",
    ################ DO ################
    "## <font color='#1881c2'>DO-Enrichment Analysis</font>\n\n",
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
    ".png)\n\n",
    ################ NCG ################
    "## <font color='#1881c2'>NCG-Enrichment Analysis</font>\n\n",
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
    ".png)\n\n",
    ################ DGN ################
    "## <font color='#1881c2'>DGN-Enrichment Analysis</font>\n\n",
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
    ".png)\n\n",
    collapse="")
}


.XYZ_ENRICH <- function(out.vecLR, i){
    patternName <- paste0("pattern", i)
    paste0(
    ################ GO_BP ################
    "## <font color='#1881c2'>GO-Enrichment Analysis (BP : Biological Process)</font>\n\n",
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
    ".png)\n\n",
    ################ GO_MF ################
    "## <font color='#1881c2'>GO-Enrichment Analysis (MF : Molecular Function)</font>\n\n",
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
    ".png)\n\n",
    ################ GO_CC ################
    "## <font color='#1881c2'>GO-Enrichment Analysis (CC : Cellular Component)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_A ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (A : Anatomy)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_B ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (B : Organisms)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_C ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (C : Diseases)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_D ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (D : Drugs)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_E ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (E : Analytical, Diagnostic and Therapeutic Techniques and Equipment)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_F ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (F : Psychiatry and Psychology)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_G ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (G : Phenomena and Processes)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_H ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (H : Disciplines and Occupations)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_I ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (I : Anthropology, Education, Sociology and Social Phenomena)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_J ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (J : Technology and Food and Beverages)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_K ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (K : Humanities)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_L ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (L : Information Science)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_M ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (M : Persons)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_N ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (N : Health Care)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_V ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (V : Publication Type)</font>\n\n",
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
    ".png)\n\n",
    ################ MeSH_Z ################
    "## <font color='#1881c2'>MeSH-Enrichment Analysis (Z : Geographical Locations)</font>\n\n",
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
    ".png)\n\n",
    ################ Reactome ################
    "## <font color='#1881c2'>Reactome Pathway-Enrichment Analysis</font>\n\n",
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
    ".png)\n\n",
    ################ DO ################
    "## <font color='#1881c2'>DO-Enrichment Analysis</font>\n\n",
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
    ".png)\n\n",
    ################ NCG ################
    "## <font color='#1881c2'>NCG-Enrichment Analysis</font>\n\n",
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
    ".png)\n\n",
    ################ DGN ################
    "## <font color='#1881c2'>DGN-Enrichment Analysis</font>\n\n",
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
    ".png)\n\n",
    collapse="")
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

.genePlot <- function(sce, assayNames, input, out.dir, GeneInfo, LR){
    GeneName <- GeneInfo$GeneName
    LigandGeneID <- unique(LR$GENEID_L)
    ReceptorGeneID <- unique(LR$GENEID_R)
    LigandGeneName <- vapply(LigandGeneID, function(x){
        GeneName[which(GeneName$ENTREZID == x)[1], "SYMBOL"]}, "")
    ReceptorGeneName <- vapply(ReceptorGeneID, function(x){
        GeneName[which(GeneName$ENTREZID == x)[1], "SYMBOL"]}, "")

    # Plot (Ligand)
    lapply(seq_along(LigandGeneID), function(x){
        Ligandfile <- paste0(out.dir, "/figures/Ligand/",
            LigandGeneID[x], ".png")
        target <- which(rownames(input) == LigandGeneID[x])
        if(length(target) != 0){
            g <- .smallTwoDplot(sce, assayNames, LigandGeneID[x],
                LigandGeneName[x], "reds")
            ggsave(Ligandfile, plot=g, dpi=200, width=6, height=6)
        }
    })

    # Plot (Receptor)
    lapply(seq_along(ReceptorGeneID), function(x){
        Receptorfile <- paste0(out.dir, "/figures/Receptor/",
            ReceptorGeneID[x], ".png")
        target <- which(rownames(input) == ReceptorGeneID[x])
        if(length(target) != 0){
            g <- .smallTwoDplot(sce, assayNames, ReceptorGeneID[x],
                ReceptorGeneName[x], "blues")
            ggsave(Receptorfile, plot=g, dpi=200, width=6, height=6)
        }
    })
}

.ligandPatternPlot <- function(numLPattern, celltypes, sce, col.ligand, ClusterL, out.dir, twoD){
    vapply(seq_len(numLPattern), function(i){
        label.ligand <- unlist(vapply(names(celltypes), function(x){
                metadata(sce)$sctensor$ligand[paste0("Dim", i), x]}, 0.0))
        label.ligand[] <- smoothPalette(label.ligand,
            palfunc=colorRampPalette(col.ligand, alpha=TRUE))
        LPatternfile <- paste0(out.dir, "/figures/Pattern_", i, "__", ".png")
        png(filename=LPatternfile, width=1000, height=1000, bg="transparent")
        par(ask=FALSE)
        par(ps=20)
        plot(twoD, col=label.ligand, pch=16, cex=2, bty="n",
            xaxt="n", yaxt="n", xlab="", ylab="",
            main="")
        dev.off()
    }, 0L)
}

.receptorPatternPlot <- function(numRPattern, celltypes, sce, col.receptor, ClusterR, out.dir, twoD){
    vapply(seq_len(numRPattern), function(i){
        label.receptor <- unlist(vapply(names(celltypes), function(x){
                metadata(sce)$sctensor$receptor[paste0("Dim", i), x]}, 0.0))
        label.receptor[] <- smoothPalette(label.receptor,
            palfunc=colorRampPalette(col.receptor, alpha=TRUE))
        RPatternfile = paste0(out.dir, "/figures/Pattern__", i, "_", ".png")
        png(filename=RPatternfile, width=1000, height=1000, bg="transparent")
        par(ask=FALSE)
        par(ps=20)
        plot(twoD, col=label.receptor, pch=16, cex=2, bty="n",
            xaxt="n", yaxt="n", xlab="", ylab="",
            main="")
        dev.off()
    }, 0L)
}

.GenerateFC <- function (x, thr){
    if (thr == "E1") {
        a <- 0.3213536
        b <- 0.1211649
    }
    else if (thr == "E2") {
        a <- 0.7019536
        b <- 0.3638012
    }
    else if (thr == "E5") {
        a <- 1.9079161
        b <- 0.6665197
    }
    else if (thr == "E10") {
        a <- 4.4291731
        b <- 0.8141489
    }
    else if (thr == "E50") {
        a <- 21.430831
        b <- 1.161132
    }
    else if (thr == "E100") {
        a <- 30.733509
        b <- 1.131347
    }
    else if (is.numeric(thr)) {
        stop("Wrong thr!!!")
    }
    10^(a * exp(-b * log10(x + 1)))
}


# Random Matrix
.matRnbinom <- function(m, disp, rn.index, num.Cell){
    t(vapply(rn.index, function(x) {
        rnbinom(n = num.Cell, mu = m[x], size = 1/disp[x])
    }, 1.0*seq_len(num.Cell)))
}

# Setting DEG
.matFC <- function(nDEG, num.Gene, num.Cell, CCI, m, row.index, rn.index){
    # Set FC
    fc.matrix <- matrix(1, nrow=num.Gene, ncol=length(num.Cell))
    for(x in seq_along(CCI)){
        lp <- which(CCI[[x]]$LPattern == 1)
        r <- row.index[[x]]
        fc.matrix[r, lp] <- .GenerateFC(m[rn.index][r], CCI[[x]]$fc)
    }
    for(x in seq_along(CCI)){
        lp <- which(CCI[[x]]$RPattern == 1)
        r <- row.index[[x]] + sum(nDEG)
        fc.matrix[r, lp] <- .GenerateFC(m[rn.index][r], CCI[[x]]$fc)
    }
    fc.matrix
}

# Assign DEG
.setDEG <- function(original.matrix, fc.matrix, num.Cell, rn.index, row.index, m, disp){
    for (i in seq_along(num.Cell)) {
        deg.index <- which(fc.matrix[, i] != 1)
        col.index <- sum(num.Cell[1:i - 1]) + 1:sum(num.Cell[i])
        if(length(deg.index) != 0){
            original.matrix[deg.index, col.index] <- t(sapply(deg.index,
                function(x){
                    rnbinom(n = num.Cell[i],
                        mu = fc.matrix[x, i] * m[rn.index[x]],
                        size = 1/disp[rn.index[x]])
            }))
        }
    }
    original.matrix
}

# Dropout Matrix
.matDrop <- function(original.matrix, lambda, nCell){
    mean.vector <- apply(original.matrix, 1, mean)
    var.vector <- apply(original.matrix, 1, var)
    droprate <- exp(-lambda * mean.vector^2)
    droprate.matrix <- vapply(seq_len(sum(nCell)), function(y){
        unlist(lapply(droprate, function(x){
            rbinom(1, 1, prob = (1 - x))
        }))
    }, seq_len(nrow(original.matrix)))
}

.simulateDropoutCounts <- function(nGene, nCell, cciInfo, lambda, seed){
    # Set Parameters
    CCI <- lapply(grep("CCI", names(cciInfo)), function(x){
        cciInfo[[x]]
    })
    LPattern <- lapply(grep("CCI", names(cciInfo)), function(x){
        cciInfo[[x]]$LPattern
    })
    RPattern <- lapply(grep("CCI", names(cciInfo)), function(x){
        cciInfo[[x]]$RPattern
    })
    nDEG <- unlist(lapply(grep("CCI", names(cciInfo)), function(x){
        cciInfo[[x]]$nGene
    }))
    fc <- unlist(lapply(grep("CCI", names(cciInfo)), function(x){
        cciInfo[[x]]$fc
    }))
    set.seed(seed)
    data("m")
    data("v")
    disp <- (v - m)/m^2
    rn.index <- sample(which(disp > 0), nGene, replace = TRUE)
    row.index <- list()
    start <- 1
    for(i in seq_along(nDEG)){
        if(i == 1){
            row.index[[i]] <- 1:nDEG[i]
        }else{
            row.index[[i]] <- start:(start+nDEG[i]-1)
        }
        start <- start + nDEG[i]
    }

    # Check
    if(sum(nDEG) > cciInfo$nPair){
        stop("Please specify larger cciInfo$nPair!")
    }
    if(length(nCell) != length(cciInfo$CCI1$LPattern)){
        stop(paste0("Please specify the length of nCell is",
            " same as cciInfo$CCI*$LPattern ",
            "and cciInfo$CCI*$RPattern!"))
    }

    # Original matrix
    original.matrix <- .matRnbinom(m, disp, rn.index, sum(nCell))

    # Setting DEG
    fc.matrix <- .matFC(nDEG, nGene, nCell, CCI, m, row.index, rn.index)
    original.matrix <- .setDEG(original.matrix, fc.matrix,
        nCell, rn.index, row.index, m, disp)

    # Setting Dropout
    droprate.matrix <- .matDrop(original.matrix, lambda, nCell)
    testdata.matrix <- original.matrix * droprate.matrix

    # Naming
    rownames(testdata.matrix) <- paste0("Gene", seq_len(nrow(testdata.matrix)))
    colnames(testdata.matrix) <- paste0("Cell", seq_len(ncol(testdata.matrix)))
    celltypes <- colnames(testdata.matrix)
    names(celltypes) <- unlist(lapply(seq_along(nCell), function(x){
        paste0("Celltype", rep(x, length=nCell[x]))
    }))

    # L-R list
    rn <- rownames(testdata.matrix)
    start1 <- sum(nDEG)+1
    end1 <- 2*sum(nDEG)
    start2 <- end1 + 1
    end2 <- end1 + cciInfo$nPair - end1/2
    start3 <- end2 + 1
    end3 <- end2 + cciInfo$nPair - end1/2
    LR <- rbind(
              cbind(rn[seq_len(sum(nDEG))], rn[start1:end1]),
              cbind(rn[start2:end2], rn[start3:end3]))
    LR <- as.data.frame(LR)
    colnames(LR) <- c("GENEID_L", "GENEID_R")

    # L-R <-> CCI relationship
    LR_CCI <- seq_len(cciInfo$nPair)
    names(LR_CCI) <- c(
        unlist(lapply(seq_along(nDEG), function(x){
            rep(paste0("CCI", x), nDEG[x])
        })),
        rep("nonDEG", cciInfo$nPair - sum(nDEG)))
    rownames(LR) <- LR_CCI

    # Set random seed as default mode
    set.seed(NULL)

    # Output
    return(
        list(simcount = testdata.matrix,
        LR=LR,
        celltypes=celltypes,
        LR_CCI=LR_CCI))
}
