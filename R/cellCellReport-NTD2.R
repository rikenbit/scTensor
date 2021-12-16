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
        lr.evidence <- metadata(sce)$lr.evidence
        LR <- .extractLR(sce, lr.evidence,
            c("GENEID_L", "GENEID_R", "SOURCEID", "SOURCEDB"))
        # SQLite connection
        con = dbConnect(SQLite(), metadata(sce)$lrbase)
        taxid <- dbGetQuery(con, "SELECT * FROM METADATA")
        taxid <- taxid[which(taxid$NAME == "TAXID"), "VALUE"]
        dbDisconnect(con)
        ###########################################
        # Taxonomy ID based information retrieval
        ###########################################
        # biomaRt Setting
        ah <- .annotationhub_taxid(taxid)
        # GeneName, Description, GO, Reactome, MeSH
        GeneInfo <- .geneInformation_taxid(sce, ah, taxid, LR)
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
        e$ah <- ah
        e$.HCLUST <- .HCLUST
        e$.OUTLIERS <- .OUTLIERS
        e$top <- top
        e$GeneInfo <- GeneInfo
        e$out.dir <- out.dir
        e$.smallTwoDplot <- .smallTwoDplot
        e$input <- input
        e$twoD <- twoD
        e$.hyperLinks <- .hyperLinks
        e$LR <- LR
        e$taxid <- taxid
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
        png(filename=paste0(out.dir, "/figures/CCIHypergraph.png"),
            width=2000, height=950)
        invisible(.CCIhyperGraphPlot_2(metadata(sce)$sctensor,
            twoDplot=twoD, label=celltypes))
        dev.off()
        # Plot（Gene-wise Hypergraph）
        invisible(g <- .geneHyperGraphPlot_2(out.vecLR, GeneInfo, out.dir))

        # Rmd（ligand, selected）
        message("ligand.Rmd is created...")
        outLg <- file(paste0(out.dir, "/ligand.Rmd"), "w")
        writeLines(.LIGAND_HEADER, outLg, sep="\n")
        writeLines(.LIGAND_BODY_2(out.vecLR, GeneInfo, index, selected), outLg, sep="\n")
        close(outLg)
        # Rmd（receptor, selected）
        message("receptor.Rmd is created...")
        outRp <- file(paste0(out.dir, "/receptor.Rmd"), "w")
        writeLines(.RECEPTOR_HEADER, outRp, sep="\n")
        writeLines(.RECEPTOR_BODY_2(out.vecLR, GeneInfo, index, selected), outRp, sep="\n")
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
        message("\nindex.Rmd is created...")
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
        message("index.Rmd is compiled to index.html...\n")
        render(paste0(out.dir, "/index.Rmd"), quiet=TRUE)
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

.geneHyperGraphPlot_2 <- function(out.vecLR, GeneInfo, out.dir){
    # Setting
    convertGeneName <- function(geneid, GeneInfo){
        if(!is.null(GeneInfo$GeneName)){
            genename <- GeneInfo$GeneName[
                which(GeneInfo$GeneName$ENTREZID == geneid),
                "SYMBOL"][1]
            if(is.na(genename[1])){
                genename = geneid
            }
            genename            
        }else{
            geneid            
        }
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
    plot.igraph(g, layout=l)
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
        plot.igraph(g,
            vertex.color=tmp_nodecolor,
            edge.color=tmp_edgecolor, layout=l)
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
            if(!is.null(pval)){                
                png(filename=paste0(out.dir, "/figures/Tagcloud/", xx,
                    "_", names(out.vecLR)[x],
                    ".png"), width=1000, height=1000)
                if(length(pval) == 1){
                    t <- sapply(t, function(x){.shrink2(x, thr=7)})
                    .SinglePlot(t)
                }else{
                    negLogPval <- -log10(pval+1E-10)
                    target <- seq_len(min(100, length(pval)))
                    if(length(pval) <= 5){
                        t <- sapply(t, function(x){.shrink2(x, thr=10)})
                    }else{
                        t <- sapply(t, function(x){.shrink2(x, thr=25)})
                    }
                    tagcloud(t[target], weights = negLogPval[target],
                    col = smoothPalette(negLogPval[target], palfunc = .palf),
                    order = "size", algorithm = "fill",
                    scale.multiplier=0.8)
                }
                dev.off()
            }
        })
    })
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
