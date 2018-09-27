.cellType <- function(input, celltypes, mergeas="mean"){
    all.celltypes = unique(names(celltypes))
    if(mergeas == "mean"){
        ct <- sapply(all.celltypes, function(x){
                rowMeans(
                    input[, which(names(celltypes) == x), drop=FALSE])})
    }else if(mergeas == "sum"){
        ct <- sapply(all.celltypes, function(x){
                rowSums(
                    input[, which(names(celltypes) == x), drop=FALSE])})
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
        out <- sapply(seq_len(ncol(A1)), function(x){
            base::outer(A1[, x], A2[, x])}, simplify=FALSE)
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
    colnames(out) <- c("Mode1", "Mode2", "Mode3", "Value", "Rank")
    veccore <- c()
    counter <- 1
    for(i in seq_len(d[1])){
        for(j in seq_len(d[2])){
            for(k in seq_len(d[3])){
                veccore <- c(veccore, core[i,j,k])
                out[counter, seq_len(3)] <- c(i, j, k)
                out[counter, 4] <- core[i,j,k]
                counter <- counter + 1
            }
        }
    }
    index <- nrow(out) - base::rank(veccore) + 1
    out[, 5] <- index
    out[order(out[,4], decreasing=TRUE), ]
}

.celltypemergedtensor <- function(input, LR, celltypes, mergeas, outer){
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
                base::outer(as.vector(L[l,]), as.vector(R[r,]), outer), along=3)
            Pair.name <- c(Pair.name,
                paste(LR$GENEID_L[i], LR$GENEID_R[i], sep="_"))
        }
    }
    list(tnsr=tnsr, pairname=Pair.name)
}

.fastPossibleCombination <- function(input, LR, celltypes, num.sampling, outer){
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
                target <- sapply(unique(names(celltypes)), function(x){
                    sample(which(names(celltypes) == x), 1)})
                pre_tnsr[,,1] <- pre_tnsr[,,1] +
                base::outer(as.vector(L[l, target]),
                as.vector(R[r, target]), outer)
            }
            tnsr <- abind(tnsr, pre_tnsr, along=3)
            Pair.name <- c(Pair.name, paste(LR$GENEID_L[i],
                LR$GENEID_R[i], sep="_"))
        }
    }
    dimnames(tnsr) <- list(unique(names(celltypes)),
        unique(names(celltypes)), Pair.name)
    list(tnsr=tnsr, pairname=Pair.name)
}

.slowPossibleCombination <- function(input, LR, celltypes, outer){
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
            pre_tnsr <- base::outer(as.vector(L[l,]), as.vector(R[r,]))
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
    dimnames(tnsr) <- list(unique(names(celltypes)),
        unique(names(celltypes)), Pair.name)
    list(tnsr=tnsr, pairname=Pair.name)
}

.cellCellDecomp.Third <- function(input, LR, celltypes, ranks, centering,
    mergeas, outer, comb, num.sampling, decomp){
    if(centering){
        fout <- .celltypemergedtensor(input, LR, celltypes, mergeas, outer)
    }else{
        if(comb == "random"){
            fout <- .fastPossibleCombination(input, LR, celltypes,
                num.sampling, outer)
        }else if(comb == "all"){
            fout <- .slowPossibleCombination(input, LR, celltypes, outer)
        }
    }
    if(dim(fout$tnsr)[1] >= 100 || dim(fout$tnsr)[2] >= 100){
        stop(paste0("Extreamly Large tensor will be generated!\n",
            "This problem will be solved by the next version of scTensor."))
    }
    tnsr <- as.tensor(fout$tnsr)
    Pair.name <- fout$pairname
    if(decomp){
        cat(paste0(paste(dim(tnsr), collapse=" * "), " Tensor is created\n"))
        out <- try(NTD(X=tnsr, rank=ranks, algorithm="KL"))
        if(is(out)[1] == "try-error"){
            out <- NTD(X=tnsr, rank=ranks, algorithm="KL")
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
        lrpairout <- .cellCellReconst(A1, A2, A3, index, type="tensor",
            names=list(colnames(A1), colnames(A2), colnames(A3)),
            collapse=TRUE)
        return(list(index=index, ligand=A1,
            receptor=A2, lrpair=A3, cellcellpattern=lrpairout))
    }else{
        tnsr_cc <- modeSum(tnsr, m=3, drop=TRUE)@data
        dimnames(tnsr_cc) <- list(unique(names(celltypes)),
            unique(names(celltypes)))
        return(list(cellcellpattern=tnsr_cc, cellcelllrpairpattern=tnsr))
    }
}

.cellCellDecomp.Second <-
function(input, LR, celltypes, rank, centering, mergeas, outer, comb,
    num.sampling, decomp){
    if(centering){
        fout <- .celltypemergedtensor(input, LR, celltypes, mergeas, outer)
    }else{
        if(comb == "random"){
            fout <- .fastPossibleCombination(input, LR, celltypes,
                num.sampling, outer)
        }else if(comb == "all"){
            fout <- .slowPossibleCombination(input, LR, celltypes, outer)
        }
    }
    tnsr <- as.tensor(fout$tnsr)
    Pair.name <- fout$pairname
    mat.tnsr <- modeSum(tnsr, m=3, drop=TRUE)@data
    dimnames(mat.tnsr) <- list(unique(names(celltypes)),
        unique(names(celltypes)))
    if(decomp){
        out <- try(NMF(mat.tnsr, J=rank, algorithm="KL"))
        if(is(out)[1] == "try-error"){
            out <- NMF(mat.tnsr, J=rank, algorithm="KL")
        }
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

.cellCellDecomp.Pearson <- function(input, celltypes){
    out <- .cellType(input, celltypes)
    out <- cor(out, method="pearson")
    out <- as.matrix(out)
    return(list(cellcellpattern=out))
}

.cellCellDecomp.Spearman <- function(input, celltypes){
    out <- .cellType(input, celltypes)
    out <- cor(out, method="spearman")
    out <- as.matrix(out)
    return(list(cellcellpattern=out))
}

.cellCellDecomp.Distance <- function(input, celltypes){
    out <- .cellType(input, celltypes)
    out <- 1 / dist(t(out))
    out <- as.matrix(out)
    diag(out) <- 1
    return(list(cellcellpattern=out))
}

.cellCellDecomp.Pearson.LR <- function(input, LR, celltypes){
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

.cellCellDecomp.Spearman.LR <- function(input, LR, celltypes){
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

.cellCellDecomp.Distance.LR <- function(input, LR, celltypes){
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

.cellCellDecomp.PossibleCombination <- function(input, LR, celltypes,
    thr1=2^5, thr2=25){
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
            pre_tnsr_bin <- base::outer(as.vector(L[l,]), as.vector(R[r,]))
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

.CCIhyperGraphPlot <- function(outobj, twoDplot=NULL, vertex.size=18, xleft=1.75,
    ybottom=-0.5, xright=1.85, ytop=0.5, label="", emph=NULL){
    # Number of Patterns
    numLPattern <- nrow(outobj$ligand)
    numRPattern <- nrow(outobj$receptor)

    #
    # Step.1 : Background Network
    #
    edgewd_L <- unlist(sapply(seq_len(numLPattern), function(x){
        rep(x, numRPattern)
        }, simplify=FALSE))
    edgewd_R <- rep(seq_len(numRPattern), numLPattern)
    edgewd_Strength <- sapply(
        seq_len(numLPattern*numRPattern), function(x){
            targetL <- which(
                outobj$index[, "Mode1"] == edgewd_L[x])
            targetR <- which(
                outobj$index[, "Mode2"] == edgewd_R[x])
            sum(outobj$index[intersect(targetL, targetR), 4])
        })
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
        if(1 <= maLR && maLR <= 12){
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
            sapply(seq_len(numLPattern), function(i){
                label.ligand <- unlist(sapply(names(label),
                    function(x){outobj$ligand[paste0("Dim", i), x]}))
                label.ligand[] <- smoothPalette(label.ligand,
                    palfunc=colorRampPalette(col.ligand, alpha=TRUE))
                par(new=TRUE)
                par(oma = c(LOMA_1, LOMA_2+(i-1)*omasize,
                    LOMA_3, oma4-omasize*i))
                plot(twoDplot, col=label.ligand, pch=16, cex=0.5, bty="n",
                    xaxt="n", yaxt="n", xlab="", ylab="",
                    main=paste0("(", i, ",*,*)"))
                })

            #
            # Step.3 : Receptor Plot
            #
            # Color
            col.receptor <- .setColor("blues")
            # Constant
            ROMA_1 = 4
            ROMA_2 = 42
            ROMA_3 = 48.85
            sapply(seq_len(numRPattern), function(i){
                label.receptor <- unlist(sapply(names(label),
                    function(x){outobj$receptor[paste0("Dim", i), x]}))
                label.receptor[] <- smoothPalette(label.receptor,
                    palfunc=colorRampPalette(col.receptor, alpha=TRUE))
                par(new=TRUE)
                par(oma = c(ROMA_1, ROMA_2+(i-1)*omasize,
                    ROMA_3, oma4-omasize*i))
                plot(twoDplot, col=label.receptor, pch=16, cex=0.5,
                    bty="n", xaxt="n", yaxt="n", xlab="", ylab="",
                    main=paste0("(*,", i, ",*)"))
                })
        }else{
            warning(paste0("LR plot can be performed when \n",
                "the maximum number of Ligand/Receptor patterns are \n",
                "higher than 1 and smaller than 8 for now."))
        }
    }
}

.omasize <- function(l, r){
    # 1 - 12
    c(44.7, 44.7, 22.4, 14.9, 11.15, 8.95,
        7.45, 6.4, 5.6,4.98, 4.47, 4.06)[max(l, r)]
}

.oma4 <- function(l, r){
    # 1 - 12
    c(131.4, 131.4, 109.2, 101.5, 97.75, 95.5,
        94.1, 93.1, 92.3, 91.7, 91.2, 90.7)[max(l, r)]
}

.smallTwoDplot <- function(input, geneid, genename, twoD, color){
    target <- which(rownames(input) == geneid)
    exp <- log10(input[target, ]+1)
    label <- smoothPalette(exp,
        palfunc=colorRampPalette(.setColor(color), alpha=TRUE))
    # Plot
    par(ps=25)
    par(oma=c(2,2,2,2))
    plot(twoD, col=label, pch=16, cex=3, bty="n", xaxt="n", yaxt="n",
        xlab="", ylab="", main=genename)
    # Gradient
    xleft = max(twoD[,1])+1
    ybottom=quantile(twoD[,2])[2]
    xright = max(twoD[,1])+2
    ytop=quantile(twoD[,2])[4]
    gradient.rect(xleft, ybottom, xright, ytop, col=smoothPalette(sort(exp),
        palfunc=colorRampPalette(.setColor(color),
        alpha=TRUE)), gradient="y")
    text(xleft-1, ybottom+(ytop-ybottom)*0/4, round(quantile(exp)[1]))
    text(xleft-1, ybottom+(ytop-ybottom)*1/4, round(quantile(exp)[2]))
    text(xleft-1, ybottom+(ytop-ybottom)*2/4, round(quantile(exp)[3]))
    text(xleft-1, ybottom+(ytop-ybottom)*3/4, round(quantile(exp)[4]))
    text(xleft-1, ybottom+(ytop-ybottom)*4/4, round(quantile(exp)[5]))
}

.ensembl <- function(species){
    if(species == "Hsa"){
        useMart("ensembl", dataset="hsapiens_gene_ensembl")
    }else if(species == "Mmu"){
        useMart("ensembl", dataset="mmusculus_gene_ensembl")
    }else if(species == "Rno"){
        useMart("ensembl", dataset="rnorvegicus_gene_ensembl")
    }else if(species == "Bta"){
        useMart("ensembl", dataset="btaurus_gene_ensembl")
    }else if(species == "Cel"){
        useMart("ensembl", dataset="celegans_gene_ensembl")
    }else if(species == "Dme"){
        useMart("ensembl", dataset="dmelanogaster_gene_ensembl")
    }else if(species == "Dre"){
        useMart("ensembl", dataset="drerio_gene_ensembl")
    }else if(species == "Gga"){
        useMart("ensembl", dataset="ggallus_gene_ensembl")
    }else if(species == "Pab"){
        useMart("ensembl", dataset="pabelii_gene_ensembl")
    }else if(species == "Xtr"){
        useMart("ensembl", dataset="xtropicalis_gene_ensembl")
    }else if(species == "Ssc"){
        useMart("ensembl", dataset="sscrofa_gene_ensembl")
    }else{
        NULL
    }
}

.geneinformation <- function(sce, ens, spc, LR){
    targetGeneID <- unique(c(LR$GENEID_L, LR$GENEID_R))
    if(!is.null(ens)){
        # Gene Name
        cat(paste0("Related gene names are retrieved from Biomart",
        	" by biomaRt package...\n"))
        GeneName <- getBM(
            attributes=c('external_gene_name', 'entrezgene'),
            filters = 'entrezgene',
            values = targetGeneID,
            mart = ens)
        # Description
        cat("Related gene descriptions are retrieved from Biomart by biomaRt package...\n")
        Description <- getBM(
            attributes=c('description', 'entrezgene'),
            filters = 'entrezgene',
            values = targetGeneID,
            mart = ens)
        # GO
        cat("Related GO IDs are retrieved from Biomart by biomaRt package...\n")
        GO <- getBM(
            attributes=c('go_id', 'entrezgene'),
            filters = 'entrezgene',
            values = targetGeneID,
            mart = ens)
    }else{
        GeneName <- NULL
        Description <- NULL
        GO <- NULL
    }
    # MeSH
    cat(paste0("Related MeSH IDs are retrieved from ",
        "MeSH.XXX.eg.db-type package...\n"))
    MeSHobj <- eval(parse(text=gsub("LRBase", "MeSH",
        LRBaseDbi::lrPackageName(metadata(sce)$lrbase))))
    MeSH <- MeSHDbi::select(MeSHobj, columns=c("MESHID", "GENEID"),
        keytype="GENEID",
        keys=targetGeneID)
    if(spc != "Pab"){
        # Reactome
        cat(paste0("Related Reactome IDs are retrieved from ",
            "reactome.db package...\n"))
        Reactome <- toTable(reactomeEXTID2PATHID)
        targetReactome <- unlist(sapply(targetGeneID,
            function(x){which(Reactome$gene_id == x)}))
        Reactome <- Reactome[targetReactome, ]
    }else{
        Reactome <- NULL
    }
    # Output
    list(GeneName=GeneName, Description=Description, GO=GO,
        Reactome=Reactome, MeSH=MeSH)
}

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

    embedBySpc <- function(spc, genename1, geneid1, description1,
                go1, reactome1, mesh1,
                genename2, geneid2, description2,
                go2, reactome2, mesh2){
        if(spc == "Ath"){
            paste0("[", geneid1,
                "](https://www.ncbi.nlm.nih.gov/gene/",
                geneid1, ")<br>",
                "Reactome: ", reactome1, "<br>",
                "MeSH: ", mesh1,
                "|",
                "[", geneid2,
                "](https://www.ncbi.nlm.nih.gov/gene/",
                geneid2, ")<br>",
                "Reactome: ", reactome2, "<br>",
                "MeSH: ", mesh2,
                "|"
            )
        }else if(spc == "Pab"){
            paste0(
                "[", genename1,
                "](https://www.ncbi.nlm.nih.gov/gene/",
                geneid1, ")<br>",
                "Description: ", description1, "<br>",
                "GO: ", go1, "<br>",
                "MeSH: ", mesh1,
                "|",
                "[", genename2,
                "](https://www.ncbi.nlm.nih.gov/gene/",
                geneid2, ")<br>",
                "Description: ", description2, "<br>",
                "GO: ", go2, "<br>",
                "MeSH: ", mesh2,
                "|"
            )
        }else{
            paste0(
                "[", genename1,
                "](https://www.ncbi.nlm.nih.gov/gene/",
                geneid1, ")<br>",
                "Description: ", description1, "<br>",
                "GO: ", go1, "<br>",
                "Reactome: ", reactome1, "<br>",
                "MeSH: ", mesh1,
                "|",
                "[", genename2,
                "](https://www.ncbi.nlm.nih.gov/gene/",
                geneid2, ")<br>",
                "Description: ", description2, "<br>",
                "GO: ", go2, "<br>",
                "Reactome: ", reactome2, "<br>",
                "MeSH: ", mesh2,
                "|"
            )
        }
    }

    convertGeneName <- function(geneid, geneInfo){
        genename <- geneInfo$GeneName[
            which(geneInfo$GeneName$entrezgene == geneid),
            "external_gene_name"][1]
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
            which(geneInfo$Description$entrezgene ==
                geneid), "description"][1]
        if(length(description) == 0 || is.na(description)){
            description = ""
        }
        description
    }

    convertGeneOntology <- function(geneid, geneInfo){
        GO <- unique(unlist(geneInfo$GO[
            which(geneInfo$GO$entrezgene == geneid),
            "go_id"]))
        GO <- GO[which(GO != "")]
        GO <- gsub(":", "%3A", GO)
        if(length(GO) != 0){
            tmp_GO <- list()
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
    # PubMed (L and R)
    PubMed <- convertPubMed(ligandGeneID, receptorGeneID, lr)
    # Embedding
    paste0(XYZ,
        embedBySpc(spc,
            GeneName_L, ligandGeneID, Description_L,
            GO_L, Reactome_L, MeSH_L,
            GeneName_R, receptorGeneID, Description_R,
            GO_R, Reactome_R, MeSH_R
            ),
        "![](figures/Ligand_", ligandGeneID, ".png)", "|",
        "![](figures/Receptor_", receptorGeneID, ".png)", "|",
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
  sapply(ranking, function(thr){
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
  })
}

.shrink <- function(x){
    block <- strsplit(x , " & ")[[1]]
    l <- length(block)
    nc <- sapply(block, nchar)
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
    "        toc_depth: 2\n",
    "---\n",
    "# <font color='#1881c2'>(", paste(index[i, seq_len(2)], collapse=","),
    ",\\*) Pattern</font>\n\n",
    "![](figures/CCIHypergraph_",
    paste(index[i, seq_len(2)], collapse="_"), ".png)\n\n")
    sub("XXXXX", TITLE, HEADER)
}

.XYZ_HEADER2 <- function(i, top){
    paste0("# <font color='#1881c2'>(\\*,\\*,", i, ") Pattern",
    " (Top", top, " LR-pairs)</font>\n\n",
    paste0("![](figures/GeneHypergraph_", i, ".png)\n\n", collapse=""),
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
        "output: ",
        "BiocStyle::html_document\nvignette: >\n ",
        "%\\VignetteIndexEntry{Vignette Title}\n ",
        "%\\VignetteEngine{knitr::rmarkdown}\n ",
        "%\\VignetteEncoding{UTF-8}\n---\n")
    sub("YYYYY", author, sub("XXXXX", title, HEADER))
}

# 1. About Dataset and scTensor Algorithm
.BODY1 <- paste0("\n\n# About scTensor Algorithm\n\n",
    "![](Workflow.jpeg)\n",
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
    "type = \"scatter\", ",
    "text = rownames(twoD), ",
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
.BODY3 <- function(numLPattern){
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
    "labCol = colnames(l)",
    ")\n",
    "```\n" # Bottom
    )
    BODY3 <- paste0(BODY3,
        paste(sapply(seq_len(numLPattern), function(i){
            paste0("![](figures/Pattern_", i, "__.png)")
        }), collapse=""))
    BODY3 <- paste0(BODY3, "\n")

}

# 4. R-Pattern
.BODY4 <- function(numRPattern){
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
        "labCol = colnames(r)",
        ")\n",
        "```\n" # Bottom
        )
    BODY4 <- paste0(BODY4,
        paste(sapply(seq_len(numRPattern), function(i){
            paste0("![](figures/Pattern__", i, "_.png)")
        }), collapse=""))
    BODY4 <- paste0(BODY4, "\n")
}

# 5. LR-Pattern
.BODY5 <- paste0("\n\n# LR-pair Patterns\n\n",
"```{r}\n", # Top
"lr <- metadata(sce)$sctensor$lrpair\n",
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
   "![](figures/GeneHypergraph.png){ width=100% }\n",
    "[Details of Ligand Gene-centric Overview](ligand.html)\n\n",
    "[Details of Receptor Gene-centric Overview](receptor.html)\n"
    )

# 8. (Ligand, Receptor, LR-pair)-Pattern
.BODY8 <- function(selected, rmdfiles, index, corevalue){
    if(length(selected) != 0){
        htmlfiles <- gsub("Rmd", "html", rmdfiles)
        BODY8 <- sapply(selected, function(x){
            i <- selected[x]
            paste0("\n\n## (", paste(index[i, seq_len(3)],
                collapse=","),
                ") Pattern : (", round(corevalue[i], 2), " %)\n",
                "[Details of (", paste(index[i, seq_len(3)],
                collapse=","),
                ") Pattern", "](", htmlfiles[i], ")\n")
            })
        BODY8 <- paste(BODY8, collapse = "\n")
        BODY8 <- paste0("# (Ligand-Cell, Receptor-Cell, LR-pair)",
            " Patterns\n\n", BODY8)
    }else{
        BODY8 <- "# (Ligand-Cell, Receptor-Cell, LR-pair) Patterns\n\n"
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

.sapply_pb <- function(X, FUN, ...)
{
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- sapply(X, wrapper, ...)
  close(pb)
  res
}

.FTT <- function(x){
    sqrt(x) + sqrt(x + 1)
}

.eachVecLR <- function(x, e){
    index <- e$index
    sce <- e$sce
    .HCLUST <- e$.HCLUST
    .OUTLIERS <- e$.OUTLIERS
    top <- e$top
    spc <- e$spc
    .sapply_pb <- e$.sapply_pb
    GeneInfo <- e$GeneInfo
    temp <- e$temp
    .smallTwoDplot <- e$.smallTwoDplot
    input <- e$input
    twoD <- e$twoD
    .hyperLinks <- e$.hyperLinks
    LR <- e$LR
    .eachVecLR <- e$.eachVecLR
    .eachRender <- e$.eachRender
    .XYZ_HEADER1 <- e$.XYZ_HEADER1
    .XYZ_HEADER2 <- e$.XYZ_HEADER2

    # Each LR-Pattern Vector
    vecLR <- metadata(sce)$sctensor$lrpair[x, ]
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
        TARGET <- TARGET[1:TOP]
    }else{
        TOP <- "full"
    }
    # Eigen Value
    Value <- metadata(sce)$sctensor$lrpair[x, TARGET]
    Percentage <- Value / sum(metadata(sce)$sctensor$lrpair[x, ]) * 100
    # Hyper Link (Too Heavy)
    cat("Hyper-links are embedded...\n")
    GeneName <- GeneInfo$GeneName
    LINKS <- .sapply_pb(seq_along(TARGET), function(xx){
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

            # Plot (Ligand/Receptor)
            Ligandfile <- paste0(temp, "/figures/Ligand_", LigandGeneID, ".png")
            Receptorfile <- paste0(temp, "/figures/Receptor_", ReceptorGeneID, ".png")
            png(filename=Ligandfile, width=1000, height=1000)
                .smallTwoDplot(input, LigandGeneID, LigandGeneName, twoD, "reds")
            dev.off()
            png(filename=Receptorfile, width=1000, height=1000)
            .smallTwoDplot(input, ReceptorGeneID, ReceptorGeneName, twoD, "blues")
            dev.off()

            # Return Hyper links
            .hyperLinks(Ranking, LigandGeneID,
            ReceptorGeneID, LR, Value[xx],
            Percentage[xx],
            spc, GeneInfo, PvalueLR[xx],
            QvalueLR[xx])
        })
    LINKS <- paste(LINKS, collapse="")

    # Output object
    list(
        ClusterLR=ClusterLR,
        PvalueLR=PvalueLR,
        QvalueLR=QvalueLR,
        TARGET=TARGET,
        TOP=TOP,
        Value=Value,
        Percentage=Percentage,
        LINKS=LINKS
    )
}

.eachRender <- function(x, e, SelectedLR){
    index <- e$index
    temp <- e$temp
    out.vecLR <- e$out.vecLR
    .XYZ_HEADER1 <- e$.XYZ_HEADER1
    .XYZ_HEADER2 <- e$.XYZ_HEADER2

    indexLR <- index[x, "Mode3"]
    TARGET <- out.vecLR[, paste0("pattern", indexLR)]$TARGET
    LINKS <- out.vecLR[, paste0("pattern", indexLR)]$LINKS

    # Bottom part of Rmarkdown
    XYZ_BOTTOM <- paste(
        c(.XYZ_HEADER2(indexLR, length(TARGET)),
        LINKS),
        collapse="")

    # Each (x,y,z)-rmdfile
    RMDFILE <- paste0(c("pattern", index[x, seq_len(3)]), collapse="_")
    RMDFILE <- paste0(RMDFILE, ".Rmd")
    cat("\n")
    cat(paste0(RMDFILE, " is created...\n"))
    sink(file = paste0(temp, "/", RMDFILE))
    cat(paste(
        c(.XYZ_HEADER1(index, x), XYZ_BOTTOM),
        collapse=""))
    sink()

    # Rendering
    cat(paste0(RMDFILE, " is compiled to ",
        gsub(".Rmd", ".html", RMDFILE), "\n"))
    render(paste0(temp, "/", RMDFILE), quiet=TRUE)
}

.palf <- colorRampPalette(c("#4b61ba", "gray", "#a87963", "red"))

.shrink2 <- function(x){
    block <- strsplit(x , " ")[[1]]
    l <- length(block)
    nc <- sapply(block, nchar)
    if(l >= 2){
        out <- block[1]
        for(i in 2:l){
            end <- length(out)
            tmp <- block[i]
            if(nchar(out[end])+nchar(tmp) <= 25){
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
        c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D")
    }else if(col == "blues"){
        c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B")
    }else if(col == "greens"){
        c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B")
    }else if(col == "many"){
        c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC")
    }else{
        stop("Wrong col is specified in .setColor")
    }
}

.geneHyperGraphPlot <- function(out.vecLR, GeneInfo, temp){
    # Setting
    convertGeneName <- function(geneid, geneInfo){
        genename <- geneInfo$GeneName[
            which(geneInfo$GeneName$entrezgene == geneid),
            "external_gene_name"][1]
        if(length(genename) == 0){
            genename = geneid
        }
        genename
    }

    # Node
    nodes <- sapply(seq_len(ncol(out.vecLR)), function(x){
        names(out.vecLR["TARGET", x][[1]])
    }, simplify=FALSE)
    Lnodes <- lapply(nodes, function(x){
        sapply(x, function(xx){
            strsplit(xx, "_")[[1]][1]
        })
    })
    Rnodes <-lapply(nodes, function(x){
        sapply(x, function(xx){
            strsplit(xx, "_")[[1]][2]
        })
    })
    Lnodes <- lapply(Lnodes, function(x){
            sapply(x, function(xx){
                convertGeneName(xx, GeneInfo)
        })
    })
    Rnodes <- lapply(Rnodes, function(x){
            sapply(x, function(xx){
                convertGeneName(xx, GeneInfo)
        })
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
    freqLnodes <- sapply(uniqueLnodes, function(x){
        length(which(unlist(Lnodes) == x))
        })
    freqRnodes <- sapply(uniqueRnodes, function(x){
        length(which(unlist(Rnodes) == x))
        })
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
    edge.cols <- unlist(sapply(seq_len(ncol(out.vecLR)), function(x){
            rep(cols[x], length(out.vecLR["TARGET", x][[1]]))
        }, simplify=FALSE))

    # Setting
    V(g)$size <- freq
    E(g)$color <- edge.cols
    E(g)$width <- 0.5
    l <- layout_with_dh(g)

    # All Pattern
    png(filename=paste0(temp, "/figures/GeneHypergraph.png"),
        width=2500, height=2500)
    plot.igraph(g, layout=l)
    legend("topleft",
        legend=c("ligand", "receptor",
            colnames(out.vecLR)),
        col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5),
            cols[seq_len(ncol(out.vecLR))]),
        pch=16)
    dev.off()

    # Each Pattern
    sapply(seq_len(ncol(out.vecLR)), function(x){
        tmp_edgecolor <- edge.cols
        tmp_edgecolor[which(tmp_edgecolor != cols[x])] <- rgb(0,0,0,0.1)
        tmp_nodecolor <- V(g)$color
        grayout <- setdiff(
            setdiff(
                names(V(g)),
                Lnodes[[x]]
                ), Rnodes[[x]]
            )
        target <- unlist(sapply(grayout, function(xx){
            which(names(V(g)) == xx)
        }))
        tmp_nodecolor[target] <- rgb(0,0,0,0.1)

        # Plot
        png(filename=paste0(
            temp, "/figures/GeneHypergraph_",
            gsub("pattern", "", colnames(out.vecLR)[x]),
            ".png"),
            width=2500, height=2500)
        plot.igraph(g,
            vertex.color=tmp_nodecolor,
            edge.color=tmp_edgecolor, layout=l)
        legend("topleft",
            legend=c("ligand", "receptor",
                colnames(out.vecLR)[x]),
            col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5),
                cols[x]),
            pch=16)
        dev.off()
    })
}

.LIGAND_HEADER <- paste0(
    "# <font color='#1881c2'>Details of Ligand Gene-centric Overview</font>\n\n",
    "![](figures/GeneHypergraph.png){ width=100% }\n\n",
    "|Rank|Frequency|Ligand Gene|Candidate Partner Receptor Genes|Related CCIs|\n",
    "|----|----|----|----|----|\n"
)

.RECEPTOR_HEADER <- paste0(
    "# <font color='#1881c2'>Details of Receptor Gene-centric Overview</font>\n\n",
    "![](figures/GeneHypergraph.png){ width=100% }\n\n",
    "|Rank|Frequency|Receptor Gene|Candidate Partner Ligand Genes|Related CCIs|\n",
    "|----|----|----|----|----|\n"
)

.LIGAND_BODY <- function(out.vecLR, GeneInfo, index, selected){
    # Setting
    convertGeneName <- function(geneid, geneInfo){
        genename <- geneInfo$GeneName[
            which(geneInfo$GeneName$entrezgene == geneid),
            "external_gene_name"][1]
        if(length(genename) == 0){
            genename = geneid
        }
        genename
    }

    # Node
    nodes <- sapply(seq_len(ncol(out.vecLR)), function(x){
        names(out.vecLR["TARGET", x][[1]])
    }, simplify=FALSE)
    Lnodes <- lapply(nodes, function(x){
        sapply(x, function(xx){
            strsplit(xx, "_")[[1]][1]
        })
    })
    Rnodes <-lapply(nodes, function(x){
        sapply(x, function(xx){
            strsplit(xx, "_")[[1]][2]
        })
    })
    Lnodes <- lapply(Lnodes, function(x){
            sapply(x, function(xx){
                convertGeneName(xx, GeneInfo)
        })
    })
    Rnodes <- lapply(Rnodes, function(x){
            sapply(x, function(xx){
                convertGeneName(xx, GeneInfo)
        })
    })
    uniqueLnodes <-  unique(unlist(Lnodes))
    uniqueRnodes <-  unique(unlist(Rnodes))

    Freq <- sapply(uniqueLnodes, function(x){
            length(which(unlist(Lnodes) == x))
        })
    Freq <- sort(Freq, decreasing=TRUE)
    Rank <- rank(-Freq)
    Ligand <- names(Freq)
    Receptor <- sapply(Ligand, function(x){
        unique(unlist(sapply(seq_len(length(Lnodes)), function(xx){
            Rnodes[[xx]][which(x == Lnodes[[xx]])]
        })))
    })
    Receptor <- unlist(lapply(Receptor, function(x){
        paste(x, collapse=", ")
        }))
    CCI <- sapply(Ligand, function(x){
        hit <- sapply(seq_len(ncol(out.vecLR)), function(xx){
            length(which(Lnodes[[xx]] == x))
            })
        vec <- seq_len(ncol(out.vecLR))[which(hit != 0)]
        target <- unlist(sapply(vec, function(xx){
                p <- gsub("pattern", "", colnames(out.vecLR)[xx])
                which(p == index[selected, "Mode3"])
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
            paste(sapply(seq_len(nrow(out)), function(xx){
                paste0("[See the details of (",
                    paste(out[xx,], collapse=","),
                    ")-Pattern](",
                    "pattern_",
                    paste(out[xx,], collapse="_"),
                    ".html)")
            }), collapse=" ")
        }
    })
    paste0(
    "|", Rank,
    "|", Freq,
    "|", Ligand,
    "|", Receptor,
    "|", CCI, "|\n", collapse="")
}

.RECEPTOR_BODY <- function(out.vecLR, GeneInfo, index, selected){
    # Setting
    convertGeneName <- function(geneid, geneInfo){
        genename <- geneInfo$GeneName[
            which(geneInfo$GeneName$entrezgene == geneid),
            "external_gene_name"][1]
        if(length(genename) == 0){
            genename = geneid
        }
        genename
    }

    # Node
    nodes <- sapply(seq_len(ncol(out.vecLR)), function(x){
        names(out.vecLR["TARGET", x][[1]])
    }, simplify=FALSE)
    Lnodes <- lapply(nodes, function(x){
        sapply(x, function(xx){
            strsplit(xx, "_")[[1]][1]
        })
    })
    Rnodes <-lapply(nodes, function(x){
        sapply(x, function(xx){
            strsplit(xx, "_")[[1]][2]
        })
    })
    Lnodes <- lapply(Lnodes, function(x){
            sapply(x, function(xx){
                convertGeneName(xx, GeneInfo)
        })
    })
    Rnodes <- lapply(Rnodes, function(x){
            sapply(x, function(xx){
                convertGeneName(xx, GeneInfo)
        })
    })
    uniqueLnodes <-  unique(unlist(Lnodes))
    uniqueRnodes <-  unique(unlist(Rnodes))

    Freq <- sapply(uniqueRnodes, function(x){
            length(which(unlist(Rnodes) == x))
        })
    Freq <- sort(Freq, decreasing=TRUE)
    Rank <- rank(-Freq)
    Receptor <- names(Freq)
    Ligand <- sapply(Receptor, function(x){
        unique(unlist(sapply(seq_len(length(Rnodes)), function(xx){
            Lnodes[[xx]][which(x == Rnodes[[xx]])]
        })))
    })
    Ligand <- unlist(lapply(Ligand, function(x){
        paste(x, collapse=", ")
        }))
    CCI <- sapply(Receptor, function(x){
        hit <- sapply(seq_len(ncol(out.vecLR)), function(xx){
            length(which(Rnodes[[xx]] == x))
            })
        vec <- seq_len(ncol(out.vecLR))[which(hit != 0)]
        target <- unlist(sapply(vec, function(xx){
                p <- gsub("pattern", "", colnames(out.vecLR)[xx])
                which(p == index[selected, "Mode3"])
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
            paste(sapply(seq_len(nrow(out)), function(xx){
                paste0("[See the details of (",
                    paste(out[xx,], collapse=","),
                    ")-Pattern](",
                    "pattern_",
                    paste(out[xx,], collapse="_"),
                    ".html)")
            }), collapse=" ")
        }
    })
    paste0(
    "|", Rank,
    "|", Freq,
    "|", Receptor,
    "|", Ligand,
    "|", CCI, "|\n", collapse="")
}