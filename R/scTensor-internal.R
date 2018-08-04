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
        ppiout <- .cellCellReconst(A1, A2, A3, index, type="tensor",
            names=list(colnames(A1), colnames(A2), colnames(A3)),
            collapse=TRUE)
        return(list(index=index, ligand=A1,
            receptor=A2, ppi=A3, cellcellpattern=ppiout))
    }else{
        tnsr_cc <- modeSum(tnsr, m=3, drop=TRUE)@data
        dimnames(tnsr_cc) <- list(unique(names(celltypes)),
            unique(names(celltypes)))
        return(list(cellcellpattern=tnsr_cc, cellcellppipattern=tnsr))
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
    return(list(cellcellpattern=tnsr_cc, cellcellppipattern=tnsr_bin))
}

.lrPlot <- function(outobj, twoDplot=NULL, vertex.size=18, xleft=1.75,
    ybottom=-0.5, xright=1.85, ytop=0.5, label="", emph=NULL){
    # Number of Patterns
    numLPattern <- nrow(outobj$ligand)
    numRPattern <- nrow(outobj$receptor)

    #
    # Step.1 : Background Network
    #
    edgewd <- matrix(0, nrow=numLPattern*numRPattern, ncol=3)
    counter <- 1
    for(i in seq_len(numLPattern)){
        for(j in seq_len(numRPattern)){
            targetL <- which(outobj$index[, "Mode1"] == i)
            targetR <- which(outobj$index[, "Mode2"] == j)
            edgewd[counter, ] <- c(i, j,
                sum(outobj$index[intersect(targetL, targetR), 4]))
            counter <- counter + 1
        }
    }
    colnames(edgewd) <- c("L", "R", "Strength")

    # Node name (Top<Ligand> and Bottom<Receptor>)
    nodesSetTop <- paste0("L", seq_len(numLPattern))
    nodesSetBottom <- paste0("R", seq_len(numRPattern))

    # Empty Graph
    g <- graph.empty()

    # Add nodes
    g <- add.vertices(g, nv=length(nodesSetTop),
        attr=list(name=nodesSetTop, type=rep(TRUE, length(nodesSetTop))))
    g <- add.vertices(g, nv=length(nodesSetBottom),
        attr=list(name=nodesSetBottom, type=rep(TRUE, length(nodesSetBottom))))

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
        palfunc=colorRampPalette(brewer.pal(9, "Greens"), alpha=TRUE))
    if(!is.null(emph)){
        target <- numLPattern * (emph[1] - 1) + emph[2]
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
        palfunc=colorRampPalette(brewer.pal(9, "Greens"), alpha=TRUE)),
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
            col.ligand <- brewer.pal(9, "Reds")
            # Constant
            LOMA_1 = 48.85
            LOMA_2 = 42
            LOMA_3 = 4
            for(i in seq_len(numLPattern)){
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
            }

            #
            # Step.3 : Receptor Plot
            #
            # Color
            col.receptor <- brewer.pal(9, "Blues")
            # Constant
            ROMA_1 = 4
            ROMA_2 = 42
            ROMA_3 = 48.85
            for(i in seq_len(numRPattern)){
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
            }
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
        palfunc=colorRampPalette(brewer.pal(9, color), alpha=TRUE))
    # Plot
    par(ps=25)
    par(oma=c(2,2,2,2))
    plot(twoD, col=label, pch=16, cex=3, bty="n", xaxt="n", yaxt="n",
        xlab="", ylab="", main=genename)
    # Gradient
    xleft=19
    ybottom=-5
    xright=20
    ytop=10
    gradient.rect(xleft, ybottom, xright, ytop, col=smoothPalette(sort(exp),
        palfunc=colorRampPalette(brewer.pal(9, color),
        alpha=TRUE)), gradient="y")
    text(17, ybottom+(ytop-ybottom)*0/4, round(quantile(exp)[1]))
    text(17, ybottom+(ytop-ybottom)*1/4, round(quantile(exp)[2]))
    text(17, ybottom+(ytop-ybottom)*2/4, round(quantile(exp)[3]))
    text(17, ybottom+(ytop-ybottom)*3/4, round(quantile(exp)[4]))
    text(17, ybottom+(ytop-ybottom)*4/4, round(quantile(exp)[5]))
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
    if(!is.null(ens)){
        targetGeneID <- unique(c(LR$GENEID_L, LR$GENEID_R))
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
    }else{
        GeneName <- NULL
        Description <- NULL
        GO <- NULL
        Reactome <- NULL
        MeSH <- NULL
    }
    # Output
    list(GeneName=GeneName, Description=Description, GO=GO,
        Reactome=Reactome, MeSH=MeSH)
}

## helper for vector dividing
.div <- function(x, d=1) {
    delta <- ceiling(length(x) / d)
    y <- lapply(seq_len(d), function(i){
        as.vector(na.omit(x[((i-1)*delta+1):(i*delta)]))
    })
    return(y)
}

.hyperLinks <- function(Ranking, XYZ, LigandGeneID, ReceptorGeneID,
    LR, Value, Percentage, j, spc, GeneInfo, pvalue, qvalue){
    # Ranking
    XYZ <- paste0(XYZ, "|", Ranking, "|")
    # Gene Name (Ligand)
    GeneName_L <- GeneInfo$GeneName[
        which(GeneInfo$GeneName$entrezgene == LigandGeneID),
        "external_gene_name"][1]
    if(length(GeneName_L) == 0){
        GeneName_L = LigandGeneID
    }
    # Gene Name (Receptor)
    GeneName_R <- GeneInfo$GeneName[
        which(GeneInfo$GeneName$entrezgene == ReceptorGeneID),
        "external_gene_name"][1]
    if(length(GeneName_R) == 0){
        GeneName_R = ReceptorGeneID
    }
    # Description (Ligand)
    Description_L <- GeneInfo$Description[
        which(GeneInfo$Description$entrezgene ==
            LigandGeneID), "description"][1]
    if(length(Description_L) == 0){
        Description_L = ""
    }
    # Description (Receptor)
    Description_R <- GeneInfo$Description[
        which(GeneInfo$Description$entrezgene ==
            ReceptorGeneID), "description"][1]
    if(length(Description_R) == 0){
        Description_R = ""
    }
    # Gene Ontology (Ligand)
    GO_L <- unique(unlist(GeneInfo$GO[
        which(GeneInfo$GO$entrezgene ==
            LigandGeneID), "go_id"]))
    GO_L <- GO_L[which(GO_L != "")]
    GO_L <- gsub(":", "%3A", GO_L)
    if(length(GO_L) != 0){
        tmp_GO_L <- list()
        GO_L_loc <- .div(seq_along(GO_L), ceiling(length(GO_L) / 100))
        for(i in seq_along(GO_L_loc)){
            mi <- min(GO_L_loc[[i]])
            ma <- max(GO_L_loc[[i]])
            tmp_GO_L[[i]] <- paste0("[", mi, "-", ma, "](",
                "http://amigo.geneontology.org/",
                "goose?query=SELECT+*+FROM+term+WHERE+acc%3D%27",
                paste0(GO_L[mi:ma], collapse="%27+OR+acc%3D%27"),
                "%27%3B&mirror=bbop)")
        }
        GO_L = paste(unlist(tmp_GO_L), collapse=" ")
    }else{
        GO_L = ""
    }
    # Gene Ontology (Receptor)
    GO_R <- unique(unlist(GeneInfo$GO[
        which(GeneInfo$GO$entrezgene ==
            ReceptorGeneID), "go_id"]))
    GO_R <- GO_R[which(GO_R != "")]
    GO_R <- gsub(":", "%3A", GO_R)
    if(length(GO_R) != 0){
        tmp_GO_R <- list()
        GO_R_loc <- .div(seq_along(GO_R),
            ceiling(length(GO_R) / 100))
        for(i in seq_along(GO_R_loc)){
            mi <- min(GO_R_loc[[i]])
            ma <- max(GO_R_loc[[i]])
            tmp_GO_R[[i]] <- paste0("[", mi, "-", ma, "](",
                "http://amigo.geneontology.org/",
                "goose?query=SELECT+*+FROM+term+WHERE+acc%3D%27",
                paste0(GO_R[mi:ma],
                collapse="%27+OR+acc%3D%27"),
                "%27%3B&mirror=bbop)")
        }
        GO_R = paste(unlist(tmp_GO_R), collapse=" ")
    }else{
        GO_R = ""
    }
    # Reactome (Ligand)
    Reactome_L <- unique(unlist(GeneInfo$Reactome[
        which(GeneInfo$Reactome$gene_id ==
        LigandGeneID), "DB_ID"]))
    Reactome_L <- Reactome_L[which(Reactome_L != "")]
    if(length(Reactome_L) != 0){
        tmp_Reactome_L <- list()
        Reactome_L_loc <- .div(seq_along(Reactome_L),
            ceiling(length(Reactome_L) / 100))
        for(i in seq_along(Reactome_L_loc)){
            mi <- min(Reactome_L_loc[[i]])
            ma <- max(Reactome_L_loc[[i]])
            tmp_Reactome_L[[i]] <- paste0("[", mi, "-", ma, "](",
                "https://reactome.org/content/query?q=",
                paste0(Reactome_L[mi:ma], collapse="+"),
                "&types=Pathway&cluster=true)")
        }
        Reactome_L = paste(unlist(tmp_Reactome_L), collapse=" ")
    }else{
        Reactome_L = ""
    }
    # Reactome (Receptor)
    Reactome_R <- unique(unlist(GeneInfo$Reactome[
        which(GeneInfo$Reactome$gene_id ==
        ReceptorGeneID), "DB_ID"]))
    Reactome_R <- Reactome_R[which(Reactome_R != "")]
    if(length(Reactome_R) != 0){
        tmp_Reactome_R <- list()
        Reactome_R_loc <- .div(seq_along(Reactome_R),
            ceiling(length(Reactome_R) / 100))
        for(i in seq_along(Reactome_R_loc)){
            mi <- min(Reactome_R_loc[[i]])
            ma <- max(Reactome_R_loc[[i]])
            tmp_Reactome_R[[i]] <- paste0("[", mi, "-", ma, "](",
                "https://reactome.org/content/query?q=",
                paste0(Reactome_R[mi:ma], collapse="+"),
                "&types=Pathway&cluster=true)")
        }
        Reactome_R = paste(unlist(tmp_Reactome_R), collapse=" ")
    }else{
        Reactome_R = ""
    }
    # MeSH (Ligand)
    MeSH_L <- GeneInfo$MeSH[which(GeneInfo$MeSH$GENEID == LigandGeneID),
        "MESHID"]
    MeSH_L <- MeSH_L[which(MeSH_L != "")]
    if(length(MeSH_L) != 0){
        tmp_MeSH_L <- list()
        MeSH_L_loc <- .div(seq_along(MeSH_L),
            ceiling(length(MeSH_L) / 100))
        for(i in seq_along(MeSH_L_loc)){
            mi <- min(MeSH_L_loc[[i]])
            ma <- max(MeSH_L_loc[[i]])
            tmp_MeSH_L[[i]] <- paste0("[", mi, "-", ma, "](",
                "https://www.ncbi.nlm.nih.gov/mesh?term=",
                paste0(MeSH_L[mi:ma], collapse="%20OR%20"), ")")
        }
        MeSH_L = paste(unlist(tmp_MeSH_L), collapse=" ")
    }else{
        MeSH_L = ""
    }
    # MeSH (Receptor)
    MeSH_R <- GeneInfo$MeSH[which(GeneInfo$MeSH$GENEID == ReceptorGeneID),
        "MESHID"]
    MeSH_R <- MeSH_R[which(MeSH_R != "")]
    if(length(MeSH_R) != 0){
        tmp_MeSH_R <- list()
        MeSH_R_loc <- .div(seq_along(MeSH_R),
            ceiling(length(MeSH_R) / 100))
        for(i in seq_along(MeSH_R_loc)){
            mi <- min(MeSH_R_loc[[i]])
            ma <- max(MeSH_R_loc[[i]])
            tmp_MeSH_R[[i]] <- paste0("[", mi, "-", ma, "](",
                "https://www.ncbi.nlm.nih.gov/mesh?term=",
                paste0(MeSH_R[mi:ma], collapse="%20OR%20"), ")")
        }
        MeSH_R = paste(unlist(tmp_MeSH_R), collapse=" ")
    }else{
        MeSH_R = ""
    }
    # PubMed (L and R)
    targetPubMed <- intersect(which(LR$GENEID_L == LigandGeneID),
        which(LR$GENEID_R == ReceptorGeneID))
    PubMed <- LR$SOURCEID[targetPubMed]
    if(length(targetPubMed) != 0 && !is.null(LR$SOURCEID[targetPubMed])){
        PubMed <- unique(strsplit(LR$SOURCEID[targetPubMed], "\\|")[[1]])
        tmp_PubMed <- list()
        PubMed_loc <- .div(seq_along(PubMed),
            ceiling(length(PubMed) / 100))
        for(i in seq_along(PubMed_loc)){
            mi <- min(PubMed_loc[[i]])
            ma <- max(PubMed_loc[[i]])
            tmp_PubMed[[i]] <- paste0("[", mi, "-", ma, "](",
                "https://www.ncbi.nlm.nih.gov/pubmed/?term=",
                paste0(PubMed[mi:ma], collapse="%20OR%20"), ")")
        }
        PubMed = paste(unlist(tmp_PubMed), collapse=" ")
    }else{
        PubMed = ""
    }

    # Embedding
    if(spc != "Pab"){
        XYZ <- paste0(XYZ,
            "[", GeneName_L,
            "](https://www.ncbi.nlm.nih.gov/gene/",
            LigandGeneID, ")<br>",
            "Description: ", Description_L, "<br>",
            "GO: ", GO_L, "<br>",
            "Reactome: ", Reactome_L, "<br>",
            "MeSH: ", MeSH_L,
            "|",
            "[", GeneName_R,
            "](https://www.ncbi.nlm.nih.gov/gene/",
            ReceptorGeneID, ")<br>",
            "Description: ", Description_R, "<br>",
            "GO: ", GO_R, "<br>",
            "Reactome: ", Reactome_R, "<br>",
            "MeSH: ", MeSH_R,
            "|"
        )
    }else if(spc %in% c("Cel", "Dme", "Pab", "Ssc")){
        XYZ <- paste0(XYZ,
            "[", GeneName_L,
            "](https://www.ncbi.nlm.nih.gov/gene/",
            LigandGeneID, ")<br>",
            "Description: ", Description_L, "<br>",
            "GO: ", GO_L, "<br>",
            "MeSH: ", MeSH_L,
            "|",
            "[", GeneName_R,
            "](https://www.ncbi.nlm.nih.gov/gene/",
            ReceptorGeneID, ")<br>",
            "Description: ", Description_R, "<br>",
            "GO: ", GO_R, "<br>",
            "MeSH: ", MeSH_R,
            "|"
        )
    }else if(spc == "Ath"){
        XYZ <- paste0(XYZ,
            "[", LigandGeneID,
            "](https://www.ncbi.nlm.nih.gov/gene/",
            LigandGeneID, ")<br>",
            "MeSH: ", MeSH_L,
            "|",
            "[", ReceptorGeneID,
            "](https://www.ncbi.nlm.nih.gov/gene/",
            ReceptorGeneID, ")<br>",
            "MeSH: ", MeSH_R,
            "|"
        )
    }else{
        XYZ <- paste0(XYZ,
            "[", LigandGeneID,
            "](https://www.ncbi.nlm.nih.gov/gene/",
            LigandGeneID, ")",
            "|",
            "[", ReceptorGeneID,
            "](https://www.ncbi.nlm.nih.gov/gene/",
            ReceptorGeneID, ")",
            "|"
        )
    }
    # Figure, Core Tensor value, PubMed
    XYZ <- paste0(XYZ,
        "![](figures/Ligand_", LigandGeneID, ".png)", "|",
        "![](figures/Receptor_", ReceptorGeneID, ".png)", "|",
        round(Value[j], 3), " (", round(Percentage[j], 3), "%)",
        "|", round(pvalue, 3),
        "|", round(qvalue, 3),
        "|", PubMed, "|\n")
    XYZ
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
