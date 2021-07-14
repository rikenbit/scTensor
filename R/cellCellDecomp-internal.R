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
    dimnames(tnsr) <- list(unique(names(celltypes)),
        unique(names(celltypes)), Pair.name)
    list(tnsr=tnsr, pairname=Pair.name)
}

.cellCellDecomp.Third_3 <- function(input, LR, celltypes, ranks, rank, centering,
        mergeas, outerfunc, comb, num.sampling, num.perm, decomp, thr1, thr2, L1_A, L2_A, verbose){
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
    if(decomp){
        message(paste0(paste(dim(tnsr), collapse=" * "), " Tensor is created"))
        out <- try(MultiCX(Y=tnsr, modes=1:2, thr=0.9))
        if(is(out)[1] == "try-error"){
            out <- MultiCX(Y=tnsr, modes=1:2, thr=0.9)
        }
        A1 <- t(out$C[[1]]$C)
        A2 <- t(out$C[[2]]$C)
        rownames(A1) <- paste0("Dim", seq_len(nrow(A1)))
        colnames(A1) <- unique(names(celltypes))
        rownames(A2) <- paste0("Dim", seq_len(nrow(A2)))
        colnames(A2) <- unique(names(celltypes))
        dimnames(out$U@data) <- list(
            Lpattern=paste0("Dim", seq_len(nrow(A1))),
            Rpattern=paste0("Dim", seq_len(nrow(A2))),
            LR=Pair.name
        )
        if (nrow(A1) * nrow(A2) == 1){
            core <- t(sum(out$U@data))
        }else{
            core <- modeSum(out$U, m=3, drop=TRUE)@data
            if(nrow(A1) == 1 || nrow(A2) == 1){
                dim(core) <- c(nrow(A1), nrow(A2))
            }
        }
        index <- .core2table_2(core)
        if(is.vector(index)){
            index <- t(index)
        }
        return(list(index=index, ligand=A1,
            receptor=A2, lrpair=out$U,
            recerror=out$RecError, relchange=out$RelChange))
    }else{
        tnsr_cc <- modeSum(tnsr, m=3, drop=TRUE)@data
        dimnames(tnsr_cc) <- list(unique(names(celltypes)),
            unique(names(celltypes)))
        return(list(cellcellpattern=tnsr_cc, cellcelllrpairpattern=tnsr))
    }
}

.cellCellDecomp.Third_2 <- function(input, LR, celltypes, ranks, rank, centering,
        mergeas, outerfunc, comb, num.sampling, num.perm, decomp, thr1, thr2,
        L1_A, L2_A, verbose){
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
        out <- try(NTD(X=tnsr, L1_A=L1_A, L2_A=L2_A,
            rank=ranks, modes=1:2, num.iter=30,
            verbose=verbose, algorithm="Frobenius"))
        if(is(out)[1] == "try-error"){
            out <- NTD(X=tnsr, L1_A=L1_A, L2_A=L2_A,
                rank=ranks, modes=1:2, num.iter=30, algorithm="Frobenius")
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
        mergeas, outerfunc, comb, num.sampling, num.perm, decomp, thr1, thr2,
        L1_A, L2_A, verbose){
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
        out <- try(NTD(X=tnsr, L1_A=L1_A, L2_A=L2_A,
            rank=ranks, num.iter=30,
            verbose=verbose, algorithm="Frobenius"))
        if(is(out)[1] == "try-error"){
            out <- try(NTD(X=tnsr, L1_A=L1_A, L2_A=L2_A,
                rank=ranks, num.iter=30,
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
    decomp, thr1, thr2, L1_A, L2_A, verbose){
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
    decomp, thr1, thr2, L1_A, L2_A, verbose){
    out <- .cellType(input, celltypes)
    out <- cor(out, method="pearson")
    out <- as.matrix(out)
    return(list(cellcellpattern=out))
}

.cellCellDecomp.Spearman <- function(input, LR, celltypes, ranks, rank,
    centering, mergeas, outerfunc, comb, num.sampling, num.perm,
    decomp, thr1, thr2, L1_A, L2_A, verbose){
    out <- .cellType(input, celltypes)
    out <- cor(out, method="spearman")
    out <- as.matrix(out)
    return(list(cellcellpattern=out))
}

.cellCellDecomp.Distance <- function(input, LR, celltypes, ranks, rank,
    centering, mergeas, outerfunc, comb, num.sampling, num.perm,
    decomp, thr1, thr2, L1_A, L2_A, verbose){
    out <- .cellType(input, celltypes)
    out <- 1 / dist(t(out))
    out <- as.matrix(out)
    diag(out) <- 1
    return(list(cellcellpattern=out))
}

.cellCellDecomp.Pearson.LR <- function(input, LR, celltypes, ranks, rank,
    centering, mergeas, outerfunc, comb, num.sampling, num.perm,
    decomp, thr1, thr2, L1_A, L2_A, verbose){
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
    decomp, thr1, thr2, L1_A, L2_A, verbose){
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
    decomp, thr1, thr2, L1_A, L2_A, verbose){
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
    decomp, thr1, thr2, L1_A, L2_A, verbose){
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
    decomp, thr1, thr2, L1_A, L2_A, verbose){
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
    decomp, thr1, thr2, L1_A, L2_A, verbose){
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
    decomp, thr1, thr2, L1_A, L2_A, verbose){
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
    # 3-Order, 2-mode decomposition : MultilinearCX
    "cx" = .cellCellDecomp.Third_3,
    # 3-Order, 2-mode decomposition : NTD2
    "ntd2" = .cellCellDecomp.Third_2,
    # 3-Order : NTD
    "ntd" = .cellCellDecomp.Third,
    # 2-Order : NMF
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