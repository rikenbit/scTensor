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
