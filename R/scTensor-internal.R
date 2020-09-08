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

.importAssays <- function(sce, assayNames){
    if(assayNames %in% names(assays(sce))){
        eval(parse(text=paste0("assays(sce)$", assayNames)))
    }else{
        stop("Please specify the valid assayNames (cf. names(assays(sce)))")
    }
}

'%ni%' <- Negate('%in%')
