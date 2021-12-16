.geneInformation_taxid <- function(sce, ah, taxid, LR){
    targetGeneID <- as.character(unique(c(LR$GENEID_L, LR$GENEID_R)))
    # Gene symbol
    if("SYMBOL" %in% AnnotationDbi::columns(ah) && !is.null(ah)){
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
        target <- as.character(unique(GeneName[, "ENTREZID"]))
        GeneName <- data.frame(
            ENTREZID=target,
            SYMBOL=as.character(g$matching[target]),
            stringsAsFactors = FALSE)
    }else{
        GeneName <- NULL
    }

    # GENENAME
    if("GENENAME" %in%  AnnotationDbi::columns(ah) && !is.null(ah)){
        message("Related gene descriptions are retrieved from AnnotationHub...")
        Description <- AnnotationDbi::select(ah, columns=c("GENENAME", "ENTREZID"),
            keytype="ENTREZID", keys=targetGeneID)
        Description <- Description[, c("ENTREZID", "GENENAME")]
    }else{
        Description <- NULL
    }

    # GO
    if("GO" %in%  AnnotationDbi::columns(ah) && !is.null(ah)){
        message("Related GO IDs are retrieved from AnnotationHub...")
        GO <- AnnotationDbi::select(ah, columns=c("GO", "ENTREZID"),
            keytype="ENTREZID", keys=targetGeneID)
        GO <- GO[, c("ENTREZID", "GO")]
    }else{
        GO <- NULL
    }

    # ENSG
    if("ENSEMBL" %in%  AnnotationDbi::columns(ah) && !is.null(ah)){
        message("Related Ensembl Gene IDs are retrieved from AnnotationHub...")
        ENSG <- AnnotationDbi::select(ah, columns=c("ENSEMBL", "ENTREZID"),
            keytype="ENTREZID", keys=targetGeneID)
        ENSG <- ENSG[, c("ENTREZID", "ENSEMBL")]
    }else{
        ENSG <- NULL
    }

    # ENSP
    if("ENSEMBLPROT" %in%  AnnotationDbi::columns(ah) && !is.null(ah)){
        message("Related Ensembl Protein IDs are retrieved from AnnotationHub...")
        ENSP <- AnnotationDbi::select(ah, columns=c("ENSEMBLPROT", "ENTREZID"),
            keytype="ENTREZID", keys=targetGeneID)
        ENSP <- ENSP[, c("ENTREZID", "ENSEMBLPROT")]
    }else{
        ENSP <- NULL
    }

    # UniProtKB
    if("UNIPROT" %in%  AnnotationDbi::columns(ah) && !is.null(ah)){
        message("Related UniProtKB IDs are retrieved from AnnotationHub...")
        UniProtKB <- AnnotationDbi::select(ah, columns=c("UNIPROT", "ENTREZID"),
            keytype="ENTREZID", keys=targetGeneID)
        UniProtKB <- UniProtKB[, c("ENTREZID", "UNIPROT")]
    }else if("UNIPROTID" %in%  AnnotationDbi::columns(ah) && !is.null(ah)){
        message("Related UniProtKB IDs are retrieved from AnnotationHub...")
        UniProtKB <- AnnotationDbi::select(ah, columns=c("UNIPROTID", "ENTREZID"),
            keytype="ENTREZID", keys=targetGeneID)
        UniProtKB <- UniProtKB[, c("ENTREZID", "UNIPROTID")]
        colnames(UniProtKB) <- c("ENTREZID", "UNIPROT")
    }else{
        UniProtKB <- NULL
    }

    # MeSH
    ah <- AnnotationHub()
    mcah <- mcols(ah)
    msg <- gsub("LRBaseDb", "MeSHDb", ah[metadata(sce)$ahid]$title)
    ahid <- mcah@rownames[which(mcah@listData$title == msg)]
    if(length(ahid) != 0){
        message(paste0("Related MeSH IDs are retrieved from ",
            "AnnotationHub..."))
        MeSHfile <- ah[[ahid]]
        MeSHobj <- MeSHDbi::MeSHDb(MeSHfile)
        MeSH <- MeSHDbi::select(MeSHobj, columns=c("MESHID", "GENEID"),
            keytype="GENEID",
            keys=targetGeneID)
    }else{
        MeSH <- NULL
    }

    # Reactome
    Reactome <- toTable(reactomeEXTID2PATHID)
    targetReactome <- unlist(lapply(targetGeneID,
        function(x){which(Reactome$gene_id == x)}))
    if(length(targetReactome) != 0){
        message(paste0("Related Reactome IDs are retrieved from ",
            "reactome.db package..."))
        Reactome <- Reactome[targetReactome, ]
    }else{
        Reactome <- NULL
    }

    # Output
    list(GeneName=GeneName, Description=Description, GO=GO,
        ENSG=ENSG, ENSP=ENSP, UniProtKB=UniProtKB,
        Reactome=Reactome, MeSH=MeSH)
}

.geneInformation <- function(sce, ah, spc, LR){
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
        target <- as.character(unique(GeneName[, "ENTREZID"]))
        GeneName <- data.frame(
            ENTREZID=target,
            SYMBOL=as.character(g$matching[target]),
            stringsAsFactors = FALSE)

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
    MeSHname <- paste0("MeSH.", gsub(".eg.db.sqlite", "",
        strsplit(metadata(sce)$lrbase, "LRBase.")[[1]][3]), ".eg.db")
    MeSH.load <- eval(parse(text=paste0("try(requireNamespace('", MeSHname, "', quietly=TRUE), silent=TRUE)")))
    if(!MeSH.load){
        eval(parse(text=paste0("try(BiocManager::install('",
            MeSHname, "', update=FALSE, ask=FALSE), silent=TRUE)")))
    }
    MeSH.load2 <- eval(parse(text=paste0("try(require('", MeSHname, "', quietly=TRUE), silent=TRUE)")))
    if(MeSH.load2){
        eval(parse(text=paste0("library(", MeSHname, ")")))
        message(paste0("Related MeSH IDs are retrieved from ",
            "MeSH.XXX.eg.db-type package..."))
        MeSHobj <- eval(parse(text=MeSHname))
        MeSH <- MeSHDbi::select(MeSHobj, columns=c("MESHID", "GENEID"),
            keytype="GENEID",
            keys=targetGeneID)
    }else{
        MeSH <- NULL
    }

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
