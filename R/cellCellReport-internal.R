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

.NULLPlot <- function(){
    plot(1, 1, col="white", ann=FALSE, xaxt="n", yaxt="n", axes=FALSE)
    par(ps=100)
    text(1, 1, "None", col=rgb(0,0,0,0.5))
}

.SinglePlot <- function(x){
    plot(1, 1, col="white", ann=FALSE, xaxt="n", yaxt="n", axes=FALSE)
    par(ps=100)
    text(1, 1, x, col="red")
}

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

.palf <- colorRampPalette(c("#4b61ba", "gray", "#a87963", "red"))

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

.ggdefault_cols <- function(n){
    hcl(h=seq(15, 375-360/n, length=n)%%360, c=100, l=65)
}

.extractLR <- function(sce, lr.evidence, cols){
    # SQLite connection
    con = dbConnect(SQLite(), metadata(sce)$lrbase)
    LR <- dbGetQuery(con, "SELECT * FROM DATA")
    dbDisconnect(con)
    # Extract corresponding rows
    if(lr.evidence[1] == "all"){
        target <- seq(nrow(LR))
    }else{
        targetKnown <- unique(unlist(lapply(.knowndbs, function(x){grep(x, LR$SOURCEDB)})))
        if(lr.evidence[1] == "known"){
            target <- targetKnown
        }else if(lr.evidence[1] == "putative"){
            target <- setdiff(seq(nrow(LR)), targetKnown)
        }else{
            target <- unique(unlist(lapply(lr.evidence, function(x){
                which(x == LR$SOURCEDB)
            })))
            if(length(target) == 0){
                cat("Please specify the valid lr.evidence!\n")
                cat("In this LRBase package, following databases are avilable as lr.evidence.\n")
                print(as.matrix(table(LR$SOURCEDB)))
                stop("\n")
            }
        }
    }
    .uniqueLR(LR[target, cols])
}

.knowndbs <- c("DLRP", "IUPHAR", "HPMR", "CELLPHONEDB", "SINGLECELLSIGNALR")

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


# For previous LRBase (v-1.0.0 - v-1.2.0)
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
    lr, taxid, value, percentage, GeneInfo, pvalue, qvalue){
    ## helper for vector dividing
    div <- function(x, d=1) {""
        delta <- ceiling(length(x) / d)
        y <- lapply(seq_len(d), function(i){
            as.vector(na.omit(x[((i-1)*delta+1):(i*delta)]))
        })
        return(y)
    }

    embedLink <- function(genename1, geneid1, description1,
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

    convertGeneName <- function(geneid, GeneInfo){
        if(!is.null(GeneInfo$GeneName)){
            genename <- GeneInfo$GeneName[
                which(GeneInfo$GeneName$ENTREZID == geneid),
                "SYMBOL"][1]
            if(is.null(genename)){
                genename = geneid
            }else{
                if(is.na(genename)){
                    genename = geneid
                }
                if(length(genename) == 0 || genename %in% c("", NA)){
                    genename = geneid
                }
            }            
        }else{
            genename = ""
        }
        genename
    }

    convertGeneDescription <- function(geneid, GeneInfo){
        if(!is.null(GeneInfo$Description)){
            description <- GeneInfo$Description[
                which(GeneInfo$Description$ENTREZID ==
                    geneid), "GENENAME"][1]
            if(length(description) == 0 || description %in% c("", NA)){
                description = ""
            }
            description
        }else{
            description = ""
        }
    }

    convertGeneOntology <- function(geneid, GeneInfo){
        if(!is.null(GeneInfo$GO)){
            GO <- unique(unlist(GeneInfo$GO[
                which(GeneInfo$GO$ENTREZID == geneid),
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
        }else{
            GO = ""
        }
        GO
    }

    convertReactome <- function(geneid, GeneInfo){
        if(!is.null(GeneInfo$Reactome)){
            Reactome <- unique(unlist(GeneInfo$Reactome[
                which(GeneInfo$Reactome$gene_id ==
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
        }else{
            Reactome = ""
        }
        Reactome
    }

    convertMeSH <- function(geneid, GeneInfo){
        if(!is.null(GeneInfo$MeSH)){
            MeSH <- GeneInfo$MeSH[which(GeneInfo$MeSH$GENEID == geneid),
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
        }else{
            MeSH = ""
        }
        MeSH
    }

    convertUniProtKB <- function(geneid, GeneInfo){
        if(!is.null(GeneInfo$UniProtKB)){
            UniProtKB <- unique(unlist(GeneInfo$UniProtKB[
                which(GeneInfo$UniProtKB$ENTREZID == geneid),
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
        }else{
            UniProtKB = ""
        }
        UniProtKB
    }

    # taxidあるかチェック
    convertSTRING <- function(geneid, GeneInfo, taxid){
        if(!is.null(GeneInfo$ENSP)){
            ENSP <- unique(unlist(GeneInfo$ENSP[
                which(GeneInfo$ENSP$ENTREZID == geneid),
                "ENSEMBLPROT"]))
            ENSP <- ENSP[which(ENSP != "")]
            if(length(ENSP) != 0 && !is.na(taxid)){
                STRING <- paste0("https://string-db.org/network/",
                    taxid, ".", ENSP)
                STRING <- paste(
                    paste0("[", seq_along(STRING), "](", STRING, ")"),
                    collapse=" ")
            }else{
                STRING <- ""
            }
        }else{
            STRING = ""
        }
        STRING
    }

    convertRefEx <- function(geneid, GeneInfo, taxid){
        if(!is.null(GeneInfo$GeneName)){
            genename <- GeneInfo$GeneName[
                which(GeneInfo$GeneName$ENTREZID == geneid),
                "SYMBOL"][1]
            if(length(genename) != 0){
                ref <- "http://refex.dbcls.jp/genelist.php?gene_name%5B%5D="
                if(taxid == "9606"){
                    paste0(ref, genename, "&lang=en&db=human")
                }else if(taxid == "10090"){
                    paste0(ref, genename, "&lang=en&db=mouse")
                }else if(taxid == "10116"){
                    paste0(ref, genename, "&lang=en&db=rat")
                }else{
                    ""
                }
            }else{
                ""
            }
        }else{
            ""
        }
    }

    convertEA <- function(geneid, GeneInfo){
        if(!is.null(GeneInfo$ENSG)){
            ENSG <- GeneInfo$ENSG[
                which(GeneInfo$GeneName$ENTREZID == geneid),
                "ENSEMBL"][1]
            if(length(ENSG) != 0){
                paste0("https://www.ebi.ac.uk/gxa/genes/", tolower(ENSG))
            }else{
                ""
            }
        }else{
            ""
        }
    }

    convertSEA <- function(geneid, GeneInfo){
        if(!is.null(GeneInfo$GeneName)){
            genename <- GeneInfo$GeneName[
                which(GeneInfo$GeneName$ENTREZID == geneid),
                "SYMBOL"][1]
            if(length(genename) != 0){
                paste0("https://www.ebi.ac.uk/gxa/sc/search?species=&q=",
                    genename)
            }else{
                ""
            }
        }else{
            ""
        }
    }

    convertSCDB <- function(geneid, GeneInfo, taxid){
        if(!is.null(GeneInfo$GeneName)){
            genename <- GeneInfo$GeneName[
                which(GeneInfo$GeneName$ENTREZID == geneid),
                "SYMBOL"][1]
            if(length(genename) != 0 && taxid == "9606"){
                paste0("https://bioinfo.uth.edu/scrnaseqdb/",
                    "index.php?r=site/rankGene&gene=",
                    genename,
                    "&check=0")
            }else{
                ""
            }
        }else{
            ""
        }
    }

    convertPANGLAO <- function(geneid, GeneInfo, taxid){
        if(!is.null(GeneInfo$GeneName)){
            genename <- GeneInfo$GeneName[
                which(GeneInfo$GeneName$ENTREZID == geneid),
                "SYMBOL"][1]
            if(length(genename) != 0){
                ref <- "https://panglaodb.se/search.html?query="
                if(taxid == "9606"){
                    paste0(ref, genename, "&species=3")
                }else if(taxid == "10090"){
                    paste0(ref, genename, "&species=2")
                }else{
                    ""
                }
            }else{
                ""
            }
        }else{
            ""
        }
    }

    convertCMAP <- function(geneid, GeneInfo, taxid){
        if(!is.null(GeneInfo$GeneName)){
            genename <- GeneInfo$GeneName[
                which(GeneInfo$GeneName$ENTREZID == geneid),
                "SYMBOL"][1]
            if(length(genename) != 0){
                ref <- "https://clue.io/command?q="
                if(taxid == "9606"){
                    paste0(ref, genename)
                }else{
                    ""
                }
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
    GeneName_L <- convertGeneName(ligandGeneID, GeneInfo)
    # Gene Name (Receptor)
    GeneName_R <- convertGeneName(receptorGeneID, GeneInfo)
    # Description (Ligand)
    Description_L <- convertGeneDescription(ligandGeneID, GeneInfo)
    # Description (Receptor)
    Description_R <- convertGeneDescription(receptorGeneID, GeneInfo)
    # Gene Ontology (Ligand)
    GO_L <- convertGeneOntology(ligandGeneID, GeneInfo)
    # Gene Ontology (Receptor)
    GO_R <- convertGeneOntology(receptorGeneID, GeneInfo)
    # Reactome (Ligand)
    Reactome_L <- convertReactome(ligandGeneID, GeneInfo)
    # Reactome (Receptor)
    Reactome_R <- convertReactome(receptorGeneID, GeneInfo)
    # MeSH (Ligand)
    MeSH_L <- convertMeSH(ligandGeneID, GeneInfo)
    # MeSH (Receptor)
    MeSH_R <- convertMeSH(receptorGeneID, GeneInfo)
    # UniProtKB（Ligand）
    UniProtKB_L <- convertUniProtKB(ligandGeneID, GeneInfo)
    # UniProtKB（Receptor）
    UniProtKB_R <- convertUniProtKB(receptorGeneID, GeneInfo)
    # STRING（Ligand）
    STRING_L <- convertSTRING(ligandGeneID, GeneInfo, taxid)
    # STRING（Receptor）
    STRING_R <- convertSTRING(receptorGeneID, GeneInfo, taxid)
    # RefEx（Ligand）
    RefEx_L <- convertRefEx(ligandGeneID, GeneInfo, taxid)
    # RefEx（Receptor）
    RefEx_R <- convertRefEx(receptorGeneID, GeneInfo, taxid)
    # EA（Ligand）
    EA_L <- convertEA(ligandGeneID, GeneInfo)
    # EA（Receptor）
    EA_R <- convertEA(receptorGeneID, GeneInfo)
    # SEA（Ligand）
    SEA_L <- convertSEA(ligandGeneID, GeneInfo)
    # SEA（Receptor）
    SEA_R <- convertSEA(receptorGeneID, GeneInfo)
    # SCDB（Ligand）
    SCDB_L <- convertSCDB(ligandGeneID, GeneInfo, taxid)
    # SCDB（Receptor）
    SCDB_R <- convertSCDB(receptorGeneID, GeneInfo, taxid)
    # PANGLAO（Ligand）
    PANGLAO_L <- convertPANGLAO(ligandGeneID, GeneInfo, taxid)
    # PANGLAO（Receptor）
    PANGLAO_R <- convertPANGLAO(receptorGeneID, GeneInfo, taxid)
    # CMAP（Ligand）
    CMAP_L <- convertCMAP(ligandGeneID, GeneInfo, taxid)
    # CMAP（Receptor）
    CMAP_R <- convertCMAP(receptorGeneID, GeneInfo, taxid)
    # PubMed (L and R)
    PubMed <- convertPubMed(ligandGeneID, receptorGeneID, lr)

    # Embedding
    paste0(XYZ,
        embedLink(GeneName_L, ligandGeneID, Description_L,
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

.eachVecLR <- function(x, e){
    p <- e$p
    index <- e$index
    sce <- e$sce
    ah <- e$ah
    .HCLUST <- e$.HCLUST
    .OUTLIERS <- e$.OUTLIERS
    top <- e$top
    GeneInfo <- e$GeneInfo
    out.dir <- e$out.dir
    .smallTwoDplot <- e$.smallTwoDplot
    input <- e$input
    twoD <- e$twoD
    .hyperLinks <- e$.hyperLinks
    LR <- e$LR
    taxid <- e$taxid
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
    spc <- gsub(".eg.db.sqlite", "",
        strsplit(metadata(sce)$lrbase, "LRBase.")[[1]][3])
    # GO-Check
    if("GO" %ni% AnnotationDbi::columns(ah)){
        goenrich <- FALSE
    }
    # MeSH-Check
    MeSHname <- paste0("MeSH.", spc, ".eg.db")
    MeSH.load <- eval(parse(text=paste0("try(requireNamespace('", MeSHname, "', quietly=TRUE), silent=TRUE)")))
    if(!MeSH.load){
        eval(parse(text=paste0("try(BiocManager::install('",
            MeSHname, "', update=FALSE, ask=FALSE), silent=TRUE)")))
    }
    MeSH.load2 <- eval(parse(text=paste0("try(require('", MeSHname, "', quietly=TRUE), silent=TRUE)")))
    if(MeSH.load2){
        eval(parse(text=paste0("library(", MeSHname, ")")))
        meshannotation <- MeSHname
    }else{
        meshenrich <- FALSE
        meshannotation <- NA
    }
    # Reactome-Check
    reactomespc <- .REACTOMESPC[[spc]]
    if(is.null(reactomespc)){
        reactomeenrich <- FALSE
    }
    # DO/NCG/DGN-Check
    if(taxid != "9606"){
        doenrich <- FALSE
        ncgenrich <- FALSE
        dgnenrich <- FALSE
    }
    Enrich <- suppressWarnings(.ENRICHMENT(all, sig,
        meshannotation, reactomespc,
        goenrich, meshenrich, reactomeenrich, doenrich, ncgenrich, dgnenrich, p, ah))

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
    if(!is.null(GeneInfo$GeneName)){
        GeneName <- GeneInfo$GeneName[, c("SYMBOL", "ENTREZID")]
    }else{
        GeneName <- rep(FALSE, length=length(TARGET))
    }
    LINKS <- vapply(seq_along(TARGET), function(xx){
        # IDs
        Ranking <- xx
        L_R <- strsplit(names(TARGET[xx]), "_")
        LigandGeneID <- L_R[[1]][1]
        ReceptorGeneID <- L_R[[1]][2]
        if(!is.null(GeneInfo$GeneName)){
            LigandGeneName <- GeneName[which(GeneName[,2] == LigandGeneID), 1]
            ReceptorGeneName <- GeneName[which(GeneName[,2] == ReceptorGeneID), 1]

            if(length(LigandGeneName) == 0 || LigandGeneName %in% c("", NA)){
                LigandGeneName <- LigandGeneID
            }
            if(length(ReceptorGeneName) == 0 || ReceptorGeneName %in% c("", NA)){
                ReceptorGeneName <- ReceptorGeneID
            }
        }else{
            LigandGeneName <- LigandGeneID
            ReceptorGeneName <- ReceptorGeneID
        }
        # Return Hyper links
        .hyperLinks(Ranking, LigandGeneID,
        ReceptorGeneID, LR, taxid, Value[xx],
        Percentage[xx],
        GeneInfo, PvalueLR[xx],
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

.genePlot <- function(sce, assayNames, input, out.dir, GeneInfo, LR){
    GeneName <- GeneInfo$GeneName
    LigandGeneID <- unique(LR$GENEID_L)
    ReceptorGeneID <- unique(LR$GENEID_R)
    if(!is.null(GeneName)){
        LigandGeneName <- vapply(LigandGeneID, function(x){
            GeneName[which(GeneName$ENTREZID == x)[1], "SYMBOL"]}, "")
        ReceptorGeneName <- vapply(ReceptorGeneID, function(x){
            GeneName[which(GeneName$ENTREZID == x)[1], "SYMBOL"]}, "")        
    }else{
        LigandGeneName <- LigandGeneID
        ReceptorGeneName <- ReceptorGeneID
    }

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
        par(ps=20)
        plot(twoD, col=label.receptor, pch=16, cex=2, bty="n",
            xaxt="n", yaxt="n", xlab="", ylab="",
            main="")
        dev.off()
    }, 0L)
}
