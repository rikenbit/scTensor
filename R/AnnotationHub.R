.annotationhub_taxid <- function(taxid){
    ###################################
    # OrgDb Search
    ###################################
    ahub <- AnnotationHub()
    target1 <- which(mcols(ahub)[, "taxonomyid"] == taxid)
    target2 <- which(mcols(ahub)[, "rdataclass"] %in% c("OrgDb"))
    target <- intersect(target1, target2)
    ah_ids <- rev(rownames(mcols(ahub))[target])
    hit <- FALSE
    for(i in seq_along(ah_ids)){
        ah_id <- ah_ids[i]
        ah <- ahub[[ah_id]]
        check <- length(which(AnnotationDbi::columns(ah) == "ENTREZID")) == 1
        if(check){
            hit <- TRUE
            break
        }
    }
    if(hit){
        ah
    }else{
        ###################################
        # EnsDb Search
        ###################################
        target3 <- which(mcols(ahub)[, "rdataclass"] %in% c("EnsDb"))
        target <- intersect(target1, target3)
        ah_ids <- rev(rownames(mcols(ahub))[target])
        hit <- FALSE
        for(i in seq_along(ah_ids)){
            ah_id <- ah_ids[i]
            ah <- ahub[[ah_id]]
            check <- length(which(AnnotationDbi::columns(ah) == "ENTREZID")) == 1
            if(check){
                hit <- TRUE
                break
            }
        }
        if(hit){
            ah
        }else{
            NULL
        }
    }
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
