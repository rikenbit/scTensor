.GOENRICHMENT <- function(all, sig, ah, category, p){
    goParams <- try(new("GOHyperGParams",
        geneIds=sig,
        universeGeneIds=all,
        annotation=ah,
        ontology=category,
        pvalueCutoff=p,
        conditional=FALSE,
        testDirection="over"), silent=TRUE)
    if(class(goParams) != "try-error"){
        # Hyper geometric p-value
        out <- try(summary(hyperGTest(goParams)), silent=TRUE)
        if(is(out)[1] == "try-error"){
            list(Term=NULL, Pvalue=NULL)
        }else{
            list(Term=out$Term, Pvalue=out$Pvalue)
        }        
    }else{
        list(Term=NULL, Pvalue=NULL)        
    }
}

.MeSHENRICHMENT <- function(all, sig, e, category, p){
    if(class(e$meshannotation) != "MeSHDb"){
        list(Term=NULL, Pvalue=NULL)
    }else{
        meshParams <- new("MeSHHyperGParams",
            geneIds = sig,
            universeGeneIds = all,
            annotation = "e$meshannotation",
            meshdb = "e$meshdb",
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

.DOENRICHMENT <- function(all, sig, p){
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

.NCGENRICHMENT <- function(all, sig, p){
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

.DGNENRICHMENT <- function(all, sig, p){
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

.NOSIG <- list(Term=NULL, PValue=NULL)

.ENRICHMENT <- function(all, sig, e,
    reactomespc, goenrich, meshenrich, reactomeenrich,
    doenrich, ncgenrich, dgnenrich, p, ah){
    # GO
    if(goenrich){
        cat("GO-Enrichment Analysis is running...(1/3)\n")
        BP <- .GOENRICHMENT(all, sig, ah, "BP", p)
        cat("GO-Enrichment Analysis is running...(2/3)\n")
        MF <- .GOENRICHMENT(all, sig, ah, "MF", p)
        cat("GO-Enrichment Analysis is running...(3/3)\n")
        CC <- .GOENRICHMENT(all, sig, ah, "CC", p)
    }else{
        BP <- .NOSIG
        MF <- .NOSIG
        CC <- .NOSIG
    }
    # MeSH
    if(meshenrich){
        cat("MeSH-Enrichment Analysis is running...(1/16)\n")
        A <- .MeSHENRICHMENT(all, sig, e, "A", p)
        cat("MeSH-Enrichment Analysis is running...(2/16)\n")
        B <- .MeSHENRICHMENT(all, sig, e, "B", p)
        cat("MeSH-Enrichment Analysis is running...(3/16)\n")
        C <- .MeSHENRICHMENT(all, sig, e, "C", p)
        cat("MeSH-Enrichment Analysis is running...(4/16)\n")
        D <- .MeSHENRICHMENT(all, sig, e, "D", p)
        cat("MeSH-Enrichment Analysis is running...(5/16)\n")
        E <- .MeSHENRICHMENT(all, sig, e, "E", p)
        cat("MeSH-Enrichment Analysis is running...(6/16)\n")
        F <- .MeSHENRICHMENT(all, sig, e, "F", p)
        cat("MeSH-Enrichment Analysis is running...(7/16)\n")
        G <- .MeSHENRICHMENT(all, sig, e, "G", p)
        cat("MeSH-Enrichment Analysis is running...(8/16)\n")
        H <- .MeSHENRICHMENT(all, sig, e, "H", p)
        cat("MeSH-Enrichment Analysis is running...(9/16)\n")
        I <- .MeSHENRICHMENT(all, sig, e, "I", p)
        cat("MeSH-Enrichment Analysis is running...(10/16)\n")
        J <- .MeSHENRICHMENT(all, sig, e, "J", p)
        cat("MeSH-Enrichment Analysis is running...(11/16)\n")
        K <- .MeSHENRICHMENT(all, sig, e, "K", p)
        cat("MeSH-Enrichment Analysis is running...(12/16)\n")
        L <- .MeSHENRICHMENT(all, sig, e, "L", p)
        cat("MeSH-Enrichment Analysis is running...(13/16)\n")
        M <- .MeSHENRICHMENT(all, sig, e, "M", p)
        cat("MeSH-Enrichment Analysis is running...(14/16)\n")
        N <- .MeSHENRICHMENT(all, sig, e, "N", p)
        cat("MeSH-Enrichment Analysis is running...(15/16)\n")
        V <- .MeSHENRICHMENT(all, sig, e, "V", p)
        cat("MeSH-Enrichment Analysis is running...(16/16)\n")
        Z <- .MeSHENRICHMENT(all, sig, e, "Z", p)
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
        DO <- .DOENRICHMENT(all, sig, p)
    }else{
        DO <- .NOSIG
    }
    # NCG
    if(ncgenrich){
        cat("NCG-Enrichment Analysis is running...(1/1)\n")
        NCG <- .NCGENRICHMENT(all, sig, p)
    }else{
        NCG <- .NOSIG
    }
    # DGN
    if(dgnenrich){
        cat("DGN-Enrichment Analysis is running...(1/1)\n")
        DGN <- .DGNENRICHMENT(all, sig, p)
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
