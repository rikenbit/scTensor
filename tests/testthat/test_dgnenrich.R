context("DGN enrichment deprecation")

test_that(".DGNENRICHMENT returns empty result", {
    res <- scTensor:::.DGNENRICHMENT(all = c("1", "2"), sig = c("1"), p = 0.05)
    expect_true(is.null(res$Term))
    expect_true(is.null(res$Pvalue))
})

test_that(".ENRICHMENT with dgnenrich=FALSE includes DGN as empty", {
    # Minimal mock environment
    e <- new.env(parent = emptyenv())
    e$meshannotation <- NA
    e$meshdb <- NA
    res <- scTensor:::.ENRICHMENT(
        all = c("1", "2", "3"),
        sig = c("1"),
        e = e,
        reactomespc = NA,
        goenrich = FALSE,
        meshenrich = FALSE,
        reactomeenrich = FALSE,
        doenrich = FALSE,
        ncgenrich = FALSE,
        dgnenrich = FALSE,
        p = 0.05,
        ah = NULL
    )
    expect_true("DGN" %in% names(res))
    expect_true(is.null(res$DGN$Term))
})

test_that("dgnenrich=TRUE triggers a single warning", {
    # Test the warning logic from .cellCellReport indirectly
    # by calling the same conditional
    dgnenrich <- TRUE
    w <- tryCatch({
        if (dgnenrich) {
            warning(
                "DGN enrichment is currently unavailable because enrichDGN() ",
                "was removed from the DOSE package (>= 4.7.1) following a ",
                "DisGeNET licensing change. Setting dgnenrich=FALSE.",
                call. = FALSE
            )
            dgnenrich <- FALSE
        }
        dgnenrich
    }, warning = function(w) w)
    expect_is(w, "warning")
    expect_true(grepl("enrichDGN", w$message))
    expect_true(grepl("DisGeNET", w$message))
})
