## code to prepare `assembly` dataset goes here

internal_data <- load("R/sysdata.rda")
internal_data
levels(hg19[["V1"]])
if (!dir.exists("inst/extdata")) {
    dir.create("inst/extdata")
}
for (i in internal_data) {
    obj <- get(i, inherits = FALSE)
    saveRDS(obj, file = file.path("inst/extdata", paste0(i, ".rds")))
}

fct_cols <- names(hg19)[vapply(hg19, is.factor, logical(1L))]
# add hg38 assembly -------------------------------------------------
if (!"hg38" %in% internal_data) {
    hg38 <- data.table::fread(
        "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz",
        header = FALSE
    )
    hg38 <- hg38[V1 %in% levels(hg19[["V1"]])]
    # check if they have the same value with hg19
    is_same <- hg38[, purrr::imap_lgl(.SD, function(x, i) {
        setequal(unique(x), levels(hg19[[i]]))
    }), .SDcols = fct_cols]
    cat(is_same)
    if (all(is_same)) {
        hg38[, (fct_cols) := purrr::imap(.SD, function(x, i) {
            factor(x, levels(hg19[[i]]))
        }), .SDcols = fct_cols]
    }
    data.table::setDF(hg38)
}

