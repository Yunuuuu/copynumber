## code to prepare `assembly` dataset goes here

internal_data <- load("R/sysdata.rda")
internal_data
levels(hg19[["V1"]])
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
usethis::use_data(
    hg16, hg17, hg18, hg19, hg38,
    mm7, mm8, mm9,
    internal = TRUE,
    overwrite = TRUE
)

