supported_assembly <- c(
  "hg16", "hg17", "hg18", "hg19", "hg38",
  "mm7", "mm8", "mm9"
)

#' @keywords internal
get_assembly <- function(x) {
  if (inherits(x, "data.frame")) {
    data <- as.data.frame(x, check.names = FALSE)[1:5]
    names(data) <- c("V1", "V2", "V3", "V4", "V5")
    for (i in c("V1", "V4", "V5")) {
      if (!is.factor(data[[i]])) {
        data[[i]] <- factor(data[[i]])
      }
    }
  } else {
    data <- system.file("extdata", paste0(x, ".rds"),
      package = "copynumber", mustWork = TRUE
    )
    data <- readRDS(data)
  }
  data
}

validate_assembly <- function(x) {
  if (is.character(x) && length(x) == 1L) {
    if (!any(x == supported_assembly)) {
      stop(sprintf(
        "Only assembly %s are be supported, or you can provide a customized data.frame with at least 5 columns",
        paste0(supported_assembly, collapse = ", ")
      ), call. = FALSE)
    }
  } else if (inherits(x, "data.frame")) {
    if (ncol(x) < 5L) {
      stop("A customized assembly must have 5 columns", call. = FALSE)
    }
    message("Cytoband format must in follows:\nV1: chromosome\nV2: start_pos\nV3: end_pos\nV4: band\nV5: gieStain")
  } else {
    stop(sprintf(
      "assembly must be a string (one of %s) or a data.frame",
      paste0(supported_assembly, collapse = ", ")
    ), call. = FALSE)
  }
}
