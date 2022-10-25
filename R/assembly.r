#' add new assembly to copynumber pakcage
#' @param name a scalar string, indicates the genome assembly cytoBand should be
#' added into copynumber package. Notes: should be the same with the object name
#' as `copynumber` package use `get` to derive this object.
#' @return a character vector of `supported_assembly` invisibly
#' @export 
add_assembly <- function(name) {
    envir <- topenv(environment(NULL))
    supported_assembly <- c(supported_assembly, name)
    if (!exists(name, envir = envir, inherits = FALSE)) {
        unlockBinding("supported_assembly", envir)
        assignInMyNamespace("supported_assembly", supported_assembly)
        lockBinding("supported_assembly", envir)
    } else {
        stop(name, " already exits in copynumber package namespace")
    }
    invisible(supported_assembly)
}

supported_assembly <- c(
    "hg16", "hg17", "hg18", "hg19", "hg38",
    "mm7", "mm8", "mm9"
)

get_assembly <- function(x) {
    data <- get(x, inherits = TRUE)
    if (!inherits(data, "data.frame")) {
        stop("assembly should be a data.frame")
    }
    data <- as.data.frame(data)
    names(data) <- c("V1", "V2", "V3", "V4", "V5")
    if (!is.factor(data$V1)) {
        data$V1 <- factor(data$V1)
    }
    if (!is.factor(data$V4)) {
        data$V4 <- factor(data$V4)
    }
    if (!is.factor(data$V5)) {
        data$V5 <- factor(data$V5)
    }
    data
}
