#' add new assembly to copynumber pakcage
#' @param assembly A single symbol, indicates the genome assembly cytoBand
#' should be added into copynumber package. the object name will be used by
#' `copynumber` package to derive this object. `copynumber` use [get] function
#' to find the object. If `NULL`, the supported_assembly will be returned
#' directly. 
#' @return a character vector of `supported_assembly` invisibly
#' @export 
add_assembly <- function(assembly = NULL) {
    if (is.null(assembly)) return(invisible(supported_assembly))
    envir <- topenv(environment(NULL))
    assembly <- substitute(assembly)
    if (!is.name(assembly)) {
        stop("assembly must be a simple symbol", call. = FALSE)
    }
    name <- deparse(assembly)
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
        stop("assembly should be a data.frame", call. = FALSE)
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
