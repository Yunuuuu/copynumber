#' 3K aCGH data
#'
#' A subset of the aCGH data set taken from the reference below.
#'
#'
#' @name lymphoma
#' @docType data
#' @format Data frame containing 3091 probes with log2-ratio copy numbers for
#' 21 samples. The first column contains the chromosome numbers, the second
#' gives the local probe positions (in base pairs), while the subsequent
#' columns contain the copy number measurements for the individual samples.
#' @source Eide et al., "Genomic alterations reveal potential for higher grade
#' transformation in follicular lymphoma and confirm parallel evolution of
#' tumor cell clones", Blood 116:1489-1497, 2010
#' @examples
#'
#' # Get data
#' data(lymphoma)
#'
NULL


#' Subset of 244K aCGH data
#'
#' A subset of the 244K MicMa data set containing copy number measurements for
#' six samples on chromosome 17.
#'
#'
#' @name micma
#' @docType data
#' @format Data frame containing 7658 probes with log2-ratio copy numbers for 6
#' samples on chromosome 17. The first column contains the chromosome numbers,
#' the second gives the local probe positions (in base pairs), while the
#' subsequent columns contain the copy number measurements for the individual
#' samples.
#' @source Mathiesen et al., "High resolution analysis of copy number changes
#' in disseminated tumor cells of patients with breast cancer", Int J Cancer
#' 131(4):E405:E415, 2011
#' @examples
#'
#' # Get data
#' data(micma)
#'
NULL


#' Artificial SNP array data
#'
#' Artificial SNP array data containing a logR track and a BAF track
#'
#'
#' @name SNPdata
#' @aliases logR BAF
#' @docType data
#' @format Two corresponding data sets containing 10000 probes with logR and
#' BAF measurements, respectively, for 2 samples. The two first columns in both
#' data sets contain chromosome numbers and local probe positions (in base
#' pairs), while the subsequent columns contain logR-values and BAF-values in
#' the two data sets, respectively.
#' @examples
#'
#' # Get data
#' data(logR)
#' data(BAF)
#'
NULL

#' @import S4Vectors
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges mcols<-
NULL
