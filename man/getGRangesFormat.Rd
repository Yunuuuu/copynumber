% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getGRangesFormat.r
\name{getGRangesFormat}
\alias{getGRangesFormat}
\title{Get segments on the GRanges format}
\usage{
getGRangesFormat(segments)
}
\arguments{
\item{segments}{a data frame containing segmentation results found by e.g.
\code{\link{pcf}}, \code{\link{multipcf}} or \code{\link{aspcf}}.}
}
\value{
The segments converted to the GRanges container class.
}
\description{
The segments data frame obtained e.g. by \code{pcf}, \code{multipcf} or
\code{aspcf} is converted to the GRanges format.
}
\details{
GRanges, in the GenomicRanges package, is the standard BioConductor
containers for range data. For some applications it may therefore be useful
to convert segmentation results to this format.
}
\examples{

# load lymphoma data
data(lymphoma)
# Run pcf
seg <- pcf(data = lymphoma, gamma = 12)
# Obtain the GRanges format
gr <- getGRangesFormat(seg)

}
\author{
Gro Nilsen
}
