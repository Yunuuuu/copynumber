####################################################################
## Author: Gro Nilsen, Knut Liest�l and Ole Christian Lingj�rde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liest�l et al. (2012), BMC Genomics
####################################################################

### GAMMAPLOT: do pcf-segmentation for 10 different gamma-values using data for one sample at chromosome 1
### make a 4 by 3 plot with data in first panel, segmentation results in the next 10 and the number of segments
### found for each gamma in the last panel.


## Input:
### data: a numeric matrix with chromosome numbers in the first column, local probe positions in the second, and copy number data for one or more samples in subsequent columns. The header the copy number column(s) should be a sample identifier
### pos.unit: the unit used to represent the probe positions. Allowed options are "mbp" (mega base pairs), "kbp" (kilo base pairs) or "bp" (base pairs). By default assumed to be "bp".
### gammaRange: a vector of length two giving the lowest and highest value of gamma to be applied. 10 values within this range is then used in the pcf-segmentation(rounding to integer values is done if necessary). Default range is c(4,40)
### dowins: logical value indicating whether data should be winsorized before running \code{pcf}. Default is TRUE
### sample: a scalar indicating which sample is to used. The number should correspond to the sample's place (in order of appearance) in the data. Default is to use the first sample present in the data input
### chrom: a scalar indicating which chromosome is to be used. Default is chromosome 1
### cv: logical; should cross-validation be done?
### K: number of folds in K-fold cv
### cex
### col
### seg.col
### ...: other optional parameters to be passed to pcf


## Required by:
###  none

## Requires:
### adjustSeg
### connectSeg
### convert.unit
### get.xticks
### get.yticks
### pcf
### subsetData
### winsorize



#' Plot segmentation results for several values of gamma
#'
#' Data for one sample on one chromosome is segmented by \code{pcf} for 10
#' values of gamma, and results are visualized in a multi-grid plot.
#'
#' Data for one sample and one chromosome is selected, and \code{pcf} is run on
#' this data subset while applying 10 different gamma-values (within the given
#' range). The output is a multi-grid plot with the data shown in the first
#' panel, the segmentation results for the various gammas in the subsequent 10
#' panels, and the number of segments found for each gamma in the last panel.
#'
#' If \code{cv = TRUE} a K-fold cross-validation is also performed. For each
#' fold, a random (100/K) per cent of the data are set to be missing, and
#' \code{pcf} is run using the different values of \code{gamma}. The missing
#' probe values are then predicted by the estimated value of their closest
#' non-missing neighbour (see \code{pcf} on this), and the prediction error for
#' this fold is then calculated as the sum of the squared difference between
#' the predicted and the observed values. The process is repeated over the K
#' folds, and the average prediction errors are finally plotted along with the
#' number of segments in the last panel of the plot. The value of gamma for
#' which the minimum prediction error is found is marked by an asterix. Note
#' that such cross-validation tends to favor small values of gamma, and the
#' suitability of the so-called optimal gamma from this procedure should be
#' critically assessed.
#'
#' @param data either a data frame or the name of a tab-separated file from
#' which copy number data can be read. The rows of the data frame or file
#' should represent the probes. Column 1 must hold numeric or character
#' chromosome numbers, column 2 the numeric local probe positions, and
#' subsequent column(s) the numeric copy number measurements for one or more
#' samples. The header of copy number columns should give sample IDs.
#' @param pos.unit the unit used to represent the probe positions. Allowed
#' options are "mbp" (mega base pairs), "kbp" (kilo base pairs) or "bp" (base
#' pairs). By default assumed to be "bp".
#' @param gammaRange a vector of length two giving the lowest and highest value
#' of gamma to be applied. 10 (approximately) equally spaced values within this
#' range are applied in the pcf-segmentation. Default range is
#' \code{c(10,100)}.
#' @param dowins logical value indicating whether data should be winsorized
#' before running \code{pcf}. Default is TRUE.
#' @param sample an integer indicating which sample is to be segmented. The
#' number should correspond to the sample's place (in order of appearance) in
#' \code{data}. Default is to use the first sample present in the data input.
#' @param chrom a number or character indicating which chromosome is to be
#' segmented. Default is chromosome 1.
#' @param cv logical value indicating whether K-fold cross-validation should be
#' done, see details.
#' @param K the number of folds to use in K-fold cross-validation, default is
#' 5.
#' @param cex size of data points, default is 2.
#' @param col color used to plot data points, default is "grey".
#' @param seg.col color used to plot segments, default is "red".
#' @param \dots other optional parameters to be passed to \code{pcf}.
#' @return If \code{cv = TRUE} a list containing: \item{gamma}{the gamma values
#' applied.} \item{pred.error}{the average prediction error for each value of
#' gamma.} \item{opt.gamma}{the gamma for which the average prediction error is
#' minimized.}
#' @note This function applies \code{graphics::par(fig)}, and is therefore not compatible
#' with other setups for arranging multiple plots in one device such as
#' \code{graphics::par(mfrow,mfcol)}.
#' @author Gro Nilsen, Knut Liestoel, Ole Christian Lingjaerde
#' @seealso \code{\link{pcf}},\code{\link{winsorize}}
#' @examples
#'
#' # Micma data
#' data(micma)
#'
#' plotGamma(micma, chrom = 17)
#' @export
plotGamma <- function(data, pos.unit = "bp", gammaRange = c(10, 100), dowins = TRUE, sample = 1, chrom = 1, cv = FALSE, K = 5, cex = 2, col = "grey", seg.col = "red", ...) {
  nGamma <- 10 # number of gamma-values to test within the range

  # Data is either data frame or a file name
  # select data for chromosome 1 and the first sample represented in data
  data <- subsetData(data, chrom = chrom[1], sample = sample[1])
  pos <- data[, 2]

  # winsorize data before segmentation:
  use.data <- data
  if (dowins) {
    use.data <- winsorize(data = use.data, pos.unit = pos.unit, verbose = FALSE)
  }


  # Set plot-parameters: limits, tickmarks,axis labels:

  # Want to scale the x-axis to fit the desired unit given in plot.unit (default is mega base pairs)
  scale.fac <- convert.unit(unit1 = "mbp", unit2 = pos.unit)
  x <- pos * scale.fac

  # limits:
  xlim <- c(0, max(x))
  # Take out the 5% most extreme observations:
  q <- 0.00
  data.ylim <- stats::quantile(use.data[, 3], probs = c(q/2, (1 - q/2)), names = FALSE, type = 4, na.rm = TRUE)

  # tickmarks:
  at.x <- get.xticks(xlim[1], xlim[2], unit = "mbp", ideal.n = 6)
  at.y <- get.yticks(data.ylim[1], data.ylim[2])

  # f <- 1-0.013*12
  mgp.x <- c(1.5, 0.1, 0)
  mgp.y <- c(2.3, 0.5, 0)
  cex.axis <- 0.9
  cex.lab <- 0.8
  cex.main <- 1.1
  tcl <- -0.25
  # get rgb components:
  q <- grDevices::col2rgb(col)
  col2 <- grDevices::rgb(q[1], q[2], q[3], maxColorValue = 255, alpha = 100)

  # empty plot
  plot.size <- c(11.6, 8.2)
  if (grDevices::dev.cur() <= 1) { # to make Sweave work
    grDevices::dev.new(width = plot.size[1], height = plot.size[2])
  }
  graphics::par(mfrow = c(4, 3), oma = c(3, 3, 1, 0), mar = c(1, 1, 2, 1))

  # Plot data or wins data
  plot(x, use.data[, 3], ylab = "", xlab = "", main = "data", pch = ".", cex = cex, cex.main = cex.main, col = col, xlim = xlim, xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", ylim = data.ylim)
  graphics::axis(side = 2, cex.axis = cex.axis, at = at.y, mgp = mgp.y, las = 1, tcl = tcl)

  # Find vector of gammas to be applied:
  gamma <- seq(from = gammaRange[1], to = gammaRange[2], length.out = nGamma)
  gamma <- round(gamma, digits = 1)

  # Run pcf on the data with the various choice of gammas:
  segments <- vector("list", 0) # empty list
  nSeg <- rep(0, 0)
  for (g in 1:nGamma) {
    seg <- pcf(data = use.data, gamma = gamma[g], pos.unit = pos.unit, verbose = FALSE, ...)
    segments <- c(segments, list(seg)) # add this segmentation to segments-list
    nSeg <- c(nSeg, nrow(seg)) # add number of segments found for this gamma

    # Plot data or wins data
    plot(x, use.data[, 3], ylab = "", xlab = "", main = paste("gamma = ", gamma[g], sep = ""), pch = ".", cex = cex, cex.main = cex.main, col = col2, xlim = xlim, xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", ylim = data.ylim)
    if (g > 8) {
      # add xaxis
      graphics::axis(side = 1, cex.axis = cex.axis, at = at.x, labels = as.integer(at.x), mgp = mgp.x, tcl = tcl)
    }
    if (g %in% c(3, 6, 9)) {
      # add yaxis
      graphics::axis(side = 2, cex.axis = cex.axis, at = at.y, mgp = mgp.y, las = 1, tcl = tcl)
    }

    # Plot segments
    use.seg <- segments[[g]]

    # Adjust start and stop such that segments are connected, and x-axis is represented as mbp
    a <- adjustSeg(chrom = use.seg[, 2], char.arms = use.seg[, 3], start = use.seg[, 4], stop = use.seg[, 5], nPos = use.seg[, 6], type = "sample", xaxis = "pos", unit = pos.unit, connect = TRUE, op = list(plot.unit = "mbp"))

    graphics::segments(x0 = a$use.start, y0 = use.seg[, 7], x1 = a$use.stop, y1 = use.seg[, 7], col = seg.col, lwd = 1, lty = 1)

    # Connect segments vertically:
    connectSeg(a$sep.arm, nSeg = nrow(use.seg), a$use.stop, seg.mean = use.seg[, 7], col = seg.col, lwd = 1, lty = 1)
  }


  # Plot number of segments for each gamma in last panel:
  if (cv) {
    graphics::par(mar = c(1, 2, 2, 3)) # ,oma=c(2,0,0,0))
    plot(gamma, nSeg, xlab = "", ylab = "", axes = FALSE, main = "", type = "b", pch = 19, col = "black", yaxs = "r")

    graphics::axis(side = 1, cex.axis = cex.axis, mgp = mgp.x, tcl = tcl, at = gamma)

    graphics::box()
    graphics::par(xpd = TRUE)
    graphics::legend("topleft", legend = "# segments", pch = 19, lty = 1, col = "black", cex = 1, inset = c(-.1, -.2), bty = "n")
    graphics::par(xpd = FALSE)
  } else {
    graphics::par(mar = c(1, 2, 2, 1))
    graphics::barplot(height = nSeg, names.arg = as.character(gamma), beside = TRUE, space = 0.5, xlab = "", ylab = "", axes = FALSE, mgp = mgp.x, col = "lightgrey", main = "# segments", cex.main = cex.main)
    graphics::abline(h = 0) # xline
  }
  # add yaxis
  at <- pretty(x = c(0, nSeg), n = 3)
  graphics::axis(side = 2, cex.axis = cex.axis, las = 1, mgp = mgp.y, tcl = tcl) # ,at=at)
  graphics::mtext(text = "gamma", side = 1, line = 1.2, cex = cex.lab, outer = TRUE, at = c(0.85, 0.85))

  if (cv) {
    # K-fold cross-validation:
    sse <- matrix(NA, ncol = K, nrow = nGamma)
    nProbe <- nrow(use.data)
    # Find cumulative bin lengths
    cum.bin <- getCumBin(nProbe, K)

    # Get random set of probes
    cv.set <- sample(1:nProbe, nProbe)

    j <- 1
    for (k in 1:K) {
      test.set <- cv.set[j:cum.bin[k]]

      # Replace data values that belong to test set by NA:
      training.data <- use.data
      training.data[test.set, 3] <- NA

      # Run pcf on the training data with the various choice of gammas, and calculate squared residual error in test set
      for (g in 1:nGamma) {
        yhat <- pcf(data = training.data, gamma = gamma[g], pos.unit = pos.unit, return.est = TRUE, verbose = FALSE)$estimates
        # Calculate sum of squared errors for this test set:
        sse[g, k] <- sum((use.data[test.set, 3] - yhat[test.set, 3])^2, na.rm = TRUE)
      }

      j <- cum.bin[k] + 1

      cat(paste("cv progress: ", (100/K) * k, "%"), "\n")
    }

    # Calculate average sse over the K runs:
    pred.error <- apply(sse, 1, mean)
    # Add residual error to last panel plot
    graphics::par(new = TRUE)
    opt <- which.min(pred.error)
    res.pch <- rep(19, length(pred.error))
    res.pch[opt] <- NA
    plot(gamma, pred.error, axes = FALSE, ylab = "", xlab = "", type = "b", pch = res.pch, col = "blue")
    graphics::points(gamma[opt], min(pred.error), pch = 42, col = "blue", cex = 3)
    at <- pretty(x = pred.error, n = 3)
    graphics::axis(4, las = 1, cex.axis = cex.axis, tcl = tcl, mgp = c(2.3, 0.5, 0), col.axis = "blue") # ,at=at)
    # graphics::mtext(side=4,line=2.3,"cv residual error",cex=cex.lab,col="forestgreen")
    graphics::par(xpd = TRUE)
    graphics::legend("topright", legend = "cv pred.error", pch = 19, lty = 1, col = "blue", cex = 1, inset = c(-.07, -.2), bty = "n")
  }

  # x- and y-lab for entire device:
  graphics::mtext(text = "Position (mbp)", side = 1, line = 1.2, cex = cex.lab, outer = TRUE, at = c(0.33, 0.33))
  graphics::mtext(text = "Log R", side = 2, line = 1.5, cex = cex.lab, outer = TRUE)

  if (cv) {
    return(list(gamma = gamma, pred.error = pred.error, optGamma = gamma[opt]))
  }
}



getCumBin <- function(N, K) {
  fl <- floor(N/K)
  ce <- ceiling(N/K)

  if (fl == ce) {
    binsize <- rep(fl, K)
  } else {
    binsize <- c(rep(ce, round((N/K - fl) * K)), rep(fl, round((1 - (N/K - fl)) * K)))
  }
  return(cumsum(binsize))
}
