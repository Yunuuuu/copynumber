####################################################################
## Author: Gro Nilsen, Knut Liest�l and Ole Christian Lingj�rde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liest�l et al. (2012), BMC Genomics
####################################################################


# Function that plots aberrations given limits - by genome og chromosomes

## Input:
### segments: segmentation results from pcf or multipcf
### upper.lim, lower.lim: a vector with limit(s) to be applied for abberration calling
### pos.unit: unit used for positions (bp,kbp,mbp)
### chrom: a vector with chromosomes to be plotted. If specified the frequencies are plotted with one panel for each chromosome
### layout: number of columns and rows in plot
### ... : other optional plot parameters


## Required by: none

## Requires:
### numericChrom
### checkChrom
### genomeHeat
### chromosomeHeat
### is.multiseg
### checkSegments
### pullOutContent

# Main function for heatmap plotting:


#' Plot areas with copy number aberrations
#'
#' Create plots reflecting the location of aberrated segments. Results may be
#' visualized over the entire genome or by chromosomes.
#'
#' For each sample, the aberrated regions are shown in the color specified in
#' \code{colors[1]} (default dodgerblue) if the segment value is below
#' \code{thres.loss} and the color specified in \code{colors[2]} (default red)
#' if the segment value is above \code{thres.gain}. Non-aberrated regions are
#' shown in white. Each row in the plot represents a sample, while probe
#' positions are reflected along the x-axis.
#'
#' @param segments a data frame containing the segmentation results found by
#' either \code{\link{pcf}} or \code{\link{multipcf}}.
#' @param thres.gain a numeric vector giving the threshold value(s) to be
#' applied for calling gains.
#' @param thres.loss a numeric vector of same length as \code{thres.gain}
#' giving the threshold value(s) to be applied for calling losses. Default is
#' to use the negative value of \code{thres.gain}.
#' @param pos.unit the unit used to represent the probe positions. Allowed
#' options are "mbp" (mega base pairs), "kbp" (kilo base pairs) or "bp" (base
#' pairs). By default assumed to be "bp".
#' @param chrom a numeric or character vector with chromosome number(s) to
#' indicate which chromosome(s) is (are) to be plotted. If unspecified the
#' whole genome is plotted.
#' @param layout the vector of length two giving the number of rows and columns
#' in the plot window. Default is \code{c(1,1)}.
#' @param \dots other optional graphical parameters. These include the plot
#' arguments \code{xlab}, \code{ylab}, \code{main}, \code{cex.main},
#' \code{mgp}, \code{cex.lab}, \code{cex.axis}, \code{mar} and \code{title}
#' (see \code{\link{par}} on these), as well as \code{plot.size},
#' \code{plot.unit}, \code{plot.ideo}, \code{ideo.frac}, \code{cyto.text},
#' \code{assembly} and \code{cex.cytotext} (see \code{\link{plotSample}} on
#' these). In addition, a range of graphical arguments specific for this plot
#' function may be specified: \describe{ \item{list("colors")}{a character
#' vector of length two where the first and second element specifies the color
#' used to represent loss and gain, respectively. Default is
#' c("dodgerblue","red").} \item{list("sample.labels")}{a logical value
#' indicating whether sample labels are to be plotted along the y-axis. Default
#' is TRUE.} \item{list("sep.samples")}{a number in the range 0 to 0.4 used to
#' create some space between samples. 0 implies that there is no space. Default
#' is 2/nsample, where nsample is the number of samples found in
#' \code{segments}.} \item{list("sample.line")}{a numeric scalar giving the
#' margin line where the sample labels should be written, starting at 0
#' counting outwards. Default is 0.2.} \item{list("sample.cex")}{the size of
#' the sample labels.} }
#' @note This function applies \code{graphics::par(fig)}, and is therefore not compatible
#' with other setups for arranging multiple plots in one device such as
#' \code{graphics::par(mfrow,mfcol)}.
#' @author Gro Nilsen
#' @examples
#'
#' # Load lymphoma data
#' data(lymphoma)
#'
#' # Run pcf to obtain estimated copy number values
#' seg <- pcf(data = lymphoma, gamma = 12)
#'
#' # Plot aberrations for the entire genome
#' plotAberration(segments = seg, thres.gain = 0.15)
#'
#' # Plot aberrations for the first 4 chromosomes:
#' plotAberration(segments = seg, thres.gain = 0.1, chrom = c(1:4), layout = c(2, 2))
#'
#' @export
plotAberration <- function(segments, thres.gain, thres.loss = -thres.gain, pos.unit = "bp", chrom = NULL, layout = c(1, 1), ...) {
  # If chrom is unspecified, the whole genome is plotted. Otherwise, selected chromosomes are plotted
  type <- ifelse(is.null(chrom), "genome", "bychrom")

  # Check input in segments:
  segments <- pullOutContent(res = segments, what = "segments")

  # Check segments (convert multi segments to unisegments etc)
  segments <- checkSegments(segments, type)

  # Check and if necessary modify chrom to be plotted
  chrom <- checkChrom(data = NULL, segments = segments, chrom)

  # Convert segments back to data frame (is returned as list in checkSegments)
  segments <- as.data.frame(segments)

  # Get sampleIds
  sampleID <- unique(segments[, 1])

  # Check pos.unit input:
  if (!pos.unit %in% c("bp", "kbp", "mbp")) {
    stop("pos.unit must be one of bp, kbp and mbp", call. = FALSE)
  }

  # Making sure number of upperlimits and lowerlimits are the same:
  upper.lim <- thres.gain
  lower.lim <- thres.loss
  nT <- min(length(upper.lim), length(lower.lim))
  upper.lim <- upper.lim[1:nT]
  lower.lim <- lower.lim[1:nT]

  # plot heatmap, either over genome or by chromosomes:
  switch(type,
    genome = genomeAberration(segments, upper.lim, lower.lim, pos.unit, sampleID, layout, ...),
    bychrom = chromosomeAberration(segments, upper.lim, lower.lim, pos.unit, sampleID, chrom, layout, ...)
  )
}



# Function that plots aberrations by genome

## Input:
### segments: segmentation data frame
### upper.lim, lower.lim: a vector with limit(s) to be applied for abberration calling
### pos.unit: unit used for positions (bp,kbp,mbp)
### sampleID: ids for the samples to be plotted
### layout: number of columns and rows in plot
### ... : other optional plot parameters


## Required by: plotHeatmap

## Requires:
### getHeatParameters
### getGlobal.xlim
### adjustSegPos
### framedim
### colorSetup
### getCol
### addChromlines
### chromPattern
### addArmlines

genomeAberration <- function(segments, upper.lim, lower.lim, pos.unit, sampleID, layout, ...) {
  nT <- length(upper.lim)
  nr <- layout[1]
  nc <- layout[2]
  rc <- nr * nc
  nSample <- length(sampleID)

  op <- getHeatParameters(type = "genome", nc = nc, nr = nr, nSample = nSample, upper.lim = upper.lim, lower.lim = lower.lim, ab = TRUE, ...)

  # Set global xlimits if not specified by user:
  if (is.null(op$xlim)) {
    op$xlim <- getGlobal.xlim(op = op, pos.unit = pos.unit, chrom = unique(segments[, 2]))
  }


  # Check if there should be more than one file/window with plot(s), and get file.name accordingly
  if (grDevices::dev.cur() <= 1) { # to make Sweave work
    grDevices::dev.new(width = op$plot.size[1], height = op$plot.size[2], record = TRUE)
  }

  # Initialize:
  row <- 1
  clm <- 1
  new <- FALSE

  # Division of plotting window:
  frames <- framedim(nr, nc)


  for (t in 1:nT) {
    # Frame dimensions for plot t:
    fig.t <- c(frames$left[clm], frames$right[clm], frames$bot[row], frames$top[row])
    graphics::par(fig = fig.t, new = new, oma = c(0, 0, 0.5, 0), mar = op$mar)

    # Empty plot with correct dimensions:
    plot(1, 1,
      type = "n", ylim = c(0, nSample), ylab = op$ylab, xlab = op$xlab, xlim = op$xlim,
      xaxs = "i", yaxt = "n", xaxt = "n", yaxs = "i", cex.lab = 0.9, mgp = op$mgp, main = ""
    )
    # Let background be white:
    graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = "white")
    # Add shifting white/grey pattern to backgroud to separate chromosomes:
    chromPattern(pos.unit, op)


    # main title for this plot
    graphics::title(main = op$main[t], line = op$main.line, cex.main = op$cex.main)


    # Plot heatmap for each sample:
    for (i in 1:nSample) {
      sample.segments <- segments[segments[, 1] == sampleID[i], ]

      # Adjust positions to be plotted along xaxis; i.e. get global positions, scale according to plot.unit, and get left and right pos for rectangles
      # to be plotted (either continuous or 1 probe long):
      x <- adjustSegPos(chrom = sample.segments[, 2], char.arms = sample.segments[, 3], start = sample.segments[, 4], stop = sample.segments[, 5], type = "genome", unit = pos.unit, op = op)
      xleft <- x$use.start
      xright <- x$use.stop

      # Find appropriate colour for each probe:
      heat.col <- rep("white", nrow(sample.segments))
      heat.col[sample.segments[, 2] %% 2 == 0] <- "grey95"
      heat.col[sample.segments[, 7] > upper.lim[t]] <- op$colors[2]
      heat.col[sample.segments[, 7] < lower.lim[t]] <- op$colors[1]
      ytop <- i - op$sep.samples
      ybottom <- i - (1 - op$sep.samples)

      # Plot rectangles with appropriate color for each probe:
      graphics::rect(xleft, ybottom, xright, ytop, col = heat.col, border = NA)

      # Add sampleid on yaxis
      if (op$sample.labels) {
        graphics::axis(side = 2, at = (ytop - (ytop - ybottom)/2), labels = sampleID[i], line = op$sample.line, tcl = 0, cex.axis = op$sample.cex, las = 1, mgp = op$mgp, tick = FALSE)
      }
    } # endfor


    # Separate chromosomes and arms by vertical lines:
    op$chrom.lty <- 1
    addChromlines(chromosomes = segments[, 2], xaxis = "pos", unit = pos.unit, cex = op$cex.chrom, op = op)
    addArmlines(chromosomes = segments[, 2], xaxis = "pos", unit = pos.unit, cex = op$cex.chrom, op = op)

    # Box:
    graphics::abline(v = op$xlim)
    graphics::abline(h = c(0, nSample))

    # Get new page, or update column/row:
    if (t %% (nr * nc) == 0) {
      # Start new plot page (prompted by user)
      grDevices::devAskNewPage(ask = TRUE)

      # Reset columns and row in layout:
      clm <- 1
      row <- 1
      new <- FALSE
    } else {
      # Update column and row index:
      if (clm < nc) {
        clm <- clm + 1
      } else {
        clm <- 1
        row <- row + 1
      } # endif
      new <- TRUE
    } # endif
  } # endfor
}

# Function that plots aberrations by chromosomes (each chrom in separate panel)

## Input:
### segments: segmentation data frame
### upper.lim, lower.lim: a vector with limits(s) to be applied for abberration calling
### pos.unit: unit used for positions (bp,kbp,mbp)
### sampleID: the IDs for samples to be plotted
### chrom: the chromosomes to be plotted
### layout: number of columns and rows in plot
### ... : other optional plot parameters


## Required by: plotHeatmap

## Requires:
### getHeatParameters
### adjustSegPos
### framedim
### plotIdeogram
### chromMax
### get.xticks


chromosomeAberration <- function(segments, upper.lim, lower.lim, pos.unit, sampleID, chrom, layout, ...) {
  nT <- length(upper.lim)
  nr <- layout[1]
  nc <- layout[2]
  nChrom <- length(chrom)
  nSample <- length(sampleID)


  # Default plot options:
  op <- getHeatParameters(type = "bychrom", nc = nc, nr = nr, nSample = nSample, upper.lim = upper.lim, lower.lim = lower.lim, chrom = chrom, ab = TRUE, ...)

  # Margins for entire plot in window:
  if (all(op$title == "")) {
    oma <- c(0, 0, 0, 0)
  } else {
    oma <- c(0, 0, 1, 0)
  }
  mar <- c(0.2, 0.2, 0.3, 0.2)


  # Divide the plotting window by the function "framedim":
  frames <- framedim(nr, nc)


  # make separate plots for each value of limits
  for (t in 1:nT) {
    # Start new window/file:
    if (grDevices::dev.cur() <= 1) { # to make Sweave work
      grDevices::dev.new(width = op$plot.size[1], height = op$plot.size[2], record = TRUE)
    }

    # Initialize row and column index:
    row <- 1
    clm <- 1
    new <- FALSE

    # Separate plots for each chromosome:
    for (c in 1:nChrom) {
      # Frame dimensions for plot c:
      fig.c <- c(frames$left[clm], frames$right[clm], frames$bot[row], frames$top[row])
      graphics::par(fig = fig.c, new = new, oma = oma, mar = mar)

      # Make list with frame dimensions:
      frame.c <- list(left = frames$left[clm], right = frames$right[clm], bot = frames$bot[row], top = frames$top[row])

      # Select relevant chromosome number
      k <- chrom[c]

      # Pick out indeces for this chromosome
      ind.c <- which(segments[, 2] == k)
      chrom.segments <- segments[ind.c, ]

      # Find maximum position for this chromosome (and scale according to plot unit)
      scale.fac <- convert.unit(unit1 = op$plot.unit, unit2 = pos.unit)
      xlim <- c(0, max(chrom.segments[, 5])) * scale.fac

      # Plot ideogram below heatmaps:
      if (op$plot.ideo) {
        # Ideogram frame:
        ideo.frame <- frame.c
        ideo.frame$top <- frame.c$bot + (frame.c$top - frame.c$bot) * op$ideo.frac

        graphics::par(fig = unlist(ideo.frame), new = new, mar = op$mar.i)
        # Plot ideogram and get maximum probe position in ideogram:
        plotIdeogram(chrom = k, cyto.text = op$cyto.text, cyto.data = op$assembly, cex = op$cex.cytotext, unit = op$plot.unit)

        # Get maximum position for this chromosome:
        xmaxI <- chromMax(chrom = k, cyto.data = op$assembly, pos.unit = op$plot.unit)
        xlim <- c(0, xmaxI)

        new <- TRUE
      }

      # Plot-dimensions:
      frame.c$bot <- frame.c$bot + (frame.c$top - frame.c$bot) * op$ideo.frac
      graphics::par(fig = unlist(frame.c), new = new, mar = op$mar)

      # Limits:
      if (!is.null(op$xlim)) {
        xlim <- op$xlim
      }

      # empty plot set up:
      plot(1, 1,
        type = "n", ylim = c(0, nSample), xlim = xlim, ylab = op$ylab, xlab = op$xlab,
        xaxs = "i", yaxt = "n", xaxt = "n", yaxs = "i", mgp = op$mgp, main = ""
      )

      # Let background be white:
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = "white")

      # main title for this plot
      graphics::title(main = op$main[c], line = op$main.line, cex.main = op$cex.main)

      # Plot heatmap for each sample:
      for (i in 1:nSample) {
        sample.segments <- chrom.segments[chrom.segments[, 1] == sampleID[i], ]
        # Adjust positions to be plotted along xaxis; i.e. scale according to plot.unit, and get left and right pos for rectangles to be plotted (either continuous
        # or 1 probe long):
        x <- adjustSegPos(chrom = sample.segments[, 2], char.arms = sample.segments[, 3], start = sample.segments[, 4], stop = sample.segments[, 5], type = "chromosome", unit = pos.unit, op = op)
        xleft <- x$use.start
        xright <- x$use.stop

        # Find appropriate colour for each probe:
        heat.col <- rep("white", nrow(sample.segments))
        # heat.col[sample.segments[,2]%%2==0] <- "grey95"
        heat.col[sample.segments[, 7] > upper.lim[t]] <- op$colors[2]
        heat.col[sample.segments[, 7] < lower.lim[t]] <- op$colors[1]
        ytop <- i - op$sep.samples
        ybottom <- i - (1 - op$sep.samples)

        graphics::rect(xleft, ybottom, xright, ytop, col = heat.col, border = NA)

        # Add sampleid on yaxis
        if (op$sample.labels) {
          graphics::axis(side = 2, at = (ytop - (ytop - ybottom)/2), labels = sampleID[i], line = op$sample.line, tcl = 0, cex.axis = op$sample.cex, las = 1, mgp = op$mgp, tick = FALSE)
        }
      } # endfor

      if (!op$plot.ideo) {
        # Add xaxis:
        at.x <- get.xticks(xlim[1], xlim[2], unit = op$plot.unit, ideal.n = 6)
        graphics::axis(side = 1, tcl = -0.2, at = at.x, cex.axis = op$cex.axis, mgp = op$mgp)
        graphics::title(xlab = op$xlab, cex.lab = op$cex.lab, line = op$mgp[1])
      }


      # Add box around plot:
      graphics::abline(v = xlim)
      graphics::abline(h = c(0, nSample))



      # If page is full; start plotting on new page
      if (c %% (nr * nc) == 0 && c != nChrom) {
        # Add main title to page:
        graphics::title(op$title[t], outer = TRUE)

        grDevices::devAskNewPage(ask = TRUE)

        # Reset columns and row in layout:
        clm <- 1
        row <- 1
        new <- FALSE
      } else {
        # Update column and row index:
        if (clm < nc) {
          clm <- clm + 1
        } else {
          clm <- 1
          row <- row + 1
        } # endif
        new <- TRUE
      } # endif
    } # endfor

    # Add main title to page:
    graphics::title(op$title[t], outer = TRUE)
    if (t != nT) {
      grDevices::devAskNewPage(ask = TRUE)
    }
  } # endfor
} # endfunction
