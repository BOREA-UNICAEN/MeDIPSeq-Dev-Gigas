# This is a modified version of the Plotranges function written by Karline and
# included in the package LIM. This is used in the selectedGeneOntologies.r
# script.
##==============================================================================
## Plots ranges and names                               ##
##==============================================================================

Plotranges <- function (min, max,  value = NULL,
    labels=NULL, log="", pch = 16, pch.col="black",
    line.col = "gray", seg.col="black",
    xlim = NULL, main = NULL, xlab = NULL, ylab = NULL,
    lab.cex = 1.0, mark = NULL, ...)  {


  ##-----------------------------------------------------------------
  ## constructing the data
  ##-----------------------------------------------------------------
  if (! is.vector(value))
    value<-as.vector(unlist(value))
  ranges   <- cbind(min,max,value)
  if (log =="x") {

	  minflow<-min(ranges[ranges!=0])     ## minimum, different from 0
    ranges[ranges==0] <- minflow
    min[min==0]       <- minflow        ## replace 0 with minimum
    max[max==0]       <- minflow
    value[value==0]   <- minflow
  }
  numflows <- length(min)

  ##-----------------------------------------------------------------
  if (is.null(labels)) labels <- names(min)
  if (is.null(labels)) labels <- names(max)
  if (is.null(labels)) labels <- as.character(1:numflows)

  labelwidth   <- max(strwidth(labels, "inch")*lab.cex, na.rm = TRUE)
  labelheight  <- strheight("M", "inch")
  plot.new()

  ##-----------------------------------------------------------------
  ## new margins
  ##-----------------------------------------------------------------

  nmar         <- nm <- par("mar")
  #nmar[2]      <- nmar[4] + (labelwidth + 0.1)/labelheight
  #par(mar = nmar)

  y            <- seq(0, 1, length = numflows)
  if (is.null (xlim))
    xlim <- range(ranges, na.rm=TRUE)
  ylim <- c(0, 1)

  plot.window(xlim = xlim, ylim = ylim, log = log)

  loffset <- (labelwidth + 0.1)/labelheight
  labs    <- labels[1:numflows]
  mtext(labs, side = 2, line = 1.1, at = y, adj = 0.5,
        las = 2, cex = par("cex") * lab.cex, ...)

  abline  (h = y, lty = "dotted", col = line.col)
  if(!is.null(value))
    points  (value, y, pch = pch, col = pch.col)
  segments(min,y,max,y,col=seg.col,lty=1)
  if (! is.null(mark)) {
    text(labels=rep("*",length(mark)),
         rep(xlim[2],length(mark)),mark)
  }
  axis(1)
  box()
  title(main = main, xlab = xlab, ylab = ylab, ...)
  invisible()
  #par("mar"=nm)

}
