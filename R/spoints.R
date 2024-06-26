#····································································
#   spoints.R (npsp package)
#····································································
#   spoints  S3 generic
#       spoints.default
#       spoints.data.grid
#       spoints.SpatialPointsDataFrame
#
#   Based on image.plot and drape.plot functions from package fields:
#   fields, Tools for spatial data
#   Copyright 2004-2013, Institute for Mathematics Applied Geosciences
#   University Corporation for Atmospheric Research
#   Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
#   (c) R. Fernandez-Casal
#   Created: Mar 2014, Modified: Apr 2023
#
#   NOTE: Press Ctrl + Shift + O to show document outline in RStudio
#····································································



#····································································
# spoints ----
#····································································
#' Scatter plot with a color scale
#'
#' \code{spoints} (generic function) draws a scatter plot with points filled with different colors
#' and (optionally) adds a legend strip with the color scale
#' (calls \code{\link{splot}} and \code{\link{plot.default}}).
#'
#' @seealso \code{\link{splot}}, \code{\link{simage}}, \code{\link{spersp}}, 
#' \code{\link{image}}, \code{\link[fields]{image.plot}}, \code{\link{data.grid}}, 
#' \code{\link{plot.default}}.
#' @return Invisibly returns a list with the following 3 components:
#' \item{bigplot}{plot coordinates of the main plot. These values may be useful for 
#' drawing a plot without the legend that is the same size as the plots with legends.}
#' \item{smallplot}{plot coordinates of the secondary plot (legend strip).}
#' \item{old.par}{previous graphical parameters (\code{par(old.par)} 
#' will reset plot parameters to the values before entering the function).}
#' @section Side Effects: After exiting, the plotting region may be changed 
#' (\code{\link{par}("plt")}) to make it possible to add more features to the plot
#' (set \code{reset = FALSE} to avoid this).
#' @author
#' Based on \code{\link[fields]{image.plot}} function from package \pkg{fields}:
#' fields, Tools for spatial data.
#' Copyright 2004-2013, Institute for Mathematics Applied Geosciences.
#' University Corporation for Atmospheric Research.
#'
#' Modified by Ruben Fernandez-Casal <rubenfcasal@@gmail.com>.
#' @keywords hplot
#' @export
#····································································
spoints <- function(x, ...) UseMethod("spoints")
# S3 generic function spoints
#····································································


#' @rdname spoints  
#' @method spoints default
#' @param x object used to select a method. In the default method, it provides the \code{x} 
#' coordinates for the plot (and optionally the \code{y} coordinates; 
#' any reasonable way of defining the coordinates is acceptable, 
#' see the function \code{\link{xy.coords}} for details).
#' @param y y coordinates. Alternatively, a single argument \code{x} can be provided.
#' @param s numerical vector containing the values used for coloring the points. 
#' @param legend logical; if \code{TRUE} (default), the plotting region is splitted into two parts,
#' drawing the main plot in one and the legend with the color scale in the other.
#' If \code{FALSE} only the (coloured) main plot is drawn and the arguments related
#' to the legend are ignored (\code{\link{splot}} is not called).
#' @param bigplot plot coordinates for main plot. If not passed, and \code{legend}
#' is TRUE, these will be determined within the function.
#' @param smallplot plot coordinates for legend strip. If not passed, and \code{legend}
#' is TRUE, these will be determined within the function.
#' @param add logical; if \code{TRUE} the scatter plot is just added 
#' to the existing plot.
#' @param reset logical; if \code{FALSE} the plotting region
#' (\code{\link{par}("plt")}) will not be reset to make it possible to add more features
#' to the plot (e.g. using functions such as points or lines). If \code{TRUE} (default) 
#' the plot parameters will be reset to the values before entering the function.
#' @param pch vector of plotting characters or symbols: see \code{\link{points}}.
#' @param cex numerical vector giving the amount by which plotting characters
#' and symbols should be scaled relative to the default. This works as a multiple
#' of \code{\link{par}("cex")}.
#' @param xlab label for the x axis, defaults to a description of \code{x}.
#' @param ylab label for the y axis, defaults to a description of \code{y}.
#' @param asp the y/x aspect ratio, see \code{\link{plot.window}}.
#' @param ... additional graphical parameters (to be passed to the main plot function
#' or \code{sxxxx.default}; e.g. \code{xlim, ylim,} ...). NOTE:
#' graphical arguments passed here will only have impact on the main plot.
#' To change the graphical defaults for the legend use the \code{\link{par}}
#' function beforehand (e.g. \code{par(cex.lab = 2)} to increase colorbar labels).
#' @inheritParams splot
#' @keywords hplot
#' @examples
#' with( aquifer, spoints(lon, lat, head, main = "Wolfcamp aquifer data"))
#' @export
#····································································
spoints.default <- function(x, y = NULL, s, slim = range(s, finite = TRUE), col = jet.colors(128),
    breaks = NULL, legend = TRUE, horizontal = FALSE, legend.shrink = 1.0,
    legend.width = 1.2, legend.mar = ifelse(horizontal, 3.1, 5.1), legend.lab = NULL,
    bigplot = NULL, smallplot = NULL, lab.breaks = NULL, axis.args = NULL,
    legend.args = NULL, add = FALSE, reset = TRUE,
    pch = 16, cex = 1.5, xlab = NULL, ylab = NULL, asp = NA, ...) {
#····································································
    if (legend)
        # image in splot checks breaks and other parameters...
        res <- splot(slim = slim, col = col, breaks = breaks, horizontal = horizontal,
            legend.shrink = legend.shrink, legend.width = legend.width,
            legend.mar = legend.mar, legend.lab = legend.lab,
            bigplot = bigplot, smallplot = smallplot, lab.breaks = lab.breaks,
            axis.args = axis.args, legend.args = legend.args, add = add)
    else {
        if (missing(bigplot)) {
          old.par <- par(no.readonly = TRUE)
          bigplot <- old.par$plt
        } else
          old.par <- par(plt = bigplot)
          # old.par <- par(plt = bigplot, no.readonly = TRUE)
        # par(xpd = FALSE)
        res <- list(bigplot = bigplot, smallplot = NA, old.par = old.par)    
    }
    if (add & !reset) {
        # Creo que realmente no haria falta...
        warning("'reset' argument ignored when 'add = TRUE'")
        reset <- TRUE
    } 
    if (reset) on.exit(par(res$old.par))
    if (is.null(breaks)) {
        # Compute breaks (in 'cut.default' style...)
        ds <- diff(slim)
        if (ds == 0) ds <- abs(slim[1L])
        breaks <- seq.int(slim[1L] - ds/1000, slim[2L] + ds/1000, length.out = length(col) + 1)
        # Only if !missing(slim) else breaks <- length(col) + 1?
    }
    icol <- cut(as.numeric(s), breaks, labels = FALSE, include.lowest = TRUE, right = FALSE) # Use .bincode instead of cut?
    if (!add) {
        # if (is.null(xlab)) xlab <- deparse(substitute(x))
        # if (is.null(ylab)) ylab <- if (!missing(y)) deparse(substitute(y)) else "Y"
        plot(x, y, type = "p", pch = pch, cex = cex, col = col[icol], xlab = xlab, ylab = ylab, asp = asp, ...)
    } else
        points(x, y, pch = pch, cex = cex, col = col[icol], ...)
    # if (reset) par(res$old.par)
    return(invisible(res))
#····································································
}   # spoints


#····································································
# spoints S3 methods ----
#····································································
#' @rdname spoints  
#' @method spoints data.grid
#' @export
spoints.data.grid <- function(x, s = x[[1]], xlab = NULL, ylab = NULL, ...) {
#····································································
    if (!inherits(x, "data.grid") | x$grid$nd != 2L)
        stop("function only works for two-dimensional gridded data ('data.grid'-class objects)")
    ns <- dimnames(x)
    if (is.null(xlab)) xlab <- ns[1]
    if (is.null(ylab)) ylab <- ns[2]
    res <- spoints.default(coords(x), s = s, xlab = xlab, ylab = ylab, ...)  
    return(invisible(res))
#····································································
} # spoints.grid.par


#' @rdname spoints  
#' @method spoints SpatialPointsDataFrame
#' @param data.ind integer (or character) with the index (or name) of the data component.
#' @details \code{spoints.SpatialPointsDataFrame} sets default values for some of the arguments 
# from attributes of the object \code{x} (if present; see e.g. \code{\link{precipitation}}). 
#' from attributes of the object \code{x} (if present; see e.g. \code{precipitation}). 
#' @param main an overall title for the plot.
#' @export
#····································································
spoints.SpatialPointsDataFrame <- function(x, data.ind = 1, main, xlab, ylab, legend.lab, ...) {
  #····································································
  # if (!requireNamespace("sp")) stop("'sp' package required")
  if(dimensions(x) != 2) 
    stop("function only works for two-dimensional 'SpatialPointsDataFrame'")
  if(!is.null(labels <- attr(x, "labels"))) {
    if(missing(main)) main <- labels$y
    if(missing(xlab)) xlab <- labels$x[1]
    if(missing(ylab)) ylab <- labels$x[2]
    if(missing(legend.lab)) legend.lab <- labels$scale
  }
  # Coordenadas de los datos
  coord <- coordinates(x)
  y <- x[[data.ind]]
  spoints(coord, s = y, main = main, xlab = xlab, ylab = ylab, legend.lab = legend.lab, ...)
  if(!is.null(border <- attr(x, "border")) && is(border, "SpatialPolygons"))
    plot(border, border = "darkgray", lwd = 2, add = TRUE)
  if(!is.null(interior <- attr(x, "interior")) && is(interior, "SpatialPolygons"))
    plot(interior, border = "lightgray", lwd = 2, add = TRUE)
  if(!is.null(interior) || !is.null(border))
    spoints(coord[,1], coord[,2], y, add = TRUE)
}




