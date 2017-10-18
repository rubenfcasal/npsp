#--------------------------------------------------------------------
#   scattersplot.R (npsp package)
#--------------------------------------------------------------------
#   scattersplot  S3 generic
#       scattersplot.default
#       scattersplot.SpatialPointsDataFrame
#
#   (c) R. Fernandez-Casal         Last revision: Jul 2017
#--------------------------------------------------------------------


#' Exploratory scatter plots
#' 
#' Draws (in a 2 by 2 layout) the following plots: 
#' a scatter plot with a color scale, the scatter plots of the response against the (first two) 
#' coordinates and the histogram of the response values.
#' @param  x 	object used to select a method.
#' @param ... additional graphical parameters (to be passed to \code{\link{spoints}}).
#' @seealso \code{\link{splot}}, \code{\link{spoints}}, \code{\link[stats]{lowess}}, 
#' \code{\link[stats]{density}}
#' @keywords hplot
#' @export
#--------------------------------------------------------------------
scattersplot <- function(x, ...) UseMethod("scattersplot")
# S3 generic function splot.scatter
#--------------------------------------------------------------------

#--------------------------------------------------------------------
#' @rdname scattersplot  
#' @method scattersplot default
#' @param z vector of data (response variable).
#' @param main an overall title for the plot.
#' @param xlab a title for the axis corresponding to the first coordinate.
#' @param ylab a title for the axis corresponding to the second coordinate.
#' @param zlab a title for the axis corresponding to the response.
#' @param col color table used to set up the color scale (see \code{\link{spoints}}).
#' @param lowess logical. If \code{TRUE}, a \code{\link[stats]{lowess}} smooth is added to 
#' the plots of the response against the coordinates. 
#' @param density logical. If \code{TRUE}, a kernel \code{\link[stats]{density}} estimate  
#' is added to the histogram. 
#' @param omd a vector of the form \code{c(x1, x2, y1, y2)} giving the region inside outer margins 
#' in normalized device coordinates (i.e. fractions of the device region).
#' @details Standard generic function with a default method, in which argument \code{x} 
#' is a matrix with the spatial coordinates (each row is a point).
#' @export
scattersplot.default <- function(x, z, main, xlab, ylab, zlab,
                                 col = hot.colors(128), lowess = TRUE, density = TRUE,
                                 omd = c(0.05, 0.95, 0.01, 0.95), ...) {
  #--------------------------------------------------------------------
  if (missing(xlab)) xlab <- paste0(deparse(substitute(x)), "[, 1]")
  if (missing(ylab)) ylab <- paste0(deparse(substitute(x)), "[, 2]")
  if (missing(zlab)) zlab <- deparse(substitute(z))
  if (missing(main)) main <- paste("Exploratory plots of", zlab)
  old.par <- par(mfrow = c(2,2), omd = omd)
  # Scatter plot with a color scale
  spoints(x[, 1], x[, 2], z, xlab = xlab, ylab = ylab, 
          col = col, ...)
  # plot(z, x[, 2], z, xlab = main, ylab = xlab))
  plot(x[, 2], z, xlab = xlab, ylab = zlab)
  if (lowess) lines(lowess(x[, 2], z), lty = 2, lwd = 2, col = 'darkgray')
  plot(x[, 1], z, xlab = ylab, ylab = zlab)
  if (lowess) lines(lowess(x[, 1], z), lty = 2, lwd = 2, col = 'darkgray')
  hist(z, xlab = zlab, main = "", freq = !density)
  if (density) lines(density(z), col = 'darkgray')
  par(old.par)
  title(main = main)
}


#--------------------------------------------------------------------
#' @rdname scattersplot  
#' @method scattersplot SpatialPointsDataFrame
#' @param data.ind integer (or character) with the index (or name) of the data component.
#' @details \code{scattersplot.SpatialPointsDataFrame} sets default values for some of the arguments 
# from attributes of the object \code{x} (if present; see e.g. \code{\link{precipitation}}). 
#' from attributes of the object \code{x} (if present; see e.g. \code{precipitation}). 
#' @export
scattersplot.SpatialPointsDataFrame <- function(x, data.ind = 1, 
                                                main, xlab, ylab, zlab, col = hot.colors(128), lowess = TRUE, density = FALSE,  
                                                omd = c(0.05, 0.95, 0.01, 0.95), ...) {
  #--------------------------------------------------------------------
  if(dimensions(x) != 2) 
    stop("function only works for two-dimensional 'SpatialPointsDataFrame'")
  if(!is.null(labels <- attr(x, "labels"))) {
    if(missing(main)) main <- labels$name
    if(missing(xlab)) xlab <- labels$x[1]
    if(missing(ylab)) ylab <- labels$x[2]
    if(missing(zlab)) zlab <- labels$y
  }
  # Coordenadas de los datos
  coord <- coordinates(x)
  scattersplot.default(coord, x[[data.ind]], main = main, xlab = xlab, ylab = ylab, zlab = zlab,
                       col = col, lowess = lowess, density = density, omd = omd, ...)
}
