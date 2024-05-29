#····································································
#   np.den.R (npsp package)
#····································································
#   bin.den  S3 class and methods
#       bin.den(x, nbin)
#   as.bin.den()
#       as.bin.den.bin.data(object, ...)
#   np.den  S3 class and methods
#   np.den()  S3 generic
#       np.den.default(x, nbin, h, degree, drv, ncv, ...)
#       np.den.bin.den(x, h, degree, drv, ncv, ...)
#       np.den.bin.data(x, h, degree, drv, ncv, ...)   
#       np.den.svar.bin(x, h, degree, drv, ncv, ...) 
#
#   (c) Ruben Fernandez-Casal
#   Created: Oct 2013
#
#   NOTE: Press Ctrl + Shift + O to show document outline in RStudio
#····································································


#····································································
# bin.den(x, nbin = NULL) {
#····································································
#' Linear binning for density estimation
#' 
#' Creates a \code{bin.den}-\code{\link{class}} (gridded binned density) object 
#' with linear binning counts.
#' @aliases bin.den-class
#' @param  x vector or matrix of covariates (e.g. spatial coordinates). 
#'    Columns correspond with dimensions and rows with observations.
#' @param  nbin vector with the number of bins on each dimension.
#' @param type character, binning method: \code{"linear"} (default) or \code{"simple"}.
#' @details If parameter \code{nbin} is not specified is set to \code{pmax(25, rule.binning(x))}.
#' @return Returns an S3 object of \code{\link{class}} \code{bin.den} (extends \code{\link{data.grid}}). 
#'    A list with the following 3 components:
#' \item{binw}{vector or array (dimension \code{nbin}) with the bin counts (weights).}
#' \item{grid}{a \code{\link{grid.par}}-\code{\link{class}} object with the grid parameters.}
#' \item{data}{a list with a component \code{$x} with argument \code{x}.}
#' @seealso \code{\link{np.den}}, \code{\link{h.cv}}, \code{\link{bin.data}}, 
#' \code{\link{locpol}}, \code{\link{rule.binning}}.
#' @examples 
#' binden <- bin.den(earthquakes[, c("lon", "lat")], nbin = c(30,30))
#' bindat <- binning(earthquakes[, c("lon", "lat")], earthquakes$mag, nbin = c(30,30))
#' all.equal(binden, as.bin.den(bindat))
#' @export
# Interface to the Fortran routine "bin_den"
bin.den <- function(x, nbin = NULL, type = c("linear", "simple")) {
#····································································
    x <- as.matrix(x)
    ny <- nrow(x)               # number of data
    # Remove missing values  
    ok <- complete.cases(x)     # observations having no missing values across x
    if (any(!ok)) {
        warning("missing values removed")
        x <- x[ok,]
        ny <- nrow(x)
    }    
    nd <- ncol(x)                                 # number of dimensions
    if(is.null(nbin)) nbin <- pmax(25L, rule.binning(x)) else   # dimensions of the binning grid
      if (nd != length(nbin)) stop("arguments 'x' and 'nbin' have incompatible dimensions.")
    nt <- prod(nbin)
    type <- match.arg(type)
    itype <- ifelse(type == "simple", 1, 0)
    # Let's go FORTRAN!
    # subroutine bin_den(nd, nbin, x, ny, bin_min, bin_max, bin_w)
    ret <-.Fortran( "bin_den", nd = as.integer(nd), nbin = as.integer(nbin),
                  xt = as.double(t(x)), ny = as.integer(ny),
                  min = double(nd), max = double(nd), 
                  binw = double(nt), itype = as.integer(itype), PACKAGE = "npsp")
    # Construir o resultado
    result <- with( ret, data.grid( binw = binw,
              grid = grid.par(n = nbin, min = min, max = max, 
              dimnames = dimnames(x)[[2]])) )
    result$data <- list(x = x)
    oldClass(result) <- c("bin.den", "data.grid")
    return(result)
#····································································
} # bin.den


#····································································
#' @rdname bin.den  
#' @param object (gridded data) used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @export
#····································································
as.bin.den <- function(object, ...) {
  UseMethod("as.bin.den")
} # S3 generic function as.bin.den


#····································································
#' @rdname bin.den 
#' @method as.bin.den data.grid
#' @param weights.ind integer or character with the index or name of the component 
#'  containing the bin counts/weights.
#' @export
as.bin.den.data.grid <- function(object, weights.ind = 1, ...) {
#····································································
  if (!inherits(object, "data.grid"))
    stop("function only works for objects of class (or extending) 'data.grid'")
  x <- coords(object$grid)   
  binw <- object[[weights.ind]] 
  index <- is.finite(binw)
  binw[!index] <- 0
  if (!all(index <- binw > 0)) x <- x[index, ]
  result <- list(binw = binw, grid = object$grid, data = list(x = x))
  oldClass(result) <- c("bin.den", "data.grid")
  return(result)
}



#····································································
#' @rdname bin.den  
#' @method as.bin.den bin.den
#' @export
as.bin.den.bin.den <- function(object, ...) {
#····································································
  if (inherits(object, "svar.bin")) {
    warning("Conversion not yet implemented; using 'as.bin.den.data.grid()'...")
    return(as.bin.den.data.grid(object, weights.ind = 'binw'))
  }  
  if (!is.null(object$mask)) object$binw[!object$mask] <- 0
  result <- with( object, list(binw = binw, grid = grid, data = list(x = data$x)))
  oldClass(result) <- c("bin.den", "data.grid")
  return(result)
}



#····································································
# np.den(x, ...)
#····································································
#' Local polynomial density estimation 
#' 
#' Estimates a multidimensional probability density function (and its first derivatives) 
#' using local polynomial kernel smoothing of linearly binned data.
#' @aliases np.den-class
#' @param  x 	a (data) object used to select a method.
#' @param  ... 	further arguments passed to or from other methods.
#' @return Returns an S3 object of class \code{np.den} (locpol den + bin den + grid par.). 
#' A \code{\link{bin.den}} object with the additional (some optional) 3 components:
#' \item{est}{vector or array (dimension \code{nbin}) with the local polynomial density estimates. }
#' \item{locpol}{a list with 6 components:
#' \itemize{
#'    \item{\code{degree} degree of the polinomial.}
#'    \item{\code{h} bandwidth matrix.}
#'    \item{\code{rm} residual mean (of the escaled bin counts).}
#'    \item{\code{rss} sum of squared residuals (of the escaled bin counts).}
#'    \item{\code{ncv} number of cells ignored (in each dimension).}
#' }}
#' \item{deriv}{(if requested) matrix of first derivatives.} 
#' @seealso \code{\link{bin.den}}, \code{\link{binning}}, \code{\link{h.cv}}, 
#' \code{\link{data.grid}}.
#' @examples
#' bin.den <- binning(earthquakes[, c("lon", "lat")], nbin = c(30,30))
#' h.den <- h.cv(bin.den) 
#' den <- np.den(bin.den, h = h.den$h)
#' plot(den, main = 'Estimated log(density)')
#' @export
#····································································
np.den <- function(x, ...) {
  UseMethod("np.den")
}


#····································································
# np.den.default(x, nbin = NULL, h = NULL, degree = 1 + as.numeric(drv), drv = FALSE, ncv = 0, ...) {    
#····································································
#' @rdname np.den
#' @method np.den default
#' @inheritParams locpol.default
#' @details Standard generic function with a default method (interface to the 
#' fortran routine \code{lp_data_grid}), in which argument \code{x} 
#' is a vector or matrix of covariates (e.g. spatial coordinates).
#' In this case, the data are binned (calls \code{\link{bin.den}}) and the local fitting
#' procedure is applied to the scaled bin counts (calls \code{\link{np.den.bin.den}}).
#' 
#' If parameter \code{nbim} is not specified is set to \code{rep(25, ncol(x))}. 
#'
#' A multiplicative triweight kernel is used to compute the weights.
#' 
#' If \code{ncv > 1}, estimates are computed by leaving out cells with indexes within 
#' the intervals \eqn{[x_i - ncv + 1, x_i + ncv - 1]}, at each dimension i, where \eqn{x} 
#' denotes the index of the estimation position. 
#' @references
#' Wand, M.P. and Jones, M.C. (1995) \emph{Kernel Smoothing}. Chapman and Hall, London.
#' @export
np.den.default <- function(x, nbin = NULL, h = NULL, degree = 1 + as.numeric(drv), 
                            drv = FALSE, ncv = 0, ...) {    
  xbin <- bin.den(x, nbin)
  return(locpol.bin.den(xbin, h = h, degree = degree, drv = drv, ncv = ncv, ...))
}


#····································································
# np.den.bin.den(x, h, degree, drv, ncv, ...) ----    
#····································································
#' @rdname np.den
#' @method np.den bin.den
#' @export
np.den.bin.den <- locpol.bin.den 



#····································································
#' @rdname np.den
#' @method np.den bin.data
#' @export
np.den.bin.data <- function(x, h = NULL, degree = 1 + as.numeric(drv), drv = FALSE, 
                              ncv = 0, ...) {    
#····································································
    result <- np.den.bin.den(x, h = h, degree = degree, drv = drv, ncv = ncv, ...)
    oldClass(result) <- c("np.den", "bin.data", "bin.den", "data.grid")
    return(result)
}


#····································································
#' @rdname np.den
#' @method np.den svar.bin
#' @export
np.den.svar.bin <- function(x, h = NULL, degree = 1 + as.numeric(drv), drv = FALSE, 
                              ncv = 0, ...) {    
#····································································
    result <- np.den.bin.den(x, h = h, degree = degree, drv = drv, ncv = ncv, ...)
    oldClass(result) <- c("np.den", "svar.bin", "bin.data", "bin.den", "data.grid")
    return(result)
}