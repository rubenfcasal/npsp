#--------------------------------------------------------------------
#   mask.R (npsp package)
#--------------------------------------------------------------------
#   npsp.tolerance(level, warn)
#   mask  S3 generic
#           mask.default(x, tol.mask = 0, ...)
#           mask.bin.den(x, mask, warn, ...)
#           mask.bin.data(x, mask, set.NA, warn, ...)
#           mask.locpol.bin(x, mask, set.NA, warn, ...)
#
#   (c) Ruben Fernandez-Casal
#   Created: Jul 2015
#--------------------------------------------------------------------
#   PENDENTE:
#       - Exemplos mask
#       - MASK DE MASKED LOCPOL?
#--------------------------------------------------------------------


#--------------------------------------------------------------------
#' npsp Tolerances
#'
#' Returns a (convergence, taper, approximation,...) tolerance.
#' Defaults to \code{.Machine$double.eps^(1/level)}, typically about \code{1e-8}.
#' @param level numerical,
#' @param warn logical; If \code{TRUE} (the default) a warning message is issued
#' when \code{level < 1}.
#' @return Returns \code{.Machine$double.eps^(1/level)} if \code{level >= 1},
#' in other case \code{1 - .Machine$double.eps}.
#' @seealso \code{\link{.Machine}}
#' @examples
#' curve(npsp.tolerance, 1, 1000)
#' abline(h = npsp.tolerance(0, FALSE), lty = 2)
#' @export
npsp.tolerance <- function(level = 2, warn = TRUE)
    if (any(index <- level < 1)) {
        if (warn) warning("'level' should be >= 1")
        ifelse(index, 1 - .Machine$double.eps, .Machine$double.eps^(1/level))
    } else return(.Machine$double.eps^(1/level))
#--------------------------------------------------------------------


#--------------------------------------------------------------------
# mask(x, ...)
#--------------------------------------------------------------------
#' Mask methods
#'
#' Filters the data that satisfy a condition.
#'
#' @param x object used to select a method (binned data, ...).
#' @param ... further arguments passed to or from other methods
#' @seealso \code{\link{locpol}}, \code{\link{locpolhcv}}, \code{\link{binning}},
#' \code{\link{np.svar}}, \code{\link{npsp.tolerance}}.
#' @export
mask <- function(x, ...) UseMethod("mask")
#--------------------------------------------------------------------


#--------------------------------------------------------------------
#' @rdname mask
#' @method mask default
#' @param tol.mask tolerance.
#' @return \code{mask.default} returns the logical vector \code{x > tol.mask}.
#' @examples
#' mask(1:10, 5)
#' @export
mask.default <- function(x, tol.mask = 0, ...)  as.numeric(x) > tol.mask
# NOTA: E importante que mask sexa vector...
#--------------------------------------------------------------------


#--------------------------------------------------------------------
#' @rdname mask
#' @method mask data.grid
#' @param  set.NA logical; If \code{TRUE}, the values corresponding
#' to masked cells are set to \code{NA}.
#' @param  window spatial window (values outside this window will be masked), currently an sp-object of class 
#' extending \code{\link[sp:SpatialPolygons-class]{SpatialPolygons}}.
#' @export
mask.data.grid <- function(x, mask = NULL, window = NULL,
                          set.NA = FALSE, warn = FALSE, ...) 
{ #--------------------------------------------------------------------
  if(!is.null(window)){
    # if (!is.null(mask)) warning("'mask' argument is ignored.")
    spp.grid <- SpatialPoints(coords(x))
    proj4string(spp.grid) <- proj4string(window) # CRS("+init=epsg:28992 +units=km")
    mask <- array(!is.na(over(spp.grid, as(window, 'SpatialPolygons'))), 
                  dim = x$grid$n)
    x$window <- window
  } else 
    if(is.null(mask)) 
      mask <- mask.default(x[[1]], npsp.tolerance(2)) else
        if (!is.logical(mask) || length(mask) != length(x$binw))
          stop("'mask' must be a logical vector of appropriate order.")
  if (!is.null(x$mask)) {
    mask <- mask & x$mask
    if (warn) warning("'x' is already masked.")
  }
  # A botch...
  if (set.NA || warn)
    index <- !sapply(x, is.list) & (sapply(x, length) == prod(dim(x)))
  if (warn && sapply(x[index], function(d) any(!mask & (d > 0))))
    warning("some data will be masked...")
  if (set.NA) 
    sapply(x[index], function(d) is.na(d) <- !mask)
  x$mask <- mask
  return(x)
} #------------------------------------------------------------------



#--------------------------------------------------------------------
#' @rdname mask
#' @method mask bin.den
#' @param mask logical; vector (or array) indicating the selected values (not masked).
#' @param warn logical; If \code{TRUE} a warning message is generated when original data is masked.
#' @return \code{mask.bin.den}, \code{mask.bin.data} and \code{mask.locpol.bin}
#' return an object of the same class as \code{x} with the additional component \code{$mask}
#' and optionally \code{$window}.
#' @examples
#' bin <- binning(aquifer[,1:2], aquifer$head, nbin = c(41,41), set.NA = TRUE)
#' str(mask(bin, mask(bin$binw), warn = TRUE))
#' str(mask(bin, mask(bin$binw, 1)))
#' @export
mask.bin.den <- function(x, mask = mask.default(x$binw, npsp.tolerance(2)), window = NULL,
                         set.NA = FALSE, warn = TRUE, ...) 
{ #------------------------------------------------------------------
# PENDENTE: establecer binw[mask] <- -1 en estimacion locpol
# PENDENTE: Ollo con binw = -1
# COIDADO: mask <- as.numeric(x$binw) > 0    # hcv.data => predict
    if (!is.null(x$mask) && warn) warning("'x' is already masked.")
     res <- mask.data.grid(x, mask = mask, window = window,
                           set.NA = FALSE, warn = FALSE)
    nomask <- which(!res$mask & x$binw > 0)
    if (set.NA) {
      is.na(res$biny) <- !mask
      res$binw[!mask] <- -1
      if (warn && length(nomask)) warning("\nsome data will be masked...")
    } else {
      if(length(nomask)) res$nomask <- nomask
      # if (warn) warning("some data will not be masked...")
    }  
    return(res)
} #------------------------------------------------------------------


#--------------------------------------------------------------------
#' @rdname mask
#' @method mask bin.data
#' @param  filter.lp logical; If \code{TRUE}, masked nodes will be leaved out
#' in local polynomial estimation.
#' @export
# mask.bin.data <- function(x, window = NULL, mask = mask.default(x$binw, npsp.tolerance(2)),
mask.bin.data <- function(x, mask = NULL, window = NULL,
                          set.NA = FALSE, warn = FALSE, filter.lp = TRUE, ...) 
{ #------------------------------------------------------------------
    # if (filter.lp) warn <- TRUE
    if (is.null(window) & is.null(mask)) mask <- mask.default(x$binw, npsp.tolerance(2))
    res <- mask.bin.den(x, mask = mask, window = window, warn = warn, set.NA = set.NA)
    # if (set.NA) is.na(res$biny) <- !res$mask
    if (filter.lp && !set.NA) {
      if (warn && !is.null(res$nomask)) warning("\nsome data will not be masked...")
      mask <- res$mask
      mask[res$nomask] <- TRUE
      res$binw[!mask] <- -1
    }  
    return(res)
} #------------------------------------------------------------------



#--------------------------------------------------------------------
#' @rdname mask
#' @method mask locpol.bin
#' @export
mask.locpol.bin <- function(x, mask = mask.default(x$binw, npsp.tolerance(2)), window = NULL,
                          set.NA = FALSE, warn = TRUE, filter.lp = TRUE, ...) 
{ #------------------------------------------------------------------
    res <- mask.bin.data(x, mask = mask, window = window, 
                         set.NA = set.NA, warn = FALSE, filter.lp = filter.lp)
    if (set.NA | filter.lp) {
      is.na(res$est) <- !res$mask
      if (warn && !is.null(res$nomask)) warning("\nsome data will be masked...")
      res$est[res$nomask] <- NA
    }  
    return(res)
} #------------------------------------------------------------------
