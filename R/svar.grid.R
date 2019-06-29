#--------------------------------------------------------------------
#   svar.grid (npsp package)
#--------------------------------------------------------------------
#   svar.grid(svar, log, ...)  S3 generic
#       svar.grid.fitsvar(svar, log, ...)
#       svar.grid.svarmod(svar, n, min, max, log, ...)
#
#   (c) R. Fernandez-Casal
#   Created: May 2017                          Last changed: 
#--------------------------------------------------------------------


#--------------------------------------------------------------------
#   svar.grid(svar, log, ...)
#--------------------------------------------------------------------
#' Discretize a (semi)variogram model
#' 
#' Discretizes a variogram model (to speed up variogram evaluation). 
#' Constructor function of the \code{svar.grid-\link{class}}.
#' 
#' @param  svar  (fitted) variogram model (a \code{\link{svarmod}} 
#' or \code{\link{fitsvar}} object).
#' @param  log logical. If \code{TRUE}, the variogram is discretized
#' in  (base 2) logarithmic scale.
#' @param  ... further arguments passed to or from other methods.
#' @export
#----------------------------------------------------------------------
svar.grid  <- function(svar, log = TRUE, ...) UseMethod("svar.grid")
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# @rdname svar.grid
# @method svar.grid fitsvar
# @export
# svar.grid.fitsvar <- function(svar, log = TRUE, ...){
#   #----------------------------------------------------------------------
#   if (!inherits(svar, "fitsvar"))
#     stop("function only works for objects of class (or extending) 'fitsvar'.")
#   # if (esv$svar$type != "isotropic")
#   if (svar$esv$grid$nd != 1)
#     stop("'svar' must be isotropic.")
#   n <- svar$esv$grid$n
#   u <- svar$fit$u
#   u <- c(10*.Machine$double.eps, u, 1.5*u[length(u)])
#   if (log) u <- log2(u)
#   w <- svar$fit$w
#   res <- binning(u, c(w[1]/2, w, 2*w[length(w)]), nbin = 2*n)
#   u <- drop(coords(res))
#   if (log) u <- 2^u
#   res$sv <- sv(svar, u, discretize = FALSE)
#   res$log = log
#   res <- c(res, svar)
#   oldClass(res) <- c("svar.grid", "bin.data", 
#                      "bin.den", "data.grid", oldClass(svar))
#   return(res)
# }


#----------------------------------------------------------------------
#' @rdname svar.grid  
#' @method svar.grid svarmod
#' @param  n number of lags. Defaults to 256. 
#' @param  min minimun lag. Defaults to \code{10*.Machine$double.eps}. 
#' @param  max maximum lag. Defaults to \code{1.1*svar$range}. 
#' @export
svar.grid.svarmod <- function(svar, log = TRUE, n = 256, min = 10*.Machine$double.eps, 
                              max = 1.1*svar$range, ...){
  #----------------------------------------------------------------------
  if (!inherits(svar, "svarmod"))
    stop("function only works for objects of class (or extending) 'svarmod'.")
  if (!inherits(svar,"isotropic"))
    stop("'svar' must be isotropic.")
  if (log) {
    min <- log2(min)
    max <- log2(max)
  } 
  res <- data.grid(binw = rep(1, prod(n)), 
                   grid = grid.par(n = n, min = min, max = max))
  u <- drop(coords(res))
  if (log) u <- 2^u
  res$sv <- sv(svar, u, discretize = FALSE)
  res$log = log
  res <- c(res, svar)
  oldClass(res) <- c("svar.grid", "bin.den", "data.grid", oldClass(svar))
  return(res)
}
