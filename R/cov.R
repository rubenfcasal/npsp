#····································································
#   cov.R (npsp package)
#····································································
#   covar(x, h, ...)  S3 generic
#       covar.svarmod(x, h, sill = x$sill, ...)
#       covar.np.svar(x, h, sill = NULL, ...)
#   varcov(x, coords, ...)  S3 generic
#       varcov.isotropic(x, coords, sill = x$sill, range.taper, discretize, ...)
#       varcov.np.svar(x, coords, sill = max(x$est), range.taper, ...) 
#
#   (c) R. Fernandez-Casal
#   Created: Apr 2013
#
#   NOTE: Press Ctrl + Shift + O to show document outline in RStudio
#····································································


#····································································
#   covar(x, h, ...) ----
#····································································
#' Covariance values 
#' 
#' Computes covariance values (or pseudo-covariances) given a variogram model
#' or covariance estimates given a semivariogram estimate.
#' 
#' @param  x  variogram model (\code{\link{svarmod}} object) or semivariogram estimate.
#' @param  h  vector (isotropic case) or matrix of lag values.
#' @param  ... further arguments passed to or from other methods.
#' @param  sill  (theoretical or estimated) variance \eqn{C(0) = \sigma^2} or pseudo-sill (unbounded variograms).
#' @return
#' A vector of (pseudo) covariance values \eqn{C(h_i) = \sigma^2 - \gamma(h_i)} or covariance estimates.
#' @seealso
#' \code{\link{sv}}, \code{\link{varcov}}.
#' @export
#····································································
covar  <- function(x, h, ...) {
  UseMethod("covar")
}


#····································································
#' @rdname covar  
#' @method covar svarmod
#' @param  discretize logical. If \code{TRUE} the variogram is previously discretized. 
#' @export
#····································································
covar.svarmod  <- function(x, h, sill = x$sill, discretize = FALSE, ...) {
    if (!inherits(x, "svarmod"))
        stop("argument 'x' must be of class (or extending) 'svarmod'.")
    result <- sv(x, h, discretize = discretize, ...)
    if (is.na(sill)) {
        warning("semivariogram model 'x' is not bounded (computing pseudo-covariances).")
        sill = max(result)
    }  
    return(sill - result)
}


#····································································
#' @rdname covar  
#' @method covar np.svar
#' @export
#····································································
covar.np.svar  <- function(x, h, sill = NULL, ...) {
    if (!inherits(x, "np.svar"))
        stop("argument 'x' must be of class (or extending) 'np.svar'.")
    result <- interp(x, data.ind ='est', newx = h)$y  # COIDADO Error in interp.grid.par 'data' has missing values
    result[h < sqrt(.Machine$double.eps)] <- 0        # COIDADO CO EFECTO NUGGET
    if (is.null(sill)) {
        warning("'sill' parameter is set to the maximum semivariance value (computing pseudo-covariances).")
        sill = max(result)
    }  
    return(sill - result)
}



#····································································
#   varcov(x, coords, ...) ----
#····································································
#' Covariance matrix 
#' 
#' Computes the covariance matrix a corresponding to a set of spatial locations 
#' given a variogram model or a semivariogram estimate.
#' 
#' @param  x  variogram model (\code{\link{svarmod}} object) or semivariogram estimate.
#' @param  coords  matrix of coordinates (columns correspond with dimensions and rows with data).
#' @param  ... further arguments passed to or from other methods.
#' @return
#' The covariance matrix of the data. 
#' @seealso
#' \code{\link{sv}}, \code{\link{covar}}.
#' @export
#····································································
varcov <- function(x, coords, ...)  {
  UseMethod("varcov")
}


#····································································
#' @rdname varcov
#' @method varcov isotropic
#' @param  sill  (theoretical or estimated) variance \eqn{C(0) = \sigma^2} or pseudo-sill (unbounded variograms).
#' @param  range.taper (optional) if provided, covariances corresponding to 
#' distances larger than this value are set to 0.
#' @param  discretize logical. If \code{TRUE} (default), the variogram is (previously) discretized. 
#' @export
varcov.isotropic <- function(x, coords, sill = x$sill, range.taper, 
                             discretize = nrow(coords) > 256, ...) {
  #····································································
  if(missing(coords)) 
    stop("argument 'coords' (spatial coordinates) must be provided")
  coords <- as.matrix(coords)
  n <- nrow(coords)
  dists <- as.vector(dist(coords))   # lower triangle of the distance matrix
  # if (discretize) x <- svar.grid(x)
  if(!missing(range.taper)) {
    index <- dists <= range.taper
    covs  <- numeric(length(dists))
    covs[index] <- covar(x, dists[index], sill = sill, discretize = discretize)        
  } else covs <- covar(x, dists, sill = sill, discretize = discretize)    
  res <- matrix(0, n, n)
  res[row(res) > col(res)] <- covs 
  res <- res + t(res)
  diag(res) <- sill
  return(res)    
}  


#····································································
#' @rdname varcov
#' @method varcov np.svar
#' @export
varcov.np.svar <- function(x, coords, sill = max(x$est), range.taper = x$grid$max, ...) {
#···································································· 
    n <- nrow(coords)
    dists <- as.vector(dist(coords))   # lower triangle of the distance matrix    
    if(range.taper > 0) {
        index <- dists <= range.taper
        covs  <- numeric(length(dists))
        covs[index] <- covar(x, dists[index], sill = sill)        
    } else covs <- covar(x, dists, sill = sill)    
    res <- matrix(0, n, n)
    res[row(res) > col(res)] <- covs 
    res <- res + t(res)
    diag(res) <- sill
    return(res)    
}  


