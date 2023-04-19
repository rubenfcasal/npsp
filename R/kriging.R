#····································································
#   kriging.R (npsp package)
#····································································
#   np.kriging(object, ...)
#       np.kriging.default(object, svm, lp.resid, ngrid, ...)
#       np.kriging.np.geo(object, ngrid, ...)
#   kriging.simple(x, y, newx, svm)
#   kriging.simple.solve(x, newx, svm)  
#
#   (c) R. Fernandez-Casal
#
#   NOTE: Press Ctrl + Shift + O to show document outline in RStudio
#····································································
# PENDENTE:
#   - spam = TRUE
#   - newx = FALSE -> CV? chol2inv.spam
#····································································

#····································································
# np.kriging(object, ...) ---- 
#····································································
#' Nonparametric (residual) kriging
#' 
#' Compute simple kriging or residual kriging predictions 
#' (and also the corresponding simple kriging standard errors  
#  *and cross-validation measures, not yet...*
#' ). Currently, only global (residual) simple kriging is implemented.
#' @aliases kriging 
#' @param  object object used to select a method: 
#' local polynomial estimate of the trend (class \code{\link{locpol.bin}}) 
#' or nonparametric geostatistical model (class extending \code{\link{np.geo}}).
#' @param ... further arguments passed to or from other methods.
#' @export
#····································································
np.kriging <- function(object, ...) {
  UseMethod("np.kriging")
} # S3 generic function np.kriging



#' @rdname np.kriging  
#' @method np.kriging default
#' @param svm semivariogram model (of class extending \code{\link{svarmod}}).
# or covariance matrix of the data.
#' @param  ngrid number of grid nodes in each dimension. 
#' @param lp.resid  residuals (defaults to \code{residuals(object)}).
#' @seealso \code{\link{np.fitgeo}}, \code{\link{locpol}}, \code{\link{np.svar}}.
#' @export
#····································································
np.kriging.default <- function(object, svm, lp.resid = NULL, ngrid = object$grid$n, ...) {
  stopifnot(inherits(object, "locpol.bin" ))
  if (!is.null(object$mask) && is.null(object$window)) # Interpolar mask?
      warning("A resized 'object' can not be masked (`ngrid != object$grid$n`).")
  if (lp.hd <- any(ngrid != object$grid$n)) {
    bin <- with(object$data, binning(x, y, nbin = ngrid, window = object$window))
    object <- locpol(bin, h = object$locpol$h)
    lp.resid <- residuals(object)
  }  
  if (is.null(lp.resid))
    lp.resid <- if(inherits(svm, "fitsvar")) svm$esv$data$y else residuals(object)
  krig.grid <- kriging.simple(x = object$data$x, y = lp.resid, newx = object, svm = svm)
  krig.grid$kpred <- object$est + krig.grid$kpred
  krig.grid$trend <- object$est
  if (!is.null(object$mask)) is.na(krig.grid$ksd) <- !object$mask
  # krig.grid$lp <- object
  # krig.grid$residuals <- lp.resid
  return(krig.grid)
}


#' @rdname np.kriging  
#' @method np.kriging np.geo
#' @export
#····································································
np.kriging.np.geo <- function(object, ngrid = object$grid$n, ...) {
  return(
    np.kriging.default(object, object$svm, object$residuals, ngrid = ngrid, ...)
  )}


#' @rdname np.kriging
#' @inheritParams locpol.default
#' @param newx vector/matrix with the (irregular) locations to predict 
#'    (columns correspond with dimensions and rows with locations) 
#'    or an object extending \code{\link{grid.par}}-\code{\link{class}}
#'    (\code{\link{data.grid}}).
#' @export
#····································································
kriging.simple <- function(x, y, newx, svm) {
  x <- as.matrix(x)
  if ( !identical(length(y), nrow(x)) )
    stop("arguments 'y' and 'x' do not have the same length.")  
  gridded <- inherits(newx, "grid.par") || inherits(newx, "data.grid")
  if (gridded) {
    grid <- if(inherits(newx, "grid.par")) newx else newx$grid
    newx <- coords(newx)
  }
  res <- kriging.simple.solve(x = x, newx = newx, svm = svm)
  kpred <- drop(y %*% res$lambda)
  ksd2 <- with(res, cov.est[1,1] - colSums( lambda * cov.pred ))
  if (gridded) {
    krig.grid <- data.grid(kpred = kpred, ksd = sqrt(ksd2), grid = grid)
    # krig.grid$kpred <- kpred
    # krig.grid$ksd <- sqrt(ksd2)
    krig.grid$kriging <- res
    oldClass(krig.grid) <- c("krig.grid", oldClass(krig.grid))
    return(krig.grid)    
  } else {
    return(list(kpred = kpred, ksd = sqrt(ksd2), 
                kriging = res)) 
  }
}


#' @rdname np.kriging
#' @export
#····································································
kriging.simple.solve <- function(x, newx, svm) {
# Evitar recalcular matriz de covarianza y factorizacion
  # spam = TRUE
  # newx = FALSE -> CV? chol2inv.spam
  if (!inherits(svm, "svarmod"))
    stop("'svm' must be a semivariogram model ('svarmod' class).")
  x <- as.matrix(x)
  newx <- as.matrix(newx)
  if(ncol(x) != ncol(newx))
    stop("arguments 'x' and 'newx' have incompatible dimensions")
  # Covariance matrix of the data
  cov.est <-  varcov(svm, coords = x)
  # Covariances between prediction grid points and data
  dist2 <- as.numeric(.DNRM2_R(newx, x))
  cov.pred <- matrix(covar(svm, dist2, discretize = TRUE), nrow = nrow(x))
  # Kriging system
  # U.est <- chol(cov.est) # Interface to the LAPACK routines DPOTRF and DPSTRF,
  ## Time of last operation: 
  ##    user  system elapsed 
  ##   20.00    1.03   22.10
  # cov.inv.est <- chol2inv(U.est) # Interface to the LAPACK routine DPOTRI.
  # lambda <- cov.inv.est %*% cov.pred
  ## Time of last operation: 
  ##    user  system elapsed 
  ##   16.03    0.12   16.25
  # L.est <- t(U.est)
  # lambda.pred2 <- backsolve(U.est,
  #                forwardsolve(L.est, cov.pred))
  ## Time of last operation: backsolve
  ##    user  system elapsed
  ##   20.57    0.08   20.65
  # lambda <- .DPOSV_R(cov.est, cov.pred)
  # # Time of last operation:
  # #  user  system elapsed
  # # 29.50    0.02   44.99
  # library(Matrix)
  # lambda <- solve(as(cov.est, "dpoMatrix"), cov.pred)
  ## Time of last operation: 
  ##  user  system elapsed 
  ## 28.72    0.14   37.75
  # library(spam)
  if (!requireNamespace("spam")) stop("'spam' package must be installed...")
  chol <- spam::chol.spam(spam::as.spam(cov.est))
  lambda <- spam::solve.spam(chol, cov.pred)
  # res <- solve.spam(chol, as.spam(cov.pred))
  return(list(lambda = lambda, cov.est = cov.est, chol = chol, cov.pred = cov.pred))
}
