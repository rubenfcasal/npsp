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
#  or covariance matrix of the data.
#' @param  ngrid number of grid nodes in each dimension. 
#' @param lp.resid  residuals (defaults to \code{residuals(object)}).
#' @param intermediate logical, determines whether the intermediate computations 
#' are included in the output (component \code{kriging}; see Value). 
#' These calculations can be reused, e.g. for bootstrap.
#' @return 
#' \code{np.kriging()}, and \code{kriging.simple()} when \code{newx} defines 
#' gridded data (extends \code{grid.par} or \code{data.grid} classes),
#' returns an S3 object of class \code{krig.grid} (kriging results + grid par.). 
#' A \code{\link{data.grid}} object with the additional (some optional) components:
#' \item{kpred}{vector or array (dimension \code{$grid$n}) with the kriging predictions. }
#' \item{ksd}{vector or array with the kriging standard deviations. }
#' \item{kriging}{(if requested) a list with 4 components:
#' \itemize{
#'    \item{\code{lambda} matrix of kriging weights (columns correspond with predictions 
#'    and rows with data)).}
#'    \item{\code{cov.est} (estimated) covariance matrix of the data.}
#'    \item{\code{chol} Cholesky factorization of \code{cov.est}.}
#'    \item{\code{cov.pred} matrix of (estimated) covariances between data (rows) 
#'    and predictions (columns).}
#' }}
#' When \code{newx} is a matrix of coordinates (where each row is a prediction location),
#' \code{kriging.simple()} returns a list with the previous components (\code{kpred}, \code{ksd} 
#' and, if requested, \code{kriging}).
#' @examples 
#' geomod <- np.fitgeo(aquifer[,1:2], aquifer$head)
#' krig.grid <- np.kriging(geomod, ngrid = c(96, 96)) # 9216 locations
#' old.par <- par(mfrow = c(1,2))
#' simage(krig.grid, 'kpred', main = 'Kriging predictions', 
#'        xlab = "Longitude", ylab = "Latitude", reset = FALSE )
#' simage(krig.grid, 'ksd', main = 'Kriging sd', xlab = "Longitude", 
#'        ylab = "Latitude" , col = hot.colors(256), reset = FALSE)
#' par(old.par)
#' @seealso \code{\link{np.fitgeo}}, \code{\link{locpol}}, \code{\link{np.svar}}.
#' @export
#····································································
np.kriging.default <- function(object, svm, lp.resid = NULL, ngrid = object$grid$n, 
                               intermediate = FALSE, ...) {
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
  krig.grid <- kriging.simple(x = object$data$x, y = lp.resid, newx = object, 
                              svm = svm, intermediate = intermediate)
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
np.kriging.np.geo <- function(object, ngrid = object$grid$n, 
                              intermediate = FALSE, ...) {
  np.kriging.default(object, object$svm, object$residuals, ngrid = ngrid, 
                       intermediate = intermediate, ...)
}


#' @rdname np.kriging
#' @inheritParams locpol.default
#' @param x vector/matrix with data locations
#'    (each component/row is an observation location). 
#' @param newx vector/matrix with the (irregular) locations to predict 
#'    (each component/row is a prediction location). 
#'    or an object extending \code{\link{grid.par}}-\code{\link{class}}
#'    (\code{\link{data.grid}}).
#' @export
#····································································
# res <- krig.grid$kriging
kriging.simple <- function(x, y, newx, svm, intermediate = FALSE) {
  x <- as.matrix(x)
  if ( !identical(length(y), nrow(x)) )
    stop("arguments 'y' and 'x' do not have the same length.")  
  gridded <- inherits(newx, "grid.par") || inherits(newx, "data.grid")
  if (gridded) {
    grid <- if(inherits(newx, "grid.par")) newx else newx$grid
    newx <- coords(newx)
  }
  res <- .kriging.simple.solve(x = x, newx = newx, svm = svm)
  kpred <- drop(y %*% res$lambda)
  ksd2 <- with(res, cov.est[1,1] - colSums( lambda * cov.pred ))
  # ksd2 <- with(res, cov.est[1,1] - spam::colSums.spam( lambda * cov.pred ))
  # ksd2 <- with(res, spam::diag(cov.est) - spam::colSums.spam( lambda * cov.pred ))
  if (gridded) {
    result <- data.grid(kpred = kpred, ksd = sqrt(ksd2), grid = grid)
    oldClass(result) <- c("krig.grid", oldClass(result))
  } else {
    result <- list(kpred = kpred, ksd = sqrt(ksd2)) 
  }
  if (intermediate) result$kriging <- res
  return(result)
}


#' @rdname npsp-internals
#' @inheritParams kriging.simple
#' @keywords internal
#····································································
.kriging.simple.solve <- function(x, newx, svm) {
# .DNRM2_R <- npsp:::.DNRM2_R  
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
  # cov.est <-  varcov(svm, coords = x)
  cov.est <-  spam::as.spam(varcov(svm, coords = x))
  # Covariances between prediction grid points and data
  dist2 <- as.numeric(.DNRM2_R(newx, x))
  cov.pred <- matrix(covar(svm, dist2, discretize = TRUE), nrow = nrow(x))
  # cov.pred <- spam::as.spam(matrix(covar(svm, dist2, discretize = TRUE), nrow = nrow(x)))
  # Kriging system
  # library(spam)
  if (!requireNamespace("spam")) stop("'spam' package must be installed...")
  # chol <- spam::chol.spam(spam::as.spam(cov.est))
  chol <- spam::chol.spam(cov.est)
  lambda <- spam::solve.spam(chol, spam::as.spam(cov.pred))
  # lambda <- spam::solve.spam(chol, cov.pred)
  # res <- solve.spam(chol, as.spam(cov.pred))
  return(list(lambda = lambda, cov.est = cov.est, chol = chol, cov.pred = cov.pred))
}


