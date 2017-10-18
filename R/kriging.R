#--------------------------------------------------------------------
#   kriging.R (npsp package)
#--------------------------------------------------------------------
#   kriging.np(lp, svm, lp.resid = residuals(lp))
#   kriging.simple(x, y, newx, svm)
#   kriging.simple.solve(x, newx, svm)  
#
#   (c) R. Fernandez-Casal         Last revision: Sep 2017
#--------------------------------------------------------------------
# PENDENTE:
#   - spam = TRUE
#   - newx = FALSE -> CV? chol2inv.spam
#--------------------------------------------------------------------

#' Nonparametric (residual) kriging
#' 
#' Compute simple kriging or residual kriging predictions 
#' (and also the corresponding simple kriging standard errors  
#  *and cross-validation measures, not yet...*
#' ). Currently, only global (residual) simple kriging is implemented.
#' @aliases kriging 
#' @param  lp local polynomial estimate of the trend function (object of class 
#'  \code{\link{locpol.bin}}).
#' @param svm semivariogram model (of class extending \code{\link{svarmod}})
#' or covariance matrix of the data.
#' @param lp.resid  residuals (defaults to \code{residuals(lp)})
#' @seealso \code{\link{locpol}}, \code{\link{np.svar}}.
#' @export
kriging.np <- function(lp, svm, lp.resid = residuals(lp)) {
  # OJO: Si lp data is masked... 
  krig.grid <- kriging.simple(x = lp$data$x, y = lp.resid, newx = lp$grid, svm = svm)
  krig.grid$kpred <- lp$est + krig.grid$kpred
  if (!is.null(lp$mask)) 
    is.na(krig.grid$ksd) <- !lp$mask
  krig.grid$lp <- lp
  krig.grid$residuals <- lp.resid
  return(krig.grid)
}

#' @rdname kriging.np
#' @inheritParams locpol.default
#' @param newx vector/matrix with the (irregular) locations to predict 
#'    (columns correspond with dimensions and rows with locations) 
#'    or an object extending \code{\link{grid.par}}-\code{\link{class}}
#'    (\code{\link{data.grid}}).
#' @export
kriging.simple <- function(x, y, newx, svm) {
  x <- as.matrix(x)
  if ( !identical(length(y), nrow(x)) )
    stop("arguments 'y' and 'x' do not have the same length.")  
  grid <- NULL
  if((inherits(newx, "grid.par")) || (inherits(newx, "data.grid"))) {
    grid <- if(inherits(newx, "grid.par")) newx else newx$grid
    newx <- coords(grid)
  }
  res <- kriging.simple.solve(x = x, newx = newx, svm = svm)
  kpred <- drop(y %*% res$lambda)
  ksd2 <- with(res, cov.est[1,1] - colSums( lambda * cov.pred ))
  if (!is.null(grid)) {
    krig.grid <- data.grid(kpred = kpred, ksd = sqrt(ksd2), grid = grid)
    krig.grid$kriging <- res
    krig.grid$svm <- svm      
    return(krig.grid)    
  } else {
    return(list(kpred = kpred, ksd = sqrt(ksd2), 
                kriging = res, svm = svm))
  }
}

#' @rdname kriging.np
#' @export
kriging.simple.solve <- function(x, newx, svm) {
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
  # cpu.time(total = FALSE)
  ## Time of last operation: 
  ##    user  system elapsed 
  ##   20.00    1.03   22.10
  # cov.inv.est <- chol2inv(U.est) # Interface to the LAPACK routine DPOTRI.
  # lambda <- cov.inv.est %*% cov.pred
  # cpu.time(total = FALSE)
  ## Time of last operation: 
  ##    user  system elapsed 
  ##   16.03    0.12   16.25
  # L.est <- t(U.est)
  # lambda.pred2 <- backsolve(U.est,
  #                forwardsolve(L.est, cov.pred))
  # cpu.time(total = FALSE)
  ## Time of last operation: backsolve
  ##    user  system elapsed
  ##   20.57    0.08   20.65
  # lambda <- .DPOSV_R(cov.est, cov.pred)
  # cpu.time(total = FALSE)
  # # Time of last operation:
  # #  user  system elapsed
  # # 29.50    0.02   44.99
  # library(Matrix)
  # lambda <- solve(as(cov.est, "dpoMatrix"), cov.pred)
  # cpu.time(total = FALSE)
  # # Time of last operation: 
  # #  user  system elapsed 
  # # 28.72    0.14   37.75
  # library(spam)
  if (!requireNamespace("spam")) stop("'spam' package must be installed...")
  chol <- spam::chol.spam(spam::as.spam(cov.est))
  lambda <- spam::solve.spam(chol, cov.pred)
  # res <- solve.spam(chol, as.spam(cov.pred))
  return(list(lambda = lambda, cov.est = cov.est, chol = chol, cov.pred = cov.pred))
}