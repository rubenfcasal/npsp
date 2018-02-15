#--------------------------------------------------------------------
#   np.geo.R (npsp package)
#--------------------------------------------------------------------
#   np.geo and fitgeo S3 classes and methods
#       plot.np.geo(x)
#
#   np.fitgeo() S3 generic
#       np.fitgeo.default(x)
#       np.fitgeo.locpol.bin(x)
#       np.fitgeo.fitgeo(x)
#
# PENDENTE:
#   - exemplos
#   - Cuidado h automático
#     puede ser recomendable emplear MCV 
#     si no se dispone de más información...
#
#   (c) R. Fernandez-Casal            Last revision: Feb 2018
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# np.geo <- function(lp, svm, svm0 = NULL, nbin = lp$grid$n) {
#--------------------------------------------------------------------
#' Nonparametric geostatistical model (S3 class "np.geo")
#' 
#' Defines a nonparametric geostatistical model 
#' (not intended to be used regularly; see \code{\link{np.fitgeo}}). 
#' Constructor function of the \code{np.geo} and \code{fitgeo} S3 \code{\link{class}}es.
#' @aliases np.geo-class, fitgeo-class
#' @param  lp local polynomial estimate of the trend function (object of class 
#'  \code{\link{locpol.bin}}).
#' @param  svm (fitted) variogram model (object of class 
#'  \code{\link{fitsvar}} or \code{\link{svarmod}}).
#' @param  svm0 (fitted) residual variogram model (object of class 
#'  \code{\link{fitsvar}} or \code{\link{svarmod}}).
#' @param  nbin number of bins on each dimension.
#' @return Returns an object of \code{\link{class}} \code{np.geo}
#' (extends \code{\link{locpol.bin}}), the \code{lp} argument with 
#' the others and the vector of residuals as additional components.
#' @seealso \code{\link{np.fitgeo}}, \code{\link{locpol}}, \code{\link{fitsvar.sb.iso}}.
#' @export
#--------------------------------------------------------------------
np.geo <- function(lp, svm, svm0 = NULL, nbin = lp$grid$n) {
  stopifnot(inherits(lp, "locpol.bin" ))
  if (lp.hd <- any(nbin != lp$grid$n))
    lp <- with(lp, locpol(data$x, data$y, nbin = nbin, h = locpol$h))
  result <- lp
  stopifnot(inherits(svm, "svarmod"))
  result$svm <- svm
  if(!is.null(svm0)) {
    stopifnot(inherits(svm0, "svarmod" ))
    result$svm0 <- svm0
  }
  result$residuals <- if(inherits(svm, "fitsvar") & !lp.hd) 
    svm$esv$data$y else residuals(lp)
  oldClass(result) <- c("np.geo", oldClass(result))
  if(inherits(svm, "fitsvar"))
    oldClass(result) <- c("fitgeo", oldClass(result))
  else
    warning("'svm' is not a fitted variogram model")
  return(result)
}


#--------------------------------------------------------------------
#' @rdname npsp-internals
#' @method residuals np.geo
#' @keywords internal
#' @export
residuals.np.geo <- function(object, ...) {
  #--------------------------------------------------------------------
  if (!inherits(object, "np.geo"))
    stop("function only works for objects of class (or extending) 'np.geo'")
  return(object$residuals)
}



#--------------------------------------------------------------------
# np.fitgeo(x, ...)
#--------------------------------------------------------------------
#' Fit a nonparametric geostatistical model
#'
#' Fits a nonparametric (isotropic) geostatistical model 
#' (joint estimation of the trend and the variogram) by calling 
#' \code{\link{locpol}},  \code{\link{np.svariso.corr}} and 
#' \code{\link{fitsvar.sb.iso}} iteratively. 
#' At each iteration, the trend estimation bandwith is updated 
#' by a call to \code{\link{h.cv}}.
#' @param x object used to select a method (binned data, binned density 
#' or binned semivariogram).
#' @param ... further arguments passed to \code{\link{h.cv}}
#' (trend bandwith selecction parameters).
#' @details Currently, only isotropic semivariogram estimation is supported.
#' @return Returns an object of \code{\link{class}} \code{fitgeo} 
#' (extends \code{\link{np.geo}}). A \code{\link{locpol.bin}} object
#' with the additional (some optional) 3 components:
#' \item{svm}{fitted variogram model (object of class 
#'  \code{\link{fitsvar}}).}
#' \item{svm0}{(if requested) fitted residual variogram model (object of class 
#'  \code{\link{fitsvar}}).}
#' \item{residuals}{model residuals.} 
#' @seealso \code{\link{locpol}}, \code{\link{fitsvar.sb.iso}}, \code{\link{np.svar}}, 
#' \code{\link{np.svariso.corr}}, \code{\link{np.geo}}.
#' @export
#--------------------------------------------------------------------
np.fitgeo <- function(x, ...) UseMethod("np.fitgeo")



#' @rdname np.fitgeo  
#' @method np.fitgeo default
#' @inheritParams locpol.default
#' @inheritParams np.svariso.corr
#' @param iter maximum number of interations (of the whole algorithm).
#  If \code{iter == 0} the (uncorrected) residual variogram is computed
#' (traditional approach in geostatistics).
#' @param  h initial bandwidth matrix for trend estimation
#' (final bandwith if \code{iter <= 1}).
#' @param tol relative convergence tolerance (semivariogram).
#' @param h.svar bandwidth matrix for variogram estimation.
#' @param dk dimension of the Shapiro-Botha variogram model (see \code{\link{fitsvar.sb.iso}}).
#' @param svm.resid logical; if \code{TRUE}, the fitted (uncorrected) residual semivariogram model
#' is computed and returned.
#' @param warn logical; sets the handling of warning messages in bandwith selection (\code{\link{h.cv}}).
#' @param plot logical; if \code{TRUE}, semivariogram estimates obtained at each iteration are plotted.
#' @details If parameter \code{h} is not specified,
#' \code{\link{h.cv}} is called to set it.
#' If parameter \code{h.svar} is not specified,
#' is set to \code{1.5*h.cv.svar.bin()$h}.
#' @export
#--------------------------------------------------------------------
np.fitgeo.default <- function(x, y, nbin = NULL, iter = 2, h = NULL, tol = 0.05, set.NA = FALSE, 
                              h.svar = NULL, maxlag = NULL, nlags = NULL, dk = 0, svm.resid = FALSE,  
                              warn = FALSE, plot = FALSE, ...) {
#--------------------------------------------------------------------
  stopifnot(!missing(x), !missing(y)) # Solo para datos geoestadisticos
  # Binning
  bin <- binning(x, y, nbin = nbin, set.NA = set.NA)
  # Trend estimation
  if(is.null(h)) 
    h <- h.cv(bin, warn = warn, ...)$h
  lp <- locpol(bin, h = h, hat.bin = TRUE)   # np.svariso.corr
  # Variogram estimation
  lp.resid <- residuals(lp)
  if(is.null(h.svar) || svm.resid || iter == 0) {
    esvar0 <- svariso(x, lp.resid, maxlag = maxlag, nlags = nlags)
    if(is.null(h.svar)) h.svar <- 1.5 * h.cv.svar.bin(esvar0, warn = warn)$h
    if(svm.resid || iter == 0) {
      esvar0 <- np.svar(esvar0, h = h.svar)
      svm0 <- fitsvar.sb.iso(esvar0, dk = dk)
    }  
  }  
  if (iter > 0) {
    esvar <- np.svariso.corr(lp, h = h.svar, maxlag = maxlag, nlags = nlags, plot = plot)
    svm <- fitsvar.sb.iso(esvar, dk = dk)
  } else svm <- svm0
  # Result
  svm$iter <- iter
  result <- if(svm.resid && iter > 0) 
    np.geo(lp, svm, svm0) else np.geo(lp, svm)
  if(iter > 1) {
    result <- np.fitgeo.fitgeo(result, iter = iter - 1, tol = tol, ... )
    result$svm$iter <- result$svm$iter + 1
  }  
  return(result)
}


#' @rdname np.fitgeo  
#' @method np.fitgeo locpol.bin
#' @param  svm (fitted) variogram model (object of class 
#'  \code{\link{fitsvar}} or \code{\link{svarmod}}).
#' @export
#--------------------------------------------------------------------
np.fitgeo.locpol.bin <- function(x, svm, iter = 1, tol = 0.05, h.svar = svm$esv$locpol$h,   
                                 dk = 0, warn = FALSE, plot = FALSE, ...) {
#--------------------------------------------------------------------
  if(!inherits(svm, "fitsvar")) stop("'svm' is not a fitted variogram model")
  return(np.fitgeo.fitgeo(np.geo(x, svm), iter = iter, tol = tol, h.svar = h.svar,   
                          dk = dk, warn = warn, plot = plot, ...))
}

#' @rdname np.fitgeo  
#' @method np.fitgeo fitgeo
#' @export
#--------------------------------------------------------------------
np.fitgeo.fitgeo <- function(x, iter = 1, tol = 0.05, h.svar = x$svm$esv$locpol$h,   
                             dk = x$svm$par$dk, warn = FALSE, plot = FALSE, ...) {
#--------------------------------------------------------------------
  lp <- x
  svm <- x$svm
  maxlag <- svm$esv$grid$max
  nlags <- svm$esv$grid$n
  for (i in 1:iter) {
    last.fit <- svm$fit$fitted.sv
    # Trend estimation
    h <- if(hasArg("h.start")) h.cv(lp, cov.bin = svm, warn = warn, ...)$h else 
      h.cv(lp, cov.bin = svm, h.start = diag(lp$locpol$h), warn = warn, ...)$h
    lp <- locpol(lp, h = h, hat.bin = TRUE)   # np.svariso.corr
    # Variogram estimation
    lp.resid <- residuals(lp)
    if(is.null(h.svar)) {
      esvar0 <- svariso(lp$data$x, lp.resid, maxlag = maxlag, nlags = nlags)
      h.svar <- 1.5 * h.cv.svar.bin(esvar0, h.start = diag(svm$esv$locpol$h), warn = warn)$h
    }  
    esvar <- np.svariso.corr(lp, h = h.svar, maxlag = maxlag, nlags = nlags, plot = plot)
    svm <- fitsvar.sb.iso(esvar, dk = dk)
    error <- sqrt(mean((last.fit/svm$fit$fitted.sv - 1)^2, na.rm = TRUE))
    if (error < tol) break
  }
  svm$last.fit <- last.fit
  svm$last.fit.error <- error
  svm$iter <- i
  # Result
  result <- np.geo(lp, svm, x$svm0)
  return(result)
}



#' Plot a nonparametric geostatistical model
#' 
#' Plots the trend estimates and the fitted variogram model.
#' @method plot fitgeo
#' @param x a nonparametric geostatistical model object. 
#' Typically an output of \code{\link{np.fitgeo}}.
#' @param y	ignored argument.
#' @param main.trend	title for the trend plot.
#' @param main.svar	title for the semivariogram plot.
#' @param ... additional graphical parameters (currently ignored).
#' @export
plot.fitgeo <- function(x, y = NULL, main.trend = 'Trend estimates', 
                        main.svar = NULL, ...) {
  old.par <- par(mfrow = c(1,2)) #, omd = c(0.01, 0.9, 0.05, 0.95))
  simage(x, main = main.trend)
  main.svar <- if (is.null(main.svar)) 
    if (x$svm$iter == 0) "Nonparametric (uncorrected) semivariogram" else 
      "Nonparametric (bias-corrected) semivariogram"
  legend <- c("estimates", "fitted model")
  lty <- c(NA, 1)
  pch <- c(1, NA)
  lwd <- c(1, 2)
  plot(x$svm, main = main.svar, lwd = lwd, legend = FALSE, col = 'darkgray') 
  if(!is.null(x$svm0$fit)) {
    with(x$svm0$fit, lines(u, fitted.sv, lty = 3, lwd = 2))
    legend <- c(legend, 'uncorrected')
    lty <- c(lty, 3)
    pch <- c(pch, NA)
    lwd <- c(lwd, 2)
  }  
  if(!is.null(x$svm$last.fit)) {
    lines(x$svm$fit$u, x$svm$last.fit, lty = 2, lwd = 2)
    legend <- c(legend, 'previous iter')
    lty <- c(lty, 2)
    pch <- c(pch, NA)
    lwd <- c(lwd, 2)
  }  
  legend("bottomright", legend = legend, lty = lty, pch = pch, lwd = lwd)
  par(old.par)
}
