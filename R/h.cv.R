#--------------------------------------------------------------------
#   h.cv.R (npsp package)
#--------------------------------------------------------------------
#   h.cv()    S3 generic
#       h.cv.bin.data(bin, objective, h.start, h.lower, h.upper, degree, 
#               ncv, cov.bin, DEalgorithm, warn, tol.mask, ...)
#                   .compute.masked(bin, cov.bin, tol.mask)
#       h.cv.bin.den(bin, h.start, h.lower, h.upper, degree, 
#               ncv, cov.bin, DEalgorithm, warn, ...)
#       h.cv.svar.bin(bin, loss, h.start, h.lower, h.upper,             
#               degree, ncv, DEalgorithm, warn, ...) 
#                   .wloss(est, teor, w, loss) 
#   hcv.data(bin, objective, h.start, h.lower, h.upper, degree, 
#               ncv, cov.dat, DEalgorithm, ...)
#
# PENDENTE:
#   - Documentar diferencias entre h.cv.bin.data e hcv.data
#     hcv.data(x, y, ...)?
#   - Documentar metodos, Incluir so referencias?
#   - h.cv e locpolhcv PODEN TER PROBLEMAS CON DATOS MISSING
#     which(is.na(lp$est)) %in% which(is.na(bin$biny))
#   - opcion en binning() para obter cov binning a partir de cov datos?
#   - optimizar calculos matriciais "GCV" e "MASE"
#
#   (c) R. Fernandez-Casal         Last revision: Jan Sep 2013
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# h.cv(bin, ...)
#--------------------------------------------------------------------
#' Cross-validation methods for bandwidth selection
#'
#' Selects the bandwidth of a local polynomial kernel (regression, density or
#' variogram) estimator using (standart or modified) CV, GCV or MASE criteria.
#'
#' @param bin object used to select a method (binned data, binned density or binned semivariogram).
#' @param ... further arguments passed to or from other methods
#'    (e.g. parameters of the optimization routine).
#' @details Currently, only diagonal bandwidths are supported.
#' @return Returns a list containing the following 3 components:
#' \item{h}{the best (diagonal) bandwidth matrix found.} 
#' \item{value}{the value of the objective function corresponding to \code{h}.}
#' \item{objective}{the criterion used.}
#' @seealso \code{\link{locpol}}, \code{\link{locpolhcv}}, \code{\link{binning}}, 
#' \code{\link{np.svar}}.
#' @export
h.cv <- function(bin, ...) UseMethod("h.cv")
#--------------------------------------------------------------------


#--------------------------------------------------------------------
# h.cv.bin.data(bin, objective = c("CV", "GCV", "MASE"),
#               h.start = NULL, h.lower = NULL, h.upper = NULL, degree = 1, 
#               ncv = ifelse(objective == "CV", 2, 0), cov.bin = NULL, 
#               DEalgorithm = FALSE, warn = FALSE, tol.mask = npsp.tolerance(2), ...)
#--------------------------------------------------------------------
#' @rdname h.cv  
#' @method h.cv bin.data
#' @inheritParams locpol.default
#' @param objective character; optimal criterion to be used ("CV", "GCV" or "MASE").
#' @param h.start vector; initial values for the parameters (diagonal elements) to be optimized over.  
#' If \code{DEalgorithm == FALSE} (otherwise not used), defaults to \code{(3 + ncv) * lag},
#' where \code{lag = bin$grid$lag}.  
#' @param h.lower vector; lower bounds on each parameter (diagonal elements) to be optimized.  
#' Defaults to \code{(1.5 + ncv) * bin$grid$lag}.
#' @param h.upper vector; upper bounds on each parameter (diagonal elements) to be optimized.  
#'  Defaults to \code{1.5 * dim(bin) * bin$grid$lag}.
#' @param DEalgorithm logical; if \code{TRUE}, the differential evolution optimization algorithm 
#' in package \pkg{DEoptim} is used.
#' @param  ncv integer; determines the number of cells leaved out in each dimension.
#' (0 to GCV considering all the data, \eqn{>0} to traditional or modified cross-validation).
#' See "Details" bellow.
#' @param cov.bin (optional) covariance matrix of the binned data or semivariogram model
#' (\code{\link{svarmod}}-class) of the (unbinned) data. Defaults to the identity matrix.
#' @param warn logical; sets the handling of warning messages
#' (normally due to the lack of data in some neighborhoods).
#' If \code{FALSE} all warnings are ignored.
#' @param tol.mask tolerance used in the aproximations. Defaults to \code{\link{npsp.tolerance}(2)}.
#' @details
#' \code{h.cv} methods use binning approximations to the objective function values 
#' (in almost all cases, an averaged squared error). 
#' If \code{ncv > 0}, estimates are computed by leaving out binning cells with indexes within 
#' the intervals \eqn{[x_i - ncv + 1, x_i + ncv - 1]}, at each dimension i, where \eqn{x} 
#' denotes the index of the estimation location. \eqn{ncv = 1} corresponds with
#' traditional cross-validation and \eqn{ncv > 1} with modified CV 
#' (it may be appropriate for dependent data; see e.g. Chu and Marron, 1991, for the one dimensional case). 
#' Setting \code{ncv >= 2} would be recommended for sparse data (as linear binning is used).
#' For standard GCV, set \code{ncv = 0} (the whole data would be used).
#' For theoretical MASE, set \code{bin = binning(x, y = trend.teor)}, \code{cov = cov.teor} and \code{ncv = 0}.
#' 
#' If \code{DEalgorithm == FALSE}, the \code{"L-BFGS-B"} method in \code{\link{optim}} is used.
#' 
#' @references
#' Chu, C.K. and Marron, J.S. (1991) Comparison of Two Bandwidth Selectors
#'   with Dependent Errors. \emph{The Annals of Statistics}, \bold{19}, 1906-1918.
#'
#' Francisco-Fernandez M. and Opsomer J.D. (2005) Smoothing parameter selection
#'  methods for nonparametric regression with spatially correlated errors. 
#'  \emph{Canadian Journal of Statistics}, \bold{33}, 539-558.
#' @examples 
#' bin <- binning(earthquakes[, c("lon", "lat")], earthquakes$mag)
#' hcv <- h.cv(bin, ncv = 2)
#' lp <- locpol(bin, h = hcv$h)
#' ## Alternatively:
#' ## lp <- locpolhcv(earthquakes[, c("lon", "lat")], earthquakes$mag, ncv = 2)
#' 
#' simage(lp, main = 'Smoothed magnitude')
#' contour(lp, add = TRUE)
#' with(earthquakes, points(lon, lat, pch = 20))
#' 
#' ## Density estimation
#' hden <- h.cv(as.bin.den(bin))
#' den <- np.den(bin, h = hden$h)
#' 
#' plot(den, main = 'Estimated log(density)')
#' @export
h.cv.bin.data <- function(bin, objective = c("CV", "GCV", "MASE"),
            h.start = NULL, h.lower = NULL, h.upper = NULL, degree = 1,
            ncv = ifelse(objective == "CV", 2, 0), cov.bin = NULL,
            DEalgorithm = FALSE, warn = TRUE, tol.mask = npsp.tolerance(2), ...) {
#--------------------------------------------------------------------
# PENDENTE: CALCULO APROXIMADO DE tr(S %*% corr.dat) PARA CGCV
# CRITERIO APROXIMADO A PARTIR DE BINNING
# OLLO: cov.bin = MATRIZ DE COVARIANZAS DATOS BINNING OU MODELO DE VARIOGRAMA
# ... parametros adicionales para la rutina de optimizacion
# Actualmente solo ventana diagonal
# h.start VECTOR con aproximacion inicial
#       por defecto h.start <- (3+ncv)*lag
# h.lower y h.upper VECTORES para rangos de busqueda
#       por defecto h.lower <- (1.5+ncv)*lag
#       por defecto h.upper <- 1.5*n*lag
# Por defecto covarianza = identidad
# Para ventana "teorica": y = trend.teor, cov.bin = cov.bin.teor, ncv = 0
# Emplear DEalgorithm para asegurar convergencia al optimo global
#--------------------------------------------------------------------
    if (!inherits(bin, "bin.data") | inherits(bin, "svar.bin"))
        stop("function only works for objects of class 'bin.data'")
    objective <- match.arg(objective)
    if  (!is.null(cov.bin) && (objective == "GCV"))
        if ( inherits(cov.bin, "svarmod")) {
            sill <- cov.bin$sill
            if (is.na(sill)) stop("semivariogram model 'cov.bin' is not bounded.")
        } else stop("'cov.bin' must be a semivariogram model when 'objective' == 'GCV'.")
    nd <- bin$grid$nd
    n <- bin$grid$n
    nt <- prod(n)
    lag <- bin$grid$lag
    if(is.null(h.lower)) h.lower <- (1.5+ncv)*lag
    if(is.null(h.upper)) h.upper <- 1.5*n*lag

    bin$locpol$hat <- NULL
    res <- .compute.masked(bin, cov.bin = cov.bin, tol.mask = tol.mask)
    mask <- res$mask
    bin$binw[!mask] <- -1   # masked nodes will be leaved out in local polynomial estimation
    w <- res$w
    sw <- res$sw
    cov.bin <- res$cov.bin
    rssc <- sum(bin$data$y^2) - sum(w * bin$biny[mask]^2) # RSS binning correction

    if(is.null(cov.bin))
        # cov.bin <- diag(nt)
        # Select objective function
        f.opt <- switch(objective,
            CV  = function(x) {
                    lp <- locpol(bin, h = diag(x, nrow = nd), degree = degree, ncv = ncv) # hat.bin = FALSE
                    rss <- lp$locpol$rss + rssc
                    return(rss/sw)
                } ,
            GCV = function(x) {
                    lp <- locpol(bin, h = diag(x, nrow = nd), degree = degree, hat.bin = TRUE, ncv = ncv)
                    rss <- lp$locpol$rss + rssc
                    hat <- lp$locpol$hat[mask, mask]
                    # tr(S) ~= sum(diag(hat))
                    return( sw * (rss / (sw - sum(diag(hat)))^2) )
                } ,
            MASE= function(x) {
                    lp <- locpol(bin, h = diag(x, nrow = nd), degree = degree, hat.bin = TRUE, ncv = ncv)
                    rss <- lp$locpol$rss + rssc
                    hat <- lp$locpol$hat[mask, mask]
                    # tr(S%*%t(S)) = sum(S^2) ~= sum(w * hat^2 * rep(1/w, each = nrow(hat)))
                    return( (rss + sum(w * hat^2 * rep(1/w, each = nrow(hat))))/sw )
                }
        ) # switch
    else {
        # correlation matrix for "GCV"
        # if  (objective == "GCV") cov.bin <- cov2cor(cov.bin)
        if  (objective == "GCV") cov.bin <- cov.bin/sill  # Solo valido para el caso homocedastico
        # Select objective function
        f.opt <- switch(objective,
            CV  = function(x) {
                    lp <- locpol(bin, h = diag(x, nrow = nd), degree = degree, hat.bin = TRUE, ncv = ncv)
                    rss <- lp$locpol$rss + rssc
                    hat <- lp$locpol$hat[mask, mask]
                    # tr(S %*% cov.dat) = sum(hat * cov.dat) ~= sum(hat * cov.bin * w)
                    return( (rss + 2 * sum(hat * cov.bin * w))/sw )  # trace(A%*%B) = sum(A*t(B))
                } ,
            GCV = function(x) {
                    lp <- locpol(bin, h = diag(x, nrow=nd), degree = degree, hat.bin = TRUE, ncv = ncv)
                    rss <- lp$locpol$rss + rssc
                    hat <- lp$locpol$hat[mask, mask]
                    return( sw * (rss / (sw - sum(hat * cov.bin * w))^2) )
                } ,
            MASE= function(x) {
                    lp <- locpol(bin, h = diag(x, nrow = nd), degree = degree, hat.bin = TRUE, ncv = ncv)
                    rss <- lp$locpol$rss + rssc
                    hat <- lp$locpol$hat[mask, mask]
                    # tr(S %*% cov.dat %*% t(S)) = sum(cov.dat * crossprod(S)) ~= sum(cov.bin * crossprod(sqrt(w) * hat))
                    return( (rss + sum(cov.bin * crossprod(sqrt(w) * hat)))/sw )
                }
        ) # switch
    } # if(is.null(cov.bin))
    # Minimization of the objective function
    if(as.logical(warn)) {
        lp.warn <- FALSE
        w.handler <- function(w){ # warning handler
            lp.warn <<- TRUE
            invokeRestart("muffleWarning")
        }
        f.opt.bak <- f.opt
        f.opt <- function(x) withCallingHandlers(f.opt.bak(x),  warning = w.handler)
    } else {
        ops <- options(warn = -1)
        on.exit(options(ops))
    }
    if(!DEalgorithm) {
        if(is.null(h.start)) h.start <- (3+ncv)*lag
        res <- stats::optim( h.start, f.opt, method = "L-BFGS-B",
                      lower = h.lower, upper = h.upper, ...)
        res <- list(h = diag(res$par, nrow = nd), value = res$value, objective = objective)
    } else {
        if (!requireNamespace("DEoptim")) stop("'DEalgorithm' requires 'DEoptim' package")
        res <- DEoptim::DEoptim( f.opt, lower = h.lower, upper = h.upper, ...)
        res <- list(h = diag(res$optim$bestmem, nrow = nd), value = res$optim$bestval, objective = objective)
    }
    if (warn && lp.warn) {
        warning("The objective function was evaluated at too-small bandwidths \n",
            "  (there were not enough data in the neighborhoods). \n",
            "  The current optimum bandwidth may be an approximate solution, \n",
            "  (try increasing h.lower = ", paste(round(h.lower,5), collapse =', '),").")
    }
    return(res)
#--------------------------------------------------------------------
} # h.cv.bin.data




#--------------------------------------------------------------------
#   .compute.masked(bin, cov.bin = NULL, tol.mask = npsp.tolerance(2))
#--------------------------------------------------------------------
#' @rdname npsp-internals
#' @param cov.bin (optional) covariance matrix of the binned data or semivariogram model
#' (\code{\link{svarmod}}-class) of the (unbinned) data.
#' @return \code{.compute.masked} returns a list with components:
#' \item{mask}{logical vector \code{bin$binw > tol.mask}.}
#' \item{w}{\code{x$binw[mask]}.}
#' \item{sw}{\code{sum(w)}.}
#' \item{hat}{(optional) \code{bin$locpol$hat[mask, mask]}.}
#' \item{cov.bin}{(optional) masked (aproximated) covariance matrix of the binned data.}
#' @keywords internal
.compute.masked <- function(bin, cov.bin = NULL, tol.mask = npsp.tolerance(2)){
#--------------------------------------------------------------------
# COIDADO: sw <- sum(w) # normalmente = length(x$data$y) salvo en semivariogramas
# NOTA: facer bin$binw[!mask] <- -1 despois de chamar  # masked nodes will be leaved out in local polynomial estimation
# PENDENTE: incluir parametros: coords = FALSE,  corr = FALSE?
    if (tol.mask <= 0) stop("'tol.mask' must be > 0")
    mask <- mask.bin.den(bin, mask.default(bin$binw, tol.mask), warn = FALSE)$mask # Coidado pode remaskear...
    w <- bin$binw[mask]
    sw <- sum(w)
    res <- list(mask = mask, w = w, sw = sw)
    # NOTA: Chequear inherits(lp, "locpol.bin")?
    if (!is.null(bin$locpol$hat)) res$hat <- bin$locpol$hat[mask, mask]
    # sum(bin$locpol$hat[,!mask]) ~= 0
    if(!is.null(cov.bin)) {
        p <- (d <- dim(cov.bin))[1L]
        if (!is.numeric(cov.bin) || length(d) != 2L || p != d[2L]){   # check if not square matrix...
            if (!inherits(cov.bin, "svarmod"))                        # semivariogram model
                stop("'cov.bin' must be a square matrix or a semivariogram model ('svarmod' class).")
            if (!inherits(cov.bin, "isotropic"))
                stop("function only works for isotropic semivariograms")
            # v.bin <- svar.teor(0.24*sqrt(sum(xlag^2)))
            # v.bin <- svar.teor(0.467*sqrt(sum(xlag^2)))             # linear binning 'lag mean'...
            v.bin <- sv(cov.bin, sqrt(sum(bin$grid$lag^2))/3)         # simple binning 'lag mean'...
            cov.bin <- varcov(cov.bin, coords = coords(bin)[mask,])
            diag(cov.bin) <- cov.bin[1, 1] + v.bin * (1/w - 1)
        } else {
           if (p == prod(dim(bin)))
                cov.bin <- cov.bin[mask, mask]
           else if (p != sum(mask))
                stop("'cov.bin' must be a square matrix of appropriate order.")
        }
        res$cov.bin <- cov.bin
    }
    return(res)
} #------------------------------------------------------------------






#--------------------------------------------------------------------
# h.cv.bin.den(bin, h.start = NULL, h.lower = NULL, h.upper = NULL,
#            degree = 1, ncv = 2, DEalgorithm = FALSE, warn = FALSE,...)
#--------------------------------------------------------------------
#' @rdname h.cv
#' @method h.cv bin.den
#' @export
h.cv.bin.den <- function(bin, h.start = NULL, h.lower = NULL, h.upper = NULL,
            degree = 1, ncv = 2, DEalgorithm = FALSE, ...) {
#--------------------------------------------------------------------
# NOTA: MASK NON POSIBLE ACTUALMENTE (np.den -> locpol.bin.den -> lp_data_grid)
    # if (!inherits(bin, "bin.den") || inherits(bin, "bin.data"))
    if (!inherits(bin, "bin.den"))
        stop("function only works for objects of class 'bin.den'")
    nd <- bin$grid$nd
    sw <- sum(bin$binw)     # normalmente = length(bin$data$y) salvo en semivariogramas
    f.opt <- function(x) np.den(bin, h = diag(x, nrow = nd), degree = degree, ncv = ncv)$locpol$rss/sw
    # Minimization of the objective function
    lag <- bin$grid$lag
    if(is.null(h.lower)) h.lower <- (1.5 + ncv) * lag else stopifnot(length(h.lower) == nd)
    if(is.null(h.upper)) h.upper <- 1.5 * bin$grid$n * lag else stopifnot(length(h.upper) == nd)
    if(!DEalgorithm) {
        if(is.null(h.start)) h.start <- (3+ncv)*lag else stopifnot(length(h.start) == nd)
        res <- optim( h.start, f.opt, method = "L-BFGS-B",
                      lower = h.lower, upper = h.upper, ...)
        return(list(h = diag(res$par, nrow = nd), value = res$value, objective = "CV"))
    } else {
        if (!requireNamespace("DEoptim")) stop("'DEalgorithm' requires 'DEoptim' package")
        res <- DEoptim::DEoptim( f.opt, lower = h.lower, upper = h.upper, ...)
        return(list(h = diag(res$optim$bestmem, nrow = nd), value = res$optim$bestval, objective = "CV"))
    }
}


#--------------------------------------------------------------------
# h.cv.svar.bin(bin, loss = c("ARSE", "ARAE", "ASE", "AAE"),
#             h.start = NULL, h.lower = NULL, h.upper = NULL,
#             degree = 1, ncv = 1, DEalgorithm = FALSE, warn = FALSE, ...) {
#--------------------------------------------------------------------
#' @rdname h.cv
#' @method h.cv svar.bin
#' @param loss character; CV error. See "Details" bellow.
#' @details
#' The different options for the argument \code{loss} in \code{h.cv.svar.bin()} define the CV error 
#' considered in semivariogram estimation:
#' \describe{
#'  \item{\code{"ASE"}}{Averaged squared error}
#'  \item{\code{"ARSE"}}{Averaged relative squared error}
#'  \item{\code{"AAE"}}{Averaged absolute error}
#'  \item{\code{"ARAE"}}{Averaged relative absolute error}
#' }
#' @export
h.cv.svar.bin <- function(bin, loss = c("ARSE", "ARAE", "ASE", "AAE"),
            h.start = NULL, h.lower = NULL, h.upper = NULL,
            degree = 1, ncv = 1, DEalgorithm = FALSE, warn = FALSE, ...) {
#--------------------------------------------------------------------
    if (!inherits(bin, "svar.bin"))
        stop("function only works for objects of class 'svar.bin'")
    nd <- bin$grid$nd
    loss <- match.arg(loss)
    # Se podria optimizar el calculo cuando loss = 'ASE'...
    f.opt <- function(x) {
        esvar <- locpol.svar.bin(bin, h = diag(x, nrow = nd), degree = degree, ncv = ncv)
        return(with(esvar, .wloss(biny, est, binw, loss)))
    }
    # Minimization of the objective function
    lag <- bin$grid$lag
    if(is.null(h.lower)) h.lower <- (1.5 + ncv) * lag else stopifnot(length(h.lower) == nd)
    if(is.null(h.upper)) h.upper <- 1.5 * bin$grid$n * lag else stopifnot(length(h.upper) == nd)
    if(!as.logical(warn)) {
        ops <- options(warn = -1)
        on.exit(options(ops))
    }
    if(!DEalgorithm) {
        if(is.null(h.start)) h.start <- (3+ncv)*lag else stopifnot(length(h.start) == nd)
        res <- optim( h.start, f.opt, method = "L-BFGS-B",
                      lower = h.lower, upper = h.upper, ...)
        return(list(h = diag(res$par, nrow = nd), value = res$value, objective = "CV", loss = loss))
    } else {
        if (!requireNamespace("DEoptim")) stop("'DEalgorithm' requires 'DEoptim' package")
        res <- DEoptim::DEoptim( f.opt, lower = h.lower, upper = h.upper, ...)
        return(list(h = diag(res$optim$bestmem, nrow = nd), value = res$optim$bestval, objective = "CV", loss = loss))
    }
}


#--------------------------------------------------------------------
#' @rdname npsp-internals
# @param est vector of estimates.
# @param teor vector of theoretical/reference values.
# @param w vector of weights (with the same length as \code{est}).
#' @keywords internal
.wloss <- function(est, teor, w, loss = c('ASE','ARSE','AAE','ARAE')) {
#--------------------------------------------------------------------
#' @rdname h.cv
    loss <- match.arg(loss)
    return(switch(loss,
        # Averaged squared error
        ASE = stats::weighted.mean((est - teor)^2, w, na.rm = TRUE),
        # Averaged relative squared error
        ARSE = stats::weighted.mean((est/pmax(teor, npsp.tolerance()) - 1)^2, w, na.rm = TRUE),
        # Averaged absolute error
        AAE = stats::weighted.mean(abs(est - teor), w, na.rm = TRUE),
        # Averaged relative absolute error
        ARAE = stats::weighted.mean(abs(est/pmax(teor, npsp.tolerance()) - 1), w, na.rm = TRUE)

    ))
}


#--------------------------------------------------------------------
# hcv.data(bin, objective = c("CV", "GCV", "MASE"),
#         h.start = NULL, h.lower = NULL, h.upper = NULL, degree = 1,
#         ncv = ifelse(objective == "GCV", 0, 1), cov = NULL,
#         DEalgorithm = FALSE, warn = FALSE, ...)
#--------------------------------------------------------------------
#' @rdname h.cv
#' @param cov.dat covariance matrix of the data or semivariogram model
#' (of class extending \code{\link{svarmod}}). Defaults to the identity matrix
#' (uncorrelated data).
#' @details
#' \code{hcv.data} evaluates the objective function at the original data
#' (combining a binning approximation to the nonparametric estimates with a linear interpolation),
#' this can be very slow (and memory demanding; consider using \code{h.cv} instead).
#' If \code{ncv > 1} (modified CV), a similar algorithm to that in \code{h.cv} is used,
#' estimates are computed by leaving out binning cells with indexes within
#' the intervals \eqn{[x_i - ncv + 1, x_i + ncv - 1]}.
#' @export
hcv.data <- function(bin, objective = c("CV", "GCV", "MASE"),
            h.start = NULL, h.lower = NULL, h.upper = NULL, degree = 1,
            ncv = ifelse(objective == "CV", 1, 0), cov.dat = NULL,
            DEalgorithm = FALSE, warn = TRUE, ...) {
# CRITERIO BINNING "EXACTO"
# cov.dat = MATRIZ DE COVARIANZAS DOS DATOS ORIXINAIS OU MODELO DE SEMIVARIOGRAMA
# COIDADO CO TEMPO DE CPU SE NUMERO DE DATOS GRANDE
#--------------------------------------------------------------------
    if (!inherits(bin, "bin.data"))
        stop("function only works for objects of class (or extending) 'bin.data'")
    if (inherits(bin, "svar.bin"))
        stop("not supported for objects of class (or extending) 'svar.bin'")
    nd <- bin$grid$nd
    n <- bin$grid$n
    ny <- length(bin$data$y)
    objective <- match.arg(objective)
    lag <- bin$grid$lag
    if(is.null(h.lower)) h.lower <- (1.5+ncv)*lag else stopifnot(length(h.lower) == nd)
    if(is.null(h.upper)) h.upper <- 1.5*n*lag else stopifnot(length(h.upper) == nd)
    bin <- mask.bin.data(bin, mask = mask.default(bin$binw, 0), set.NA = FALSE,
                            warn = FALSE, filter.lp = TRUE) # Coidado pode remaskear...
    dmax <- sqrt(.Machine$double.xmax)
    tol <- npsp.tolerance()
    if(is.null(cov.dat))
        # cov.dat <- diag(ny)
        f.opt <- switch(objective,
            CV  = if (ncv==1) function(x) {
                    lp <- locpol(bin, h = diag(x, nrow = nd), degree = degree, hat.bin = TRUE, ncv = 0)
                    if(!is.null(lp$locpol$nrl0)) return(dmax)
                    lpdat <- predict(lp, hat.data = TRUE)
                    res <- 1 - diag(lpdat$y.hat)
                    if (any(abs(res) < tol)) {
                        warning("Not enough data in neighborhood")
                        return(dmax)
                    }
                    return(mean(((lpdat$y.est - lp$data$y) / res)^2))
                } else function(x) {
                    lp <- locpol(bin, h = diag(x, nrow = nd), degree = degree, ncv = ncv)
                    if(!is.null(lp$locpol$nrl0)) return(dmax)
                    y.est <- predict(lp)  # hat.data = FALSE
                    return(mean((y.est - lp$data$y)^2))
                  },
            GCV = function(x) {
                    lp <- locpol(bin, h = diag(x, nrow=nd), degree = degree, hat.bin = TRUE, ncv = ncv)
                     if(!is.null(lp$locpol$nrl0)) return(dmax)
                    lpdat <- predict(lp, hat.data = TRUE)
                    return( mean((lpdat$y.est - lp$data$y)^2) /
                        (1 - mean(diag(lpdat$y.hat)))^2 )
                  },
            MASE= function(x) {
                    lp <- locpol(bin, h = diag(x, nrow = nd), degree = degree, hat.bin = TRUE, ncv = ncv)
                    if(!is.null(lp$locpol$nrl0)) return(dmax)
                    lpdat <- predict(lp, hat.data = TRUE)
                    # tr(S%*%t(S)) = sum(S^2)
                    return( (sum((lpdat$y.est - lp$data$y)^2) +  sum(lpdat$y.hat^2))/ny )
                  }
        ) # switch
    else {
        p <- (d <- dim(cov.dat))[1L]
        if (!is.numeric(cov.dat) || length(d) != 2L || p != d[2L] || p != ny){
            if (!inherits(cov.dat, "svarmod"))
                stop("'cov.dat' must be a square matrix or a semivariogram model ('svarmod' class).")
            if (!inherits(cov.dat, "isotropic"))
                stop("function only works for isotropic semivariograms")
            cov.dat <- varcov(cov.dat, coords = bin$data$x)
        }
        # if  (objective == "GCV") cov.dat <- cov2cor(cov.dat)  # correlation matrix for "GCV"
        if  (objective == "GCV") cov.dat <- cov.dat/cov.dat[1,1]  # Solo valido para el caso homocedastico
        # Select objective function
        f.opt <- switch(objective,
            CV  = if (ncv==1) function(x) {
                    lp <- locpol(bin, h = diag(x, nrow = nd), degree = degree, hat.bin = TRUE, ncv = 0)
                    if(!is.null(lp$locpol$nrl0)) return(dmax)
                    lpdat <- predict(lp, hat.data = TRUE)
                    res <- 1 - diag(lpdat$y.hat)
                    if (any(abs(res) < tol)) {
                        warning("Not enough data in neighborhood")
                        return(dmax)
                    }
                    # tr(S %*% cov.dat) = sum(hat * cov.dat)
                    return((sum( (( lpdat$y.est - lp$data$y) / res )^2)
                            + 2 * sum(lpdat$y.hat * cov.dat ))/ny)  # OJO: Matriz de suavizado con todos los datos
                } else function(x) {
                    lp <- locpol(bin, h = diag(x, nrow = nd), degree = degree, hat.bin = TRUE, ncv = ncv)
                    if(!is.null(lp$locpol$nrl0)) return(dmax)
                    lpdat <- predict(lp, hat.data = TRUE)
                    return( (sum((lpdat$y.est - lp$data$y)^2) + 2 * sum(lpdat$y.hat * cov.dat ))/ny )
                  },
            GCV = function(x) {
                    lp <- locpol(bin, h = diag(x, nrow=nd), degree = degree, hat.bin = TRUE, ncv = ncv)
                    if(!is.null(lp$locpol$nrl0)) return(dmax)
                    lpdat <- predict(lp, hat.data = TRUE)
                    return( ny * ( sum((lpdat$y.est - lp$data$y)^2) /
                            (ny - sum(lpdat$y.hat * cov.dat))^2) )
                  },
            MASE= function(x) {
                    lp <- locpol(bin, h = diag(x, nrow = nd), degree = degree, hat.bin = TRUE, ncv = ncv)
                    if(!is.null(lp$locpol$nrl0)) return(dmax)
                    lpdat <- predict(lp, hat.data = TRUE)
                    # tr(S %*% cov.dat %*% t(S)) = sum(cov.dat * crossprod(S))
                    return( (sum((lpdat$y.est - lp$data$y)^2) +
                        sum(cov.dat * crossprod(lpdat$y.hat)))/ny )
                  }
        ) # switch
    } # if(is.null(cov.dat))

    # Minimization of the objective function
    if(as.logical(warn)) {
        lp.warn <- FALSE
        w.handler <- function(w){ # warning handler
            lp.warn <<- TRUE
            invokeRestart("muffleWarning")
        }
        f.opt.bak <- f.opt
        f.opt <- function(x) withCallingHandlers(f.opt.bak(x),  warning = w.handler)
    } else {
        ops <- options(warn = -1)
        on.exit(options(ops))
    }
    # dmax <- 1.1 * f.opt(h.upper)
    if(!DEalgorithm) {
        if(is.null(h.start)) h.start <- (3+ncv)*lag
            else stopifnot(length(h.start) == nd) # Verificar dentro de limites? en optim?
        res <- stats::optim( h.start, f.opt, method = "L-BFGS-B",
                      lower = h.lower, upper = h.upper, ...)
        res <- list(h = diag(res$par, nrow = nd), value = res$value, objective = objective)
    } else {
        if (!requireNamespace("DEoptim")) stop("'DEalgorithm' requires 'DEoptim' package")
        res <- DEoptim::DEoptim( f.opt, lower = h.lower, upper = h.upper, ...)
        res <- list(h = diag(res$optim$bestmem, nrow = nd), value = res$optim$bestval, objective = objective)
    }
    if (warn && lp.warn) {
        warning("The objective function was evaluated at too-small bandwidths \n",
            "  (there were not enough data in the neighborhoods). \n",
            "  The current optimum bandwidth may be an approximate solution, \n",
            "  (try increasing h.lower = ", paste(round(h.lower,5), collapse =', '),").")
    }
    return(res)
#--------------------------------------------------------------------
} # hcv.data







