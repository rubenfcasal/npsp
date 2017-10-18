#--------------------------------------------------------------------
#   svarmod.sb.R (npsp package)
#--------------------------------------------------------------------
#   kappasb(x, dk)
#   disc.sb(nx, dk, rmax)
#   fitsvar.sb.iso(esv, dk, nx, rmax, min.contrib, method, iter, tol)
#
#   (c) R. Fernandez-Casal
#   Created: Apr 2013                          Last changed: May 2017
#--------------------------------------------------------------------
# PENDENTE:
#   - documentacion
#   - @examples
#--------------------------------------------------------------------



#--------------------------------------------------------------------
#   kappasb(x, dk = 0)
#--------------------------------------------------------------------
#' Coefficients of an extended Shapiro-Botha variogram model 
#' 
#' Computes the coefficients of an extended Shapiro-Botha variogram model. 
#' @details If \code{dk >= 1}, the coefficients are computed as: 
#' \deqn{\kappa_d(x) = (2/x)^{(d-2)/2} \Gamma(d/2) J_{(d-2)/2}(x)} 
#' where \eqn{J_p} is the Bessel function of order \eqn{p}. 
#' \cr If \code{dk == 0}, the coefficients are computed as: 
#' \deqn{\kappa _\infty(x) = e^{-x^2}}
#' (corresponding to a model valid in any spatial dimension). 
#' \cr NOTE: some authors denote these functions as \eqn{\Omega_d}.
#' 
#' @param  x   numeric vector (on which the kappa function will be evaluated).
#' @param  dk  dimension of the kappa function.
#' @return
#' A vector with the coefficients of an extended Shapiro-Botha variogram model. 
#' @references
#' Shapiro, A. and Botha, J.D. (1991) Variogram fitting with a general class of 
#'   conditionally non-negative definite functions. \emph{Computational Statistics 
#'   and Data Analysis}, \bold{11}, 87-96. 
#' @seealso
#' \code{\link{svarmod.sb.iso}}, \code{\link{besselJ}}.
#' @examples
#' kappasb(seq(0, 6*pi, len = 10), 2)
#'   
#' curve(kappasb(x/5, 0), xlim = c(0, 6*pi), ylim = c(-1, 1), lty = 2)
#' for (i in 1:10) curve(kappasb(x, i), col = gray((i-1)/10), add = TRUE)
#' abline(h = 0, lty = 3)
#' @export
kappasb <- function(x, dk = 0) {
#--------------------------------------------------------------------
  if ( zeros <- any(index <- x < sqrt(.Machine$double.eps)) ) {
      dx <- dim(x)
      x <- x[!index]
      # Alternativamente poderiase facer pmax(x, .Machine$double.eps^0.5)
  }    
  res <- switch( min(dk + 1, 5),
      exp(-x*x),            # dk = 0
      cos(x),               # dk = 1
      besselJ(x, nu = 0),   # dk = 2
      sin(x)/x,             # dk = 3
      gamma(dk/2) * (2/x)^((dk-2)/2) * besselJ(x, nu = (dk-2)/2) # dk >=4    
    ) # switch
  if (zeros) { 
      res[!index] <- res
      res[index] <- 1
      dim(res) <- dx
  }     
  return(res)
} 



#--------------------------------------------------------------------
#   disc.sb(nx, dk = 0, rmax = 1)
#--------------------------------------------------------------------
#' Discretization nodes of a Shapiro-Botha variogram model
#' 
#' Computes the discretization nodes of a `nonparametric' extended Shapiro-Botha 
#' variogram model, following Gorsich and Genton (2004), as the scaled roots of 
#' Bessel functions.
#' @param  nx  number of discretization nodes.
#' @param  dk  dimension of the kappa function (\code{dk >= 1}, see Details below).
#' @param  rmax  maximum lag considered.
#' @details
#' If \code{dk >= 1}, the nodes are computed as: 
#' \deqn{x_i = q_i/rmax; i = 1,\ldots, nx,} where 
#' \eqn{q_i} are the first \eqn{n} roots of \eqn{J_{(d-2)/2}}, \eqn{J_p} 
#' is the Bessel function of order \eqn{p} and \eqn{rmax} 
#' is the maximum lag considered. The computation of the zeros of the Bessel  
#' function is done using the efficient algorithm developed by Ball (2000).
#' 
#' If \code{dk == 0} (corresponding to a model valid in any spatial dimension), 
#' the nodes are computed so the gaussian variogram models involved have
#' practical ranges: 
#    \deqn{r_i = ( 1 + 1.2(i-1))rmax/nx; i = 1,\ldots, nx.}
#'    \deqn{r_i = 2 ( 1 + (i-1) ) rmax/nx; i = 1,\ldots, nx.}
#' @references
#' Ball, J.S. (2000) Automatic computation of zeros of Bessel functions and other
#'   special functions. \emph{SIAM Journal on Scientific Computing}, \bold{21}, 
#'   1458-1464.  
#'
#' Gorsich, D.J. and Genton, M.G. (2004) On the discretization of nonparametric 
#'   covariogram estimators. \emph{Statistics and Computing}, \bold{14}, 99-108.
#' @seealso
#' \code{\link{kappasb}}, \code{\link{fitsvar.sb.iso}}.
#' @examples
#' disc.sb( 12, 1, 1.0)
#' 
#' nx <- 1
#' dk <- 0
#' x <- disc.sb(nx, dk, 1.0)
#' h <- seq(0, 1, length = 100)
#' plot(h, kappasb(x * h, 0), type="l", ylim = c(0, 1))
#' abline(h = 0.05, lty = 2)
#' @export
disc.sb <- function(nx, dk = 0, rmax = 1) {   
#--------------------------------------------------------------------
    if (dk == 0) 
        # OLLO: equiv. modelos gausianos, poden aparecer inestabilidades
        # Nodos de discretizacion "xeometricos"
        # return( 1.732/(seq(1, by = -1/nx, length = nx) * rmax) )
        # return( sqrt(3)/(seq(1/nx, by = 1.2/nx, length = nx) * rmax) )
        # return( sqrt(3)/(seq(1/nx, 1, length = nx) * rmax) )
        return( sqrt(3)/(seq(1/nx, 1, length = nx) * 2 * rmax) )
    # dk > 1
    # Let's go FORTRAN!
    #   subroutine disc_sbv(nx, x, dim, range)
    ret <-.Fortran( "disc_sbv", nx = as.integer(nx), x = double(nx), 
                  dk = as.integer(dk), as.double(rmax))
    return(ret$x)
}                  



#--------------------------------------------------------------------
#   fitsvar.sb.iso(esv, dk = ncol(esv$data$x), nx = NULL, rmax = esv$grid$max, 
#        min.contrib = 10, method = c("cressie", "equal", "npairs", "linear"), 
#        iter = 10, tol = sqrt(.Machine$double.eps)) {
#--------------------------------------------------------------------
#' Fit an isotropic Shapiro-Botha variogram model
#' 
#' Fits a `nonparametric' isotropic Shapiro-Botha variogram model by WLS through
#' quadratic programming.
#' Following Gorsich and Genton (2004), the nodes are selected as the scaled 
#' roots of Bessel functions (see \code{\link{disc.sb}}).
#' @aliases fitsvar
#' @aliases fitsvar-class
#' @param  esv pilot semivariogram estimate, a \code{\link{np.svar}}-\code{\link{class}} 
#'   (or \code{\link{svar.bin}}) object. Typically an output of the function
#'   \code{\link{np.svariso}}. 
#' @param  dk  dimension of the kappa function (\code{dk == 0} corresponds to a model 
#'   valid in any dimension; if \code{dk > 0}, it should be greater than 
#'   or equal to the dimension of the spatial process \code{ncol(esv$data$x)}).
#' @param  nx  number of discretization nodes. Defaults to \code{min(nesv - 1, 50)},
#' where \code{nesv} is the number of semivariogram estimates.
#' @param  rmax  maximum lag considered in the discretization
#'   (range of the fitted variogram on output).
#' @param  min.contrib  minimum number of (equivalent) contributing pairs 
#' (pilot estimates with a lower number are ignored, with a warning).
#' @param  method  string indicating the WLS fitting method to be used
#'   (e.g. \code{method = "cressie"}). See "Details" below.
#' @param  iter  maximum number of interations of the WLS algorithm (used only  
#'   if \code{method == "cressie"}).
#' @param  tol  absolute convergence tolerance (used only  
#'   if \code{method == "cressie"}).
#' @details
#' The fit is done using a (possibly iterated) weighted least squares criterion, minimizing: 
#' \deqn{WLS(\theta) = \sum_i w_i[(\hat{\gamma}(h_i)) -	\gamma(\theta; h_i)]^2.}
#' The different options for the argument \code{method} define the WLS algorithm used:
#' \describe{
#'  \item{\code{"cressie"}}{The default method. The procedure is
#'  iterative, with \eqn{w_i = 1} (OLS) used for the first step
#'  and with the weights recalculated at each iteration,
#'  following Cressie (1985), until convergence: \deqn{w_i =
#'  N(h_i)/\gamma(\hat{\theta}; h_i)^2,} where \eqn{N(h_i)}
#'  is the (equivalent) number of contributing pairs in the
#'  estimation at lag \eqn{h_i}.}
#'  \item{\code{"equal"}}{Ordinary least squares: \eqn{w_i = 1}.}
#'  \item{\code{"npairs"}}{\eqn{w_i = N(h_i).}} 
#'  \item{\code{"linear"}}{\eqn{w_i = N(h_i)/h_i^2} 
#'  (default fitting method in \pkg{gstat} package).} 
#' } 
#' Function \code{\link[quadprog]{solve.QP}} of \pkg{quadprog} package is used
#' to solve a strictly convex quadratic program. To avoid problems, the Choleski decomposition
#' of the matrix corresponding to the original problem is computed using \code{\link{chol}} with \code{pivot = TRUE}.
#' If this matrix is only positive semi-definite (non-strictly convex QP),
#' the number of discretization nodes will be less than \code{nx}.
#' @return
#' Returns the fitted variogram model, an object of \code{\link{class}} \code{fitsvar}.
#' A \code{\link{svarmod}} object
# (extending \code{\link{sb.iso}}: a \code{\link{svarmod}}) object
#' with additional components \code{esv} (pilot semivariogram estimate) and \code{fit} containing:
#' \item{u}{vector of lags/distances.}
#' \item{sv}{vector of pilot semivariogram estimates.}
#' \item{fitted.sv}{vector of fitted semivariances.}
#' \item{w}{vector of (least squares) weights.}
#' \item{wls}{value of the objective function.}
#' \item{method}{string indicating the WLS fitting method used.}
#' \item{iter}{number of WLS iterations (if \code{method == "cressie"}).}   
#' 
#' @references
#' Ball, J.S. (2000) Automatic computation of zeros of Bessel functions and other
#'   special functions. \emph{SIAM Journal on Scientific Computing}, \bold{21}, 
#'   1458-1464.
#' 
#' Cressie, N. (1985) Fitting variogram models by weighted least squares.
#'   \emph{Mathematical Geology}, \bold{17}, 563-586. 
#' 
#' Cressie, N. (1993) \emph{Statistics for Spatial Data}. New York. Wiley.
#' 
#' Fernandez Casal R., Gonzalez Manteiga W. and  Febrero Bande M. (2003) 
#' Flexible Spatio-Temporal Stationary Variogram Models, 
#' \emph{Statistics and Computing}, \bold{13}, 127-136.
#'
#' Gorsich, D.J. and Genton, M.G. (2004) On the discretization of nonparametric 
#'   covariogram estimators. \emph{Statistics and Computing}, \bold{14}, 99-108.
#' 
#' Shapiro, A. and Botha, J.D. (1991) Variogram fitting with a general class of 
#'   conditionally non-negative definite functions. \emph{Computational Statistics 
#'   and Data Analysis}, \bold{11}, 87-96. 
#' @seealso
#' \code{\link{svarmod.sb.iso}}, \code{\link{disc.sb}}, \code{\link{plot.fitsvar}}.
#' @examples
#' # Trend estimation
#' lp <- locpol(aquifer[,1:2], aquifer$head, nbin = c(41,41),
#'              h = diag(100, 2), hat.bin = TRUE)
#'                                # 'np.svariso.corr()' requires a 'lp$locpol$hat' component
#'
#' # Variogram estimation
#' esvar <- np.svariso.corr(lp, maxlag = 150, nlags = 60, h = 60, plot = FALSE)
#'
#' # Variogram fitting
#' svm2 <- fitsvar.sb.iso(esvar)  # dk = 2
#' svm3 <- fitsvar.sb.iso(esvar, dk = 0) # To avoid negative covariances...
#' svm4 <- fitsvar.sb.iso(esvar, dk = 10) # To improve fit...
#'
#' plot(svm4, main = "Nonparametric bias-corrected semivariogram and fitted models", legend = FALSE)
#' plot(svm3, add = TRUE)
#' plot(svm2, add = TRUE, lty = 3)
#' legend("bottomright", legend = c("NP estimates", "fitted model (dk = 10)", "dk = 0", "dk = 2"),
#'             lty = c(NA, 1, 1, 3), pch = c(1, NA, NA, NA), lwd = c(1, 2, 1, 1))
#' @export
#--------------------------------------------------------------------
fitsvar.sb.iso <- function(esv, dk = 4*ncol(esv$data$x), nx = NULL, rmax = esv$grid$max, 
                            min.contrib = 10, method = c("cressie", "equal", "npairs", "linear"), 
                            iter = 10, tol = sqrt(.Machine$double.eps)) {
  #   PENDENTE:
  #     - Rounding errors?  w <- w/sum(w)
  #     - rematar documentacion: details, examples, ...
  #     - Version preliminar, final 'fit.svar.sb' valida para modelos anisotropicos
  #     - verificar missing values
  #     - nodes un vector de ptos de discretizacion
  #--------------------------------------------------------------------
  if (!inherits(esv, "svar.bin"))
    stop("function only works for objects of class (or extending) 'svar.bin'.")
  # if (esv$svar$type != "isotropic")
  if (esv$grid$nd != 1)
    stop("pilot variogram estimates 'esv' must be isotropic.")
  method <- match.arg(method)
  if (method != "cressie") iter <- 1
  # if (!requireNamespace(quadprog)) stop("'quadprog' package is required.")
  # Let's go...
  u <- as.numeric(coords(esv))
  if (inherits(esv, "np.svar")) {
    # "np.svar" class
    v <- esv$est
    if (!is.null(esv$locpol$hat)) {
      # Aproximacion varianza estilo Cressie (suponiendo independencia) para ajuste wls
      # PENDIENTE: ESCRIBIR/REVISAR ESTAS CUENTAS
      n <- 1 / with(esv, rowSums(locpol$hat^2 /
                                   pmax(matrix(binw, nrow=grid$n, ncol=grid$n, byrow=TRUE), 1))) # num equivalente de aportacions
    } else {
      n <- esv$binw # num de aportacions
    }
  } else {
    # "svar.bin" class
    v <- esv$biny
    n <- esv$binw # num de aportacions
  }
  if (!all(index <- n >= min.contrib)) {
    warning("some pilot semivariogram estimates will be ignored (contrib < min.contrib)")
    u <- u[index]
    v <- v[index]
    n <- n[index]
  }
  n.esv <- length(u)
  # Discretization points
  if (is.null(nx)) {
    nx <- min(n.esv - 1, 50)
  } else {
    if (nx >= length(u)) {
      warning("'nx' must be less than the number of variogram estimates (corrected)")
      nx <- n.esv - 1
    }
  }
  x <- disc.sb( nx, dk, rmax)
  # x <- sqrt(3)/u[1:nx]
  # x <- sqrt(3)/u[2:(nx + 1)]
  # M(1:nesv,1:npar)
  # M <- cbind(-outer(u, x, function(u, x) kappasb(u*x, dk)), 1)
  M <- cbind( -kappasb(outer(u, x), dk), 1)
  # Reescalar pesos
  w <- switch(method,
              npairs =  n,
              linear =   n/pmax(u^2, sqrt(.Machine$double.eps)),
              1 # default: "cressie", "equal"
  )
  w <- w/sum(w)
  # (iterative) quadratic programming
  n.par <- nx + 1
  Amat <- diag(n.par)
  Amat[1:nx, n.par] <- -1
  # wls loop
  i <- 0
  conver <- FALSE
  while( (i < iter) & !conver) {
    i <- i + 1
    # D = t(M) %*% diag(w) %*% M
    Dmat <- suppressWarnings(chol(crossprod(sqrt(w) * M), pivot = TRUE))
    n.par2 <- attr(Dmat, "rank")
    pivot <- attr(Dmat, "pivot")[1:n.par2]
    inu <- match(n.par, pivot, nomatch = 0)
    if (!inu) stop('the algorithm has failed to achieve an optimal solution.')
    # max(abs(crossprod(Dmat[1:n.par2, 1:n.par2])  - crossprod(sqrt(w) * M[, pivot])))
    # Invert the upper triangular matrix (is there a more efficient way?)
    Dmat <- backsolve(Dmat[1:n.par2, 1:n.par2], diag(n.par2))
    # d = t(M) %*% diag(w) %*% v
    dvec <- drop((w * v) %*% M[, pivot]) # crossprod?
    # solve min(-d^T b + 1/2 b^T D b) with the constraints A^T b >= b_0
    res <- quadprog::solve.QP(Dmat, dvec, Amat[pivot, pivot], factorized = TRUE) # bvec <- rep(0, n.par2)
    fit <- drop(M[, pivot] %*% res$solution)
    if (i > 1) {
      # fitted values absolute difference convergence criteria
      conver <-  max(abs(fit - last.fit)) < tol
      w <- n/pmax(fit^2, sqrt(.Machine$double.eps))
      w <- w/sum(w)
    }
    wls <- sum(w*(v-fit)^2)
    last.fit <- fit
  }
  if (method == "cressie") {
    iter <- i
    if (!conver) warning("the wls algorithm did not converge. \n",
                         "   Maximum number of iterations exceeded; \n",
                         "   the current values may be an approximate solution, \n",
                         "   or 'tol' is too big or 'iter' is too small.")
  }
  # Return the fitted Shapiro-Botha model
  nu <- res$solution[inu]
  z <- res$solution[-inu]
  x <- x[pivot[-inu]]
  # Rounding errors in solve.QP, the constraints A^T b >= b_0 might not (exactly) hold
  if (any(pivot <- z < 100*.Machine$double.eps)) {
    if (!all(pivot)) {
      z <- z[!pivot]
      x <- x[!pivot]
    } else {
      z <- 0
      x <- x[1]
    }
  }
  nu <- max(nu, sum(z))
  result <- svarmod.sb.iso( dk = dk, x = x, z = z, nu = nu, range = rmax)
  # Add fitting details
  result$fit <- list(u = u, sv = v, fitted.sv = fit, w = w, wls = wls,
                     method = method, iter = iter)
  result$esv <- esv
  oldClass(result) <- c("fitsvar", oldClass(result))
  return(result)
  #--------------------------------------------------------------------
} # fitsvar.sb.iso




