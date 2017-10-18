#--------------------------------------------------------------------
#   npsp-rules.R (library(npsp))
#--------------------------------------------------------------------
#   rule
#   rule.binning              S3 methods
#       rule.binning.default
#
#--------------------------------------------------------------------
#    npsp es un software libre y viene sin GARANTIA ALGUNA.
#    Usted puede redistribuirlo bajo ciertas circunstancias.
#    Escriba 'license()' para detalles de distribucion.
#
#   (c) R. Fernandez-Casal
#   Creation: Jul 2015                       Last revision: Aug 2015
#--------------------------------------------------------------------

#--------------------------------------------------------------------
#' npsp Rules
#'
#' Compute the number of classes for a histogram,
#' the number of nodes of a binning grid, etc.
#' @param  x  data vector or object used to select a method.
#' @param  d  (spatial) dimension.
#' @param  rule character; rule to be used.
#' @param  ... further arguments passed to or from other methods.
#' @details
#' The Rice Rule, \eqn{m = \lceil 2 n^{1/3} \rceil,}
#' is a simple alternative to Sturges's rule (\code{\link{nclass.Sturges}}).
#' @seealso
#' \code{\link{hist}}, \code{\link{nclass.Sturges}}, \code{\link{nclass.scott}},
#' \code{\link{nclass.FD}}, 
# \code{\link{np.svariso}}, \code{\link{svar.bin}},
#' \code{\link{binning}}, \code{\link{np.den}}, \code{\link{bin.den}}.
#' @return
#' The rule values (vector or scalar).
#' @importFrom grDevices nclass.FD nclass.Sturges nclass.scott
#' @export
#--------------------------------------------------------------------
rule <- function(x, d = 1, rule = c("Rice", "Sturges", "scott", "FD"), ...){
  rule <- match.arg(rule)
    res <- switch(rule,
        Rice  = .rice.rule(x, ...),
        Sturges = nclass.Sturges(x),
        scott = nclass.scott(x),
        FD = nclass.FD(x) )
    return(res ^(1/d))
}

#--------------------------------------------------------------------
#' @rdname npsp-internals
#' @param  a  scale values.
#' @param  b  exponent values.
#' @keywords internal
# @export
#--------------------------------------------------------------------
.rice.rule <- function(x, a = 2, b = 3, ...) ceiling(a * x ^ (1 / b))
# PENDENTE:
# ["https://en.wikipedia.org/wiki/Histogram#Mathematical_definition"]
# nb <- round(log(ny)/log(2) + 1) # ter en conta o num de dimensions?
#     <- round(log(xxx)/log(2) + 1) # <- round((log(xxx)/log(1.1) + 1)^(1/2))^2
# Regla de Sturges k=1+log2(n)
# Regla de Scott K=(2n)^(1/3)
# (4*n)^(2/5)

#--------------------------------------------------------------------
#' @rdname rule
# @param  x  object used to select a method.
# @param  ... further arguments passed to or from other methods.
#' @return
#' \code{rule.binning} returns a vector with the suggested number of bins
#' on each dimension.
#' @export
#--------------------------------------------------------------------
rule.binning <- function(x, ...) UseMethod("rule.binning")

#--------------------------------------------------------------------
#' @rdname rule
#' @method rule.binning default
#' @param  a  scale values.
#' @param  b  exponent values.
#' @export
#' @return
#' \code{rule.binning.default} returns \code{rep(ceiling(a * nrow(x) ^ (1 / b)), d)}.
#--------------------------------------------------------------------
rule.binning.default <- function(x, d = ncol(x), a = 2, b = d + 1, ...) {
    x <- as.matrix(x)
    ny <- nrow(x)
    if (!missing(d)) d <- as.integer(d)
    # return( rep( round(rule(ny, rule = "Rice", b = b)), d) )
    return( rep( .rice.rule(ny, a = a, b = b), d) )
}


#--------------------------------------------------------------------
#' @rdname rule
# @param  x  object used to select a method.
# @param  ... further arguments passed to or from other methods.
#' @return
#' \code{rule.svar} returns the suggested number of bins
#' for variogram estimation.
#' @export
#--------------------------------------------------------------------
rule.svar <- function(x, ...) UseMethod("rule.svar")

#--------------------------------------------------------------------
#' @rdname rule
#' @method rule.svar default
#' @return
#' \code{rule.svar.default} returns \code{ceiling(a * (nrow(x)^2 / 4) ^ (1 / b))}.
#' @export
#--------------------------------------------------------------------
rule.svar.default <- function(x, d = ncol(x), a = 2, b = d + 1, ...) {
    x <- as.matrix(x)
    ny <- nrow(x)
    if (!missing(d)) d <- as.integer(d)
    return( .rice.rule(ny^2/4, a = a, b = b))
}

#--------------------------------------------------------------------
#' @rdname rule
#' @method rule.svar bin.den
#' @export
#--------------------------------------------------------------------
rule.svar.bin.den <- function(x, ...) {
    # nlags <- round(2 * mean(dim(x)))
    # maxlag <- with(x$grid, 0.45*sqrt(sum(max - min)^2))     # 45% of largest lag
    # h <- 2*maxlag/nlags
    # return( list(nlags = nlags, maxlag = maxlag, h = h, dk = dk ))
    return(round(2 * mean(dim(x))))
}
