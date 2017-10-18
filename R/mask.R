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
#   Created: Jul 2015                          Last changed: Jan 2016
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
#' @method mask bin.den
#' @param mask logical; vector indicating the selected values (not masked).
#' @param  set.NA logical; If \code{TRUE}, the values corresponding
#' to masked cells are set to \code{NA}.
#' @param warn logical; If \code{TRUE} a warning message is generated when original data is masked.
#' @return \code{mask.bin.den}, \code{mask.bin.data} and \code{mask.locpol.bin}
#' return an object of the same class as \code{x} with the additional component \code{$mask}.
#' @examples
#' bin <- binning(aquifer[,1:2], aquifer$head, nbin = c(41,41), set.NA = TRUE)
#' str(mask(bin, mask(bin$binw), warn = TRUE))
#' str(mask(bin, mask(bin$binw, 1)))
#' @export
mask.bin.den <- function(x, mask = mask.default(x$binw, npsp.tolerance(2)),
                         warn = TRUE, ...) {
#--------------------------------------------------------------------
# PENDENTE: establecer binw[mask] <- -1 en estimacion locpol
# PENDENTE: Ollo con binw = -1
# COIDADO: mask <- as.numeric(x$binw) > 0    # hcv.data => predict
    if (!is.logical(mask) || length(mask) != length(x$binw))
           stop("'mask' must be a logical vector of appropriate order.")
    if (!is.null(x$mask)) {
        mask <- mask & x$mask
        if (warn) warning("'x' is already masked.")
    }
    if (warn && any(!mask & (x$binw > 0)))
        warning("some data will be masked...")
    # if (filter.lp) res$binw[!mask] <- 0        
    x$mask <- mask
    return(x)
} #------------------------------------------------------------------




#--------------------------------------------------------------------
#' @rdname mask
#' @method mask bin.data
#' @param  filter.lp logical; If \code{TRUE}, masked nodes will be leaved out
#' in local polynomial estimation.
#' @export
mask.bin.data <- function(x, mask = mask.default(x$binw, npsp.tolerance(2)),
                          set.NA = TRUE, warn = TRUE, filter.lp = set.NA, ...) {
#--------------------------------------------------------------------
    # if (filter.lp) warn <- TRUE
    res <- mask.bin.den(x, mask = mask, warn = warn)
    if (set.NA) is.na(res$biny) <- !res$mask
    if (filter.lp) res$binw[!mask] <- -1
    return(res)}
#--------------------------------------------------------------------



#--------------------------------------------------------------------
#' @rdname mask
#' @method mask locpol.bin
#' @export
mask.locpol.bin <- function(x, mask = mask.default(x$binw, npsp.tolerance(2)),
                          set.NA = TRUE, warn = TRUE, filter.lp = set.NA, ...) {
#--------------------------------------------------------------------
    res <- mask.bin.data(x, mask = mask, set.NA = set.NA, warn = warn, filter.lp = filter.lp)
    if (set.NA) is.na(res$est) <- !res$mask  # Non seria preciso se os datos binning son enmascarados...
    return(res)}
