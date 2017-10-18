#--------------------------------------------------------------------
#   coords.R (npsp package)
#--------------------------------------------------------------------
#   coords  S3 generic
#           coords.grid.par(x, ...)
#           coords.data.grid(x, ...)
#   coordvalues  S3 generic
#           coordvalues.grid.par(x, ...)
#           coordvalues.data.grid(x, ...)
#
#   (c) Ruben Fernandez-Casal
#   Created: Aug 2012                          Last changed: Dec 2015
#--------------------------------------------------------------------


#--------------------------------------------------------------------
#' (spatial) coordinates
#'
#' Retrieves the (spatial) coordinates of the object.
#'
#' @param  x 	a (spatial) object used to select a method.
#' @param  ... 	further arguments passed to or from other methods.
#' @return A matrix of coordinates (columns correspond with dimensions and rows with data).
#' @seealso \code{\link{coordvalues}}.
#' @export
coords <- function(x, ...) UseMethod("coords")
# S3 generic function coordinates
# PENDENTE: renomear, deberia funcionar como 'coordinates' dentro do NAMESPACE
#--------------------------------------------------------------------

#--------------------------------------------------------------------
#' Coordinate values
#'
#' Returns the coordinate values in each dimension.
#' @param  x 	a (spatial) object used to select a method.
#' @param  ... 	further arguments passed to or from other methods.
#' @return A list with the (unique) coordinates along each axis.
#' @seealso \code{\link{coords}}.
#' @export
coordvalues <- function(x, ...) UseMethod("coordvalues")
# S3 generic function coordsvalues
# PENDENTE: renombrar 'coordinatevalues'?
#--------------------------------------------------------------------


#--------------------------------------------------------------------
# grid.par
#--------------------------------------------------------------------

#--------------------------------------------------------------------
#' @rdname coords
#' @method coords grid.par
#' @export
coords.grid.par <- function(x, ...) {
#--------------------------------------------------------------------
    if (!inherits(x, "grid.par"))
      stop("function only works for objects of class (or extending) 'grid.par'")
    co <- coordvalues.grid.par(x)
    # if (grid$nd == 1) return(co[[1]])
    cc <- drop(do.call("expand.grid", co))
    return(do.call("cbind", lapply(cc, as.numeric)))
}


#--------------------------------------------------------------------
#' @rdname coordvalues
#' @method coordvalues grid.par
#' @export
coordvalues.grid.par <- function(x, ...) {
#--------------------------------------------------------------------
    if (!inherits(x, "grid.par"))
      stop("function only works for objects of class (or extending) 'grid.par'")
    result <- list()
    for (i in seq_along(x$n))  result[[i]] <- with(x, min[i]+(0:(n[i]-1))*lag[i])
    names(result) <- dimnames(x)
    return(result)
}


#--------------------------------------------------------------------
# data.grid
#--------------------------------------------------------------------

#--------------------------------------------------------------------
#' @rdname coords
#' @method coords data.grid
#' @param  masked logical; If \code{TRUE}, only the coordinates corresponding
#' to unmasked cells are returned (see \code{\link{mask}}).
#' @export
coords.data.grid <- function(x, masked = FALSE, ...) {
#--------------------------------------------------------------------
# COIDADO: variable mask en data grid
    if (!is(x, "data.grid"))
      stop("function only works for objects of class (or extending) 'data.grid'")
    res <- coords.grid.par(x$grid)
    if (masked && !is.null(x$mask)) res <- res[x$mask, ]
    return(res)
}

#--------------------------------------------------------------------
#' @rdname coordvalues
#' @method coordvalues data.grid
#' @export
coordvalues.data.grid <- function(x, ...) {
#--------------------------------------------------------------------
    if (!is(x, "data.grid"))
      stop("function only works for objects of class (or extending) 'data.grid'")
    return( coordvalues.grid.par(x$grid) )
}
