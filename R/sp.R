#' Convert npsp object to sp object
#'
#' Converts a npsp object to a \link[sp]{sp} object.
#'
#' @param obj a \code{\link{npsp}} object.
#' @param ... further arguments passed to or from other methods.
#' @seealso \code{\link{as.data.frame.data.grid}}
#' @export
as.sp <- function(obj, ...) UseMethod("as.sp")


#' @rdname as.sp
#' @method as.sp grid.par
#' @export
as.sp.grid.par <- function(obj, ...)
  return(with(obj, GridTopology(cellcentre.offset = min, cellsize = lag, cells.dim = n)))


#' @rdname as.sp
#' @method as.sp data.grid
#' @param data.ind integer or character; vector with indexes or names of the data components.
#' @param proj4string an object of class \code{\link[sp]{CRS-class}}.
as.sp.data.grid <- function(obj, data.ind = NULL, proj4string = CRS(as.character(NA)), ...){
  return(SpatialGridDataFrame(as.sp(obj$grid), 
                              as.data.frame(obj, data.ind = data.ind, sp = TRUE), 
                              proj4string = proj4string))
}



