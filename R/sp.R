#····································································
#   sp.R (npsp package)
#····································································
#   (c) R. Fernandez-Casal
#
#   NOTE: Press Ctrl + Shift + O to show document outline in RStudio
#····································································

#' Convert npsp object to sp object
#'
#' Converts a npsp object to a \link[sp:00sp]{sp} object.
#'
#' @param obj a \code{\link{npsp}} object.
#' @param ... further arguments passed to or from other methods.
#' @seealso \code{\link{as.data.frame.data.grid}}
#' @export
#····································································
as.sp <- function(obj, ...)  { 
  UseMethod("as.sp")
}


#' @rdname as.sp
#' @method as.sp grid.par
#' @return \code{as.sp.grid.par} returns a \code{\link[sp]{GridTopology-class}} object.
#' @export
#····································································
as.sp.grid.par <- function(obj, ...) { 
  with(obj, GridTopology(cellcentre.offset = min, cellsize = lag, cells.dim = n))
}  


#' @rdname as.sp
#' @method as.sp data.grid
#' @param data.ind integer or character; vector with indexes or names of the data components.
#' @param proj4string a \code{\link[sp]{CRS-class}} object.
#' @return \code{as.sp.data.grid} returns a \code{\link[sp]{SpatialGridDataFrame-class}} object. 
#' @export
#····································································
as.sp.data.grid <- function(obj, data.ind = NULL, proj4string = CRS(as.character(NA)), ...){
  SpatialGridDataFrame(as.sp(obj$grid), 
                              as.data.frame(obj, data.ind = data.ind, sp = TRUE), 
                              proj4string = proj4string)
}



