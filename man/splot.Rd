% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splot.R
\name{splot}
\alias{splot}
\alias{scolor}
\alias{jet.colors}
\alias{hot.colors}
\title{Utilities for plotting with a color scale}
\usage{
splot(
  slim = c(0, 1),
  col = jet.colors(128),
  breaks = NULL,
  horizontal = FALSE,
  legend.shrink = 0.9,
  legend.width = 1.2,
  legend.mar = ifelse(horizontal, 3.1, 5.1),
  legend.lab = NULL,
  bigplot = NULL,
  smallplot = NULL,
  lab.breaks = NULL,
  axis.args = NULL,
  legend.args = NULL,
  add = FALSE
)

scolor(s, col = jet.colors(128), slim = range(s, finite = TRUE))

jet.colors(n)

hot.colors(n, rev = TRUE)
}
\arguments{
\item{slim}{limits used to set up the color scale.}

\item{col}{color table used to set up the color scale (see \code{\link{image}} for
details).}

\item{breaks}{(optional) numeric vector with the breakpoints for the color scale: 
must have one more breakpoint than \code{col} and be in increasing order.}

\item{horizontal}{logical; if \code{FALSE} (default) legend will be a vertical strip on the
right side. If \code{TRUE} the legend strip will be along the bottom.}

\item{legend.shrink}{amount to shrink the size of legend relative to the
full height or width of the plot.}

\item{legend.width}{width in characters of the legend strip. Default is
1.2, a little bigger that the width of a character.}

\item{legend.mar}{width in characters of legend margin that has the axis.
Default is 5.1 for a vertical legend and 3.1 for a horizontal legend.}

\item{legend.lab}{label for the axis of the color legend. Default is no
label as this is usual evident from the plot title.}

\item{bigplot}{plot coordinates for main plot. If not passed these will be 
determined within the function.}

\item{smallplot}{plot coordinates for legend strip. If not passed these
will be determined within the function.}

\item{lab.breaks}{if breaks are supplied these are text string labels to
put at each break value. This is intended to label axis on a transformed
scale such as logs.}

\item{axis.args}{additional arguments for the axis function used to create
the legend axis (see \code{\link[fields]{image.plot}} for details).}

\item{legend.args}{arguments for a complete specification of the legend
label. This is in the form of list and is just passed to the \code{\link{mtext}}
function. Usually this will not be needed (see \code{\link[fields]{image.plot}} 
for details).}

\item{add}{logical; if \code{TRUE} the legend strip is just added 
to the existing plot (the graphical parameters are not changed).}

\item{s}{values to be converted to the color scale.}

\item{n}{number of colors (\code{>= 1}) to be in the palette.}

\item{rev}{logical; if \code{TRUE}, the palette is reversed (decreasing overall luminosity).}
}
\value{
\code{splot} invisibly returns a list with the following 3 components:
\item{bigplot}{plot coordinates of the main plot. These values may be useful for 
drawing a plot without the legend that is the same size as the plots with legends.}
\item{smallplot}{plot coordinates of the secondary plot (legend strip).}
\item{old.par}{previous graphical parameters (\code{par(old.par)} 
will reset plot parameters to the values before entering the function).}

\code{jet.colors} and \code{hot.colors} return a character vector of colors (similar to 
\code{\link{heat.colors}} or \code{\link{terrain.colors}}; see \code{\link{rgb}}).
}
\description{
\code{splot} is designed to combine a standard R plot with 
a legend representing a (continuous) color scale. This is done by splitting 
the plotting region into two parts. Keeping one for the main chart and 
putting the legend in the other.  
For instance, \code{sxxxx} functions (\code{\link{spoints}}, \code{\link{simage}} 
and \code{\link{spersp}}) draw the corresponding high-level plot (\code{xxxx}),
after calling \code{splot}, to include a legend strip for the color scale. 

These functions are based on function \code{\link[fields]{image.plot}} of package
\pkg{fields}, see its documentation for additional information.

\code{jet.colors} and \code{hot.colors} create a color table useful for contiguous 
color scales and \code{scolor} assigns colors to a numerical vector.
}
\details{
\code{scolor} converts a real valued vector to a color scale. The range 
\code{slim} is divided into \code{length(col) + 1} pieces of equal length.
Values which fall outside the range of the scale are coded as \code{NA}.

\code{jet.colors} generates a rainbow style color table similar to the MATLAB (TM) 
jet color scheme. It may be appropriate to distinguish between values above and 
below a central value (e.g. between positive and negative values).

\code{hot.colors} generates a color table similar to the MATLAB (TM) 
hot color scheme (reversed by default). It may be appropriate to represent values  
ranging from 0 to some maximum level (e.g. density estimation). 
The default value \code{rev = TRUE} may be adecuate to grayscale convertion.
}
\section{Side Effects}{
 After exiting \code{splot}, the plotting region may be changed 
(\code{\link{par}("plt")}) to make it possible to add more features to the plot.
}

\examples{
# Plot equivalent to spoints():
scale.range <- range(aquifer$head)
res <- splot(slim = scale.range)
with( aquifer, plot(lon, lat, col = scolor(head, slim = scale.range), 
       pch = 16, cex = 1.5, main = "Wolfcamp aquifer data"))
par(res$old.par) # restore graphical parameters       
# Multiple plots with a common legend:
# regularly spaced 2D data...
set.seed(1)
nx <- c(40, 40) # ndata =  prod(nx)
x1 <- seq(-1, 1, length.out = nx[1])
x2 <- seq(-1, 1, length.out = nx[2])
trend <- outer(x1, x2, function(x,y) x^2 - y^2)
y <- trend + rnorm(prod(nx), 0, 0.1)
scale.range <- c(-1.2, 1.2)
scale.color <- heat.colors(64)
# 1x2 plot with some room for the legend...
old.par <- par(mfrow = c(1,2), omd = c(0.05, 0.85, 0.05, 0.95))
image( x1, x2, trend, zlim = scale.range, main = 'Trend', col = scale.color)
image( x1, x2, y, zlim = scale.range, main = 'Data', col = scale.color)
par(old.par)
# the legend can be added to any plot...
splot(slim = scale.range, col = scale.color, add = TRUE)
## note that argument 'zlim' in 'image' corresponds with 'slim' in 'sxxxx' functions. 
}
\seealso{
\code{\link{spoints}}, \code{\link{simage}}, \code{\link{spersp}}, 
\code{\link{image}}, \code{\link[fields]{image.plot}}.
}
\author{
Based on \code{\link[fields]{image.plot}} function from package \pkg{fields}:
fields, Tools for spatial data. 
Copyright 2004-2013, Institute for Mathematics Applied Geosciences. 
University Corporation for Atmospheric Research.

Modified by Ruben Fernandez-Casal <rubenfcasal@gmail.com>.
}
\keyword{hplot}
