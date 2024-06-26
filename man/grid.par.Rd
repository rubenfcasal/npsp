% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grid.par.R
\name{grid.par}
\alias{grid.par}
\alias{grid.par-class}
\title{Grid parameters (S3 class "grid.par")}
\usage{
grid.par(
  n,
  min,
  max = min + (n - 1) * lag,
  lag = (max - min)/(n - 1),
  dimnames = names(min)
)
}
\arguments{
\item{n}{integer vector; number of nodes in each dimension.}

\item{min}{vector; minimum values of the coordinates.}

\item{max}{vector; maximum values of the coordinates (optional).}

\item{lag}{vector; lag in each dimension (optional).}

\item{dimnames}{character vector; names used to label the dimensions.}
}
\value{
Returns an object of class \code{grid.par}, a list with the arguments as components 
   and an additional component \code{$nd = length(n)}.
}
\description{
Defines a full regular (spatial) grid. 
Constructor function of the \code{grid.par}-\code{\link{class}}.
}
\details{
All parameters must have the same length. 
Only one of the arguments \code{max} or \code{lag} must be specified. 

S3 'version' of the \code{\link[sp:GridTopology-class]{GridTopology}}-\code{\link{class}} of the \pkg{sp} package.
}
\examples{
grid.par(n = c(100, 100), min = c(-10, 42), max = c(-7.5, 44))
grid.par(n = c(100, 100), min = c(-10, 42), lag = c(0.03, 0.02))
}
\seealso{
\code{\link{data.grid}}.
}
