% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npsp-geoR.R
\name{npsp-geoR}
\alias{npsp-geoR}
\alias{as.variogram}
\alias{variogram}
\alias{as.variogram.svar.bin}
\alias{as.variogram.np.svar}
\alias{as.variomodel}
\alias{variomodel}
\alias{as.variomodel.svarmod}
\title{Interface to package "geoR"}
\usage{
as.variogram(x, ...)

\method{as.variogram}{svar.bin}(x, ...)

\method{as.variogram}{np.svar}(x, ...)

as.variomodel(m, ...)

\method{as.variomodel}{svarmod}(m, ...)
}
\arguments{
\item{x}{semivariogram estimate (e.g. \code{\link{svar.bin}} or \code{\link{np.svar}} object).}

\item{...}{further arguments passed to or from other methods.}

\item{m}{variogram model (e.g. \code{\link{svarmod}} object).}
}
\value{
\code{as.variogram()} returns an object of the (not fully documented) \pkg{geoR}-class 
\code{variogram}.

\code{as.variomodel()} returns an object of the \pkg{geoR}-class \code{variomodel}.
}
\description{
Utilities to interact with the \pkg{geoR} package.
}
\details{
\code{as.variogram} tries to convert a semivariogram estimate \eqn{\hat{\gamma}(h_i)} 
to an object of the (not fully documented) \pkg{geoR}-class \code{variogram} 
(see e.g. \code{\link[geoR]{variog}}).

\code{as.variomodel} tries to convert a semivariogram model \eqn{\gamma(pars; h)} 
to an object of the \pkg{geoR}-class \code{variomodel} 
(see e.g. \code{\link[geoR]{variofit}}).
}
\seealso{
\code{\link[geoR]{variog}}, \code{\link[geoR]{variofit}}, \code{\link{variomodel}}, 
\code{\link{svar.bin}}, \code{\link{np.svar}}.
}
