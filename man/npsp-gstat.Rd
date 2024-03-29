% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npsp-gstat.R
\name{npsp-gstat}
\alias{npsp-gstat}
\alias{as.vgm}
\alias{as.vgm.variomodel}
\alias{as.vgm.svarmod}
\alias{vgm.tab.svarmod}
\alias{as.vgm.sb.iso}
\title{Interface to package "gstat"}
\usage{
as.vgm(x, ...)

\method{as.vgm}{variomodel}(x, ...)

\method{as.vgm}{svarmod}(x, ...)

vgm.tab.svarmod(x, h = seq(0, x$range, length = 1000), sill = x$sill, ...)

\method{as.vgm}{sb.iso}(x, h = seq(0, x$range, length = 1000), sill = x$sill, ...)
}
\arguments{
\item{x}{variogram model object (used to select a method).}

\item{...}{further arguments passed to or from other methods.}

\item{h}{vector of lags at which the covariogram is evaluated.}

\item{sill}{sill of the covariogram (or pseudo-sill).}
}
\value{
A \code{variogramModel}-\code{\link{class}} object of the \pkg{gstat} package.
}
\description{
Utilities to interact with the \pkg{gstat} package.
}
\details{
Tries to convert a variogram object to \code{\link[gstat]{vgm}} 
(\code{variogramModel}-\code{\link{class}} of \pkg{gstat} package).
S3 generic function. 

\code{as.vgm.variomodel} tries to convert an object of class \code{variomodel}
defined in \pkg{geoR} (interface to \code{\link[gstat:vgm]{as.vgm.variomodel}}
defined in \pkg{gstat}).

\code{vgm.tab.svarmod} converts a \code{svarmod} object to a 
\code{variogramModel}-\code{\link{class}} object of type \code{"Tab"} 
(one-dimensional covariance table).

\code{as.vgm.sb.iso} is an alias of \code{vgm.tab.svarmod}.
}
\seealso{
\code{\link[gstat]{vgm}}, \code{\link{svarmod}}.
}
