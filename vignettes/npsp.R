## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(cache = FALSE, fig.height=4, fig.width=7, 
                      fig.align = 'center', out.width = '100%')
# knitr::spin("npsp_intro2.R", knit = FALSE)

## -----------------------------------------------------------------------------
library(npsp)

## -----------------------------------------------------------------------------
library(sp)
summary(precipitation)

## -----------------------------------------------------------------------------
spoints(precipitation)

## ----scattersplot, fig.height=9, fig.width=9----------------------------------
scattersplot(precipitation)

## -----------------------------------------------------------------------------
x <- coordinates(precipitation)
y <- precipitation$y
lp <- locpol(x, y, nbin = c(120, 120), h = diag(c(5, 5)))

## -----------------------------------------------------------------------------
attr <- attributes(precipitation)
border <- attr$border 
# Mask grid points outside of continental part of USA
lp <- mask(lp, window = border)
interior <- attr$interior 
slim <- range(y)
col <- jet.colors(256)
simage(lp,  slim = slim, main = "Trend estimates",
       col = col, xlab = "Longitude", ylab = "Latitude", asp = 1)
plot(border, border = "darkgray", lwd = 2, add = TRUE)
plot(interior, border = "lightgray", lwd = 1, add = TRUE)

## -----------------------------------------------------------------------------
bin <- binning(x, y, nbin = c(120, 120), window = border)
# lp <- locpol(bin, h = diag(c(5, 5)))
lp2 <- locpol(bin, h = diag(c(2, 2)))
simage(lp2, slim = slim, main = "Trend estimates",
       col = col, xlab = "Longitude", ylab = "Latitude", asp = 1)
plot(border, border = "darkgray", lwd = 2, add = TRUE)
plot(interior, border = "lightgray", lwd = 1, add = TRUE)

## -----------------------------------------------------------------------------
bin <- binning(x, y, window = border)
lp0.h <- h.cv(bin)$h
lp0 <- locpol(bin, h = lp0.h, hat.bin = TRUE) 

## -----------------------------------------------------------------------------
svar.bin <- svariso(x, residuals(lp0), nlags = 60, maxlag = 20)  
svar.h <- h.cv(svar.bin)$h

svar.np <- np.svar(svar.bin, h = svar.h)
svar.np2 <- np.svariso.corr(lp0, nlags = 60, maxlag = 20, 
                            h = svar.h, plot = FALSE) 

## -----------------------------------------------------------------------------
svm0 <- fitsvar.sb.iso(svar.np, dk = 0) 
svm1 <- fitsvar.sb.iso(svar.np2, dk = 0) 
# plot...
plot(svm1, main = "Nonparametric bias-corrected semivariogram\nand fitted models", 
     legend = FALSE, xlim = c(0,max(coords(svar.np2))), 
     ylim = c(0,max(svar.np2$biny, na.rm = TRUE)))
plot(svm0, lwd = c(1, 1), add = TRUE)
plot(svar.np, type = "p", pch = 2, add = TRUE)
# abline(h = c(svm1$nugget, svm1$sill), lty = 3)
# abline(v = 0, lty = 3)
legend("bottomright", legend = c("corrected", 'biased'),
            lty = c(1, 1), pch = c(1, 2), lwd = c(1, 1))

## ----geomod, fig.height=4, fig.width=12---------------------------------------
geomod <- np.fitgeo(x, y, nbin = c(30, 30), maxlag = 20, svm.resid = TRUE, window = border)
plot(geomod)

## ----kriging, fig.height=4, fig.width=12--------------------------------------
krig.grid <- np.kriging(geomod, ngrid = c(120, 120))
# Plot kriging predictions and kriging standard deviations
old.par <- par(mfrow = c(1,2), omd = c(0.0, 0.98, 0.0, 1))
simage(krig.grid, 'kpred', main = 'Kriging predictions', slim = slim,
      xlab = "Longitude", ylab = "Latitude" , col = col, asp = 1, reset = FALSE)
plot(border, border = "darkgray", lwd = 2, add = TRUE)
plot(interior, border = "lightgray", lwd = 1, add = TRUE)
simage(krig.grid, 'ksd', main = 'Kriging sd', xlab = "Longitude", 
       ylab = "Latitude" , col = hot.colors(256), asp = 1, reset = FALSE)
plot(border, border = "darkgray", lwd = 2, add = TRUE)
plot(interior, border = "lightgray", lwd = 1, add = TRUE)
par(old.par)

