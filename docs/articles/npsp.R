## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(fig.height=5, fig.width=7, fig.align = 'center')

## ------------------------------------------------------------------------
library(npsp)

## ---- fig.height=5, fig.width=7------------------------------------------
#   ?aquifer
str(aquifer)
#   Scatter plot with a color scale
with(aquifer, spoints(lon, lat, head, main = "Wolfcamp aquifer data"))

## ---- fig.height=7, fig.width=7------------------------------------------
summary(aquifer)
scattersplot(aquifer[,1:2], aquifer$head/100, 
             main = "Wolfcamp aquifer data", xlab = 'lon', 
             ylab='lat', zlab = 'piezometric-head  (hundreds)',
             col = jet.colors(128))

## ------------------------------------------------------------------------
#   fig.height=5, fig.width=5}
cpu.time(reset=TRUE)
bin <- binning(aquifer[,1:2], aquifer$head, nbin = c(41,41), set.NA = TRUE) 
simage(bin, main = 'Binning averages')
points(bin$data$x, col = 'darkgray')
cpu.time(total = FALSE)

## ------------------------------------------------------------------------
nbin.hd <- c(128, 128)
bin.nna <- binning(aquifer[,1:2], aquifer$head, nbin = nbin.hd)
with(bin.nna$data, summary((y - interp(bin.nna, newx = x)$y)/y ))

## ------------------------------------------------------------------------
#   fig.height=7, fig.width=7}
h.den <- diag(c(55,55))#    alternatively: h.cv(as.bin.den(bin))$h
den <- np.den(bin, h = h.den, degree = 0)
plot(den, main = 'Estimated log(density)')
#   Index with grid nodes far from data
mask <- log(den$est) > -15
bin <- mask(bin, mask = mask)
cpu.time(total = FALSE)

## ------------------------------------------------------------------------
lp <- locpol(bin, h = diag(75, 2), hat.bin = TRUE)  
#                                  np.svariso.corr: 'lp' must have a '$locpol$hat' component
# Perspective plot with a color scale
spersp(lp, main = 'Trend estimates', zlab = 'piezometric-head levels', 
       theta = 120)   
cpu.time(total = FALSE)

## ------------------------------------------------------------------------
lp.resid <- lp$data$y - predict(lp)
  maxlag <- 0.55*sqrt(sum(diff(apply(aquifer[,1:2], 2, range))^2))
esvar <- np.svariso(aquifer[,1:2], lp.resid, maxlag = 150, nlags = 60, h = 60)
svm <- fitsvar.sb.iso(esvar)  #   dk = 2
plot(svm, main = "Nonparametric semivariogram and fitted model")
cpu.time(total = FALSE)

## ------------------------------------------------------------------------
esvar2 <- np.svariso.corr(lp, maxlag = 150, nlags = 60, h = 60, plot = TRUE)
svm2 <- fitsvar.sb.iso(esvar2)  #   dk = 2
plot(svm2, main = "Nonparametric bias-corrected semivariogram and fitted models", 
     lwd = 2) 
with(svm$fit, lines(u, fitted.sv, lty = 2))
cpu.time(total = FALSE)

## ----eval = TRUE, warning=FALSE------------------------------------------
#   Example (speeding computations...):
bin2 <- binning(aquifer[,1:2], aquifer$head, nbin = c(21,21))
# Warning: There is not enough data in some neighborhoods ('NRL < NINDRL'):
h.cv(bin2, h.start = c(50, 25), objective = "GCV", ncv = 0)
# cov.bin <-  varcov(svm2, coords = coords(bin2))
lp.h <- h.cv(bin2, h.start = c(50, 25), objective = "GCV", ncv = 0, cov.bin = svm2)
lp.h
cpu.time(total = FALSE)

## ------------------------------------------------------------------------
lp <- locpol(lp, h = lp.h$h, hat.bin = TRUE)   # np.svariso.corr
# Perspective plot with a color scale
spersp(lp, main = 'Trend estimates', zlab = 'piezometric-head levels', theta = 120)   
cpu.time(total = FALSE)


## ------------------------------------------------------------------------
lp.resid <- residuals(lp)
esvar <- np.svariso(aquifer[,1:2], lp.resid, maxlag = 150, nlags = 60, h = 60)
svm <- fitsvar.sb.iso(esvar)  # dk = 2
esvar3 <- np.svariso.corr(lp, maxlag = 150, nlags = 60, h = 60, plot = FALSE)
svm3 <- fitsvar.sb.iso(esvar3, dk = 0)

plot(svm3, main = "Nonparametric bias-corrected semivariogram and fitted models", 
     lwd = 2) 
# with(svm$fit, lines(u, fitted.sv, lty = 2))
with(svm2$fit, lines(u, fitted.sv, lty = 2, lwd = 2))
cpu.time(total = FALSE)

## ------------------------------------------------------------------------
lp <- locpol(aquifer[,1:2], aquifer$head, nbin = nbin.hd, 
             h = lp$locpol$h, hat.bin = FALSE)
# Perspective plot with a color scale
simage(lp, main = 'Final trend estimates\n(piezometric-head levels)')   
cpu.time(total = FALSE)

## ----kriging-------------------------------------------------------------
krig.grid <- kriging.np(lp, svm3)
cpu.time(total = FALSE)

## ----kriging.maps--------------------------------------------------------
simage(krig.grid, 'kpred', main = 'Kriging predictions', 
       col =  jet.colors(256))
simage(krig.grid, 'ksd', main = 'Kriging sd', col =  hot.colors(256))
with(aquifer, points(lon, lat, cex = 0.75))
cpu.time()

