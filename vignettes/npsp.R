## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(cache = FALSE, fig.height=5, fig.width=7, fig.align = 'center')
# knitr::opts_chunk$set(comment=NA, prompt=TRUE, dev='svg', out.width=1024, fig.height=6, fig.width=8)

## ----library, warning=TRUE-----------------------------------------------
library(npsp)
library(sp)

## ----Datos---------------------------------------------------------------
# Total Monthly Precipitation
spoints(precipitation)

n <- length(precipitation$y)
coord <- coordinates(precipitation)
attributes <- attributes(precipitation)
labels <- attributes$labels 
border <- attributes$border 
interior <- attributes$interior 
slim <- range(precipitation$y)
col <- jet.colors(256)
col2 <- hot.colors(256)
cpu.time(reset=TRUE)

## ----summary-------------------------------------------------------------
y.summary <- summary(precipitation$y) # Resumen de datos
y.summary

scattersplot(precipitation)
cpu.time(total = FALSE)

## ----np.fitgeo-----------------------------------------------------------
nbin <- c(30, 30)
geom <- np.fitgeo(coord, precipitation$y, nbin = nbin, svm.resid = TRUE)
cpu.time(total = FALSE)

## ------------------------------------------------------------------------
spp.grid <- SpatialPoints(coords(geom))
proj4string(spp.grid) <- proj4string(border) # CRS("+init=epsg:28992 +units=km")
mask.sp <- !is.na(over(spp.grid, as(border, 'SpatialPolygons')))
geom <- mask(geom, mask = mask.sp | (geom$binw > 0))
cpu.time(total = FALSE)

## ---- fig.height=5, fig.width=9------------------------------------------
plot(geom)

## ----binning-------------------------------------------------------------
cpu.time(reset=TRUE)
bin <- binning(coord, precipitation$y, nbin = nbin, set.NA = TRUE)
cpu.time(total = FALSE)

simage(bin, main = 'Binning averages and data points', slim = slim,
       col = col2, xlab = labels$x[1], ylab = labels$x[2],
       sub = paste0("(", paste(dim(bin), collapse = "x"), ")"))
plot(border, border = "darkgray", lwd = 2, add = TRUE)
plot(interior, border = "lightgray", lwd = 2, add = TRUE)

points(bin$data$x, pch ='+')
coordvs <- coordvalues(bin)
abline(v = coordvs[[1]], lty = 3)
abline(h = coordvs[[2]], lty = 3)

## ----lp0-----------------------------------------------------------------
lp0.h <- h.cv(bin)$h
lp0.h    # <- diag(c(6.85061, 2.946348))   

# Initial linear Local trend estimation
lp0 <- locpol(bin, h = lp0.h, hat.bin = TRUE) 

cpu.time(total = FALSE)

simage(lp0, main = "Initial trend estimates", slim = slim,
       col = col2, xlab = labels$x[1], ylab = labels$x[2])
plot(border, border = "darkgray", lwd = 2, add = TRUE)
plot(interior, border = "lightgray", lwd = 1, add = TRUE)
points(coord[,1], coord[,2], col="darkgray")

# Residuals
lp0.pred <- predict(lp0)
lp0.resid <- lp0$data$y - lp0.pred
lp0.r2 <- cor(lp0.pred, lp0$data$y)^2 # assuming independence...

old.par <- par(mfrow=c(1,2))
hist(lp0.resid)
plot(lp0.pred, lp0.resid, main = "Residuals vs Fitted",
     sub = paste("R^2 =", round(lp0.r2, 2)),
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, lty = 2, col = "darkgray")
par(old.par)

## ----svm0----------------------------------------------------------------
# rule.svar(coord)
nlags <- 60 
maxlag <- 10
# Empirical variogram with linear binning:
svar.bin <- svariso(coord, lp0.resid,  nlags = nlags, maxlag = maxlag)  
sv.lags <- coords(svar.bin)
svar.np.h <- h.cv(svar.bin)$h
# Local linear variogram estimation
svar.np <- np.svar(svar.bin, h = svar.np.h)  # biased
# Fitted Shapiro-Botha variogram model
svm0 <- fitsvar.sb.iso(svar.np, dk = 0)       
# Bias-corrected variogram estimation
svar.np2 <- np.svariso.corr(lp0, nlags = nlags, maxlag = maxlag, 
                            h=svar.np.h, plot = TRUE) 
# Fitted Shapiro-Botha variogram model
svm02 <- fitsvar.sb.iso(svar.np2, dk = 0) 
cpu.time(total = FALSE)

plot(svm02, main = "Nonparametric bias-corrected semivariogram\nand fitted models", 
     legend = FALSE, xlim = c(0,max(coords(svar.np2))), 
     ylim = c(0,max(svar.np2$biny, na.rm = TRUE)))
plot(svm0, add = TRUE)
plot(svar.np, type = "p", pch = 2, add = TRUE)
abline(h = c(svm02$nugget, svm02$sill), lty = 3)
abline(v = 0, lty = 3)
legend("bottomright", legend = c("corrected", 'biased'),
            lty = c(1, 1), pch = c(NA, NA), lwd = c(2, 1))

## ----hcv.data------------------------------------------------------------
# GCV criterion (Francisco-Fernandez and Opsomer, 2005)
lp.h <- h.cv(bin, objective = "GCV", cov = svm02, DEalgorithm = FALSE)$h
lp.h    # <- diag(c(11.11463, 18.59536))
cpu.time(total = FALSE)
# lp.h <- hcv.data(bin, objective = "GCV", cov = svm02, DEalgorithm = FALSE)$h
# lp.h    # <- diag(c(10.0871, 18.3956))
## Time of last operation: hcv.data
##    user  system elapsed 
##   10.39    1.39   11.87

## ----lp------------------------------------------------------------------

## ------------------------------------------------------------------------
lp <- locpol(bin, h = lp.h, hat.bin = TRUE)
cpu.time(total = FALSE)
simage(lp, main = "Trend Estimation", slim = slim, 
       xlab = labels$x[1], ylab = labels$x[2], col = col2)
plot(border, border = "darkgray", lwd = 2, add = TRUE)
plot(interior, border = "lightgray", lwd = 1, add = TRUE)
# spoints(coord[,1], coord[,2], precipitation$y, add = TRUE)
# New Residuals
lp.pred <- predict(lp)
lp.resid <- lp$data$y - lp.pred
lp.r2 <- cor(lp.pred, lp$data$y)^2 # assuming independence...
old.par <- par(mfrow=c(1,2))
hist(lp.resid)
plot(lp.pred, lp.resid, main = "Residuals vs Fitted", 
     sub = paste("R^2 =", round(lp0.r2, 2)),
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, lty = 2, col = "darkgray")
par(old.par)

## ----lp.svm2-------------------------------------------------------------

svar.bin <- svariso(coord, lp.resid,  nlags = nlags, maxlag = maxlag)  
svar.np.h <- h.cv(svar.bin)$h  # svar.np.h <- 2.257358
svar.np <- np.svar(svar.bin, h = svar.np.h)
svm <- fitsvar.sb.iso(svar.np, dk = 0)
svar.np2 <- np.svariso.corr(lp, nlags = nlags, maxlag = maxlag, 
                            h = svar.np.h, plot = TRUE, ylim = c(0, 0.3))

## ------------------------------------------------------------------------
svm2 <- fitsvar.sb.iso(svar.np2, dk = 0) # Shapiro-Botha model

cpu.time(total = FALSE)

plot(svm2, main = "Nonparametric bias-corrected semivariogram\nand fitted models", 
     legend = FALSE)
plot(svm, add = TRUE)
plot(svar.np, type = "p", pch = 2, add = TRUE)
abline(h = c(svm2$nugget, svm2$sill), lty = 3)
abline(v = 0, lty = 3)
plot(svm0, lty = 2, lwd = 1, add = TRUE)
plot(svm02, lty = 2, lwd = 2, add = TRUE)
legend("bottomright", legend = c("corrected", 'biased', "corrected initial",
       'biased initial'), lty = c(1, 1, 2, 2), pch = c(1, 2, NA, NA), 
       lwd = c(2, 1))

cpu.time(total = FALSE)

## ----bin.hd--------------------------------------------------------------

## ------------------------------------------------------------------------
bin.hd <- binning(coord, precipitation$y,
               nbin = c(120,120), set.NA = TRUE)
simage(bin.hd, main = 'Binning averages',
       xlab = labels$x[1], ylab = labels$x[2],
       sub = paste0("(", paste(dim(bin.hd), collapse = "x"), ")"),
       slim = slim)
plot(border, border = "darkgray", lwd = 2, add = TRUE)
plot(interior, border = "lightgray", lwd = 1, add = TRUE)

## ----h.den---------------------------------------------------------------
h.den <- h.cv(as.bin.den(bin), ncv = 2)$h         # diag(0.4061458, 0.6349045)
den <- np.den(bin.hd, h = h.den, degree = 0)
plot(den, main = 'Estimated log(density)',
          xlab = labels$x[1], ylab = labels$x[2],
          sub = paste0("(", paste(dim(bin.hd), collapse = "x"), ")") )
plot(border, border = "darkgray", lwd = 2, add = TRUE)
plot(interior, border = "lightgray", lwd = 1, add = TRUE)

## ----mask----------------------------------------------------------------
spp.grid <- SpatialPoints(coords(bin.hd))
proj4string(spp.grid) <- proj4string(border) # CRS("+init=epsg:28992 +units=km")
mask.sp <- !is.na(over(spp.grid, as(border, 'SpatialPolygons')))
mask <- mask.sp | (bin.hd$binw > 0)  #indice
bin.hd <- mask(bin.hd, mask = mask)

## ------------------------------------------------------------------------
lp.hd <- locpol(bin.hd, h = lp.h)

simage(lp.hd, main = "Final trend estimates", slim = slim, 
       xlab = labels$x[1], ylab = labels$x[2], col = col2)
plot(border, border = "darkgray", lwd = 2, add = TRUE)
plot(interior, border = "lightgray", lwd = 1, add = TRUE)

## ----lp.def--------------------------------------------------------------
lp <- lp.hd
trend.est <- predict(lp)
lp.resid <- lp$data$y - trend.est
lp$mask <- NULL
lp <- mask(lp, mask = mask.sp, warn = FALSE)  # warning=FALSE

cpu.time(total = FALSE)

## ----kriging-------------------------------------------------------------

krig.grid <- kriging.np(lp, svm2, lp.resid)
cpu.time(total = FALSE)

## ----kriging.maps--------------------------------------------------------

simage(krig.grid, 'kpred', main = 'Kriging predictions', slim = slim,
                  xlab = labels$x[1], ylab = labels$x[2])
plot(border, border = "darkgray", lwd = 2, add = TRUE)
plot(interior, border = "lightgray", lwd = 1, add = TRUE)

simage(krig.grid, 'ksd', main = 'Kriging sd',
                  xlab = labels$x[1], ylab = labels$x[2], col = col2)
plot(border, border = "darkgray", lwd = 2, add = TRUE)
plot(interior, border = "lightgray", lwd = 1, add = TRUE)

cpu.time()
# save.image(".RData")

