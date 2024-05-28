---
title: 'Kriging with gstat'
author: 'Ruben Fernandez-Casal (ruben.fcasal@udc.es)'
date: 'npsp 0.7.13'
output: rmarkdown::html_vignette
vignette: >
  %\VignetteEngine{knitr::knitr}
  %\VignetteIndexEntry{Kriging with gstat}
  %\usepackage[UTF-8]{inputenc}
---



# Kriging with ***gstat***

You can use the `krige` (or `krige.cv`) utilities in `gstat` package together with `as.vgm` for (global or local) kriging...


## Wolfcamp aquifer data:


```r
library(npsp)
# ?aquifer; str(aquifer); summary(aquifer)
# Scatter plot with a color scale
with(aquifer, spoints(lon, lat, head, main = "Wolfcamp aquifer data"))
```

<div class="figure" style="text-align: center">
<img src="figure/unnamed-chunk-1-1.png" alt="plot of chunk unnamed-chunk-1"  />
<p class="caption">plot of chunk unnamed-chunk-1</p>
</div>


## Kriging


### Trend estimation


```r
lp <- locpol(aquifer[,1:2], aquifer$head, h = diag(75, 2), 
             hat.bin = TRUE)  
            # np.svariso.corr: 'lp' must have a '$locpol$hat' component

# Mask grid nodes far from data
mask <- log(np.den(lp, h = diag(c(55,55)), degree = 0)$est) > -15
lp <- mask(lp, mask = mask)

spersp(lp, main = 'Trend estimates', 
       zlab = 'piezometric-head levels', theta = 120)   
```

<div class="figure" style="text-align: center">
<img src="figure/unnamed-chunk-2-1.png" alt="plot of chunk unnamed-chunk-2"  />
<p class="caption">plot of chunk unnamed-chunk-2</p>
</div>

```r
cpu.time(total = FALSE)
```

```
## Time of last operation: 
##    user  system elapsed 
##   36.09    1.57   39.01
```
    
    
### Variogram estimation


```r
lp.resid <- residuals(lp)
esvar <- np.svariso(aquifer[,1:2], lp.resid, maxlag = 150, nlags = 60, h = 60)
svm <- fitsvar.sb.iso(esvar)  # dk = 2
esvar2 <- np.svariso.corr(lp, maxlag = 150, nlags = 60, h = 60)
svm2 <- fitsvar.sb.iso(esvar2, dk = 0)  # dk = Inf
plot(svm2, main = "Nonparametric bias-corrected semivariogram and fitted models", 
     lwd = 2) 
with(svm$fit, lines(u, fitted.sv, lty = 2))
```

<div class="figure" style="text-align: center">
<img src="figure/unnamed-chunk-3-1.png" alt="plot of chunk unnamed-chunk-3"  />
<p class="caption">plot of chunk unnamed-chunk-3</p>
</div>

```r
cpu.time(total = FALSE)
```

```
## Time of last operation: 
##    user  system elapsed 
##    0.06    0.03    0.10
```

## Residual Kriging


```r
library(sp)
library(gstat)
```

```
## 
## Attaching package: 'gstat'
```

```
## The following object is masked from 'package:npsp':
## 
##     as.vgm.variomodel
```

```r
spdf <- SpatialPointsDataFrame(aquifer[,1:2], 
            data.frame(y = aquifer$head, r = lp.resid))
newdata <- SpatialPoints(coords(lp))
krig <- krige(r ~ 1, locations = spdf, newdata = newdata, 
              model = as.vgm(svm), beta = 0)
```

```
## [using simple kriging]
```

```r
krig.grid <- data.grid(kpred = lp$est + krig@data$var1.pred, 
                       ksd = sqrt(krig@data$var1.var), 
        grid = lp$grid)
krig2 <- krige(r ~ 1, locations = spdf, newdata = newdata, 
               model = as.vgm(svm2), beta = 0)
```

```
## [using simple kriging]
```

```r
krig2.grid <- data.grid(kpred = lp$est + krig2@data$var1.pred, 
                        ksd = sqrt(krig2@data$var1.var), 
        grid = lp$grid)
scale.color <- jet.colors(64)
scale.range <- c(1100, 4100)
# 1x2 plot with some room for the legend...
old.par <- par(mfrow = c(1,2), omd = c(0.01, 0.9, 0.05, 0.95),
               plt= c(0.08, 0.94, 0.1, 0.8))
spersp(krig.grid, main = 'Kriging predictions', col = scale.color, 
       legend = FALSE, theta = 120, reset = FALSE)
spersp(krig2.grid, main = 'Kriging predictions \n (bias-corrected)', 
       col = scale.color, legend = FALSE, theta = 120, reset = FALSE)
par(old.par)
splot(slim = scale.range, col = scale.color, legend.shrink = 0.6, add = TRUE)
```

<div class="figure" style="text-align: center">
<img src="figure/unnamed-chunk-4-1.png" alt="plot of chunk unnamed-chunk-4" width="100%" />
<p class="caption">plot of chunk unnamed-chunk-4</p>
</div>

```r
old.par <- par(mfrow = c(1,2), omd = c(0.05, 0.85, 0.05, 0.95))
scale.range <- c(125, 200)
scale.range <- range(krig.grid$ksd, krig2.grid$ksd, finite = TRUE)
image( krig.grid, 'ksd', zlim = scale.range, # asp = 1, 
       main = 'Kriging sd', col = scale.color)
with(aquifer, points(lon, lat, cex = 0.75))
image( krig2.grid, 'ksd', zlim = scale.range, # asp = 1, 
       main = 'Kriging sd (bias-corrected)', col = scale.color)
with(aquifer, points(lon, lat, cex = 0.75))
par(old.par)
splot(slim = scale.range, col = scale.color, add = TRUE)
```

<div class="figure" style="text-align: center">
<img src="figure/unnamed-chunk-5-1.png" alt="plot of chunk unnamed-chunk-5" width="100%" />
<p class="caption">plot of chunk unnamed-chunk-5</p>
</div>

```r
cpu.time()
```

```
## Time of last operation: 
##    user  system elapsed 
##    0.42    0.11    0.53 
## Total time:
##    user  system elapsed 
##   36.57    1.71   39.64
```

### Notes: 
* To reproduce results in SERRA paper use `data(wolfcamp)` in package `geoR`.
* Results obtained with `aquifer` data set are comparable with those in Cressie (1993, section 4.1). 



