# npsp 0.7-1

* Changes in `kriging.np`: 

  - S3 generic function. 

  - `kriging.np` renamed as `kriging.np.default`.

  - Added `kriging.np.np.geo` S3 method.

  - Added `ngrid` parameter.


# npsp 0.7-0

* Added `np.geo` S3 class (nonparametric geostatistical model),  
  constructor function and methods. 

  - Added `plot()` S3 method for `np.geo` class.

  - Added `residuals()` S3 method for `np.geo` class.

* Added `np.fitgeo()` S3 generic function and methods
 (`np.fitgeo.default`, `np.fitgeo.locpol.bin` and `np.fitgeo.np.geo`).
 
* Fixed bug in `np.svariso.hcv()` (thanks to Tomas Cotos-Yañez).


# npsp 0.6-2 

*  Added a website for the package (with pkgdown).

*  Added 'NEWS.md' and 'index.Rmd'.

*  Added some vignettes (pkgdown articles): 
   "npsp.Rmd", "precipitation.Rmd", "krigstat.Rmd", "docs/aquifer.Rmd", "docs/Introduccion.Rmd".


# npsp 0.6-1 

*  Added `scattersplot()` S3 generic function (and methods).

*  Added (some) support for `sp` classes.
 
    - Added `spoints()` and `scattersplot()` methods for objects of 
      class `SpatialPointsDataFrame`.
    
    - Added `as.sp()` generic function.
    
    - Added `as.sp.grid.par()` and `as.sp.data.grid()` methods.
 
*  Added `precipitation` data set.
 
*  Added `as.data.frame.data.grid` 

 
# npsp 0.6-0 

*  Added `kriging.np()` and `kriging.simple()` functions.
 
    - Added `kriging.simple.solve()` internal function.
 
    - Methods `as.spam()`, `chol.spam()` and `solve.spam()` imported from package `spam`.
    
    - Added `.DPOSV_R()` interface to LAPACK routine `DPOSV`.
 
*  Minor changes in `plot.fitsvar()` 
   (`lwd` parameter is passed to `lines()` when `add = TRUE`).

*   Added the registration of 'native routines' (`.Fortran` calls). 
 
    - Added 'src/init.c' and `@useDynLib npsp, .registration = TRUE`.

    - Fortran routine `binning` renamed as `binning_r`.
 
*  Updated 'README.md'.
 

# npsp 0.5-5 

*  Added `svar.grid()` S3 class (discretized semivariogram),  
   generic function (constructor) and methods. 
   
    - Added `svar.grid.fitsvar()` and `svar.grid.svarmod()` methods.

    - Added `sv()` and `plot()` S3 methods for `svar.grid` class.

*  Major changes in `varcov.isotropic()`. 

    - Returns 0 if `h < .Machine$double.eps`.
   
    - Added `discretize` parameter 
      (if `TRUE`, the default value, the variogram is previously discretized). 

*  Minor changes in `covar.svarmod()`. Added `discretize` parameter 
   (if `TRUE` the variogram is discretized as a first step). 

*  Added `plot()` S3 method for `svarmod` class.

   
# npsp 0.5-4 

*  Minor changes in `covar.svarmod()`
   (argument `...` is passed to `sv()`).

*  `fitsvar.sb.iso` returns additional components
   (`$fit$w` and `$esv`).

*  Added  `.DNRM2_R()` internal function (interface to BLAS routine DNRM2).

*  Minor changes in FORTRAN code
   (to avoid warnings and obsolescent features: `tql2.f`, `lp_module.f90`).

*  Fixed bug in `mask.bin.data()` and `mask.locpol.bin()` 
   (when `warn = FALSE`, now it is not changed by `filter.lp`).

*  Fixed bug in FORTRAN function KTW(u)
   (in the normalizing constant; thanks to Tomas Cotos-Yañez).


# npsp 0.5-3 or older 

See [npsp ChangeLog](https://github.com/rubenfcasal/npsp/blob/master/ChangeLog) (also for more info). 
