# Development version



# npsp 0.7-13 (2024-02-17)

* Added `useRaster = all(dim(x) > dev.size("px"))` argument to `image()` 
  (and `simage()`) methods for gridded data.

* Small changes in *scr/tql2.f90* (preliminary translation to Fortran 90 of the 
  former *scr/tql2.f*).


# npsp 0.7-12 (2023-06-20)

* Minor changes in `locpol()` S3 methods (for `bin.data`, `bin.den` and `svar.bin` 
  classes) so that the result extends the class of its main argument 
  (previously assumed fixed).
  
* Small changes in FORTRAN routine `besselzeros()` 
  (DFLOAT replaced by the standard DBLE; CRAN requirement). 


# npsp 0.7-11 (2023-05-01)

* Added `intermediate` argument to `np.svariso.corr()` which allows
  to return intermediate computations in `$kriging` output component 
  (these calculations can be reused, e.g. for bootstrap).

* Added `verbose` argument to `np.svariso.corr()` to avoid writing info messages 
  to the console (it can be disabled even if `plot = TRUE`).
  
* Improved documentation of `splot()`.  
  

# npsp 0.7-10 (2023-04-22)

* Changes in `h.cv.bin.data()` when `objective == "GCV"` to adapt it to the 
  heteroscedastic case.
  Warning: there may be differences with selected bandwidths in older versions.

* Changes in `simage()`, `spersp()` and `spoints()`: former argument `graphics.reset` 
  renamed as `reset`, and changed the default value to `TRUE` (to restore user's 
  graphical options).


# npsp 0.7-9 (2021-05-17)

* Added some references in the description field of 'DESCRIPTION' file.

* Avoided the use of `options(warn=-1)` in `h.cv()` methods (CRAN requirement).

* Added `on.exit(par(old.par))` in `plot.fitgeo()` and `scattersplot.default()` 
  to make sure that you the user's options are not changed (CRAN suggestion).

* Improvements in documentation (added return values, added examples in `npden()`, 
  removed `\dontrun{}` use and commented code lines in examples...).


# npsp 0.7-8 (2021-05-10)

* Renamed the admissible values of the `lost` parameter in `h.cv.svar.bin()`
  and `np.svariso.hcv()`.
  
* Changed 'NEWS.md' formatting and suppressed the default addition of 
  CRAN release dates (pkgdown).   
 
* Updated 'npsp.Rmd' vignette.


# npsp 0.7-7 (2020-12-12)

* Added `mask.window` component to `data.grid` class.

* Added new `window` parameter to `data.grid()`, `bining()`,
  `np.fitgeo.default()` and `mask()` methods.

* Minor changes in FORTRAN code (to avoid rank mismatch in 'dsytrfi.f90',
  flagged with an error in gfortran 10; CRAN policy requirement).
 
 
# npsp 0.7-6 (2019-07-10)

* Added parameter `xlim = NULL` in variogram plot methods
  (`plot.fitsvar()`, `plot.svar.bin()` and `plot.np.svar()`).

* `np.kriging()` methods now recompute residuals when 
  `any(ngrid != object$grid$n))`.


# npsp 0.7-5 (2019-06-29)

* Updated `np.fitgeo()` S3 methods.
 
* Changes in 'Makevars' to remove module files created by the Fortran compilation. 
 
* Fixed bug in `h.cv.bin.data()` (`match.arg(objective)`).
 
* 'README.md' is now generated from 'README.Rmd'.
 
* Updated roxygen documentation to avoid warnings.
  

# npsp 0.7-4 (2019-06-15)

* `kriging.np()` methods renamed as `np.kriging()` (for consistency).
 
* Minor changes in FORTRAN code (related to `error(i, label)` function, 
  to avoid LTO warnings from gcc9, which does not detect Fortran optional arguments).

* Changes in pkgdown documentation ('NEWS.Rmd', 'README.Rmd'...).
 
* Changes in 'npsp.Rmd' vignette.
 

# npsp 0.7-3 (2018-11-10)

* Added `as.data.grid()` S3 generic and `as.data.grid.SpatialGridDataFrame()` method.
 
* Added `as.bin.data.SpatialGridDataFrame()` S3 method.
 
* Added parameter `corr.svar` in `np.fitgeo()` S3 methods.

* Added parameter `asp = NA` in `spoints()` and `simage()` methods.
 

# npsp 0.7-2 (2018-02-15)

* Changed the default value of `discretize = nrow(coords) > 256` in `varcov.isotropic()`. 
 
* Added parameter `optional` in `as.data.frame.data.grid()` (for S3 compatibility).

* Changes in `svar.grid()` S3-methods:
 
    - Removed S3-method `svar.grid.fitsvar()`.
   
    - Changed the default value of `n = 256` in `svar.grid.svarmod()`.

* Changes in `sv()` S3-methods: `sv.svarmod()` and `sv.sb.iso()`.
 
* Updated roxygen2 documentation in "data.grid.R", "np.svar.R" and "svar.grid".
 

# npsp 0.7-1 (2018-02-15)

* Changes in `kriging.np`: 
 
    - S3 generic function. 

    - `kriging.np` renamed as `kriging.np.default`.

    - Added `kriging.np.np.geo` S3 method.

    - Added `ngrid` parameter.

* Added `residuals.np.geo()` and `plot.fitgeo()` S3 methods.
 
* Changes in `np.fitgeo()` S3 methods.

* Updated (roxygen2) documentation in "npsp-package.R".
 
* "npsp-plot.R" renamed as "svar.plot.R ".
 

# npsp 0.7-0 (2018-02-10)

* Added `np.geo` S3 class (nonparametric geostatistical model),  
  constructor function and methods. 

    - Added `plot()` S3 method for `np.geo` class.

    - Added `residuals()` S3 method for `np.geo` class.

* Added `np.fitgeo()` S3 generic function and methods
  (`np.fitgeo.default`, `np.fitgeo.locpol.bin` and `np.fitgeo.np.geo`).

* Fixed bug in `np.svariso.hcv()`
  (calls `h.cv.svar.bin()` instead of `h.cv.bin.data()`; thanks to Tomas Cotos-Yañez).   

   
# npsp 0.6-2 (2017-09-24)

* Added a website for the package (with pkgdown).
 
* Added 'NEWS.md' and '_pkgdown.yml'.
 
* Added some vignettes (pkgdown articles): 
  "npsp.Rmd", "precipitation.Rmd", "krigstat.Rmd", "docs/aquifer.Rmd", "docs/Introduccion.Rmd".


# npsp 0.6-1 (2017-06-16)

* Added `scattersplot()` S3 generic function (and methods).

* Added (some) support for `sp` classes.
 
    - Added `spoints()` and `scattersplot()` methods for objects of 
      class `SpatialPointsDataFrame`.
    
    - Added `as.sp()` generic function.
    
    - Added `as.sp.grid.par()` and `as.sp.data.grid()` methods.
 
* Added `precipitation` data set.
 
* Added `as.data.frame.data.grid()` S3 method.

 
# npsp 0.6-0 (2017-05-29)

* Added `kriging.np()` and `kriging.simple()` functions.
 
    - Added `kriging.simple.solve()` internal function.
 
    - Methods `as.spam()`, `chol.spam()` and `solve.spam()` imported from package `spam`.
    
    - Added `.DPOSV_R()` interface to LAPACK routine `DPOSV`.
 
* Minor changes in `plot.fitsvar()` 
  (`lwd` parameter is passed to `lines()` when `add = TRUE`).

*  Added the registration of 'native routines' (`.Fortran` calls). 
 
    - Added 'src/init.c' and `@useDynLib npsp, .registration = TRUE`.

    - Fortran routine `binning` renamed as `binning_r`.
 
* Updated 'README.md'.

   
# npsp 0.5-5 (2017-05-21)

* Added `svar.grid()` S3 class (discretized semivariogram),  
  generic function (constructor) and methods. 
   
    - Added `svar.grid.fitsvar()` and `svar.grid.svarmod()` methods.

    - Added `sv()` and `plot()` S3 methods for `svar.grid` class.

* Major changes in `varcov.isotropic()`. 

    - Returns 0 if `h < .Machine$double.eps`.
   
    - Added `discretize` parameter 
      (if `TRUE`, the default value, the variogram is previously discretized). 

* Minor changes in `covar.svarmod()`. Added `discretize` parameter 
  (if `TRUE` the variogram is discretized as a first step). 

* Added `plot()` S3 method for `svarmod` class.

   
# npsp 0.5-4 (2017-05-19)

* Minor changes in `covar.svarmod()`
  (argument `...` is passed to `sv()`).

* `fitsvar.sb.iso` returns additional components
  (`$fit$w` and `$esv`).

* Added  `.DNRM2_R()` internal function (interface to BLAS routine DNRM2).

* Minor changes in FORTRAN code
  (to avoid warnings and obsolescent features: `tql2.f`, `lp_module.f90`).

* Fixed bug in `mask.bin.data()` and `mask.locpol.bin()` 
  (when `warn = FALSE`, now it is not changed by `filter.lp`).

* Fixed bug in FORTRAN function KTW(u)
  (in the normalizing constant; thanks to Tomas Cotos-Yañez).
   

# npsp 0.5-3 (2016-09-28)

* Added 'README.md'
   
* Changes in FORTRAN code to avoid warnings compiling with -Wall -pedantic

* Changed the default value of `legend.shrink` to 1.0 in `simage.default()` and `spoints.default`.
   
* Fixed bug in `spoints.default` (when `add = TRUE`).
   

# npsp 0.5-2 (2016-07-01)

* Minor changes in `as.bin.data.data.grid()`.

* Added `as.bin.data.bin.data()`, `as.bin.den.bin.den()` and 
  `as.bin.den.data.grid()` methods.

* Removed `as.bin.den.bin.data()` 
  (`bin.den` method is now used).


# npsp 0.5-1 (2016-02-12)

* Added `h.cv.svar.bin()` (and `.wloss()` internal function).

* Minor changes in `h.cv.bin.den()`.

* Minor changes in `h.cv.bin.data()` and `hcv.data()` related to warning handling
  (the default value of `warn` parameter was also changed to TRUE).

* Changed the default value of `ncv` parameter in `h.cv.bin.data()` and
  `h.cv.bin.den()` (`ncv = 2` when `objective == `CV``).
   
* Fixed bug in `np.svariso.corr` 
  (due to extrapolations with `approx()`).

* Changes to conform to the new CRAN policy (IMPORTS).


# npsp 0.5-0 (2015-12-01)

* Major changes in `h.cv.bin.data()`.

    - Improved binning approximations of auxiliary quantities.

    - Argument `cov.bin` also admits a semivariogram model

    - Approximate computation of the covariance matrix of the binned data
      (added `.compute.masked()` internal function).

* Changes in `hcv.data()`.

    - Argument `cov` renamed as `cov.dat`, also admits a semivariogram model.

    - Improved computations when `objective = 'MASE'`.

* Minor changes in fortran code 
  (routine `lp` in `lp_module.f90` masks binning nodes with bin%w(i) < 0).

* File 'inst/CHANGES' renamed as 'ChangeLog'.


# npsp 0.4-1 (2015-07-08)

* Added `npsp.tolerance()`.

* Added `mask()` S3 generic function and methods
  (`mask.default`, `mask.bin.den`, `mask.bin.data` and `mask.locpol.bin`).

* Minor changes on `coords.data.grid()`
  (new parameter `masked`, defaults to FALSE).

* Minor changes on `predict.locpol.bin()`
  (when `!is.null(object$mask)` ...).

* Added `predict.np.den()`.


# npsp 0.4-0 (2015-03-24)
   
* Changes on `fitsvar.sb.iso()` to solve non-strictly convex quadratic programs
  (and to avoid rounding errors in `solve.QP`, 
  the constraints might not hold exactly...).
  
* Minor changes on `disc.sb()`
  (computation of the discretization nodes when `dk = 0`).
  
* Added `rule()` and `rule.binning()` default S3 method (`.rice.rule()`).

* Changed the default value of `nbin` parameter in `binning()`, `bin.den()` and
  `locpol.default()`.

* Added `rule.svar()` S3 methods.
  Changed the default value of `nlags` parameter in `svar.bin`.


# npsp 0.3-7 (2014-11-28)

* Minor changes on `h.cv.bin.data` and `hcv.data`
  (improved computations).

* Minor changes on `spoints.default`
  (`xlab` and `ylab` default values).
  
* Minor changes on `spersp.default`
  (to allow for non matrix argument `s` of appropriate length).

* Changes on `as.variogram.np.svar`
  (equivalent number of contributions).
  

# npsp 0.3-6 (2014-10-16)

* Changes on `h.cv.bin.data` and `hcv.data`
  (improved computations, `warn` parameter added, ...).

* Changed the default value of `hat.bin` argument to TRUE in `locpol.svar.bin`,
  `np.svar`, `np.svariso` and `np.svariso.corr`
  (to allow for the computation of approximated estimation variances - fitsvar.sb.iso).

* Changes on `fitsvar.sb.iso()`
  (`min.contrib`, `gstat` -> `linear` method, ...).
  

# npsp 0.3-5 (2014-08-10)

* Added `plot()` S3 methods for `svar.bin` and `np.svar` classes.

* Updated demos `aquifer` and `variogram`.
  
* `simage.default()` calls `box()` to avoid overplotting of the axis lines.

* Minor changes on FORTRAN routines `set_bin_den`, set_grid_bin ('grid_module.f90')
  and `lp` ('lp_module.f90') to avoid problems with large covariate/coordinate values.
  Warning: there may be differences with estimates computed with older versions.

* Minor changes on FORTRAN routine `predict_locpol_bin` (in 'lp_module.f90')
  to allow for extrapolations (e.g. near the grid border).
  
  
# npsp 0.3-4 (2014-05-05)

* Minor changes on `fitsvar.sb.iso`, now returns an object of class `fitsvar`
  (and inherits `sb.iso`).

* Added `plot()` S3 methods for `np.den` and `fitsvar` classes.

* Added `hot.colors()` (and `.rev.colorRampPalette()`).

* Added `cpu.time()` and `.cpu.time.ini()` (`npsp-internals`).

* Updated demo `aquifer`.
  

# npsp 0.3-3 (2014-04-04)

* Added `as.bin.data` generic function.
  
* Minor changes on `spersp.default` and `simage.default` to allow for
  non matrix arguments (of appropriate length) `z` and `s` respectively.  
  
* Minor changes on FORTRAN routine `lp` (in 'lp_module.f90') to avoid potential
  problems with memory allocation in case of error ("there is not enough data 
  in neighborhoods").
  

# npsp 0.3-2 (2014-03-24)

* Updated demo `aquifer` (illustrating the use of `np.svariso.corr()`).

* Minor changes on `bin.data`, `locpol.bin`, `svar.bin` and `np.svar` to
  allow for a dim attribute in argument `y`.  

* Minor changes on `disc.sb()` (computation of the discretization nodes when `dk = 0`).

* Fixed bug (when `degree = 0`) in fortran subroutine `lp` (`lp_module.f90`).


# npsp 0.3-1 (2014-03-16)

* Added `splot()`, `scolor()` and `jet.colors()`
  (utilities for plotting with a color scale). 

* Added `spoints()`, `spersp()` and `simage()` S3 generic functions (and methods).

* Added `persp()` and `image()` S3 methods for class `data.grid`.
  
* Changes on package demos (to not depend on package `fields` for graphic display).

* Updated documentation (`aquifer`, `earthquakes`, `locpol`, `binning`, `h.cv`, ...).
  

# npsp 0.3-0 (2014-03-03)

* Added `np.svariso.corr()` function (nonparametric bias-corrected variogram
  estimation under non-constant trend). 
  
* Renamed `svarisonp()` and `svarisohcv()` to `np.svariso()` and `np.svariso.hcv()`
  respectively.  
  
* Minor changes on `np.svariso.hcv()`.
  
* Improvements in the computation of the optimal bandwidth with the `GCV` criterion
  for dependent data (`h.cv()` and `hcv.data()` functions). 

* Updated documentation of `np.svar`, `locpol` and `binning`.


# npsp 0.2-5 (2014-01-20)

* Added `varcov()` S3 generic function (and methods).

* Changes on `covar()` (it is now an S3 generic function).

* Changes on `svarmod()` and `svarmod.sb.iso()`
  (`type` specifies a subclass of `svarmod`).
  
* Minor bug fixes in `h.cv()` and `hcv.data()` (to ensure binning/data hat
  matrix computation when needed).


# npsp 0.2-4 (2013-11-01)

* Added `aquifer` package demo.
  
* `NAMESPACE` file is now automatically generated by `roxygen2`.

* Changed dependency on package `quadprog` from `Depends` to `Imports`.

* Changes on Fortran code to conform to the Fortran 90/95 standard
  (CRAN policy requirement). The implementation of additional grid types is
  postponed until Fortran compilers used at CRAN (specially in the case of 
  Mac OS X) support the required Fortran 2003 features (mainly type-bound 
  procedures).
  
* Added explicit dependencies to src/Makevars to allow parallel make.


# npsp 0.2-3 (2013-10-19)

* Changed the default value of `hat.bin` argument to FALSE in `locpol` and
  `predict.locpol.bin`.
  
* Added a default value for `maxlag` argument in `svar.bin`, `svariso`,
  `np.svar`, `svarisonp` and `svarisohcv`.
  
* Changed the default value of `nx` argument in `fitsvar.sb.iso`
  to avoid "Error in solve.QP(Dmat, dvec, Amat, bvec) :
  matrix D in quadratic function is not positive definite!".
  
* Updated documentation.
   

# npsp 0.2-2 (2013-10-15)

* Added `bin.den` S3 class and methods and `as.bin.den` generic function.

* Added `np.den` S3 class and generic function.
  
* Added `locpol.bin.den` (alias of `np.den.bin.den`) and `h.cv.bin.den` methods.
  
* Added `earthquakes` and `aquifer` data sets.
  
* Added a `dimnames` argument to `grid.par` function
  (constructor of the class of the same name).
  
* Minor changes on `binning()` and `interp.data.grid()`.


# npsp 0.2-1 (2013-09-12)

* Completed the package documentation (using `roxygen2`).

* Changes on `svarmod`, added `svarmodels` (`svarmod.R`).

* Changes on `as.variomodel` and `as.vgm` (`npsp-geoR.R` and `npsp-gstat.R`,
  interoperability with geoR and gstat, respectively).

* Some minor changes to pass (for the first time) R CMD check without notes or warnings
  (`.onLoad` -> `.onAttach`, ...).


# npsp 0.2-0 (2013-06-25)

* Major changes in R functions `locpol.default`, `locpol.bin.data` and `svarisonp`:
  
    - Added a new option to set the degree of the local polynomial used.
   
    - Added a new option to compute (partial) derivative estimates.
   
    - Added an option to enable/disable binning hat matrix computation.   

* Added an option to enable/disable data hat matrix computation in `predict.locpol.bin`.

* Major changes in fortran code (`lp_module.f90`, `svar_module.f90`, `linreg_module.f90`):
   
    - Weighted linear regression allows for rank-deficient matrices 
      (`DGELSYR` fortran routine, a modification of LAPACK `DGELSY`). 

    - New functionalities (degree, derivatives...)
    
    - Changes on Fortran-R interfaces.
   

# npsp 0.1-8 (2013-04-01)

* Added `svarmod` and `sb.iso` (extends `svarmod`)
  S3 classes and methods.  

* Added `svar()` S3 generic and `covar()` functions.

* Added `fitsvar.sb.iso()`, `kappa.sb()` and `disc.sb()` functions.

* Added `as.vgm()` S3 generic (interoperability with gstat).
  

# npsp 0.1-7 (2013-02-19)

* Added `h.cv.bin()` and `h.cv()` functions (EXPERIMENTAL).
  
* Changes on `locpolhcv()` and `svarisohcv()`.

* Changes on `predict.locpol.bin()` (new fortran code).
  
  
# npsp 0.1-6 (2012-10-24)

* Added `predict.locpol.bin()` and S3 generic function `interp()` (`interp.grid.par()`
  and `interp.data.grid()` methods).  

* Changes on R-Fortran interfaces (parameters for type(grid_bin) :: bin,
  `lp_raw` replaces `est_bin`, ...).
  

# npsp 0.1-5 (2012-09-19)

* Added `hopt.cv()`, `locpolhcv()` and `svarisohcv()` functions (EXPERIMENTAL).

* Changes on FORTRAN and R code to handle missing values (EXPERIMENTAL)
  (NAs on input & output).
  
* `binning()` is again a standard function (interface to the fortran routine `binning`).
   

# npsp 0.1-4 (2012-09-01)

* Added `svar.bin` (extends `bin.data`), and `np.svar` (extends `svar.bin`)
  S3 classes and methods.  

* Changes on `svariso()` (returns an object of class `svar.bin`, ...).

* Added `svarisonp()` and `as.variogram()` S3 generic functions.
  `svarisonp.default` replaces `svarisonp` (interface to the fortran routine `svar_iso_np`).  

* Added `variogram` package demo.


# npsp 0.1-3 (2012-08-25)

* Added `grid.par`, `data.grid`, `bin.data` and `locpol.bin`
  S3 classes and methods.  

* Added `coords()` and `coordvalues()` S3 generic functions.

* Added `binning()` and `locpol()` S3 generic functions.
  `binning.default` replaces `set_bin` (interface to the fortran routine `set_bin`).  
  `locpol.default` replaces `locpolbin` (interface to the fortran routine `est_bin`).  

* Added `locpol.bin.data()` (interface to the fortran routine `lp_bin`).

* Added `datagrid`, `binning` and `locpol` package demos.

  
# npsp 0.1-2 (2012-08-20)

* Added `svarisonp()` (interface to the fortran routine `svar_iso_np`).

* Added `svariso()` (interface to the fortran routine `svar_iso_bin`).

* Added `as.variogram.svariso()` (interoperability with geoR).
  
  
# npsp 0.1-1 (2012-06-01)

* Added a package demo.
  
* Added `npsp-geoR.R` (interoperability with geoR).

* Some minor changes (output `set_bin` and `locpolbin`, `onLoad`...).


# npsp 0.1-0 (2012-04-17) 

* Initial version in package form.
  