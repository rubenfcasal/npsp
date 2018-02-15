npsp: Nonparametric spatial (geo)statistics
============================================

This package implements nonparametric methods 
for inference on multidimensional spatial (or spatio-temporal) processes,
which may be (especially) useful in (automatic) geostatistical modeling and interpolation.

Main functions
--------------

Nonparametric methods for inference on both spatial trend 
and variogram functions:

*  `np.fitgeo` (automatically) fits an isotropic nonparametric geostatistical model 
   by estimating the trend and the variogram (using a bias-corrected estimator) iteratively 
   (by calling `h.cv`, `locpol`, `np.svariso.corr` and `fitsvar.sb.iso` at each iteration)

*  `locpol`, `np.den` and `np.svar` use local polynomial kernel smoothing to compute
   nonparametric estimates of a multidimensional regression function (e.g. a spatial trend),
   a probability density function or a semivariogram (or their first derivatives), respectively. 
   Estimates of these functions can be constructed for any dimension 
   (depending on the amount of available memory). 
   To speed up computations, linear binning is used to discretize the data. 
   A full bandwidth matrix and a multiplicative triweight kernel is used to compute the weights. 
   Main calculations are performed in FORTRAN using the LAPACK library.

*  `np.svariso.corr` computes a bias-corrected nonparametric semivariogram estimate 
   using an iterative algorithm similar to that described in 
   Fernandez-Casal and Francisco-Fernandez (2014). 
   This procedure tries to correct the bias due to the direct use of residuals 
   (obtained, in this case, from a nonparametric estimation of the trend function) 
   in semivariogram estimation.

*  `fitsvar.sb.iso` fits a ‘nonparametric’ isotropic Shapiro-Botha variogram model by WLS. 
   Currently, only isotropic semivariogram estimation is supported.


Nonparametric residual kriging (sometimes called external drift kriging):

*  `kriging.np` computes residual kriging predictions  
   (and the corresponding simple kriging standard errors).

*  `kriging.simple` computes simple kriging predictions and standard errors.

*  Currently, only global (residual) simple kriging is implemented.  
   Users are encouraged to use `krige` (or `krige.cv`) utilities in gstat package
   together with `as.vgm` for local kriging.


Other functions
---------------

Among the other functions intended for direct access by the user, the following 
(methods for multidimensional linear binning, local polynomial kernel regression, 
density or variogram estimation) could be emphasized: 
`binning`, `bin.den`, `svar.bin`, `h.cv` and `interp`. 

There are functions for plotting data joint with a legend representing a continuous color scale
(based on `image.plot` of package `fields`):

*  `splot` allows to combine a standard R plot with a legend. 

*  `spoints`, `simage` and `spersp` draw the corresponding high-level plot 
   with a legend strip for the color scale.

There are also some functions which can be used to interact with other packages. For instance, `as.variogram` (geoR) or `as.vgm` (gstat).


News
----

 * Version 0.7-1 (2018-02-15)

    - Changes in `kriging.np` 
      (S3 generic, `kriging.np.default`, `kriging.np.np.geo`, `ngrid` parameter).


 * Version 0.7-0 (2018-02-10)

    - Added `np.geo` S3 class (nonparametric geostatistical model),  
      constructor function and methods (`plot`, `residuals`).

    - Added `np.fitgeo()` S3 generic function and methods
     (`np.fitgeo.default`, `np.fitgeo.locpol.bin` and `np.fitgeo.np.geo`).

    - Fixed bug in `np.svariso.hcv()` (thanks to Tomas Cotos-Yañez).
    

 * Version 0.6-2 (2017-09-24)
 
    - Added a website for the package (with pkgdown).
    
    - Added 'NEWS.md' and 'index.Rmd'.

    - Added some vignettes (pkgdown articles): 
      "npsp.Rmd", "precipitation.Rmd", "krigstat.Rmd", 
      "docs/aquifer.Rmd", "docs/Introduccion.Rmd".
    
 
 * Version 0.6-1 (2017-06-16)
 
    - Added `scattersplot()` S3 generic function (and methods).

    - Added (some) support for `sp` classes.
 
    - Added `precipitation` data set.
 
    - Added `as.data.frame.data.grid` 
    
 
 * Version 0.6-0 (2017-05-29)

    - Added `kriging.np()` and `kriging.simple()` functions.
 
    - Minor changes in `plot.fitsvar()`

    - Added the registration of 'native routines' (`.Fortran` calls). 
 
    - Updated 'README.md'.


 * Version 0.5-5 (2017-05-21)

     - Added `svar.grid()` S3 class (discretized semivariogram) and methods. 
   
     - Major changes in `varcov.isotropic()` and `sv.sb.iso()`: 
       added `discretize` parameter (if `TRUE`the variogram is previously discretized).  
   
     - Added `plot()` S3 method for `svarmod` class.


 * Older versions (or more info): [npsp ChangeLog](./ChangeLog)


Author
------

Ruben Fernandez-Casal (Dep. Mathematics, University of A Coruña, Spain). 
Please send comments, error reports or suggestions to [rubenfcasal@gmail.com](mailto:rubenfcasal@gmail.com).


Acknowledgments
---------------

Important suggestions and contributions to some techniques included here were
made by Sergio Castillo-Páez (Universidad de las Fuerzas Armadas ESPE, Ecuador) 
and Tomas Cotos-Yañez (Dep. Statistics, University of Vigo, Spain).


References
----------

* Fernández-Casal R., Castillo-Páez S. and Francisco-Fernández M. (2018), Nonparametric geostatistical risk mapping, *Stoch. Environ. Res. Ris. Assess.*, [DOI](https://doi.org/10.1007/s00477-017-1407-y).

* Fernández-Casal R., Castillo-Páez S. and García-Soidán P. (2017), Nonparametric estimation of the small-scale variability of heteroscedastic spatial processes, *Spa. Sta.*, **22**, 358-370, [DOI](https://doi.org/10.1016/j.spasta.2017.04.001).

* Fernandez-Casal R. and Francisco-Fernandez M. (2014) Nonparametric bias-corrected variogram estimation under non-constant trend, *Stoch. Environ. Res. Ris. Assess.*, **28**, 1247-1259.

* Fernandez-Casal R., Gonzalez-Manteiga W. and Febrero-Bande M. (2003) Flexible Spatio-Temporal Stationary Variogram Models, *Statistics and Computing*, **13**, 127-136.

* Rupert D. and Wand M.P. (1994) Multivariate locally weighted least squares regression. *The Annals of Statistics*, **22**, 1346-1370.

* Shapiro A. and Botha J.D. (1991) Variogram fitting with a general class of conditionally non-negative definite functions. *Computational Statistics and Data Analysis*, **11**, 87-96.

* Wand M.P. (1994) Fast Computation of Multivariate Kernel Estimators. *Journal of Computational and Graphical Statistics*, **3**, 433-445.

* Wand M.P. and Jones M.C. (1995) *Kernel Smoothing*. Chapman and Hall, London.