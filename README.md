npsp: Nonparametric spatial (geo)statistics
============================================

This package implements nonparametric methods which may be useful in geostatistical practice.

Main functions
--------------

`locpol`, `np.den` and `np.svar` use local polynomial kernel smoothing to compute nonparametric estimates of a multidimensional regression function (e.g. a spatial trend), a probability density function or a semivariogram (or their first derivatives), respectively. Estimates of these functions can be constructed for any dimension (the amount of available memory is the only limitation). To speed up computations, linear binning is used to discretize the data. A full bandwidth matrix and a multiplicative triweight kernel is used to compute the weights. Main calculations are performed in FORTRAN using the LAPACK library.

`np.svariso.corr` computes a bias-corrected nonparametric semivariogram estimate using an iterative algorithm similar to that described in Fernandez-Casal and Francisco-Fernandez (2014). This procedure tries to correct the bias due to the direct use of residuals, obtained from a nonparametric estimation of the trend function, in semivariogram estimation.

`fitsvar.sb.iso` fits a ‘nonparametric’ isotropic Shapiro-Botha variogram model by WLS. Currently, only isotropic semivariogram estimation is supported.

There are also functions for plotting data joint with a legend representing a continuous color scale. `splot` allows to combine a standard R plot with a legend. `spoints`, `simage` and `spersp` draw the corresponding high-level plot with a legend strip for the color scale. These functions are based on `image.plot` of package `fields`.

Among the other functions intended for direct access by the user, the following could be emphasized: `binning`, `bin.den`, `svar.bin`, `h.cv` and `interp`. There are also some functions which can be used to interact with other packages. For instance, `as.variogram` (geoR) or `as.vgm` (gstat).

Kriging is not yet implemented in this package. Users are encouraged to use `krige` (or `krige.cv`) utilities in gstat package together with `as.vgm`.


Author
------

Ruben Fernandez-Casal (Dep. Mathematics, University of A Coru\~na, Spain). Please send comments, error reports or suggestions to rubenfcasal@gmail.com.

References
----------

* Fernandez-Casal R. and Francisco-Fernandez M. (2014) Nonparametric bias-corrected variogram estimation under non-constant trend, *Stoch. Environ. Res. Ris. Assess*, **28**, 1247-1259.

* Fernandez-Casal R., Gonzalez-Manteiga W. and Febrero-Bande M. (2003) Flexible Spatio-Temporal Stationary Variogram Models, *Statistics and Computing*, **13**, 127-136.

* Rupert D. and Wand M.P. (1994) Multivariate locally weighted least squares regression. *The Annals of Statistics*, **22**, 1346-1370.

* Shapiro A. and Botha J.D. (1991) Variogram fitting with a general class of conditionally non-negative definite functions. *Computational Statistics and Data Analysis*, **11**, 87-96.

* Wand M.P. (1994) Fast Computation of Multivariate Kernel Estimators. *Journal of Computational and Graphical Statistics*, **3**, 433-445.

* Wand M.P. and Jones M.C. (1995) *Kernel Smoothing*. Chapman and Hall, London.