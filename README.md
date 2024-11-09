
<!-- 
README.md is generated from README.Rmd. 
Please edit that file 
-->

# npsp: Nonparametric spatial (geo)statistics

### Version 0.7.14

This package implements nonparametric methods for inference on
multidimensional spatial (or spatio-temporal) processes, which may be
(especially) useful in (automatic) geostatistical modeling and
interpolation.

## Main functions

Nonparametric methods for inference on both spatial trend and variogram
functions:

- `np.fitgeo()` (automatically) fits an isotropic nonparametric
  geostatistical model by estimating the trend and the variogram (using
  a bias-corrected estimator) iteratively (by calling `h.cv()`,
  `locpol()`, `np.svariso.corr()` and `fitsvar.sb.iso()` at each
  iteration).

- `locpol()`, `np.den()` and `np.svar()` use local polynomial kernel
  smoothing to compute nonparametric estimates of a multidimensional
  regression function (e.g. a spatial trend), a probability density
  function or a semivariogram (or their first derivatives),
  respectively. Estimates of these functions can be constructed for any
  dimension (depending on the amount of available memory).

- `np.svariso.corr()` computes a bias-corrected nonparametric
  semivariogram estimate using an iterative algorithm similar to that
  described in Fernandez-Casal and Francisco-Fernandez (2014). This
  procedure tries to correct the bias due to the direct use of residuals
  (obtained, in this case, from a nonparametric estimation of the trend
  function) in semivariogram estimation.

- `fitsvar.sb.iso()` fits a ‘nonparametric’ isotropic Shapiro-Botha
  variogram model by WLS. Currently, only isotropic semivariogram
  estimation is supported.

Nonparametric residual kriging (sometimes called external drift
kriging):

- `np.kriging()` computes residual kriging predictions  
  (and the corresponding simple kriging standard errors).

- `kriging.simple()` computes simple kriging predictions and standard
  errors.

- Currently, only global (residual) simple kriging is implemented.  
  Users are encouraged to use `gstat::krige()` (or `gstat::krige.cv()`)
  together with `as.vgm()` for local kriging.

## Other functions

Among the other functions intended for direct access by the user, the
following (methods for multidimensional linear binning, local polynomial
kernel regression, density or variogram estimation) could be emphasized:
`binning()`, `bin.den()`, `svar.bin()`, `h.cv()` and `interp()`. There
are functions for plotting data joint with a legend representing a
continuous color scale (based on `fields::image.plot()`):

- `splot()` allows to combine a standard R plot with a legend.

- `spoints()`, `simage()` and `spersp()` draw the corresponding
  high-level plot with a legend strip for the color scale.

There are also some functions which can be used to interact with other
packages. For instance, `as.variogram()` (geoR) or `as.vgm()` (gstat).

See the
[Reference](https://rubenfcasal.github.io/npsp/reference/index.html) for
the complete list of functions.

## Installation

`npsp` is available from CRAN, but you can install the development
version from github with:

``` r
# install.packages("devtools")
devtools::install_github("rubenfcasal/npsp")
```

Note also that, as this package requires compilation, Windows users need
to have previously installed the appropriate version of
[Rtools](https://cran.r-project.org/bin/windows/Rtools/), and OS X users
need to have installed
[Xcode](https://apps.apple.com/us/app/xcode/id497799835).

Alternatively, Windows users may install the corresponding
*npsp_X.Y.Z.zip* file in the [releases
section](https://github.com/rubenfcasal/npsp/releases/latest) of the
github repository.

For R versions 4.4.x under Windows:

``` r
install.packages('https://github.com/rubenfcasal/npsp/releases/download/v0.7-14/npsp_0.7-14.zip',
                 repos = NULL)
```

## Author

[Ruben Fernandez-Casal](https://rubenfcasal.github.io) (Dep.
Mathematics, University of A Coruña, Spain). Please send comments, error
reports or suggestions to <rubenfcasal@gmail.com>.

## Acknowledgments

Important suggestions and contributions to some techniques included here
were made by Sergio Castillo-Páez (Universidad de las Fuerzas Armadas
ESPE, Ecuador) and Tomas Cotos-Yañez (Dep. Statistics, University of
Vigo, Spain).

This research has been supported by MINECO grant MTM2017-82724-R, and by
the Xunta de Galicia (Grupos de Referencia Competitiva ED431C-2020-14
and Centro de Investigación del Sistema universitario de Galicia ED431G
2019/01), all of them through the ERDF.

## References

- Castillo-Páez S., Fernández-Casal R. and García-Soidán P. (2019). [A
  nonparametric bootstrap method for spatial
  data](https://doi.org/10.1016/j.csda.2019.01.017), **137**, *Comput.
  Stat. Data Anal.*, 1-15.

- Fernández-Casal, R., Castillo-Páez, S. and Francisco-Fernandez, M.
  (2024). [Nonparametric Conditional Risk Mapping Under
  Heteroscedasticity](https://doi.org/10.1007/s13253-023-00555-0),
  *JABES*, **29**, 56-72.

- Fernández-Casal R., Castillo-Páez S. and Francisco-Fernández M.
  (2018). [Nonparametric geostatistical risk
  mapping](https://doi.org/10.1007/s00477-017-1407-y), *Stoch. Environ.
  Res. Ris. Assess.*, **32**, 675-684.

- Fernández-Casal R., Castillo-Páez S. and García-Soidán P. (2017).
  [Nonparametric estimation of the small-scale variability of
  heteroscedastic spatial
  processes](https://doi.org/10.1016/j.spasta.2017.04.001), *Spa. Sta.*,
  **22**, 358-370.

- Fernandez-Casal R. and Francisco-Fernandez M. (2014). [Nonparametric
  bias-corrected variogram estimation under non-constant
  trend](https://doi.org/10.1007/s00477-013-0817-8), *Stoch. Environ.
  Res. Ris. Assess.*, **28**, 1247-1259.

- Fernandez-Casal R., Gonzalez-Manteiga W. and Febrero-Bande M. (2003).
  [Flexible Spatio-Temporal Stationary Variogram
  Models](https://doi.org/10.1023/A:1023204525046), *Statistics and
  Computing*, **13**, 127-136.

- Rupert D. and Wand M.P. (1994). Multivariate locally weighted least
  squares regression. *The Annals of Statistics*, **22**, 1346-1370.

- Shapiro A. and Botha J.D. (1991). Variogram fitting with a general
  class of conditionally non-negative definite functions. *Computational
  Statistics and Data Analysis*, **11**, 87-96.

- Wand M.P. (1994). Fast Computation of Multivariate Kernel Estimators.
  *Journal of Computational and Graphical Statistics*, **3**, 433-445.

- Wand M.P. and Jones M.C. (1995). *Kernel Smoothing*. Chapman and Hall,
  London.
