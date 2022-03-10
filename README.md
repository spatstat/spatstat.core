# spatstat.core

## Statistical analysis and modelling of spatial data for the spatstat family

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/spatstat.core)](http://cran.r-project.org/web/packages/spatstat.core) 
[![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.core)](https://github.com/spatstat/spatstat.core)

The original `spatstat` package has been split into
several sub-packages (See [spatstat/spatstat](https://github.com/spatstat/spatstat))

This package `spatstat.core` is one of the
sub-packages. It contains all the main user-level functions that perform
statistical analysis and modelling of spatial data,
with the exception of data on linear networks.

### Overview 

`spatstat.core` supports

- exploratory analysis (quadrat counting, kernel smoothing, K-function, pair correlation function)
- nonparametric modelling (resource selection function, prospectivity)
- parametric modelling (fitting models to point pattern data, model selection, model prediction)
- formal inference (hypothesis tests, confidence intervals)
- informal validation (model diagnostics)

### Detailed contents

- kernel estimation of intensity of a point pattern
- kernel smoothing of mark values attached to point locations
- kernel estimation of relative risk
- kernel smoothing of a line segment pattern
- image blurring
- Choi-Hall data sharpening of point locations
- bandwidth selection
- transects of an image along a line or curve
- exploratory tools (Clark-Evans index, Fry plot, Morisita plot, scan statistic)
- cluster detection (Allard-Fraley cluster set, Byers-Raftery cleaning)
- nonparametric estimation of intensity as a function of a covariate
- ROC curve, AUC
- summary functions (K-function, pair correlation function,
empty space function, nearest neighbour distance function, J-function, etc)
and multi-type versions of these functions
- mark correlation function, mark independence diagnostoc
- local summary functions (LISA)
- simulation envelopes of summary functions
- manipulation of summary functions (plot, evaluate, differentiate, smooth etc)
- spatial bootstrap
- simulation of fitted point process models
- programming tools
- fitting Poisson point process models to point pattern data (`ppm`)
- fitting spatial logistic regression models to point pattern data (`slrm`)
- fitting Cox point process models to point pattern data (`kppm`)
- fitting Neyman-Scott cluster process models to point pattern data (`kppm`)
- fitting Gibbs point process models to point pattern data (`ppm`)
- class support for fitted models (update, summary, predict, plot, coef, vcov)
- minimum contrast estimation
- hypothesis tests (quadrat test, Clark-Evans test, Berman test, Diggle-Cressie-Loosmore-Ford test, scan test, studentised permutation test, segregation test, ANOVA tests of fitted models, Dao-Genton test, balanced independent two-stage test)
- model validation (leverage, influence, dfbeta, partial residuals, added
variable plot, diagnostic plots, pseudoscore residual plots, model compensators,
Q-Q plots)

