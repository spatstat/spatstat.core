\name{spatstat.core-package}
\alias{spatstat.core-package} 
\alias{spatstat.core} 
\docType{package}
\title{The spatstat.core Package}
\description{
  The \pkg{spatstat.core} package
  belongs to the \pkg{spatstat} family of packages.
  It contains the core functionality
  for statistical analysis and modelling of spatial data.
}
\details{
  \pkg{spatstat} is
  a family of \R packages
  for the statistical analysis of spatial data.
  Its main focus is the analysis of
  spatial patterns of points in two-dimensional space.

  The original \pkg{spatstat} package
  has now been split into several
  sub-packages.


  This sub-package \pkg{spatstat.core} contains all the main
  user-level functions
  that perform statistical analysis and modelling of spatial data.

  (The main exception is that functions for linear networks
  are in the separate sub-package \pkg{spatstat.linnet}.)
}
\section{Structure of the spatstat family}{  
  The orginal \pkg{spatstat} package grew to be very large.
  It has now been divided into several \bold{sub-packages}:
  \itemize{
    \item \pkg{spatstat.utils} containing basic utilities
    \item \pkg{spatstat.sparse} containing linear algebra utilities
    \item \pkg{spatstat.data} containing datasets
    \item \pkg{spatstat.geom} containing geometrical objects
    and geometrical operations
    \item \pkg{spatstat.core} containing the main functionality
    for statistical analysis and modelling of spatial data
    \item \pkg{spatstat.linnet} containing functions for
    spatial data on a linear network
    \item \pkg{spatstat}, which simply loads the other sub-packages
    listed above, and provides documentation.
  }
  When you install \pkg{spatstat}, these sub-packages are also
  installed. Then if you load the \pkg{spatstat} package by typing
  \code{library(spatstat)}, the other sub-packages listed above will
  automatically be loaded or imported.

  For an overview of all the functions available in
  the sub-packages of \pkg{spatstat},
  see the help file for \code{"spatstat-package"}
  in the \pkg{spatstat} package.

  Additionally there are several \bold{extension packages:}
  \itemize{
    \item \pkg{spatstat.gui} for interactive graphics
    \item \pkg{spatstat.local} for local likelihood
    (including geographically weighted regression)
    \item \pkg{spatstat.Knet} for additional, computationally efficient code
    for linear networks
    \item \pkg{spatstat.sphere} (under development) for spatial data
    on a sphere, including spatial data on the earth's surface
  }
  The extension packages must be installed separately
  and loaded explicitly if needed. They also have separate documentation.
}
\section{Overview of Functionality in \pkg{spatstat.core}}{
  The \pkg{spatstat} family of packages
  is designed to support a complete statistical analysis
  of spatial data. It supports
  \itemize{
    \item creation, manipulation and plotting of point patterns;
    \item exploratory data analysis;
    \item spatial random sampling;
    \item simulation of point process models;
    \item parametric model-fitting;
    \item non-parametric smoothing and regression;
    \item formal inference (hypothesis tests, confidence intervals);
    \item model diagnostics.
  }
  For an overview, see the help file for \code{"spatstat-package"}
  in the \pkg{spatstat} package. 

  Following is a list of the functionality provided in the
  \pkg{spatstat.core} package only.
  
  \bold{To simulate a random point pattern:}
  
  \tabular{ll}{
	\code{\link[spatstat.core]{runifpoint}} \tab
	generate \eqn{n} independent uniform random points \cr
	\code{\link[spatstat.core]{rpoint}} \tab
	generate \eqn{n} independent random points \cr
	\code{\link[spatstat.core]{rmpoint}} \tab
	generate \eqn{n} independent multitype random points \cr
	\code{\link[spatstat.core]{rpoispp}} \tab
	simulate the (in)homogeneous Poisson point process \cr
	\code{\link[spatstat.core]{rmpoispp}} \tab
	simulate the (in)homogeneous multitype Poisson point process \cr
	\code{\link[spatstat.core]{runifdisc}} \tab
	generate \eqn{n} independent uniform random points in disc\cr
	\code{\link[spatstat.core]{rstrat}} \tab
	stratified random sample of points \cr
	\code{\link[spatstat.core]{rMaternI}}  \tab
	simulate the \Matern Model I inhibition process\cr
	\code{\link[spatstat.core]{rMaternII}} \tab
	simulate the \Matern Model II inhibition process\cr
	\code{\link[spatstat.core]{rSSI}} \tab
	simulate Simple Sequential Inhibition process\cr
	\code{\link[spatstat.core]{rStrauss}} \tab
	simulate Strauss process (perfect simulation)\cr
	\code{\link[spatstat.core]{rHardcore}} \tab
	simulate Hard Core process (perfect simulation)\cr
	\code{\link[spatstat.core]{rStraussHard}} \tab
	simulate Strauss-hard core process (perfect simulation)\cr
	\code{\link[spatstat.core]{rDiggleGratton}} \tab
	simulate Diggle-Gratton process (perfect simulation)\cr
	\code{\link[spatstat.core]{rDGS}} \tab
	simulate Diggle-Gates-Stibbard process (perfect simulation)\cr
	\code{\link[spatstat.core]{rPenttinen}} \tab
	simulate Penttinen process (perfect simulation)\cr
	\code{\link[spatstat.core]{rNeymanScott}} \tab
	simulate a general Neyman-Scott process\cr
	\code{\link[spatstat.core]{rPoissonCluster}} \tab
	simulate a general Poisson cluster process\cr
	\code{\link[spatstat.core]{rMatClust}} \tab
	simulate the \Matern Cluster process\cr
	\code{\link[spatstat.core]{rThomas}} \tab
	simulate the Thomas process  \cr
	\code{\link[spatstat.core]{rGaussPoisson}}  \tab
	simulate the Gauss-Poisson cluster process\cr
	\code{\link[spatstat.core]{rCauchy}} \tab
	simulate Neyman-Scott Cauchy cluster process \cr
	\code{\link[spatstat.core]{rVarGamma}} \tab
	simulate Neyman-Scott Variance Gamma cluster process \cr
	\code{\link[spatstat.core]{rthin}} \tab  random thinning  \cr
	\code{\link[spatstat.core]{rcell}} \tab
	simulate the Baddeley-Silverman cell process  \cr
	\code{\link[spatstat.core]{rmh}} \tab
	simulate Gibbs point process using Metropolis-Hastings \cr
	\code{\link[spatstat.core]{simulate.ppm}} \tab
	simulate Gibbs point process using Metropolis-Hastings \cr
	\code{\link[spatstat.core]{runifpointOnLines}} \tab
	generate \eqn{n} random points along specified line segments \cr
	\code{\link[spatstat.core]{rpoisppOnLines}} \tab
	generate Poisson random points along specified line segments 
      }

      \bold{To randomly change an existing point pattern:}
      
      \tabular{ll}{
	\code{\link[spatstat.core]{rshift}} \tab random shifting of points \cr
	\code{\link[spatstat.core]{rthin}} \tab  random thinning \cr
	\code{\link[spatstat.core]{rlabel}} \tab random (re)labelling of a multitype
	point pattern \cr
        \code{\link[spatstat.core]{quadratresample}} \tab block resampling 
      }

      \bold{To interrogate a point pattern:}

      \tabular{ll}{
	\code{\link[spatstat.core]{density.ppp}} \tab
	kernel estimation of point pattern intensity\cr
	\code{\link[spatstat.core]{densityHeat.ppp}} \tab
	diffusion kernel estimation of point pattern intensity\cr
	\code{\link[spatstat.core]{Smooth.ppp}} \tab
	kernel smoothing of marks of point pattern\cr
	\code{\link[spatstat.core]{sharpen.ppp}} \tab
	data sharpening\cr
      }

      \bold{Manipulation of pixel images:}

      An object of class \code{"im"} represents a pixel image.

      \tabular{ll}{
	\code{\link[spatstat.core]{blur}} \tab apply Gaussian blur to image\cr
	\code{\link[spatstat.core]{Smooth.im}} \tab apply Gaussian blur to image\cr
	\code{\link[spatstat.core]{transect.im}} \tab line transect of image \cr
	\code{\link[spatstat.geom]{pixelcentres}} \tab extract centres of pixels \cr
	\code{\link[spatstat.core]{rnoise}} \tab random pixel noise
     }

    \bold{Line segment patterns}

    An object of class \code{"psp"} represents a pattern of straight line
    segments.

    \tabular{ll}{
      \code{\link[spatstat.core]{density.psp}} \tab kernel smoothing of line segments\cr
      \code{\link[spatstat.core]{rpoisline}} \tab generate a realisation of the
      Poisson line process inside a window
    }

    \bold{Tessellations}

    An object of class \code{"tess"} represents a tessellation.

    \tabular{ll}{
      \code{\link[spatstat.core]{rpoislinetess}} \tab generate tessellation using Poisson line
      process 
    }

    \bold{Three-dimensional point patterns}

    An object of class \code{"pp3"} represents a three-dimensional
    point pattern in a rectangular box. The box is represented by
    an object of class \code{"box3"}.

    \tabular{ll}{
      \code{\link[spatstat.core]{runifpoint3}} \tab generate uniform random points in 3-D \cr
      \code{\link[spatstat.core]{rpoispp3}} \tab generate Poisson random points in 3-D \cr
      \code{\link[spatstat.core]{envelope.pp3}} \tab generate simulation envelopes for
      3-D pattern \cr
    }

    \bold{Multi-dimensional space-time point patterns}

    An object of class \code{"ppx"} represents a 
    point pattern in multi-dimensional space and/or time.

    \tabular{ll}{
      \code{\link[spatstat.core]{runifpointx}} \tab generate uniform random points \cr
      \code{\link[spatstat.core]{rpoisppx}} \tab generate Poisson random points 
    }
    

  \bold{Classical exploratory tools:}
  \tabular{ll}{
    \code{\link[spatstat.core]{clarkevans}} \tab Clark and Evans aggregation index \cr
    \code{\link[spatstat.core]{fryplot}} \tab Fry plot \cr
    \code{\link[spatstat.core]{miplot}} \tab Morisita Index plot
  }

  \bold{Smoothing:}
  \tabular{ll}{
    \code{\link[spatstat.core]{density.ppp}} \tab kernel smoothed density/intensity\cr
    \code{\link[spatstat.core]{relrisk}} \tab kernel estimate of relative risk\cr
    \code{\link[spatstat.core]{Smooth.ppp}} \tab spatial interpolation of marks  \cr
    \code{\link[spatstat.core]{bw.diggle}} \tab cross-validated bandwidth selection
                               for \code{\link[spatstat.core]{density.ppp}}\cr
    \code{\link[spatstat.core]{bw.ppl}} \tab likelihood cross-validated bandwidth selection
                               for \code{\link[spatstat.core]{density.ppp}}\cr
    \code{\link[spatstat.core]{bw.CvL}} \tab Cronie-Van Lieshout bandwidth selection
                                for density estimation\cr
    \code{\link[spatstat.core]{bw.scott}} \tab Scott's rule of thumb 
                                for density estimation\cr
    \code{\link[spatstat.core]{bw.abram}} \tab Abramson's rule for adaptive bandwidths\cr
    \code{\link[spatstat.core]{bw.relrisk}} \tab cross-validated bandwidth selection
                               for \code{\link[spatstat.core]{relrisk}} \cr
    \code{\link[spatstat.core]{bw.smoothppp}} \tab cross-validated bandwidth selection
                               for \code{\link[spatstat.core]{Smooth.ppp}} \cr
    \code{\link[spatstat.core]{bw.frac}} \tab bandwidth selection using window geometry\cr
    \code{\link[spatstat.core]{bw.stoyan}} \tab Stoyan's rule of thumb for bandwidth
                               for \code{\link[spatstat.core]{pcf}}
  }

  \bold{Modern exploratory tools:}
  \tabular{ll}{
    \code{\link[spatstat.core]{clusterset}} \tab Allard-Fraley feature detection  \cr
    \code{\link[spatstat.core]{nnclean}} \tab Byers-Raftery feature detection  \cr
    \code{\link[spatstat.core]{sharpen.ppp}} \tab Choi-Hall data sharpening \cr
    \code{\link[spatstat.core]{rhohat}} \tab Kernel estimate of covariate effect\cr
    \code{\link[spatstat.core]{rho2hat}} \tab Kernel estimate of effect of two covariates\cr
    \code{\link[spatstat.core]{spatialcdf}} \tab Spatial cumulative distribution function\cr
    \code{\link[spatstat.core]{roc}} \tab Receiver operating characteristic curve
  }

  \bold{Summary statistics for a point pattern:}

  \tabular{ll}{
    \code{\link[spatstat.core]{Fest}} \tab empty space function \eqn{F} \cr
    \code{\link[spatstat.core]{Gest}} \tab nearest neighbour distribution function \eqn{G} \cr
    \code{\link[spatstat.core]{Jest}} \tab \eqn{J}-function \eqn{J = (1-G)/(1-F)} \cr
    \code{\link[spatstat.core]{Kest}} \tab Ripley's \eqn{K}-function\cr
    \code{\link[spatstat.core]{Lest}} \tab Besag \eqn{L}-function\cr
    \code{\link[spatstat.core]{Tstat}} \tab Third order \eqn{T}-function \cr
    \code{\link[spatstat.core]{allstats}} \tab all four functions \eqn{F}, \eqn{G}, \eqn{J}, \eqn{K} \cr
    \code{\link[spatstat.core]{pcf}} \tab 	pair correlation function \cr
    \code{\link[spatstat.core]{Kinhom}} \tab \eqn{K} for inhomogeneous point patterns \cr
    \code{\link[spatstat.core]{Linhom}} \tab \eqn{L} for inhomogeneous point patterns \cr
    \code{\link[spatstat.core]{pcfinhom}} \tab pair correlation for inhomogeneous patterns\cr
    \code{\link[spatstat.core]{Finhom}} \tab \eqn{F} for inhomogeneous point patterns \cr
    \code{\link[spatstat.core]{Ginhom}} \tab \eqn{G} for inhomogeneous point patterns \cr
    \code{\link[spatstat.core]{Jinhom}} \tab \eqn{J} for inhomogeneous point patterns \cr
    \code{\link[spatstat.core]{localL}} \tab Getis-Franklin neighbourhood density function\cr
    \code{\link[spatstat.core]{localK}} \tab neighbourhood K-function\cr
    \code{\link[spatstat.core]{localpcf}} \tab local pair correlation function\cr
    \code{\link[spatstat.core]{localKinhom}} \tab local \eqn{K} for inhomogeneous point patterns \cr
    \code{\link[spatstat.core]{localLinhom}} \tab local \eqn{L} for inhomogeneous point patterns \cr
    \code{\link[spatstat.core]{localpcfinhom}} \tab local pair correlation for inhomogeneous patterns\cr
    \code{\link[spatstat.core]{Ksector}} \tab Directional \eqn{K}-function\cr
    \code{\link[spatstat.core]{Kscaled}} \tab locally scaled \eqn{K}-function \cr
    \code{\link[spatstat.core]{Kest.fft}} \tab fast \eqn{K}-function using FFT for large datasets \cr
    \code{\link[spatstat.core]{Kmeasure}} \tab reduced second moment measure \cr
    \code{\link[spatstat.core]{envelope}} \tab simulation envelopes for a summary
    function \cr
    \code{\link[spatstat.core]{varblock}} \tab variances and confidence intervals\cr
                   \tab for a summary function \cr
    \code{\link[spatstat.core]{lohboot}} \tab bootstrap for a summary function 
  }

  Related facilities:
  \tabular{ll}{
    \code{\link[spatstat.core]{plot.fv}} \tab plot a summary function\cr
    \code{\link[spatstat.core]{eval.fv}} \tab evaluate any expression involving
    summary functions\cr
    \code{\link[spatstat.core]{harmonise.fv}} \tab make functions compatible \cr
    \code{\link[spatstat.core]{eval.fasp}} \tab evaluate any expression involving
    an array of functions\cr
    \code{\link[spatstat.core]{with.fv}} \tab evaluate an expression for a 
    summary function\cr
    \code{\link[spatstat.core]{Smooth.fv}} \tab apply smoothing to a summary function\cr
    \code{\link[spatstat.core]{deriv.fv}} \tab calculate derivative of a summary function\cr
    \code{\link[spatstat.core]{pool.fv}} \tab pool several estimates of a summary function\cr
    \code{\link[spatstat.core]{density.ppp}} \tab kernel smoothed density\cr
    \code{\link[spatstat.core]{densityHeat.ppp}} \tab diffusion kernel smoothed density\cr
    \code{\link[spatstat.core]{Smooth.ppp}} \tab spatial interpolation of marks  \cr
    \code{\link[spatstat.core]{relrisk}} \tab kernel estimate of relative risk\cr
    \code{\link[spatstat.core]{sharpen.ppp}} \tab data sharpening  \cr
    \code{\link[spatstat.core]{rknn}} \tab theoretical distribution of nearest
    neighbour distance
 }

  \bold{Summary statistics for a multitype point pattern:}
  A multitype point pattern is represented by an object \code{X}
  of class \code{"ppp"} such that \code{marks(X)} is a factor. 
  \tabular{ll}{
    \code{\link[spatstat.core]{relrisk}} \tab kernel estimation of relative risk  \cr
    \code{\link[spatstat.core]{scan.test}} \tab spatial scan test of elevated risk  \cr
    \code{\link[spatstat.core]{Gcross},\link[spatstat.core]{Gdot},\link[spatstat.core]{Gmulti}} \tab
    multitype nearest neighbour distributions 
    \eqn{G_{ij}, G_{i\bullet}}{G[i,j], G[i.]} \cr
    \code{\link[spatstat.core]{Kcross},\link[spatstat.core]{Kdot}, \link[spatstat.core]{Kmulti}} \tab
    multitype \eqn{K}-functions 
    \eqn{K_{ij}, K_{i\bullet}}{K[i,j], K[i.]} \cr
    \code{\link[spatstat.core]{Lcross},\link[spatstat.core]{Ldot}} \tab
    multitype \eqn{L}-functions 
    \eqn{L_{ij}, L_{i\bullet}}{L[i,j], L[i.]} \cr
    \code{\link[spatstat.core]{Jcross},\link[spatstat.core]{Jdot},\link[spatstat.core]{Jmulti}} \tab
    multitype \eqn{J}-functions
    \eqn{J_{ij}, J_{i\bullet}}{J[i,j],J[i.]} \cr
    \code{\link[spatstat.core]{pcfcross}} \tab
    multitype pair correlation function \eqn{g_{ij}}{g[i,j]} \cr
    \code{\link[spatstat.core]{pcfdot}} \tab
    multitype pair correlation function \eqn{g_{i\bullet}}{g[i.]} \cr
    \code{\link[spatstat.core]{pcfmulti}} \tab
    general pair correlation function \cr
    \code{\link[spatstat.core]{markconnect}} \tab
    marked connection function \eqn{p_{ij}}{p[i,j]} \cr
    \code{\link[spatstat.core]{alltypes}} \tab  estimates of the above
    for all \eqn{i,j} pairs \cr
    \code{\link[spatstat.core]{Iest}} \tab  multitype \eqn{I}-function\cr
    \code{\link[spatstat.core]{Kcross.inhom},\link[spatstat.core]{Kdot.inhom}} \tab
    inhomogeneous counterparts of \code{Kcross}, \code{Kdot} \cr
    \code{\link[spatstat.core]{Lcross.inhom},\link[spatstat.core]{Ldot.inhom}} \tab
    inhomogeneous counterparts of \code{Lcross}, \code{Ldot} \cr
    \code{\link[spatstat.core]{pcfcross.inhom},\link[spatstat.core]{pcfdot.inhom}} \tab
    inhomogeneous counterparts of \code{pcfcross}, \code{pcfdot} \cr
    \code{\link[spatstat.core]{localKcross},\link[spatstat.core]{localKdot}} \tab
    local counterparts of \code{Kcross}, \code{Kdot} \cr
    \code{\link[spatstat.core]{localLcross},\link[spatstat.core]{localLdot}} \tab
    local counterparts of \code{Lcross}, \code{Ldot} \cr
    \code{\link[spatstat.core]{localKcross.inhom},\link[spatstat.core]{localLcross.inhom}} \tab
    local counterparts of \code{Kcross.inhom}, \code{Lcross.inhom} 
  }

  \bold{Summary statistics for a marked point pattern:}
  A marked point pattern is represented by an object \code{X}
  of class \code{"ppp"} with a component \code{X$marks}.
  The entries in the vector \code{X$marks} may be numeric, complex,
  string or any other atomic type. For numeric marks, there are the
  following functions:
  \tabular{ll}{
    \code{\link[spatstat.core]{markmean}} \tab smoothed local average of marks \cr
    \code{\link[spatstat.core]{markvar}} \tab smoothed local variance of marks \cr
    \code{\link[spatstat.core]{markcorr}} \tab mark correlation function \cr
    \code{\link[spatstat.core]{markcrosscorr}} \tab mark cross-correlation function \cr
    \code{\link[spatstat.core]{markvario}} \tab mark variogram \cr
    \code{\link[spatstat.core]{markmarkscatter}} \tab mark-mark scatterplot \cr
    \code{\link[spatstat.core]{Kmark}} \tab mark-weighted \eqn{K} function \cr
    \code{\link[spatstat.core]{Emark}} \tab mark independence diagnostic \eqn{E(r)} \cr
    \code{\link[spatstat.core]{Vmark}} \tab mark independence diagnostic \eqn{V(r)} \cr
    \code{\link[spatstat.core]{nnmean}} \tab nearest neighbour mean index \cr
    \code{\link[spatstat.core]{nnvario}} \tab nearest neighbour mark variance index 
  }
  For marks of any type, there are the following:
  \tabular{ll}{
    \code{\link[spatstat.core]{Gmulti}} \tab multitype nearest neighbour distribution \cr
    \code{\link[spatstat.core]{Kmulti}} \tab multitype \eqn{K}-function \cr
    \code{\link[spatstat.core]{Jmulti}} \tab multitype \eqn{J}-function 
  }
  Alternatively use \code{\link[spatstat.geom]{cut.ppp}} to convert a marked point pattern
  to a multitype point pattern.

  \bold{Programming tools:}
  \tabular{ll}{
    \code{\link[spatstat.core]{marktable}} \tab tabulate the marks of neighbours
    in a point pattern 
  }

  \bold{Summary statistics for a three-dimensional point pattern:}

  These are for 3-dimensional point pattern objects (class \code{pp3}).

  \tabular{ll}{
    \code{\link[spatstat.core]{F3est}} \tab empty space function \eqn{F} \cr
    \code{\link[spatstat.core]{G3est}} \tab nearest neighbour function \eqn{G} \cr
    \code{\link[spatstat.core]{K3est}} \tab \eqn{K}-function \cr
    \code{\link[spatstat.core]{pcf3est}} \tab pair correlation function
  }

  Related facilities:
  \tabular{ll}{
    \code{\link[spatstat.core]{envelope.pp3}} \tab simulation envelopes 
  }

  \bold{Summary statistics for random sets:}
  
  These work for point patterns (class \code{ppp}),
  line segment patterns (class \code{psp})
  or windows (class \code{owin}).
  
  \tabular{ll}{
    \code{\link[spatstat.core]{Hest}} \tab spherical contact distribution \eqn{H} \cr
    \code{\link[spatstat.core]{Gfox}} \tab Foxall \eqn{G}-function \cr
    \code{\link[spatstat.core]{Jfox}} \tab Foxall \eqn{J}-function
  }
  
 \bold{Model fitting (Cox and cluster models)}
  
  Cluster process models (with homogeneous or inhomogeneous intensity)
  and Cox processes can be fitted by the function \code{\link[spatstat.core]{kppm}}.
  Its result is an object of class \code{"kppm"}.
  The fitted model can be printed, plotted, predicted, simulated
  and updated.

  \tabular{ll}{
    \code{\link[spatstat.core]{kppm}} \tab  Fit model\cr
    \code{\link[spatstat.core]{plot.kppm}} \tab  Plot the fitted model\cr
    \code{\link[spatstat.core]{summary.kppm}} \tab  Summarise the fitted model\cr
    \code{\link[spatstat.core]{fitted.kppm}} \tab Compute fitted intensity \cr
    \code{\link[spatstat.core]{predict.kppm}} \tab Compute fitted intensity \cr
    \code{\link[spatstat.core]{update.kppm}} \tab Update the model \cr
    \code{\link[spatstat.core]{improve.kppm}} \tab Refine the estimate of trend \cr
    \code{\link[spatstat.core]{simulate.kppm}} \tab Generate simulated realisations \cr
    \code{\link[spatstat.core]{vcov.kppm}} \tab Variance-covariance matrix of coefficients \cr
    \code{\link[spatstat.core:methods.kppm]{coef.kppm}}
    \tab Extract trend coefficients \cr
    \code{\link[spatstat.core:methods.kppm]{formula.kppm}}
    \tab Extract trend formula \cr
    \code{\link[spatstat.core]{parameters}} \tab Extract all model parameters \cr
    \code{\link[spatstat.core]{clusterfield}} \tab Compute offspring density \cr
    \code{\link[spatstat.core]{clusterradius}} \tab Radius of support of offspring density \cr
    \code{\link[spatstat.core]{Kmodel.kppm}} \tab \eqn{K} function of fitted model \cr
    \code{\link[spatstat.core]{pcfmodel.kppm}} \tab Pair correlation of fitted model 
  }
  
  For model selection, you can also use 
  the generic functions \code{\link[stats]{step}}, \code{\link[stats]{drop1}} 
  and \code{\link[stats]{AIC}} on fitted point process models.
  For variable selection, see \code{\link[spatstat.core]{sdr}}.
  
  The theoretical models can also be simulated,
  for any choice of parameter values,
  using \code{\link[spatstat.core]{rThomas}}, \code{\link[spatstat.core]{rMatClust}},
  \code{\link[spatstat.core]{rCauchy}}, \code{\link[spatstat.core]{rVarGamma}},
  and \code{\link[spatstat.core]{rLGCP}}.
  
  Lower-level fitting functions include:

  \tabular{ll}{
    \code{\link[spatstat.core]{lgcp.estK}} \tab fit a log-Gaussian Cox process model\cr
    \code{\link[spatstat.core]{lgcp.estpcf}} \tab fit a log-Gaussian Cox process model\cr
    \code{\link[spatstat.core]{thomas.estK}} \tab fit the Thomas process model \cr
    \code{\link[spatstat.core]{thomas.estpcf}} \tab fit the Thomas process model \cr
    \code{\link[spatstat.core]{matclust.estK}} \tab fit the \Matern Cluster process model \cr
    \code{\link[spatstat.core]{matclust.estpcf}} \tab fit the \Matern Cluster process model \cr
    \code{\link[spatstat.core]{cauchy.estK}} \tab fit a Neyman-Scott Cauchy cluster process \cr
   \code{\link[spatstat.core]{cauchy.estpcf}} \tab fit a Neyman-Scott Cauchy cluster process\cr
   \code{\link[spatstat.core]{vargamma.estK}} \tab fit a Neyman-Scott Variance Gamma process\cr
   \code{\link[spatstat.core]{vargamma.estpcf}} \tab fit a Neyman-Scott Variance Gamma process\cr
    \code{\link[spatstat.core]{mincontrast}} \tab low-level algorithm for fitting models
    \cr \tab by the method of minimum contrast 
  }

  \bold{Model fitting (Poisson and Gibbs models)}
  
  Poisson point processes are the simplest models for point patterns.
  A Poisson model assumes that the points are stochastically
  independent. It may allow the points to have a non-uniform spatial
  density. The special case of a Poisson process with a uniform
  spatial density is often called Complete Spatial Randomness.
  
  Poisson point processes are included in the more general class of Gibbs point
  process models. In a Gibbs model, there is \emph{interaction}
  or dependence between points. Many different types of interaction
  can be specified.
  
  For a detailed explanation of how to fit Poisson or
  Gibbs point process models to point pattern data using \pkg{spatstat},
  see Baddeley and Turner (2005b) or Baddeley (2008).
  
  \bold{To fit a Poisson or Gibbs point process model:}

  Model fitting in \pkg{spatstat} is performed mainly by the function
  \code{\link[spatstat.core]{ppm}}. Its result is an object of class \code{"ppm"}.
  
  Here are some examples, where \code{X} is a point pattern (class
  \code{"ppp"}):
  
  \tabular{ll}{
    \emph{command} \tab \emph{model} \cr
    \code{ppm(X)} \tab Complete Spatial Randomness \cr
    \code{ppm(X ~ 1)} \tab Complete Spatial Randomness \cr
    \code{ppm(X ~ x)} \tab Poisson process with \cr
                             \tab intensity loglinear in \eqn{x} coordinate \cr
    \code{ppm(X ~ 1, Strauss(0.1))} \tab Stationary Strauss process \cr
    \code{ppm(X ~ x, Strauss(0.1))} \tab Strauss process with \cr
                              \tab conditional intensity loglinear in \eqn{x}
  }
  It is also possible to fit models that depend on
  other covariates.

  \bold{Manipulating the fitted model:}

  \tabular{ll}{
    \code{\link[spatstat.core]{plot.ppm}} \tab 		Plot the fitted model\cr
    \code{\link[spatstat.core]{predict.ppm}}
    \tab   Compute the spatial trend and conditional intensity\cr
    \tab   of the fitted point process model \cr
    \code{\link[spatstat.core]{coef.ppm}} \tab Extract the fitted model coefficients\cr
    \code{\link[spatstat.core]{parameters}} \tab Extract all model parameters\cr
    \code{\link[spatstat.core]{formula.ppm}} \tab Extract the trend formula\cr
    \code{\link[spatstat.core]{intensity.ppm}} \tab Compute fitted intensity \cr
    \code{\link[spatstat.core]{Kmodel.ppm}} \tab \eqn{K} function of fitted model \cr
    \code{\link[spatstat.core]{pcfmodel.ppm}} \tab pair correlation of fitted model \cr
    \code{\link[spatstat.core]{fitted.ppm}} \tab Compute fitted conditional intensity at quadrature points \cr
    \code{\link[spatstat.core]{residuals.ppm}} \tab Compute point process residuals at quadrature points \cr
    \code{\link[spatstat.core]{update.ppm}} \tab Update the fit \cr
    \code{\link[spatstat.core]{vcov.ppm}} \tab Variance-covariance matrix of estimates\cr
    \code{\link[spatstat.core]{rmh.ppm}} \tab Simulate from fitted model  \cr
    \code{\link[spatstat.core]{simulate.ppm}} \tab Simulate from fitted model  \cr
    \code{\link[spatstat.core]{print.ppm}} \tab Print basic information about a fitted model\cr
    \code{\link[spatstat.core]{summary.ppm}} \tab Summarise a fitted model\cr
    \code{\link[spatstat.core]{effectfun}} \tab Compute the fitted effect of one covariate\cr
    \code{\link[spatstat.core]{logLik.ppm}} \tab log-likelihood or log-pseudolikelihood\cr
    \code{\link[spatstat.core]{anova.ppm}} \tab Analysis of deviance \cr
    \code{\link[spatstat.core]{model.frame.ppm}} \tab Extract data frame used to fit model  \cr
    \code{\link[spatstat.core]{model.images}} \tab Extract spatial data used to fit model  \cr
    \code{\link[spatstat.core]{model.depends}} \tab Identify variables in the model \cr
    \code{\link[spatstat.core]{as.interact}} \tab Interpoint interaction component of model \cr
    \code{\link[spatstat.core]{fitin}} \tab Extract fitted interpoint interaction \cr
    \code{\link[spatstat.core]{is.hybrid}} \tab Determine whether the model is a hybrid \cr
    \code{\link[spatstat.core]{valid.ppm}} \tab Check the model is a valid point process \cr
    \code{\link[spatstat.core]{project.ppm}} \tab Ensure the model is a valid point process 
  }
  For model selection, you can also use 
  the generic functions \code{\link[stats]{step}}, \code{\link[stats]{drop1}} 
  and \code{\link[stats]{AIC}} on fitted point process models.
  For variable selection, see \code{\link[spatstat.core]{sdr}}.
  
  See \code{\link[spatstat.geom]{spatstat.options}} to control plotting of fitted model.
  
  \bold{To specify a point process model:}
  
  The first order ``trend'' of the model is determined by an \R 
  language formula. The formula specifies the form of the
  \emph{logarithm} of the trend.
  
  \tabular{ll}{
    \code{X ~ 1}  \tab No trend (stationary) \cr
    \code{X ~ x}   \tab Loglinear trend
      \eqn{\lambda(x,y) = \exp(\alpha + \beta x)}{lambda(x,y) =	exp(alpha + beta * x)} \cr
    \tab where \eqn{x,y} are Cartesian coordinates \cr
    \code{X ~ polynom(x,y,3)}  \tab Log-cubic polynomial trend  \cr
    \code{X ~ harmonic(x,y,2)}  \tab Log-harmonic polynomial trend \cr
    \code{X ~ Z}   \tab Loglinear function of covariate \code{Z} \cr
    \tab  \eqn{\lambda(x,y) = \exp(\alpha + \beta Z(x,y))}{lambda(x,y) =	exp(alpha + beta * Z(x,y))} 
  }

  The higher order (``interaction'') components are described by
  an object of class \code{"interact"}. Such objects are created by:
  \tabular{ll}{
    \code{\link[spatstat.core]{Poisson}()} \tab		the Poisson point process\cr
    \code{\link[spatstat.core]{AreaInter}()}	 \tab Area-interaction process\cr
    \code{\link[spatstat.core]{BadGey}()} \tab	multiscale Geyer process\cr
    \code{\link[spatstat.core]{Concom}()} \tab	connected component interaction\cr
    \code{\link[spatstat.core]{DiggleGratton}() } \tab Diggle-Gratton potential \cr
    \code{\link[spatstat.core]{DiggleGatesStibbard}() } \tab Diggle-Gates-Stibbard potential \cr
    \code{\link[spatstat.core]{Fiksel}()}	 \tab Fiksel pairwise interaction process\cr
    \code{\link[spatstat.core]{Geyer}()}	 \tab Geyer's saturation process\cr
    \code{\link[spatstat.core]{Hardcore}()}	 \tab Hard core process\cr
    \code{\link[spatstat.core]{HierHard}()}	 \tab Hierarchical multiype hard core process\cr
    \code{\link[spatstat.core]{HierStrauss}()}	 \tab Hierarchical multiype Strauss process\cr
    \code{\link[spatstat.core]{HierStraussHard}()}	 \tab
    Hierarchical multiype Strauss-hard core process\cr
    \code{\link[spatstat.core]{Hybrid}()}	 \tab Hybrid of several interactions\cr
    \code{\link[spatstat.core]{LennardJones}() } \tab Lennard-Jones potential \cr
    \code{\link[spatstat.core]{MultiHard}()} \tab 		multitype hard core process \cr
    \code{\link[spatstat.core]{MultiStrauss}()} \tab 		multitype Strauss process \cr
    \code{\link[spatstat.core]{MultiStraussHard}()} \tab 	multitype Strauss/hard core process \cr
    \code{\link[spatstat.core]{OrdThresh}()} \tab		Ord process, threshold potential\cr
    \code{\link[spatstat.core]{Ord}()} \tab 		        Ord model, user-supplied potential \cr
    \code{\link[spatstat.core]{PairPiece}()} \tab		pairwise interaction, piecewise constant \cr
    \code{\link[spatstat.core]{Pairwise}()} \tab	pairwise interaction, user-supplied potential\cr
    \code{\link[spatstat.core]{Penttinen}()} \tab	Penttinen pairwise interaction\cr
    \code{\link[spatstat.core]{SatPiece}()} \tab	Saturated pair model, piecewise  constant potential\cr
    \code{\link[spatstat.core]{Saturated}()} \tab	Saturated pair model, user-supplied potential\cr
    \code{\link[spatstat.core]{Softcore}()} \tab pairwise interaction, soft core potential\cr
    \code{\link[spatstat.core]{Strauss}()} \tab  Strauss process \cr
    \code{\link[spatstat.core]{StraussHard}()} \tab Strauss/hard core point process \cr
    \code{\link[spatstat.core]{Triplets}()} \tab Geyer triplets process
  }
  Note that it is also possible to combine several such interactions
  using \code{\link[spatstat.core]{Hybrid}}.
  
  \bold{Simulation and goodness-of-fit for fitted models:}
  
  \tabular{ll}{
    \code{\link[spatstat.core]{rmh.ppm}} \tab simulate realisations of a fitted model \cr
    \code{\link[spatstat.core]{simulate.ppm}} \tab simulate realisations of a fitted model \cr
    \code{\link[spatstat.core]{envelope}} \tab compute simulation envelopes for a
    fitted model 
  }

  \bold{Model fitting (determinantal point process models)}

  Code for fitting \emph{determinantal point process models} has 
  recently been added to \pkg{spatstat}.

  For information, see the help file for \code{\link[spatstat.core]{dppm}}.

  \bold{Model fitting (spatial logistic regression)}
  
  Pixel-based spatial logistic regression is an alternative
  technique for analysing spatial point patterns
  that is widely used in Geographical Information Systems.
  It is approximately equivalent to fitting a Poisson point process
  model.
  
  In pixel-based logistic regression, the spatial domain is
  divided into small pixels, the presence or absence of a
  data point in each pixel is recorded, and logistic regression
  is used to model the presence/absence indicators as a function
  of any covariates.
  
  Facilities for performing spatial logistic regression are
  provided in \pkg{spatstat} for comparison purposes.
  
  \bold{Fitting a spatial logistic regression}
  
  Spatial logistic regression is performed by the function
  \code{\link[spatstat.core]{slrm}}. Its result is an object of class \code{"slrm"}.
  There are many methods for this class, including methods for
  \code{print}, \code{fitted}, \code{predict}, \code{simulate},
  \code{anova}, \code{coef}, \code{logLik}, \code{terms},
  \code{update}, \code{formula} and \code{vcov}. 
  
  For example, if \code{X} is a point pattern (class
  \code{"ppp"}):
  
  \tabular{ll}{
    \emph{command} \tab \emph{model} \cr
    \code{slrm(X ~ 1)} \tab Complete Spatial Randomness \cr
    \code{slrm(X ~ x)} \tab Poisson process with \cr
                \tab intensity loglinear in \eqn{x} coordinate \cr
    \code{slrm(X ~ Z)} \tab Poisson process with \cr
                \tab intensity loglinear in covariate \code{Z}
  }

  \bold{Manipulating a fitted spatial logistic regression}
  
  \tabular{ll}{
    \code{\link[spatstat.core]{anova.slrm}} \tab Analysis of deviance \cr
    \code{\link[spatstat.core]{coef.slrm}}  \tab Extract fitted coefficients \cr
    \code{\link[spatstat.core]{vcov.slrm}}  \tab Variance-covariance matrix of fitted coefficients \cr
    \code{\link[spatstat.core]{fitted.slrm}} \tab Compute fitted probabilities or
    intensity \cr
    \code{\link[spatstat.core]{logLik.slrm}}   \tab Evaluate loglikelihood of fitted
    model \cr
    \code{\link[spatstat.core]{plot.slrm}}    \tab Plot fitted probabilities or
    intensity \cr
    \code{\link[spatstat.core]{predict.slrm}} \tab Compute predicted probabilities or
    intensity with new data \cr
    \code{\link[spatstat.core]{simulate.slrm}} \tab Simulate model
  }
  
  There are many other undocumented methods for this class,
  including methods for \code{print}, \code{update}, \code{formula}
  and \code{terms}. Stepwise model selection is
  possible using \code{step} or \code{stepAIC}.
  For variable selection, see \code{\link[spatstat.core]{sdr}}.

  \bold{Simulation}

  There are many ways to generate a random point pattern,
  line segment pattern, pixel image or tessellation
  in \pkg{spatstat}. 

  \bold{Random point patterns:}

  \tabular{ll}{
    \code{\link[spatstat.core]{runifpoint}} \tab
    generate \eqn{n} independent uniform random points \cr
    \code{\link[spatstat.core]{rpoint}} \tab
    generate \eqn{n} independent random points \cr
    \code{\link[spatstat.core]{rmpoint}} \tab
    generate \eqn{n} independent multitype random points \cr
    \code{\link[spatstat.core]{rpoispp}} \tab
    simulate the (in)homogeneous Poisson point process \cr
    \code{\link[spatstat.core]{rmpoispp}} \tab
    simulate the (in)homogeneous multitype Poisson point process \cr
    \code{\link[spatstat.core]{runifdisc}} \tab
    generate \eqn{n} independent uniform random points in disc\cr
    \code{\link[spatstat.core]{rstrat}} \tab
    stratified random sample of points \cr
    \code{\link[spatstat.core]{rMaternI}}  \tab
    simulate the \Matern Model I inhibition process\cr
    \code{\link[spatstat.core]{rMaternII}} \tab
    simulate the \Matern Model II inhibition process\cr
    \code{\link[spatstat.core]{rSSI}} \tab
    simulate Simple Sequential Inhibition process\cr
    \code{\link[spatstat.core]{rHardcore}} \tab
    simulate hard core process (perfect simulation)\cr
    \code{\link[spatstat.core]{rStrauss}} \tab
    simulate Strauss process (perfect simulation)\cr
    \code{\link[spatstat.core]{rStraussHard}} \tab
    simulate Strauss-hard core process (perfect simulation)\cr
    \code{\link[spatstat.core]{rDiggleGratton}} \tab
    simulate Diggle-Gratton process (perfect simulation)\cr
    \code{\link[spatstat.core]{rDGS}} \tab
    simulate Diggle-Gates-Stibbard process (perfect simulation)\cr
    \code{\link[spatstat.core]{rPenttinen}} \tab
    simulate Penttinen process (perfect simulation)\cr
    \code{\link[spatstat.core]{rNeymanScott}} \tab
    simulate a general Neyman-Scott process\cr
    \code{\link[spatstat.core]{rMatClust}} \tab
    simulate the \Matern Cluster process\cr
    \code{\link[spatstat.core]{rThomas}} \tab
    simulate the Thomas process  \cr
    \code{\link[spatstat.core]{rLGCP}} \tab
    simulate the log-Gaussian Cox process  \cr
    \code{\link[spatstat.core]{rGaussPoisson}}  \tab
    simulate the Gauss-Poisson cluster process\cr
    \code{\link[spatstat.core]{rCauchy}} \tab
    simulate Neyman-Scott process with Cauchy clusters \cr
    \code{\link[spatstat.core]{rVarGamma}} \tab
    simulate Neyman-Scott process with Variance Gamma clusters \cr
    \code{\link[spatstat.core]{rcell}} \tab
    simulate the Baddeley-Silverman cell process  \cr
    \code{\link[spatstat.core]{runifpointOnLines}} \tab
    generate \eqn{n} random points along specified line segments \cr
    \code{\link[spatstat.core]{rpoisppOnLines}} \tab
    generate Poisson random points along specified line segments 
  }
      
  \bold{Resampling a point pattern:}

  \tabular{ll}{
    \code{\link[spatstat.core]{quadratresample}} \tab block resampling \cr
    \code{\link[spatstat.core]{rshift}} \tab random shifting of (subsets of) points\cr
    \code{\link[spatstat.core]{rthin}} \tab  random thinning 
  }
  
  See also \code{\link[spatstat.core]{varblock}} for estimating the variance
  of a summary statistic by block resampling, and
  \code{\link[spatstat.core]{lohboot}} for another bootstrap technique.
  
  \bold{Fitted point process models:}

   If you have fitted a point process model to a point pattern dataset,
   the fitted model can be simulated.

   Cluster process models 
   are fitted by the function \code{\link[spatstat.core]{kppm}} yielding an
   object of class \code{"kppm"}. To generate one or more simulated
   realisations of this fitted model, use 
   \code{\link[spatstat.core]{simulate.kppm}}.

   Gibbs point process models 
   are fitted by the function \code{\link[spatstat.core]{ppm}} yielding an
   object of class \code{"ppm"}. To generate a simulated
   realisation of this fitted model, use \code{\link[spatstat.core]{rmh}}.
   To generate one or more simulated realisations of the fitted model,
   use \code{\link[spatstat.core]{simulate.ppm}}.

   \bold{Other random patterns:}

   \tabular{ll}{
     \code{\link[spatstat.core]{rpoisline}} \tab
     simulate the Poisson line process within a window \cr
     \code{\link[spatstat.core]{rpoislinetess}} \tab
     generate random tessellation using Poisson line process \cr
     \code{\link[spatstat.core]{rMosaicSet}} \tab
     generate random set by selecting some tiles of a tessellation \cr
     \code{\link[spatstat.core]{rMosaicField}} \tab
     generate random pixel image by assigning random values
     in each tile of a tessellation
   }

   \bold{Simulation-based inference}

   \tabular{ll}{
    \code{\link[spatstat.core]{envelope}} \tab critical envelope for Monte Carlo
    test of goodness-of-fit \cr
    \code{\link[spatstat.core]{bits.envelope}} \tab critical envelope for balanced
    two-stage Monte Carlo test \cr
    \code{\link[spatstat.core]{qqplot.ppm}} \tab diagnostic plot for interpoint
    interaction \cr
    \code{\link[spatstat.core]{scan.test}} \tab spatial scan statistic/test \cr
    \code{\link[spatstat.core]{studpermu.test}} \tab studentised permutation test\cr
    \code{\link[spatstat.core]{segregation.test}} \tab test of segregation of types 
  }

  \bold{Hypothesis tests:}

  \tabular{ll}{
    \code{\link[spatstat.core]{quadrat.test}} \tab \eqn{\chi^2}{chi^2} goodness-of-fit
    test on quadrat counts \cr
    \code{\link[spatstat.core]{clarkevans.test}} \tab Clark and Evans test \cr
    \code{\link[spatstat.core]{cdf.test}} \tab Spatial distribution goodness-of-fit test\cr
    \code{\link[spatstat.core]{berman.test}} \tab Berman's goodness-of-fit tests\cr
    \code{\link[spatstat.core]{envelope}} \tab critical envelope for Monte Carlo
    test of goodness-of-fit \cr
    \code{\link[spatstat.core]{scan.test}} \tab spatial scan statistic/test \cr
    \code{\link[spatstat.core]{dclf.test}} \tab Diggle-Cressie-Loosmore-Ford test \cr
    \code{\link[spatstat.core]{mad.test}} \tab Mean Absolute Deviation test \cr
    \code{\link[spatstat.core]{anova.ppm}} \tab Analysis of Deviance for
    point process models 
  }

  \bold{More recently-developed tests:}
  
  \tabular{ll}{
    \code{\link[spatstat.core]{dg.test}} \tab Dao-Genton test \cr
    \code{\link[spatstat.core]{bits.test}} \tab Balanced independent two-stage test \cr
    \code{\link[spatstat.core]{dclf.progress}} \tab Progress plot for DCLF test \cr
    \code{\link[spatstat.core]{mad.progress}} \tab Progress plot for MAD test \cr
  }

  \bold{Sensitivity diagnostics:}

  Classical measures of model sensitivity such as leverage and influence
  have been adapted to point process models.
  
    \tabular{ll}{
      \code{\link[spatstat.core]{leverage.ppm}} \tab Leverage for point process model\cr
      \code{\link[spatstat.core]{influence.ppm}} \tab Influence for point process model\cr
      \code{\link[spatstat.core]{dfbetas.ppm}} \tab Parameter influence\cr
      \code{\link[spatstat.core]{dffit.ppm}} \tab Effect change diagnostic
    }
  
  \bold{Diagnostics for covariate effect:}

  Classical diagnostics for covariate effects have been adapted to
  point process models.

  \tabular{ll}{
    \code{\link[spatstat.core]{parres}} \tab Partial residual plot\cr
    \code{\link[spatstat.core]{addvar}} \tab Added variable plot \cr
    \code{\link[spatstat.core]{rhohat}} \tab Kernel estimate of covariate effect\cr
    \code{\link[spatstat.core]{rho2hat}} \tab Kernel estimate of covariate effect
    (bivariate) 
  }
  
  \bold{Residual diagnostics:}
  
  Residuals for a fitted point process model, and diagnostic plots
  based on the residuals, were introduced in Baddeley et al (2005) and
  Baddeley, Rubak and \Moller (2011).
  
  Type \code{demo(diagnose)}
  for a demonstration of the diagnostics features.

  \tabular{ll}{
    \code{\link[spatstat.core]{diagnose.ppm}} \tab diagnostic plots for spatial trend\cr
    \code{\link[spatstat.core]{qqplot.ppm}} \tab diagnostic Q-Q plot for interpoint interaction\cr
    \code{\link[spatstat.data]{residualspaper}} \tab examples from Baddeley et al (2005) \cr
    \code{\link[spatstat.core]{Kcom}} \tab model compensator of \eqn{K} function \cr
    \code{\link[spatstat.core]{Gcom}} \tab model compensator of \eqn{G} function \cr
    \code{\link[spatstat.core]{Kres}} \tab score residual of \eqn{K} function \cr
    \code{\link[spatstat.core]{Gres}} \tab score residual of \eqn{G} function \cr
    \code{\link[spatstat.core]{psst}} \tab pseudoscore residual of summary function \cr
    \code{\link[spatstat.core]{psstA}} \tab pseudoscore residual of empty space function \cr
    \code{\link[spatstat.core]{psstG}} \tab pseudoscore residual of \eqn{G} function \cr
    \code{\link[spatstat.core]{compareFit}} \tab compare compensators of several fitted
    models
  }


  \bold{Resampling and randomisation procedures}

  You can build your own tests based on randomisation
  and resampling using the following capabilities:
  
  \tabular{ll}{
    \code{\link[spatstat.core]{quadratresample}} \tab block resampling \cr
    \code{\link[spatstat.core]{rshift}} \tab random shifting of (subsets of) points\cr
    \code{\link[spatstat.core]{rthin}} \tab  random thinning  
  }
}
\section{Licence}{
  This library and its documentation are usable under the terms of the "GNU 
  General Public License", a copy of which is distributed with the package.
}
\author{
  \spatstatAuthors.
}
\section{Acknowledgements}{
  Kasper Klitgaard Berthelsen,
  Ottmar Cronie,
  Tilman Davies,
  Julian Gilbey, 
  Yongtao Guan,
  Ute Hahn,
  Kassel Hingee,
  Abdollah Jalilian,
  Marie-Colette van Lieshout,
  Greg McSwiggan,
  Tuomas Rajala,
  Suman Rakshit,
  Dominic Schuhmacher,
  Rasmus Waagepetersen
  and
  Hangsheng Wang
  made substantial contributions of code.

  For comments, corrections, bug alerts and suggestions, we thank
  Monsuru Adepeju,
  Corey Anderson,
  Ang Qi Wei,
  Ryan Arellano,
  Jens \ifelse{latex}{\out{{\AA}str{\" o}m}}{Astrom},
  Robert Aue,
  Marcel Austenfeld,
  Sandro Azaele,
  Malissa Baddeley,
  Guy Bayegnak,
  Colin Beale,
  Melanie Bell,
  Thomas Bendtsen,
  Ricardo Bernhardt,
  Andrew Bevan,
  Brad Biggerstaff,
  Anders Bilgrau,
  Leanne Bischof,
  Christophe Biscio,
  Roger Bivand,
  Jose M. Blanco Moreno,
  Florent Bonneu,
  Jordan Brown,
  Ian Buller,
  Julian Burgos,
  Simon Byers,
  Ya-Mei Chang,
  Jianbao Chen,
  Igor Chernayavsky,
  Y.C. Chin,
  Bjarke Christensen,
  \ifelse{latex}{\out{Luc\'{\i{}}a}}{Lucia} Cobo Sanchez,
  Jean-Francois Coeurjolly,
  Kim Colyvas,
  Hadrien Commenges,
  Rochelle Constantine,
  Robin Corria Ainslie,
  Richard Cotton,
  Marcelino de la Cruz,
  Peter Dalgaard,
  Mario D'Antuono,
  Sourav Das,
  Peter Diggle,
  Patrick Donnelly,
  Ian Dryden,
  Stephen Eglen,
  Ahmed El-Gabbas,
  Belarmain Fandohan,
  Olivier Flores,
  David Ford,
  Peter Forbes,
  Shane Frank,
  Janet Franklin,
  Funwi-Gabga Neba,
  Oscar Garcia,
  Agnes Gault,
  Jonas Geldmann,
  Marc Genton,
  Shaaban Ghalandarayeshi,
  Jason Goldstick,
  Pavel Grabarnik,
  C. Graf,
  Ute Hahn,
  Andrew Hardegen,
  Martin \Bogsted Hansen,
  Martin Hazelton,
  Juha Heikkinen,
  Mandy Hering,
  Markus Herrmann,
  Maximilian Hesselbarth,
  Paul Hewson,
  Hamidreza Heydarian,
  Kurt Hornik,
  Philipp Hunziker,
  Jack Hywood,
  Ross Ihaka,
  \ifelse{latex}{\out{\u{C}enk I\c{c}\"{o}s}}{Cenk Icos},
  Aruna Jammalamadaka,
  Robert John-Chandran,
  Devin Johnson,
  Mahdieh Khanmohammadi,
  Bob Klaver,
  Lily Kozmian-Ledward,
  Peter Kovesi,
  Mike Kuhn,
  Jeff Laake,
  Robert Lamb,
  \ifelse{latex}{\out{Fr\'{e}d\'{e}ric}}{Frederic} Lavancier,
  Tom Lawrence,
  Tomas Lazauskas,
  Jonathan Lee,
  George Leser,
  Angela Li,
  Li Haitao,
  George Limitsios,
  Andrew Lister,
  Nestor Luambua,
  Ben Madin,
  Martin Maechler,
  Kiran Marchikanti,
  Jeff Marcus,
  Robert Mark,
  Peter McCullagh,
  Monia Mahling,
  Jorge Mateu Mahiques,
  Ulf Mehlig,
  Frederico Mestre,
  Sebastian Wastl Meyer,
  Mi Xiangcheng,
  Lore De Middeleer,
  Robin Milne,
  Enrique Miranda,
  Jesper \Moller,
  Annie \ifelse{latex}{\out{Molli{\'e}}}{Mollie},
  Ines Moncada,
  Mehdi Moradi,
  Virginia Morera Pujol,
  Erika Mudrak,
  Gopalan Nair,
  Nader Najari,
  Nicoletta Nava,
  Linda Stougaard Nielsen,
  Felipe Nunes,
  Jens Randel Nyengaard,
  Jens \Oehlschlaegel,
  Thierry Onkelinx,
  Sean O'Riordan,
  Evgeni Parilov,
  Jeff Picka,
  Nicolas Picard,
  Tim Pollington,
  Mike Porter,
  Sergiy Protsiv,
  Adrian Raftery,
  Ben Ramage,
  Pablo Ramon,
  Xavier Raynaud,
  Nicholas Read,
  Matt Reiter,
  Ian Renner,
  Tom Richardson,
  Brian Ripley,
  Ted Rosenbaum,
  Barry Rowlingson,
  Jason Rudokas,
  Tyler Rudolph,
  John Rudge,
  Christopher Ryan,
  Farzaneh Safavimanesh,
  Aila \Sarkka,
  Cody Schank,
  Katja Schladitz,
  Sebastian Schutte,
  Bryan Scott,
  Olivia Semboli,
  \ifelse{latex}{\out{Fran\c{c}ois S\'{e}m\'{e}curbe}}{Francois Semecurbe},
  Vadim Shcherbakov,
  Shen Guochun,
  Shi Peijian,
  Harold-Jeffrey Ship,
  Tammy L Silva,
  Ida-Maria Sintorn,
  Yong Song, 
  Malte Spiess,
  Mark Stevenson,
  Kaspar Stucki,
  Jan Sulavik,
  Michael Sumner,
  P. Surovy,
  Ben Taylor,
  Thordis Linda Thorarinsdottir,
  Leigh Torres,
  Berwin Turlach,
  Torben Tvedebrink,
  Kevin Ummer,
  Medha Uppala,
  Andrew van Burgel,
  Tobias Verbeke,
  Mikko Vihtakari,
  Alexendre Villers,
  Fabrice Vinatier,
  Maximilian Vogtland,
  Sasha Voss,
  Sven Wagner,
  Hao Wang,
  H. Wendrock,
  Jan Wild,
  Carl G. Witthoft,
  Selene Wong,
  Maxime Woringer,
  Luke Yates,
  Mike Zamboni
  and
  Achim Zeileis.
}
\keyword{spatial}
\keyword{package}

