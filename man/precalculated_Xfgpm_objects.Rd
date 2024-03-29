% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/8_precalculated_Xfgpm_objects.R
\docType{data}
\name{precalculated_Xfgpm_objects}
\alias{precalculated_Xfgpm_objects}
\alias{xm}
\alias{xm25}
\alias{xmc}
\alias{xmh}
\alias{xms}
\title{Precalculated Xfgpm objects}
\format{
Five objects of class \code{"Xfgpm"}:
\describe{
  \item{xm}{With 32 training points and default parameters.}
  \item{xm25}{With 32 training points and 25 iterations of the algorithm.}
  \item{xmc}{With 32 training points and customized solution space.}
  \item{xmh}{With 32 training points and customized heuristic parameters.}
  \item{xms}{With 32 training points and a time budget constraint and large number of iterations.}
}
}
\description{
A dataset containing the results of the application of
    \code{fgpm_factory} to \code{fgp_BB7} analytic black-box
    function. See \bold{Examples} for details.
}
\examples{

\dontrun{

##################################################################
## Construction of xm object with default parameters (~12 seconds)
##################################################################
set.seed(100)
n.tr <- 32
x1 <- x2 <- x3 <- x4 <- x5 <- seq(0,1,length = n.tr^(1/5))
sIn <- expand.grid(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5)
fIn <- list(f1 = matrix(runif(n.tr * 10), ncol = 10),
            f2 = matrix(runif(n.tr * 22), ncol = 22))
sOut <- fgp_BB7(sIn, fIn, n.tr)
xm <- fgpm_factory(sIn = sIn, fIn = fIn, sOut = sOut)

##################################################################
## Construction of xm25 object with 25 iterations (~20 seconds)
##################################################################
xm25 <- fgpm_factory(sIn = sIn, fIn = fIn, sOut = sOut,
                     setup = list(n.iter = 25))

##################################################################
## Construction of xmc object with customized solution space (~12 seconds)
##################################################################
myctr <- list(s_keepOn = c(1,2), # keep both scalar inputs always on
f_keepOn = c(2), # keep f2 always active
f_disTypes = list("2" = c("L2_byindex")), # only use L2_byindex distance for f2
f_fixDims = matrix(c(2,4), ncol = 1), # f2 projected in dimension 4
f_maxDims = matrix(c(1,5), ncol = 1), # f1 projected in dimension max 5
f_basTypes = list("1" = c("B-splines")), # only use B-splines projection for f1
kerTypes = c("matern5_2", "gauss")) # test only Matern 5/2 and Gaussian kernels
xmc <- fgpm_factory(sIn = sIn, fIn = fIn, sOut = sOut, ctraints = myctr)

##################################################################
## Construction of xmc object with customized heuristic parameters (~15 seconds)
##################################################################
mysup <- list(n.iter = 30, n.pop = 12, tao0 = .15, dop.s = 1.2,
              dop.f = 1.3, delta.f = 4, dispr.f = 1.1, q0 = .85,
              rho.l = .2, u.gbest = TRUE, n.ibest = 2, rho.g = .08)
xmh <- fgpm_factory(sIn = sIn, fIn = fIn, sOut = sOut, setup = mysup)

##################################################################
## Construction of xmc object with time budget constraint (~60 seconds)
##################################################################
mysup <- list(n.iter = 2000)
mytlim <- 60
xms <- fgpm_factory(sIn = sIn, fIn = fIn, sOut = sOut,
                    setup = mysup, time.lim = mytlim)

}

}
\keyword{data}
