% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MLE.GFGM.spline.R
\name{MLE.GFGM.spline}
\alias{MLE.GFGM.spline}
\title{Maximum likelihood estimation for bivariate dependent competing risks data under the generalized FGM copula with the marginal distributions approximated by splines}
\usage{
MLE.GFGM.spline(t.event, event1, event2, p, q, theta, h.plot = TRUE)
}
\arguments{
\item{t.event}{Vector of the observed failure times.}

\item{event1}{Vector of the indicators for the failure cause 1.}

\item{event2}{Vector of the indicators for the failure cause 2.}

\item{p}{Copula parameter that greater than or equal to 1.}

\item{q}{Copula parameter that greater than 1 (integer).}

\item{theta}{Copula parameter with restricted range.}

\item{h.plot}{Plot hazard functions if \code{TRUE}.}
}
\value{
\item{n}{Sample size.}
\item{g1}{Maximum likelihood estimator of the splines coefficients for the failure cause 1.}
\item{g2}{Maximum likelihood estimator of the splines coefficients for the failure cause 2.}
\item{g1.var}{Covariance matrix of splines coefficients estimates for the failure cause 1.}
\item{g2.var}{Covariance matrix of splines coefficients estimates for the failure cause 2.}
}
\description{
Maximum likelihood estimation for bivariate dependent competing risks data under the generalized FGM copula with the marginal distributions approximated by splines.
}
\details{
The copula parameter \code{q} is restricted to be a integer due to the binominal theorem.
The admissible range of \code{theta} is given in \code{Dependence.GFGM}.

To adapt our functions to dependent censoring models in Emura and Chen (2018), one can simply set \code{event2} = \code{1-event1}.
}
\examples{
con   = c(16,224,16,80,128,168,144,176,176,568,392,576,128,56,112,160,384,600,40,416,
          408,384,256,246,184,440,64,104,168,408,304,16,72,8,88,160,48,168,80,512,
          208,194,136,224,32,504,40,120,320,48,256,216,168,184,144,224,488,304,40,160,
          488,120,208,32,112,288,336,256,40,296,60,208,440,104,528,384,264,360,80,96,
          360,232,40,112,120,32,56,280,104,168,56,72,64,40,480,152,48,56,328,192,
          168,168,114,280,128,416,392,160,144,208,96,536,400,80,40,112,160,104,224,336,
          616,224,40,32,192,126,392,288,248,120,328,464,448,616,168,112,448,296,328,56,
          80,72,56,608,144,408,16,560,144,612,80,16,424,264,256,528,56,256,112,544,
          552,72,184,240,128,40,600,96,24,184,272,152,328,480,96,296,592,400,8,280,
          72,168,40,152,488,480,40,576,392,552,112,288,168,352,160,272,320,80,296,248,
          184,264,96,224,592,176,256,344,360,184,152,208,160,176,72,584,144,176)
uncon = c(368,136,512,136,472,96,144,112,104,104,344,246,72,80,312,24,128,304,16,320,
          560,168,120,616,24,176,16,24,32,232,32,112,56,184,40,256,160,456,48,24,
          200,72,168,288,112,80,584,368,272,208,144,208,114,480,114,392,120,48,104,272,
          64,112,96,64,360,136,168,176,256,112,104,272,320,8,440,224,280,8,56,216,
          120,256,104,104,8,304,240,88,248,472,304,88,200,392,168,72,40,88,176,216,
          152,184,400,424,88,152,184)
cen   = rep(630,44)

t.event = c(con,uncon,cen)
event1  = c(rep(1,length(con)),rep(0,length(uncon)),rep(0,length(cen)))
event2  = c(rep(0,length(con)),rep(1,length(uncon)),rep(0,length(cen)))

library(GFGM.copula)
MLE.GFGM.spline(t.event,event1,event2,3,2,0.75)
}
\references{
Emura T, Chen Y-H (2018) Analysis of Survival Data with Dependent Censoring, Copula-Based Approaches, JSS Research Series in Statistics, Springer, Singapore.

Shih J-H, Emura T (2018) Likelihood-based inference for bivariate latent failure time models with competing risks udner the generalized FGM copula, Computational Statistics, 33:1293-1323.

Shih J-H, Emura T (2019) Bivariate dependence measures and bivariate competing risks models under the generalized FGM copula, Statistical Papers, 60:1101-1118.
}
\seealso{
\code{\link{Dependence.GFGM}}
}
