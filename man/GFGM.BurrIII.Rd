% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GFGM.BurrIII.R
\name{GFGM.BurrIII}
\alias{GFGM.BurrIII}
\title{Generate samples from the generalized FGM copula with the Burr III margins}
\usage{
GFGM.BurrIII(n, p, q, theta, Alpha, Beta, Gamma)
}
\arguments{
\item{n}{Sample size.}

\item{p}{Copula parameter that greater than or equal to 1.}

\item{q}{Copula parameter that greater than 1.}

\item{theta}{Copula parameter with restricted range.}

\item{Alpha}{Positive shape parameter for the Burr III margin.}

\item{Beta}{Positive shape parameter for the Burr III margin.}

\item{Gamma}{Common positive shape parameter for the Burr III margins.}
}
\value{
\item{X}{\code{X} is asscoiated with the parameter \code{Alpha}.}
\item{Y}{\code{Y} is asscoiated with the parameter \code{Beta}.}
}
\description{
Generate samples from the generalized FGM copula with the Burr III margins.
}
\details{
The admissible range of \code{theta} is given in \code{Dependence.GFGM}.
}
\examples{
library(GFGM.copula)
GFGM.BurrIII(5,3,2,0.75,1,1,1)
}
\references{
Shih J-H, Emura T (2018) Likelihood-based inference for bivariate latent failure time models with competing risks udner the generalized FGM copula, Computational Statistics, 33:1293-1323.

Shih J-H, Emura T (2019) Bivariate dependence measures and bivariate competing risks models under the generalized FGM copula, Statistical Papers, 60:1101-1118.
}
\seealso{
\code{\link{Dependence.GFGM}}
}
