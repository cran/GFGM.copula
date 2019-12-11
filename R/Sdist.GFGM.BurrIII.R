#' Sub-distribution functions under the generalized FGM copula with the Burr III margins
#'
#' @param time Vector of times.
#' @param Alpha Positive shape parameter for the Burr III margin (failure cause 1).
#' @param Beta Positive shape parameter for the Burr III margin (failure cause 2).
#' @param Gamma Common positive shape parameter for the Burr III margins.
#' @param p Copula parameter that greater than or equal to 1.
#' @param q Copula parameter that greater than 1 (integer).
#' @param theta Copula parameter with restricted range.
#' @param eta Location parameter with default value 0.
#' @description Sub-distribution functions under the generalized FGM copula with the Burr III margins.
#' @details The copula parameter \code{q} is restricted to be a integer due to the binominal theorem.
#' The admissible range of \code{theta} is given in \code{Dependence.GFGM}.
#'
#' @return \item{time}{Failure times}
#' \item{Sdist.1}{Probability of an object fails due to the failure cause 1.}
#' \item{Sdist.2}{Probability of an object fails due to the failure cause 2.}
#'
#' @references Shih J-H, Emura T (2018) Likelihood-based inference for bivariate latent failure time models with competing risks udner the generalized FGM copula, Computational Statistics, 33:1293-1323.
#' @references Shih J-H, Emura T (2019) Bivariate dependence measures and bivariate competing risks models under the generalized FGM copula, Statistical Papers, 60:1101-1118.
#' @seealso \code{\link{MLE.GFGM.BurrIII}}, \code{\link{Dependence.GFGM}}
#' @export
#'
#' @examples
#' library(GFGM.copula)
#' Sdist.GFGM.BurrIII(c(1:5),1,1,1,3,2,0.75,eta = 1)

Sdist.GFGM.BurrIII = function(time,Alpha,Beta,Gamma,p,q,theta,eta = 0) {

  ### checking inputs ###
  if (eta > min(time)) {stop("time cannot be smaller than eta.")}
  if (Alpha <= 0) {stop("Alpha must be positive")}
  if (Beta <= 0) {stop("Beta must be positive")}
  if (Gamma <= 0) {stop("Gamma must be positive")}
  if (p < 1) {stop("p must be greater than or equal to 1")}
  if (q <= 1 | q != round(q)) {stop("q must be greater than 1 (integer)")}
  theta.UB = ((1+p*q)/(q-1))^(q-1)/p^q
  theta.LB = -min(1,(((1+p*q)/(q-1))^(q-1)/p^q)^2)
  if (theta > theta.UB | theta < theta.LB) {stop("theta is invalid")}

  H  = (1+(time-eta)^(-Gamma))^-1
  K = 0
  for ( i in 0:q ) {

    for ( j in 0:(q-1) ) {

      k1 = H^(Alpha*(p*j+1)+Beta*(p*i+1))/(Alpha*(p*j+1)+Beta*(p*i+1))
      k2 = H^(Alpha*(p*j+p+1)+Beta*(p*i+1))/(Alpha*(p*j+p+1)+Beta*(p*i+1))
      K  = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*(k1-(1+p*q)*k2)

    }

  }

  Sdist.1 = H^Alpha-Alpha/(Alpha+Beta)*H^(Alpha+Beta)-theta*Alpha*K

  H  = (1+(time-eta)^(-Gamma))^-1
  K = 0
  for ( i in 0:q ) {

    for ( j in 0:(q-1) ) {

      k1 = H^(Beta*(p*j+1)+Alpha*(p*i+1))/(Beta*(p*j+1)+Alpha*(p*i+1))
      k2 = H^(Beta*(p*j+p+1)+Alpha*(p*i+1))/(Beta*(p*j+p+1)+Alpha*(p*i+1))
      K  = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*(k1-(1+p*q)*k2)

    }

  }

  Sdist.2 = H^Beta-Beta/(Alpha+Beta)*H^(Alpha+Beta)-theta*Beta*K

  return(cbind(time,Sdist.1,Sdist.2))

}
