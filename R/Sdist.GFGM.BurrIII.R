#' Sub-distribution function under the generalized FGM copula with the Burr III margins
#'
#' @param t.time Time.
#' @param Alpha Positive shape parameter for the Burr III margin (failure cause 1).
#' @param Beta Positive shape parameter for the Burr III margin (failure cause 2).
#' @param Gamma Common positive shape parameter for the Burr III margins.
#' @param p Copula parameter that greater than or equal to 1.
#' @param q Copula parameter that greater than or equal to 1 (integer).
#' @param theta Copula parameter with restricted range.
#' @param cause Failure cause.
#' @param eta Location parameter with default value 0.
#' @description Sub-distribution function under the generalized FGM copula with the Burr III margins.
#' @details The original paper is submitted for review.
#'
#' The copula parameter \code{q} is restricted to be integer due to the binominal theorem.
#' The restricted range of \eqn{\theta} is given in \code{Dependence.GFGM}
#'
#' @return Probability of an object fails at \code{t.time} due to the failure cause 1 or 2.
#'
#' @references Shih and Emura (2016) Bivariate dependence measures and bivariate competing risks models under the generalized FGM copula, Statistical Papers, doi: 10.1007/s00362-016-0865-5.
#' @references Shih and Emura (2017) Likelihood inference for bivariate latent failure time models with competing risks udner the generalized FGM copula (in re-submission, Computational Statistics).
#' @seealso \code{\link{Dependence.GFGM}}
#' @export
#'
#' @examples
#' library(GFGM.copula)
#' Sdist.GFGM.BurrIII(2,1,1,1,3,2,0.5,1,1)
#' Sdist.GFGM.BurrIII(2,1,1,1,3,2,0.5,2,1)

Sdist.GFGM.BurrIII = function(t.time,Alpha,Beta,Gamma,p,q,theta,cause,eta = 0) {

  t.time = t.time
  Alpha = Alpha
  Beta = Beta
  Gamma = Gamma
  p = p
  q = q
  theta = theta

  ### checking inputs ###
  if (cause != 1 & cause != 2) {stop("cause must be 1 or 2")}
  if (eta > min(t.time)) {stop("t.time cannot be smaller than eta.")}
  if (Alpha <= 0) {stop("Alpha must be positive")}
  if (Beta <= 0) {stop("Beta must be positive")}
  if (Gamma <= 0) {stop("Gamma must be positive")}
  if (p < 1) {stop("p must be greater than or equal to 1")}
  if (q < 1 | q != round(q)) {stop("q must be greater than or equal to 1 (integer)")}
  theta.UB = ((1+p*q)/(q-1))^(q-1)/p^q
  theta.LB = -min(1,(((1+p*q)/(q-1))^(q-1)/p^q)^2)
  if (theta > theta.UB | theta < theta.LB) {stop("theta is invalid")}

  if (cause == 2) {

    temp.A = Alpha
    temp.B = Beta
    Alpha = temp.B
    Beta = temp.A

  }

  H  = (1+(t.time-eta)^(-Gamma))^-1
  K = 0
  for ( i in 0:q ) {

    for ( j in 0:(q-1) ) {

      k1 = H^(Alpha*(p*j+1)+Beta*(p*i+1))/(Alpha*(p*j+1)+Beta*(p*i+1))
      k2 = H^(Alpha*(p*j+p+1)+Beta*(p*i+1))/(Alpha*(p*j+p+1)+Beta*(p*i+1))
      K  = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*(k1-(1+p*q)*k2)

    }

  }

  return(H^Alpha-Alpha/(Alpha+Beta)*H^(Alpha+Beta)-theta*Alpha*K)

}
