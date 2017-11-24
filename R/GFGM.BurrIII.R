#' Generate samples from the generalized FGM copula with the Burr III margins
#'
#' @param n Sample size.
#' @param p Copula parameter that greater than or equal to 1.
#' @param q Copula parameter that greater than 1.
#' @param theta Copula parameter with restricted range.
#' @param Alpha Positive shape parameter for the Burr III margin.
#' @param Beta Positive shape parameter for the Burr III margin.
#' @param Gamma Common positive shape parameter for the Burr III margins.
#' @description Generate samples from the generalized FGM copula with the Burr III margins.
#' @details The admissiable range of \code{theta} is given in \code{Dependence.GFGM}
#' @return \item{X}{\code{X} is asscoiated with the parameter \code{Alpha}.}
#' \item{Y}{\code{Y} is asscoiated with the parameter \code{Beta}.}
#'
#' @references Shih and Emura (2016) Bivariate dependence measures and bivariate competing risks models under the generalized FGM copula, Statistical Papers, doi: 10.1007/s00362-016-0865-5.
#' @references Shih and Emura (2017) Likelihood inference for bivariate latent failure time models with competing risks udner the generalized FGM copula (in re-submission, Computational Statistics).
#' @seealso \code{\link{Dependence.GFGM}}
#' @importFrom stats runif integrate uniroot
#' @export
#'
#' @examples
#' library(GFGM.copula)
#' GFGM.BurrIII(5,3,2,0.75,1,1,1)

GFGM.BurrIII = function(n,p,q,theta,Alpha,Beta,Gamma) {

  U = rep(0,n)
  V = rep(0,n)
  for (i in 1:n) {

    v = runif(1)
    a = runif(1)
    K = (1-v^p)^(q-1)*(1-(1+p*q)*v^p)
    f = function(u) {u+theta*K*u*(1-u^p)^q-a}
    u = uniroot(f,c(0,1),tol = 1e-9)

    U[i] = u$root
    V[i] = v

  }
  X = (U^(-1/Alpha)-1)^(-1/Gamma)
  Y = (V^(-1/Beta)-1)^(-1/Gamma)

  return(cbind(X,Y))

}
