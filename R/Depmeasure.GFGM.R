#' Bivariate dependence measures under the generalized FGM copula
#'
#' @param p Copula parameter that greater than or equal to 1.
#' @param q Copula parameter that greater than 1.
#' @param theta Copula parameter with restricted range.
#' @description Compute Kendall's tau and Spearman's rho with their boundaries under the generalized FGM copula.
#' @details The admissible range of \code{theta} (\eqn{\theta}) is
#' \deqn{-\min\bigg\{1,\frac{1}{p^{2q}}\bigg(\frac{1+pq}{q-1}\bigg)^{2q-2}\bigg\} \leq \theta \leq \frac{1}{p^{q}}\bigg(\frac{1+pq}{q-1}\bigg)^{q-1}.}
#' See also Shih and Emura (2019) for details.
#' @return \item{theta}{Dependence parameter.}
#' \item{tau}{Kendall's tau.}
#' \item{rho}{Spearman's rho.}
#'
#' @references Shih J-H, Emura T (2019) Bivariate dependence measures and bivariate competing risks models under the generalized FGM copula, Statistical Papers, 60:1101-1118.
#' @export
#'
#' @examples
#' library(GFGM.copula)
#' Dependence.GFGM(3,2,0.75)

Dependence.GFGM = function(p,q,theta) {

  ### checking inputs ###
  if (p < 1) {stop("p must be greater than or equal to 1")}
  if (q <= 1) {stop("q must be greater than 1")}
  theta.UB = ((1+p*q)/(q-1))^(q-1)/p^q
  theta.LB = -min(1,(((1+p*q)/(q-1))^(q-1)/p^q)^2)
  if (theta > theta.UB | theta < theta.LB) {stop("theta is invalid")}

  theta.res = c(theta,theta.LB,theta.UB)
  tau.res   = c(8*(q/(2+p*q)*beta(2/p,q))^2*theta,
                8*(q/(2+p*q)*beta(2/p,q))^2*theta.LB,
                8*(q/(2+p*q)*beta(2/p,q))^2*theta.UB)
  rho.res   = 3*tau.res/2
  return(list(theta = theta.res,tau = tau.res,rho = rho.res))

}
