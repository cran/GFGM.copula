#' Sub-distribution functions under the generalized FGM copula with the marginal distributions approximated by splines
#'
#' @param time Vector of times.
#' @param g1 Splines coefficients for the failure cause 1.
#' @param g2 Splines coefficients for the failure cause 2.
#' @param tmin Lower bound of times.
#' @param tmax upper bound of times.
#' @param p Copula parameter that greater than or equal to 1.
#' @param q Copula parameter that greater than 1 (integer).
#' @param theta Copula parameter with restricted range.
#' @description Sub-distribution functions under the generalized FGM copula with the marginal distributions approximated by splines.
#' @details The original paper is submitted for review.
#'
#' The splines coefficients \code{g1} and \code{g2} are usually computed by \code{MLE.GFGM.spline}.
#' The copula parameter \code{q} is restricted to be a integer due to the binominal theorem.
#' The admissible range of \code{theta} is given in \code{Dependence.GFGM}
#'
#' @return \item{time}{Failure times}
#' \item{Sdist.1}{Probability of an object fails due to the failure cause 1.}
#' \item{Sdist.2}{Probability of an object fails due to the failure cause 2.}
#'
#' @references Shih and Emura (2016) Bivariate dependence measures and bivariate competing risks models under the generalized FGM copula, Statistical Papers, doi: 10.1007/s00362-016-0865-5.
#' @references Shih and Emura (2018) Likelihood inference for bivariate latent failure time models with competing risks udner the generalized FGM copula (in re-submission, Computational Statistics).
#' @seealso \code{\link{MLE.GFGM.spline}}, \code{\link{Dependence.GFGM}}
#' @importFrom stats nlm
#' @importFrom joint.Cox M.spline I.spline
#' @export
#'
#' @examples
#' library(GFGM.copula)
#' Sdist.GFGM.spline(seq(1,5,1),rep(0.1,5),rep(0.1,5),1,5,3,2,0.75)

Sdist.GFGM.spline = function(time,g1,g2,tmin,tmax,p,q,theta) {

  ### checking inputs ###
  if (length(time[time <= 0]) != 0) {stop("time must be positive")}
  if (length(time[time < tmin]) != 0 | length(time[time > tmax]) != 0) {stop("time must be between tmin and tmax")}
  if (length(g1) != 5) {stop("g1 is invalid")}
  if (length(g2) != 5) {stop("g1 is invalid")}
  if (p < 1) {stop("p must be greater than or equal to 1")}
  if (q <= 1 | q != round(q)) {stop("q must be greater than 1 (integer)")}
  theta.UB = ((1+p*q)/(q-1))^(q-1)/p^q
  theta.LB = -min(1,(((1+p*q)/(q-1))^(q-1)/p^q)^2)
  if (theta > theta.UB | theta < theta.LB) {stop("theta is invalid")}

  ### functions ###
  F1.spline = function(time) {return(1-exp(-I.spline(time,tmin,tmax)%*%g1))}
  F2.spline = function(time) {return(1-exp(-I.spline(time,tmin,tmax)%*%g2))}
  f1.spline = function(time) {return(M.spline(time,tmin,tmax)%*%g1*exp(-I.spline(time,tmin,tmax)%*%g1))}
  f2.spline = function(time) {return(M.spline(time,tmin,tmax)%*%g2*exp(-I.spline(time,tmin,tmax)%*%g2))}

  Sdist.F1.spline = function(time) {

    M1 = F1.spline(time)
    m2 = function(time) {return(F2.spline(time)*f1.spline(time))}
    M2 = integrate(m2,tmin,time)$value
    K = 0
    for (i in 0:q) {

      for (j in 0:(q-1)) {

        k.func = function (time) {f1.spline(time)*F2.spline(time)^(p*i+1)*F1.spline(time)^(p*j)*(1-(1+p*q)*F1.spline(time)^p)}
        k = integrate(k.func,tmin,time)$value
        K = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*k

      }

    }

    return(M1-M2-theta*K)

  }
  Sdist.F2.spline = function(time) {

    M1 = F2.spline(time)
    m2 = function(time) {return(F1.spline(time)*f2.spline(time))}
    M2 = integrate(m2,tmin,time)$value
    K = 0
    for (i in 0:q) {

      for (j in 0:(q-1)) {

        k.func = function (time) {return(f2.spline(time)*F1.spline(time)^(p*i+1)*F2.spline(time)^(p*j)*(1-(1+p*q)*F2.spline(time)^p))}
        k = integrate(k.func,tmin,time)$value
        K = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*k

      }

    }

    return(M1-M2-theta*K)

  }

  Sdist.1 = sapply(time,Sdist.F1.spline)
  Sdist.2 = sapply(time,Sdist.F2.spline)

  return(cbind(time,Sdist.1,Sdist.2))

}
