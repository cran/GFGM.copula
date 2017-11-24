#' Maximum likelihood estimation for bivariate dependent competing risks data under the generalized FGM copula with the marginal distributions approximated by splines
#'
#' @param t.event Vector of the observed failure times.
#' @param event1 Vector of the indicators for the failure cause 1.
#' @param event2 Vector of the indicators for the failure cause 2.
#' @param p Copula parameter that greater than or equal to 1.
#' @param q Copula parameter that greater than 1 (integer).
#' @param theta Copula parameter with restricted range.
#' @description Maximum likelihood estimation for bivariate dependent competing risks data under the generalized FGM copula with the marginal distributions approximated by splines.
#' @details The original paper is submitted for review.
#'
#' The copula parameter \code{q} is restricted to be a integer due to the binominal theorem.
#' The admissiable range of \code{theta} is given in \code{Dependence.GFGM}
#'
#' @return \item{g1}{Maximum likelihood estimator of the splines coefficients for the failure cause 1.}
#' \item{g2}{Maximum likelihood estimator of the splines coefficients for the failure cause 2.}
#'
#' @references Shih and Emura (2016) Bivariate dependence measures and bivariate competing risks models under the generalized FGM copula, Statistical Papers, doi: 10.1007/s00362-016-0865-5.
#' @references Shih and Emura (2017) Likelihood inference for bivariate latent failure time models with competing risks udner the generalized FGM copula (in re-submission, Computational Statistics).
#' @seealso \code{\link{Dependence.GFGM}}
#' @importFrom joint.Cox M.spline I.spline
#' @export
#'
#' @examples
#' con   = c(16,224,16,80,128,168,144,176,176,568,392,576,128,56,112,160,384,600,40,416,
#'           408,384,256,246,184,440,64,104,168,408,304,16,72,8,88,160,48,168,80,512,
#'           208,194,136,224,32,504,40,120,320,48,256,216,168,184,144,224,488,304,40,160,
#'           488,120,208,32,112,288,336,256,40,296,60,208,440,104,528,384,264,360,80,96,
#'           360,232,40,112,120,32,56,280,104,168,56,72,64,40,480,152,48,56,328,192,
#'           168,168,114,280,128,416,392,160,144,208,96,536,400,80,40,112,160,104,224,336,
#'           616,224,40,32,192,126,392,288,248,120,328,464,448,616,168,112,448,296,328,56,
#'           80,72,56,608,144,408,16,560,144,612,80,16,424,264,256,528,56,256,112,544,
#'           552,72,184,240,128,40,600,96,24,184,272,152,328,480,96,296,592,400,8,280,
#'           72,168,40,152,488,480,40,576,392,552,112,288,168,352,160,272,320,80,296,248,
#'           184,264,96,224,592,176,256,344,360,184,152,208,160,176,72,584,144,176)
#' uncon = c(368,136,512,136,472,96,144,112,104,104,344,246,72,80,312,24,128,304,16,320,
#'           560,168,120,616,24,176,16,24,32,232,32,112,56,184,40,256,160,456,48,24,
#'           200,72,168,288,112,80,584,368,272,208,144,208,114,480,114,392,120,48,104,272,
#'           64,112,96,64,360,136,168,176,256,112,104,272,320,8,440,224,280,8,56,216,
#'           120,256,104,104,8,304,240,88,248,472,304,88,200,392,168,72,40,88,176,216,
#'           152,184,400,424,88,152,184)
#' cen   = rep(630,44)
#'
#' t.event = c(con,uncon,cen)
#' event1  = c(rep(1,length(con)),rep(0,length(uncon)),rep(0,length(cen)))
#' event2  = c(rep(0,length(con)),rep(1,length(uncon)),rep(0,length(cen)))
#'
#' library(GFGM.copula)
#' MLE.GFGM.spline(t.event,event1,event2,3,2,0.75)

MLE.GFGM.spline = function(t.event,event1,event2,p,q,theta) {

  ### checking inputs ###
  n = length(t.event)
  if (length(t.event[t.event <= 0]) != 0) {stop("t.event must be positive")}
  if (length(event1) != n) {stop("the length of event1 is different from t.event")}
  if (length(event2) != n) {stop("the length of event2 is different from t.event")}
  if (length(event1[event1 == 0 | event1 == 1]) != n) {stop("elements in event1 must be either 0 or 1")}
  if (length(event2[event2 == 0 | event2 == 1]) != n) {stop("elements in event2 must be either 0 or 1")}
  temp.event = event1+event2
  if (length(temp.event[temp.event == 2]) != 0) {stop("event1 and event2 cannot be 1 simultaneously")}
  if (p < 1) {stop("p must be greater than or equal to 1")}
  if (q <= 1 | q != round(q)) {stop("q must be greater than 1 (integer)")}
  theta.UB = ((1+p*q)/(q-1))^(q-1)/p^q
  theta.LB = -min(1,(((1+p*q)/(q-1))^(q-1)/p^q)^2)
  if (theta > theta.UB | theta < theta.LB) {stop("theta is invalid")}

  ### functions ###
  nlm.logL.spline = function (par) {

    g11 = exp(par[1])
    g12 = exp(par[2])
    g13 = exp(par[3])
    g14 = exp(par[4])
    g15 = exp(par[5])
    g21 = exp(par[6])
    g22 = exp(par[7])
    g23 = exp(par[8])
    g24 = exp(par[9])
    g25 = exp(par[10])

    h1 = M.spline(t.event,tmin,tmax)%*%c(g11,g12,g13,g14,g15)
    H1 = I.spline(t.event,tmin,tmax)%*%c(g11,g12,g13,g14,g15)
    h2 = M.spline(t.event,tmin,tmax)%*%c(g21,g22,g23,g24,g25)
    H2 = I.spline(t.event,tmin,tmax)%*%c(g21,g22,g23,g24,g25)

    S1 = exp(-H1); F1 = 1-S1; f1 = S1*h1
    S2 = exp(-H2); F2 = 1-S2; f2 = S2*h2

    K1 = 0
    for (i in 0:q) {

      for (j in 0:(q-1)) {

        k1 = f1*F2^(p*i+1)*F1^(p*j)*(1-(1+p*q)*F1^p)
        K1 = K1+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*k1

      }

    }
    Sub.f1 = f1-F2*f1-theta*K1

    K2 = 0
    for (i in 0:q) {

      for (j in 0:(q-1)) {

        k2 = f2*F1^(p*i+1)*F2^(p*j)*(1-(1+p*q)*F2^p)
        K2 = K2+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*k2

      }

    }
    Sub.f2 = f2-F1*f2-theta*K2

    F.bar = 1-F1-F2+F1*F2*(1+theta*(1-F1^p)^q*(1-F2^p)^q)

    return(-sum(event1*log(Sub.f1))-sum(event2*log(Sub.f2))-sum((1-event1-event2)*log(F.bar)))

  }
  F1.spline = function(t.event) {return(1-exp(-I.spline(t.event,tmin,tmax)%*%g1))}
  F2.spline = function(t.event) {return(1-exp(-I.spline(t.event,tmin,tmax)%*%g2))}
  f1.spline = function(t.event) {return(M.spline(t.event,tmin,tmax)%*%g1*exp(-I.spline(t.event,tmin,tmax)%*%g1))}
  f2.spline = function(t.event) {return(M.spline(t.event,tmin,tmax)%*%g2*exp(-I.spline(t.event,tmin,tmax)%*%g2))}

  tmin = min(t.event)
  tmax = max(t.event)

  initial.par = rep(0,10)
  try.par = initial.par
  repeat {

    res.spline = try(nlm(nlm.logL.spline,try.par,iterlim = 5000),silent = TRUE)
    if (class(res.spline) == "try-error") {

      try.par = initial.par + runif(10,-10,0)
      next

    }
    if (res.spline$code == 1) {break}

    try.par = initial.par + runif(10,-10,0)

  }
  par = res.spline$estimate
  g1 = exp(par[1:5])
  g2 = exp(par[6:10])

  return(list(n = n,g1 = g1,g2 = g2))

}
