#' The Cramer-von Mises type statistics under the generalized FGM copula
#'
#' @param t.event Vector of the observed failure times.
#' @param event1 Vector of the indicators for the failure cause 1.
#' @param event2 Vector of the indicators for the failure cause 2.
#' @param Alpha Positive shape parameter for the Burr III margin (failure cause 1).
#' @param Beta Positive shape parameter for the Burr III margin (failure cause 2).
#' @param Gamma Common positive shape parameter for the Burr III margins.
#' @param g1 Splines coefficients for the failure cause 1.
#' @param g2 Splines coefficients for the failure cause 2.
#' @param p Copula parameter that greater than or equal to 1.
#' @param q Copula parameter that greater than 1 (integer).
#' @param theta Copula parameter with restricted range.
#' @param eta Location parameter with default value 0.
#' @param Sdist.plot Plot sub-distribution functions if \code{TRUE}.
#' @description Compute the Cramer-von Mises type statistics under the generalized FGM copula.
#' @details The copula parameter \code{q} is restricted to be a integer due to the binominal theorem.
#' The admissible range of \code{theta} is given in \code{Dependence.GFGM}.
#'
#' @return \item{S.overall}{Cramer-von Mises type statistic based on parametric and non-parametric estimators of sub-distribution functions for testing overall model.}
#' \item{S.GFGM}{Cramer-von Mises type statistic based on semi-parametric and non-parametric estimators of sub-distribution functions for testing the generalized FGM copula.}
#'
#' @references Shih J-H, Emura T (2018) Likelihood-based inference for bivariate latent failure time models with competing risks udner the generalized FGM copula, Computational Statistics, 33:1293-1323.
#' @seealso \code{\link{Dependence.GFGM}}, \code{\link{MLE.GFGM.BurrIII}}, \code{\link{MLE.GFGM.spline}}
#' @importFrom cmprsk cuminc timepoints
#' @importFrom utils globalVariables
#' @importFrom graphics legend lines plot
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
#' #res.BurrIII = MLE.GFGM.BurrIII(t.event,event1,event2,5000,3,2,0.75,eta = -71)
#' #Alpha = res.BurrIII$Alpha[1]
#' #Beta = res.BurrIII$Beta[1]
#' #Gamma = res.BurrIII$Gamma[1]
#' #res.spline = MLE.GFGM.spline(t.event,event1,event2,3,2,0.75)
#' #g1 = res.spline$g1
#' #g2 = res.spline$g2
#' #CvM.GFGM.BurrIII(t.event,event1,event2,Alpha,Beta,Gamma,g1,g2,3,2,0.75,eta = -71)

CvM.GFGM.BurrIII = function(t.event,event1,event2,Alpha,Beta,Gamma,g1,g2,p,q,
                            theta,eta = 0,Sdist.plot = TRUE) {

  event = event1+2*event2
  t.data = data.frame(t.event,event)
  t1.event = unique(sort(t.data$t.event[event == 1]))
  t2.event = unique(sort(t.data$t.event[event == 2]))
  res.cuminc = cuminc(t.event,event)

  tmin = min(t.event)
  tmax = max(t.event)

  Sdist.BIII1 = Sdist.GFGM.BurrIII(t1.event,Alpha,Beta,Gamma,p,q,theta,eta = eta)[,2]
  Sdist.BIII2 = Sdist.GFGM.BurrIII(t2.event,Alpha,Beta,Gamma,p,q,theta,eta = eta)[,3]
  Sdist.spline1 = Sdist.GFGM.spline(t1.event,g1,g2,tmin,tmax,p,q,theta)[,2]
  Sdist.spline2 = Sdist.GFGM.spline(t2.event,g1,g2,tmin,tmax,p,q,theta)[,3]
  Sdist.non1 = as.vector(timepoints(res.cuminc,t1.event)$est[1,])
  Sdist.non2 = as.vector(timepoints(res.cuminc,t2.event)$est[2,])

  S.overall = sum((Sdist.BIII1-Sdist.non1)^2)+sum((Sdist.BIII2-Sdist.non2)^2)
  S.GFGM    = sum((Sdist.spline1-Sdist.non1)^2)+sum((Sdist.spline2-Sdist.non2)^2)

  if (Sdist.plot) {

    t.grid = seq(tmin,tmax,length.out = 500)
    Sdist.BurrIII = Sdist.GFGM.BurrIII(t.grid,Alpha,Beta,Gamma,p,q,theta,eta)
    Sdist.spline  = Sdist.GFGM.spline(t.grid,g1,g2,tmin,tmax,p,q,theta)

    plot(res.cuminc$`1 1`$time,res.cuminc$`1 1`$est,xlab = "t",ylim = c(0,1),
         xlim = c(tmin,tmax),ylab = "Sub-distribution",type = "s")
    lines(Sdist.BurrIII[,1],Sdist.BurrIII[,2],type = "l",col = "blue")
    lines(Sdist.spline[,1],Sdist.spline[,2],type = "l",col = "red")
    legend("topleft",c("para 1","semi-para 1","non-para 1"),
           lty = 1,col = c("blue","red","black"))

    plot(res.cuminc$`1 2`$time,res.cuminc$`1 2`$est,xlab = "t",ylim = c(0,1),
         xlim = c(tmin,tmax),ylab = "Sub-distribution",type = "s")
    lines(Sdist.BurrIII[,1],Sdist.BurrIII[,3],type = "l",col = "blue",lty = 2)
    lines(Sdist.spline[,1],Sdist.spline[,3],type = "l",col = "red",lty = 2)
    legend("topleft",c("para 2","semi-para 2","non-para 2"),
           lty = 2,col = c("blue","red","black"))

  }
  return(list(S.overall = S.overall,S.GFGM = S.GFGM))

}
