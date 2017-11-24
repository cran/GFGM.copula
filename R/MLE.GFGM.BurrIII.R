#' Maximum likelihood estimation for bivariate dependent competing risks data under the generalized FGM copula with the Burr III margins
#'
#' @param t.event Vector of the observed failure times.
#' @param event1 Vector of the indicators for the failure cause 1.
#' @param event2 Vector of the indicators for the failure cause 2.
#' @param D Positive tunning parameter in the NR algorithm.
#' @param p Copula parameter that greater than or equal to 1.
#' @param q Copula parameter that greater than 1 (integer).
#' @param theta Copula parameter with restricted range.
#' @param eta Location parameter with default value 0.
#' @param Gamma.0 Initial guess for the common shape parameter gamma with default value 1.
#' @param epsilon.0 Positive tunning parameter in the NR algorithm with default value 1e-5.
#' @param epsilon.1 Positive tunning parameter in the NR algorithm with default value 1e-10.
#' @param epsilon.2 Positive tunning parameter in the NR algorithm with default value 1e-300.
#' @param r.1 Positive tunning parameter in the NR algorithm with default value 1.
#' @param r.2 Positive tunning parameter in the NR algorithm with default value 1.
#' @param r.3 Positive tunning parameter in the NR algorithm with default value 1.
#' @description Maximum likelihood estimation for bivariate dependent competing risks data under the generalized FGM copula with the Burr III margins.
#' @details The original paper is submitted for review.
#'
#' The copula parameter \code{q} is restricted to be a integer due to the binominal theorem.
#' The admissiable range of \code{theta} is given in \code{Dependence.GFGM}
#'
#' @return \item{n}{Sample size.}
#' \item{count}{Iteration number.}
#' \item{random}{Randomization number.}
#' \item{Alpha}{Positive shape parameter for the Burr III margin (failure cause 1).}
#' \item{Beta}{Positive shape parameter for the Burr III margin (failure cause 2).}
#' \item{Gamma}{Common shape parameter for the Burr III margins.}
#' \item{MeanX}{Mean lifetime due to failure cause 1.}
#' \item{MeanY}{Mean lifetime due to failure cause 2.}
#' \item{logL}{Log-likelihood value under the fitted model.}
#'
#' @references Shih and Emura (2016) Bivariate dependence measures and bivariate competing risks models under the generalized FGM copula, Statistical Papers, doi: 10.1007/s00362-016-0865-5.
#' @references Shih and Emura (2017) Likelihood inference for bivariate latent failure time models with competing risks udner the generalized FGM copula (in re-submission, Computational Statistics).
#' @seealso \code{\link{Dependence.GFGM}}
#' @importFrom stats qnorm runif
#' @importFrom utils globalVariables
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
#' MLE.GFGM.BurrIII(t.event,event1,event2,5000,3,2,0.75,eta = -71)

MLE.GFGM.BurrIII = function(t.event,event1,event2,D,p,q,theta,eta = 0,Gamma.0 = 1,
                            epsilon.0 = 1e-5,epsilon.1 = 1e-10,epsilon.2 = 1e-300,
                            r.1 = 1,r.2 = 1,r.3 = 1) {

  ### checking inputs ###
  n = length(t.event)
  if (length(t.event[t.event <= 0]) != 0) {stop("t.event must be positive")}
  if (length(event1) != n) {stop("the length of event1 is different from t.event")}
  if (length(event2) != n) {stop("the length of event2 is different from t.event")}
  if (length(event1[event1 == 0 | event1 == 1]) != n) {stop("elements in event1 must be either 0 or 1")}
  if (length(event2[event2 == 0 | event2 == 1]) != n) {stop("elements in event2 must be either 0 or 1")}
  temp.event = event1+event2
  if (length(temp.event[temp.event == 2]) != 0) {stop("event1 and event2 cannot be 1 simultaneously")}
  if (D <= 0) {stop("D must be positive")}
  if (epsilon.0 <= 0) {stop("epsilon.0 must be positive")}
  if (epsilon.1 <= 0) {stop("epsilon.1 must be positive")}
  if (epsilon.2 <= 0) {stop("epsilon.2 must be positive")}
  if (r.1 <= 0) {stop("r.1 must be positive")}
  if (r.2 <= 0) {stop("r.2 must be positive")}
  if (r.3 <= 0) {stop("r.3 must be positive")}
  if (Gamma.0 <= 0) {stop("Gamma.0 must be positive")}
  if (p < 1) {stop("p must be greater than or equal to 1")}
  if (q <= 1 | q != round(q)) {stop("q must be greater than 1 (integer)")}
  theta.UB = ((1+p*q)/(q-1))^(q-1)/p^q
  theta.LB = -min(1,(((1+p*q)/(q-1))^(q-1)/p^q)^2)
  if (theta > theta.UB | theta < theta.LB) {stop("theta is invalid")}

  ### functions ###
  f1_func    = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    M1 = Alpha*Gamma*t.event^(-Gamma-1)*H^(Alpha+1)
    M2 = Alpha*Gamma*t.event^(-Gamma-1)*H^(Alpha+Beta+1)
    K = 0
    for (i in 0:q) {

      for (j in 0:(q-1)) {

        K1 = H^(Alpha*(p*j+1)+Beta*(p*i+1)+1)
        K  = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*K1*(1-(1+p*q)*H^(Alpha*p))*t.event^(-Gamma-1)

      }

    }

    return(M1-M2-theta*Alpha*Gamma*K)

  }
  f1_A1_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    M1 = Gamma*t.event^(-Gamma-1)*H^(Alpha+1)*(1+Alpha*log(H))
    M2 = Gamma*t.event^(-Gamma-1)*H^(Alpha+Beta+1)*(1+Alpha*log(H))
    K = 0
    for (i in 0:q) {

      for (j in 0:(q-1)) {

        K1 = Gamma*t.event^(-Gamma-1)*H^(Alpha*(p*j+1)+Beta*(p*i+1)+1)*(1+Alpha*(p*j+1)*log(H))
        K2 = Gamma*t.event^(-Gamma-1)*H^(Alpha*(p*j+p+1)+Beta*(p*i+1)+1)*(1+Alpha*(p*j+p+1)*log(H))
        K  = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*(K1-(1+p*q)*K2)

      }

    }

    return(M1-M2-theta*K)

  }
  f1_A2_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    M1 = Gamma*t.event^(-Gamma-1)*H^(Alpha+1)*log(H)*(2+Alpha*log(H))
    M2 = Gamma*t.event^(-Gamma-1)*H^(Alpha+Beta+1)*log(H)*(2+Alpha*log(H))
    K = 0
    for (i in 0:q) {

      for (j in 0:(q-1)) {

        K1 = Gamma*t.event^(-Gamma-1)*H^(Alpha*(p*j+1)+Beta*(p*i+1)+1)*(p*j+1)*log(H)*(2+Alpha*(p*j+1)*log(H))
        K2 = Gamma*t.event^(-Gamma-1)*H^(Alpha*(p*j+p+1)+Beta*(p*i+1)+1)*(p*j+p+1)*log(H)*(2+Alpha*(p*j+p+1)*log(H))
        K  = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*(K1-(1+p*q)*K2)

      }

    }

    return(M1-M2-theta*K)

  }
  f1_B1_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    M2 = Alpha*Gamma*t.event^(-Gamma-1)*H^(Alpha+Beta+1)*log(H)
    K = 0
    for (i in 0:q) {

      for (j in 0:(q-1)) {

        K1 = Gamma*t.event^(-Gamma-1)*H^(Alpha*(p*j+1)+Beta*(p*i+1)+1)*(p*i+1)*log(H)
        K2 = Gamma*t.event^(-Gamma-1)*H^(Alpha*(p*j+p+1)+Beta*(p*i+1)+1)*(p*i+1)*log(H)
        K  = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*(K1-(1+p*q)*K2)

      }

    }

    return(-M2-theta*Alpha*K)

  }
  f1_B2_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    M2 = Alpha*Gamma*t.event^(-Gamma-1)*H^(Alpha+Beta+1)*log(H)^2
    K = 0
    for (i in 0:q) {

      for (j in 0:(q-1)) {

        K1 = Gamma*t.event^(-Gamma-1)*H^(Alpha*(p*j+1)+Beta*(p*i+1)+1)*(p*i+1)^2*log(H)^2
        K2 = Gamma*t.event^(-Gamma-1)*H^(Alpha*(p*j+p+1)+Beta*(p*i+1)+1)*(p*i+1)^2*log(H)^2
        K  = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*(K1-(1+p*q)*K2)

      }

    }

    return(-M2-theta*Alpha*K)

  }
  f1_G1_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    M1 = Alpha*t.event^(-Gamma-1)*H^(Alpha+1)*(1-Gamma*log(t.event)+(Alpha+1)*Gamma*t.event^(-Gamma)*log(t.event)*H)
    M2 = Alpha*t.event^(-Gamma-1)*H^(Alpha+Beta+1)*(1-Gamma*log(t.event)+(Alpha+Beta+1)*Gamma*t.event^(-Gamma)*log(t.event)*H)
    K = 0
    for (i in 0:q) {

      for (j in 0:(q-1)) {

        k1 = (Alpha*(p*j+1)+Beta*(p*i+1)+1)*Gamma*t.event^(-Gamma)*log(t.event)*H
        k2 = (Alpha*(p*j+p+1)+Beta*(p*i+1)+1)*Gamma*t.event^(-Gamma)*log(t.event)*H
        K1 = t.event^(-Gamma-1)*H^(Alpha*(p*j+1)+Beta*(p*i+1)+1)*(1-Gamma*log(t.event)+k1)
        K2 = t.event^(-Gamma-1)*H^(Alpha*(p*j+p+1)+Beta*(p*i+1)+1)*(1-Gamma*log(t.event)+k2)
        K  = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*(K1-(1+p*q)*K2)

      }

    }

    return(M1-M2-theta*Alpha*K)

  }
  f1_G2_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    m11 = (1-Gamma*log(t.event)+(Alpha+1)*Gamma*t.event^(-Gamma)*log(t.event)*H)*(-log(t.event)+(Alpha+1)*t.event^(-Gamma)*log(t.event)*H)
    m21 = (1-Gamma*log(t.event)+(Alpha+Beta+1)*Gamma*t.event^(-Gamma)*log(t.event)*H)*(-log(t.event)+(Alpha+Beta+1)*t.event^(-Gamma)*log(t.event)*H)
    m12 = (Alpha+1)*log(t.event)*t.event^(-Gamma)*H*(1-Gamma*log(t.event)+log(t.event)*Gamma*t.event^(-Gamma)*H)-log(t.event)
    m22 = (Alpha+Beta+1)*log(t.event)*t.event^(-Gamma)*H*(1-Gamma*log(t.event)+log(t.event)*Gamma*t.event^(-Gamma)*H)-log(t.event)
    M1  = Alpha*t.event^(-Gamma-1)*H^(Alpha+1)*(m11+m12)
    M2  = Alpha*t.event^(-Gamma-1)*H^(Alpha+Beta+1)*(m21+m22)
    K = 0
    for (i in 0:q) {

      for (j in 0:(q-1)) {

        k11 = 1-Gamma*log(t.event)+(Alpha*(p*j+1)+Beta*(p*i+1)+1)*Gamma*t.event^(-Gamma)*log(t.event)*H
        k21 = 1-Gamma*log(t.event)+(Alpha*(p*j+p+1)+Beta*(p*i+1)+1)*Gamma*t.event^(-Gamma)*log(t.event)*H
        k12 = -log(t.event)+(Alpha*(p*j+1)+Beta*(p*i+1)+1)*t.event^(-Gamma)*log(t.event)*H
        k22 = -log(t.event)+(Alpha*(p*j+p+1)+Beta*(p*i+1)+1)*t.event^(-Gamma)*log(t.event)*H
        k13 = (Alpha*(p*j+1)+Beta*(p*i+1)+1)*log(t.event)*t.event^(-Gamma)*H*(1-Gamma*log(t.event)+log(t.event)*Gamma*t.event^(-Gamma)*H)-log(t.event)
        k23 = (Alpha*(p*j+p+1)+Beta*(p*i+1)+1)*log(t.event)*t.event^(-Gamma)*H*(1-Gamma*log(t.event)+log(t.event)*Gamma*t.event^(-Gamma)*H)-log(t.event)
        K1  = t.event^(-Gamma-1)*H^(Alpha*(p*j+1)+Beta*(p*i+1)+1)*(k11*k12+k13)
        K2  = t.event^(-Gamma-1)*H^(Alpha*(p*j+p+1)+Beta*(p*i+1)+1)*(k21*k22+k23)
        K   = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*(K1-(1+p*q)*K2)

      }

    }

    return(M1-M2-theta*Alpha*K)

  }
  f1_AB_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    M2 = Gamma*t.event^(-Gamma-1)*H^(Alpha+Beta+1)*(1+Alpha*log(H))*log(H)
    K = 0
    for (i in 0:q) {

      for (j in 0:(q-1)) {

        k1 = (1+Alpha*(p*j+1)*log(H))*(p*i+1)*log(H)
        k2 = (1+Alpha*(p*j+p+1)*log(H))*(p*i+1)*log(H)
        K1 = Gamma*t.event^(-Gamma-1)*H^(Alpha*(p*j+1)+Beta*(p*i+1)+1)*k1
        K2 = Gamma*t.event^(-Gamma-1)*H^(Alpha*(p*j+p+1)+Beta*(p*i+1)+1)*k2
        K  = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*(K1-(1+p*q)*K2)

      }

    }

    return(-M2-theta*K)

  }
  f1_AG_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    m11 = m21 = (Alpha*Gamma*t.event^(-Gamma)*log(t.event)*H)/(1+Alpha*log(H))
    m12 = 1-Gamma*log(t.event)+(Alpha+1)*Gamma*t.event^(-Gamma)*log(t.event)*H
    m22 = 1-Gamma*log(t.event)+(Alpha+Beta+1)*Gamma*t.event^(-Gamma)*log(t.event)*H
    M1 = t.event^(-Gamma-1)*H^(Alpha+1)*(1+Alpha*log(H))*(m11+m12)
    M2 = t.event^(-Gamma-1)*H^(Alpha+Beta+1)*(1+Alpha*log(H))*(m21+m22)
    K = 0
    for (i in 0:q) {

      for (j in 0:(q-1)) {

        k11 = (Alpha*(p*j+1)*Gamma*t.event^(-Gamma)*log(t.event)*H)/(1+Alpha*(p*j+1)*log(H))
        k21 = (Alpha*(p*j+p+1)*Gamma*t.event^(-Gamma)*log(t.event)*H)/(1+Alpha*(p*j+p+1)*log(H))
        k12 = 1-Gamma*log(t.event)+(Alpha*(p*j+1)+Beta*(p*i+1)+1)*Gamma*t.event^(-Gamma)*log(t.event)*H
        k22 = 1-Gamma*log(t.event)+(Alpha*(p*j+p+1)+Beta*(p*i+1)+1)*Gamma*t.event^(-Gamma)*log(t.event)*H
        K1  = t.event^(-Gamma-1)*H^(Alpha*(p*j+1)+Beta*(p*i+1)+1)*(1+Alpha*(p*j+1)*log(H))*(k11+k12)
        K2  = t.event^(-Gamma-1)*H^(Alpha*(p*j+p+1)+Beta*(p*i+1)+1)*(1+Alpha*(p*j+p+1)*log(H))*(k21+k22)
        K   = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*(K1-(1+p*q)*K2)

      }

    }

    return(M1-M2-theta*K)

  }
  f1_BG_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    m21 = 1-Gamma*log(t.event)+(Alpha+Beta+1)*Gamma*t.event^(-Gamma)*log(t.event)*H
    m22 = k12 = k22 = Gamma*t.event^(-Gamma)*log(t.event)*H/log(H)
    M2  = Alpha*t.event^(-Gamma-1)*H^(Alpha+Beta+1)*log(H)*(m21+m22)
    K = 0
    for (i in 0:q) {

      for (j in 0:(q-1)) {

        k11 = 1-Gamma*log(t.event)+(Alpha*(p*j+1)+Beta*(p*i+1)+1)*Gamma*t.event^(-Gamma)*log(t.event)*H
        k21 = 1-Gamma*log(t.event)+(Alpha*(p*j+p+1)+Beta*(p*i+1)+1)*Gamma*t.event^(-Gamma)*log(t.event)*H
        K1  = t.event^(-Gamma-1)*H^(Alpha*(p*j+1)+Beta*(p*i+1)+1)*(p*i+1)*log(H)*(k11+k12)
        K2  = t.event^(-Gamma-1)*H^(Alpha*(p*j+p+1)+Beta*(p*i+1)+1)*(p*i+1)*log(H)*(k21+k22)
        K   = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*(K1-(1+p*q)*K2)

      }

    }

    return(-M2-theta*Alpha*K)

  }
  f2_func    = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    M1 = Beta*Gamma*t.event^(-Gamma-1)*H^(Beta+1)
    M2 = Beta*Gamma*t.event^(-Gamma-1)*H^(Alpha+Beta+1)
    K = 0
    for (i in 0:q) {

      for (j in 0:(q-1)) {

        K1 = H^(Alpha*(p*i+1)+Beta*(p*j+1)+1)
        K  = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*K1*(1-(1+p*q)*H^(Beta*p))*t.event^(-Gamma-1)

      }

    }

    return(M1-M2-theta*Beta*Gamma*K)

  }
  f2_A1_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    M2 = Beta*Gamma*t.event^(-Gamma-1)*H^(Alpha+Beta+1)*log(H)
    K = 0
    for (i in 0:q) {

      for (j in 0:(q-1)) {

        K1 = Gamma*t.event^(-Gamma-1)*H^(Alpha*(p*i+1)+Beta*(p*j+1)+1)*(p*i+1)*log(H)
        K2 = Gamma*t.event^(-Gamma-1)*H^(Alpha*(p*i+1)+Beta*(p*j+p+1)+1)*(p*i+1)*log(H)
        K  = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*(K1-(1+p*q)*K2)

      }

    }

    return(-M2-theta*Beta*K)

  }
  f2_A2_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    M2 = Beta*Gamma*t.event^(-Gamma-1)*H^(Alpha+Beta+1)*log(H)^2
    K = 0
    for (i in 0:q) {

      for (j in 0:(q-1)) {

        K1 = Gamma*t.event^(-Gamma-1)*H^(Alpha*(p*i+1)+Beta*(p*j+1)+1)*(p*i+1)^2*log(H)^2
        K2 = Gamma*t.event^(-Gamma-1)*H^(Alpha*(p*i+1)+Beta*(p*j+p+1)+1)*(p*i+1)^2*log(H)^2
        K  = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*(K1-(1+p*q)*K2)

      }

    }

    return(-M2-theta*Beta*K)

  }
  f2_B1_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    M1 = Gamma*t.event^(-Gamma-1)*H^(Beta+1)*(1+Beta*log(H))
    M2 = Gamma*t.event^(-Gamma-1)*H^(Alpha+Beta+1)*(1+Beta*log(H))
    K = 0
    for (i in 0:q) {

      for (j in 0:(q-1)) {

        K1 = Gamma*t.event^(-Gamma-1)*H^(Alpha*(p*i+1)+Beta*(p*j+1)+1)*(1+Beta*(p*j+1)*log(H))
        K2 = Gamma*t.event^(-Gamma-1)*H^(Alpha*(p*i+1)+Beta*(p*j+p+1)+1)*(1+Beta*(p*j+p+1)*log(H))
        K  = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*(K1-(1+p*q)*K2)

      }

    }

    return(M1-M2-theta*K)

  }
  f2_B2_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    M1 = Gamma*t.event^(-Gamma-1)*H^(Beta+1)*log(H)*(2+Beta*log(H))
    M2 = Gamma*t.event^(-Gamma-1)*H^(Alpha+Beta+1)*log(H)*(2+Beta*log(H))
    K = 0
    for (i in 0:q) {

      for (j in 0:(q-1)) {

        K1 = Gamma*t.event^(-Gamma-1)*H^(Alpha*(p*i+1)+Beta*(p*j+1)+1)*(p*j+1)*log(H)*(2+Beta*(p*j+1)*log(H))
        K2 = Gamma*t.event^(-Gamma-1)*H^(Alpha*(p*i+1)+Beta*(p*j+p+1)+1)*(p*j+p+1)*log(H)*(2+Beta*(p*j+p+1)*log(H))
        K  = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*(K1-(1+p*q)*K2)

      }

    }

    return(M1-M2-theta*K)

  }
  f2_G1_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    M1 = Beta*t.event^(-Gamma-1)*H^(Beta+1)*(1-Gamma*log(t.event)+(Beta+1)*Gamma*t.event^(-Gamma)*log(t.event)*H)
    M2 = Beta*t.event^(-Gamma-1)*H^(Alpha+Beta+1)*(1-Gamma*log(t.event)+(Alpha+Beta+1)*Gamma*t.event^(-Gamma)*log(t.event)*H)
    K = 0
    for (i in 0:q) {

      for (j in 0:(q-1)) {

        k1 = (Alpha*(p*i+1)+Beta*(p*j+1)+1)*Gamma*t.event^(-Gamma)*log(t.event)*H
        k2 = (Alpha*(p*i+1)+Beta*(p*j+p+1)+1)*Gamma*t.event^(-Gamma)*log(t.event)*H
        K1 = t.event^(-Gamma-1)*H^(Alpha*(p*i+1)+Beta*(p*j+1)+1)*(1-Gamma*log(t.event)+k1)
        K2 = t.event^(-Gamma-1)*H^(Alpha*(p*i+1)+Beta*(p*j+p+1)+1)*(1-Gamma*log(t.event)+k2)
        K  = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*(K1-(1+p*q)*K2)

      }

    }

    return(M1-M2-theta*Beta*K)

  }
  f2_G2_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    m11 = (1-Gamma*log(t.event)+(Beta+1)*Gamma*t.event^(-Gamma)*log(t.event)*H)*(-log(t.event)+(Beta+1)*t.event^(-Gamma)*log(t.event)*H)
    m21 = (1-Gamma*log(t.event)+(Alpha+Beta+1)*Gamma*t.event^(-Gamma)*log(t.event)*H)*(-log(t.event)+(Alpha+Beta+1)*t.event^(-Gamma)*log(t.event)*H)
    m12 = (Beta+1)*log(t.event)*t.event^(-Gamma)*H*(1-Gamma*log(t.event)+log(t.event)*Gamma*t.event^(-Gamma)*H)-log(t.event)
    m22 = (Alpha+Beta+1)*log(t.event)*t.event^(-Gamma)*H*(1-Gamma*log(t.event)+log(t.event)*Gamma*t.event^(-Gamma)*H)-log(t.event)
    M1  = Beta*t.event^(-Gamma-1)*H^(Beta+1)*(m11+m12)
    M2  = Beta*t.event^(-Gamma-1)*H^(Alpha+Beta+1)*(m21+m22)
    K = 0
    for (i in 0:q) {

      for (j in 0:(q-1)) {

        k11 = 1-Gamma*log(t.event)+(Alpha*(p*i+1)+Beta*(p*j+1)+1)*Gamma*t.event^(-Gamma)*log(t.event)*H
        k21 = 1-Gamma*log(t.event)+(Alpha*(p*i+1)+Beta*(p*j+p+1)+1)*Gamma*t.event^(-Gamma)*log(t.event)*H
        k12 = -log(t.event)+(Alpha*(p*i+1)+Beta*(p*j+1)+1)*t.event^(-Gamma)*log(t.event)*H
        k22 = -log(t.event)+(Alpha*(p*i+1)+Beta*(p*j+p+1)+1)*t.event^(-Gamma)*log(t.event)*H
        k13 = (Alpha*(p*i+1)+Beta*(p*j+1)+1)*log(t.event)*t.event^(-Gamma)*H*(1-Gamma*log(t.event)+log(t.event)*Gamma*t.event^(-Gamma)*H)-log(t.event)
        k23 = (Alpha*(p*i+1)+Beta*(p*j+p+1)+1)*log(t.event)*t.event^(-Gamma)*H*(1-Gamma*log(t.event)+log(t.event)*Gamma*t.event^(-Gamma)*H)-log(t.event)
        K1  = t.event^(-Gamma-1)*H^(Alpha*(p*i+1)+Beta*(p*j+1)+1)*(k11*k12+k13)
        K2  = t.event^(-Gamma-1)*H^(Alpha*(p*i+1)+Beta*(p*j+p+1)+1)*(k21*k22+k23)
        K   = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*(K1-(1+p*q)*K2)

      }

    }

    return(M1-M2-theta*Beta*K)

  }
  f2_AB_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    M2 = Gamma*t.event^(-Gamma-1)*H^(Alpha+Beta+1)*(1+Beta*log(H))*log(H)
    K = 0
    for (i in 0:q) {

      for (j in 0:(q-1)) {

        k1 = (1+Beta*(p*j+1)*log(H))*(p*i+1)*log(H)
        k2 = (1+Beta*(p*j+p+1)*log(H))*(p*i+1)*log(H)
        K1 = Gamma*t.event^(-Gamma-1)*H^(Alpha*(p*i+1)+Beta*(p*j+1)+1)*k1
        K2 = Gamma*t.event^(-Gamma-1)*H^(Alpha*(p*i+1)+Beta*(p*j+p+1)+1)*k2
        K  = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*(K1-(1+p*q)*K2)

      }

    }

    return(-M2-theta*K)

  }
  f2_AG_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    m21 = 1-Gamma*log(t.event)+(Alpha+Beta+1)*Gamma*t.event^(-Gamma)*log(t.event)*H
    m22 = k12 = k22 = Gamma*t.event^(-Gamma)*log(t.event)*H/log(H)
    M2  = Beta*t.event^(-Gamma-1)*H^(Alpha+Beta+1)*log(H)*(m21+m22)
    K = 0
    for (i in 0:q) {

      for (j in 0:(q-1)) {

        k11 = 1-Gamma*log(t.event)+(Alpha*(p*i+1)+Beta*(p*j+1)+1)*Gamma*t.event^(-Gamma)*log(t.event)*H
        k21 = 1-Gamma*log(t.event)+(Alpha*(p*i+1)+Beta*(p*j+p+1)+1)*Gamma*t.event^(-Gamma)*log(t.event)*H
        K1  = t.event^(-Gamma-1)*H^(Alpha*(p*i+1)+Beta*(p*j+1)+1)*(p*i+1)*log(H)*(k11+k12)
        K2  = t.event^(-Gamma-1)*H^(Alpha*(p*i+1)+Beta*(p*j+p+1)+1)*(p*i+1)*log(H)*(k21+k22)
        K   = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*(K1-(1+p*q)*K2)

      }

    }

    return(-M2-theta*Beta*K)

  }
  f2_BG_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    m11 = m21 = (Beta*Gamma*t.event^(-Gamma)*log(t.event)*H)/(1+Beta*log(H))
    m12 = 1-Gamma*log(t.event)+(Beta+1)*Gamma*t.event^(-Gamma)*log(t.event)*H
    m22 = 1-Gamma*log(t.event)+(Alpha+Beta+1)*Gamma*t.event^(-Gamma)*log(t.event)*H
    M1 = t.event^(-Gamma-1)*H^(Beta+1)*(1+Beta*log(H))*(m11+m12)
    M2 = t.event^(-Gamma-1)*H^(Alpha+Beta+1)*(1+Beta*log(H))*(m21+m22)
    K = 0
    for (i in 0:q) {

      for (j in 0:(q-1)) {

        k11 = (Beta*(p*j+1)*Gamma*t.event^(-Gamma)*log(t.event)*H)/(1+Beta*(p*j+1)*log(H))
        k21 = (Beta*(p*j+p+1)*Gamma*t.event^(-Gamma)*log(t.event)*H)/(1+Beta*(p*j+p+1)*log(H))
        k12 = 1-Gamma*log(t.event)+(Alpha*(p*i+1)+Beta*(p*j+1)+1)*Gamma*t.event^(-Gamma)*log(t.event)*H
        k22 = 1-Gamma*log(t.event)+(Alpha*(p*i+1)+Beta*(p*j+p+1)+1)*Gamma*t.event^(-Gamma)*log(t.event)*H
        K1  = t.event^(-Gamma-1)*H^(Alpha*(p*i+1)+Beta*(p*j+1)+1)*(1+Beta*(p*j+1)*log(H))*(k11+k12)
        K2  = t.event^(-Gamma-1)*H^(Alpha*(p*i+1)+Beta*(p*j+p+1)+1)*(1+Beta*(p*j+p+1)*log(H))*(k21+k22)
        K   = K+choose(q,i)*choose(q-1,j)*(-1)^(i+j)*(K1-(1+p*q)*K2)

      }

    }

    M1-M2-theta*K

  }
  Fb_func    = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    K = 0
    for (i in 0:q) {

      for (j in 0:q) {

        K = K+choose(q,i)*choose(q,j)*(-1)^(i+j)*H^(Alpha*(p*j+1)+Beta*(p*i+1))

      }

    }

    return(1-H^Alpha-H^Beta+H^(Alpha+Beta)+theta*K)

  }
  Fb_A1_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    K = 0
    for (i in 0:q) {

      for (j in 0:q) {

        K = K+choose(q,i)*choose(q,j)*(-1)^(i+j)*(p*j+1)*log(H)*H^(Alpha*(p*j+1)+Beta*(p*i+1))

      }

    }

    return(-log(H)*H^Alpha+log(H)*H^(Alpha+Beta)+theta*K)

  }
  Fb_A2_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    K = 0
    for (i in 0:q) {

      for (j in 0:q) {

        k = H^(Alpha*(p*j+1)+Beta*(p*i+1))
        K = K+choose(q,i)*choose(q,j)*(-1)^(i+j)*(p*j+1)^2*log(H)^2*H^(Alpha*(p*j+1)+Beta*(p*i+1))

      }

    }

    return(-log(H)^2*H^Alpha+log(H)^2*H^(Alpha+Beta)+theta*K)

  }
  Fb_B1_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    K = 0
    for (i in 0:q) {

      for (j in 0:q) {

        K = K+choose(q,i)*choose(q,j)*(-1)^(i+j)*(p*i+1)*log(H)*H^(Alpha*(p*j+1)+Beta*(p*i+1))

      }

    }

    return(-log(H)*H^Beta+log(H)*H^(Alpha+Beta)+theta*K)

  }
  Fb_B2_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    K = 0
    for (i in 0:q) {

      for (j in 0:q) {

        k = H^(Alpha*(p*j+1)+Beta*(p*i+1))
        K = K+choose(q,i)*choose(q,j)*(-1)^(i+j)*(p*i+1)^2*log(H)^2*H^(Alpha*(p*j+1)+Beta*(p*i+1))

      }

    }

    return(-log(H)^2*H^Beta+log(H)^2*H^(Alpha+Beta)+theta*K)

  }
  Fb_G1_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    M1 = Alpha*log(t.event)*t.event^(-Gamma)*H^(Alpha+1)
    M2 = Beta*log(t.event)*t.event^(-Gamma)*H^(Beta+1)
    M3 = (Alpha+Beta)*log(t.event)*t.event^(-Gamma)*H^(Alpha+Beta+1)
    K = 0
    for (i in 0:q) {

      for (j in 0:q) {

        k = H^(Alpha*(p*j+1)+Beta*(p*i+1)+1)
        K = K+choose(q,i)*choose(q,j)*(-1)^(i+j)*(Alpha*(p*j+1)+Beta*(p*i+1))*log(t.event)*t.event^(-Gamma)*k

      }

    }

    return(-M1-M2+M3+theta*K)

  }
  Fb_G2_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    M1 = Alpha*t.event^(-Gamma)*log(t.event)*H^(Alpha+1)*(-log(t.event)+(Alpha+1)*t.event^(-Gamma)*log(t.event)*H)
    M2 = Beta*t.event^(-Gamma)*log(t.event)*H^(Beta+1)*(-log(t.event)+(Beta+1)*t.event^(-Gamma)*log(t.event)*H)
    M3 = (Alpha+Beta)*t.event^(-Gamma)*log(t.event)*H^(Alpha+Beta+1)*(-log(t.event)+(Alpha+Beta+1)*t.event^(-Gamma)*log(t.event)*H)
    K = 0
    for (i in 0:q) {

      for (j in 0:q) {

        k = (Alpha*(p*j+1)+Beta*(p*i+1))*t.event^(-Gamma)*log(t.event)*H^(Alpha*(p*j+1)+Beta*(p*i+1)+1)
        K = K+choose(q,i)*choose(q,j)*(-1)^(i+j)*k*(-log(t.event)+(Alpha*(p*j+1)+Beta*(p*i+1)+1)*t.event^(-Gamma)*log(t.event)*H)

      }

    }

    return(-M1-M2+M3+theta*K)

  }
  Fb_AB_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    K = 0
    for (i in 0:q) {

      for (j in 0:q) {

        K = K+choose(q,i)*choose(q,j)*(-1)^(i+j)*(p*i+1)*(p*j+1)*log(H)^2*H^(Alpha*(p*j+1)+Beta*(p*i+1))

      }

    }

    return(log(H)^2*H^(Alpha+Beta)+theta*K)

  }
  Fb_AG_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    M1 = t.event^(-Gamma)*log(t.event)*H^(Alpha+1)*(1+Alpha*log(H))
    M2 = t.event^(-Gamma)*log(t.event)*H^(Alpha+Beta+1)*(1+(Alpha+Beta)*log(H))
    K = 0
    for (i in 0:q) {

      for (j in 0:q) {

        k = (p*j+1)*(1+(Alpha*(p*j+1)+Beta*(p*i+1))*log(H))
        K = K+choose(q,i)*choose(q,j)*(-1)^(i+j)*t.event^(-Gamma)*log(t.event)*H^(Alpha*(p*j+1)+Beta*(p*i+1)+1)*k

      }

    }

    return(-M1+M2+theta*K)

  }
  Fb_BG_func = function(par) {

    Alpha = exp(par[1])
    Beta  = exp(par[2])
    Gamma = exp(par[3])

    H = (1+t.event^(-Gamma))^(-1)
    M1 = t.event^(-Gamma)*log(t.event)*H^(Beta+1)*(1+Beta*log(H))
    M2 = t.event^(-Gamma)*log(t.event)*H^(Alpha+Beta+1)*(1+(Alpha+Beta)*log(H))
    K = 0
    for (i in 0:q) {

      for (j in 0:q) {

        k = (p*i+1)*(1+(Alpha*(p*j+1)+Beta*(p*i+1))*log(H))
        K = K+choose(q,i)*choose(q,j)*(-1)^(i+j)*t.event^(-Gamma)*log(t.event)*H^(Alpha*(p*j+1)+Beta*(p*i+1)+1)*k

      }

    }

    return(-M1+M2+theta*K)

  }
  logL_func  = function(par) {

    return(sum(event1*log(f1_func(par)))+sum(event2*log(f2_func(par)))+sum((1-event1-event2)*log(Fb_func(par))))

  }
  score_func = function(par) {

    log_f1_A1 = f1_A1_func(par)*exp(par[1])/f1_func(par)
    log_f2_A1 = f2_A1_func(par)*exp(par[1])/f2_func(par)
    log_Fb_A1 = Fb_A1_func(par)*exp(par[1])/Fb_func(par)

    s1 = sum(event1*log_f1_A1)+sum(event2*log_f2_A1)+sum((1-event1-event2)*log_Fb_A1)

    log_f1_B1 = f1_B1_func(par)*exp(par[2])/f1_func(par)
    log_f2_B1 = f2_B1_func(par)*exp(par[2])/f2_func(par)
    log_Fb_B1 = Fb_B1_func(par)*exp(par[2])/Fb_func(par)

    s2 = sum(event1*log_f1_B1)+sum(event2*log_f2_B1)+sum((1-event1-event2)*log_Fb_B1)

    log_f1_G1 = f1_G1_func(par)*exp(par[3])/f1_func(par)
    log_f2_G1 = f2_G1_func(par)*exp(par[3])/f2_func(par)
    log_Fb_G1 = Fb_G1_func(par)*exp(par[3])/Fb_func(par)

    s3 = sum(event1*log_f1_G1)+sum(event2*log_f2_G1)+sum((1-event1-event2)*log_Fb_G1)

    return(c(s1,s2,s3))

  }
  Hessian_func = function(par) {

    A1 = (f1_A2_func(par)*exp(2*par[1])+f1_A1_func(par)*exp(par[1]))/f1_func(par)
    B1 = (f2_A2_func(par)*exp(2*par[1])+f2_A1_func(par)*exp(par[1]))/f2_func(par)
    C1 = (Fb_A2_func(par)*exp(2*par[1])+Fb_A1_func(par)*exp(par[1]))/Fb_func(par)
    A2 = ((f1_A1_func(par)*exp(par[1]))/f1_func(par))^2
    B2 = ((f2_A1_func(par)*exp(par[1]))/f2_func(par))^2
    C2 = ((Fb_A1_func(par)*exp(par[1]))/Fb_func(par))^2
    log_f1_A2 = A1-A2
    log_f2_A2 = B1-B2
    log_Fb_A2 = C1-C2

    h11 = sum(event1*log_f1_A2)+sum(event2*log_f2_A2)+sum((1-event1-event2)*log_Fb_A2)

    A1 = (f1_B2_func(par)*exp(2*par[2])+f1_B1_func(par)*exp(par[2]))/f1_func(par)
    B1 = (f2_B2_func(par)*exp(2*par[2])+f2_B1_func(par)*exp(par[2]))/f2_func(par)
    C1 = (Fb_B2_func(par)*exp(2*par[2])+Fb_B1_func(par)*exp(par[2]))/Fb_func(par)
    A2 = ((f1_B1_func(par)*exp(par[2]))/f1_func(par))^2
    B2 = ((f2_B1_func(par)*exp(par[2]))/f2_func(par))^2
    C2 = ((Fb_B1_func(par)*exp(par[2]))/Fb_func(par))^2
    log_f1_B2 = A1-A2
    log_f2_B2 = B1-B2
    log_Fb_B2 = C1-C2

    h22 = sum(event1*log_f1_B2)+sum(event2*log_f2_B2)+sum((1-event1-event2)*log_Fb_B2)

    A1 = (f1_G2_func(par)*exp(2*par[3])+f1_G1_func(par)*exp(par[3]))/f1_func(par)
    B1 = (f2_G2_func(par)*exp(2*par[3])+f2_G1_func(par)*exp(par[3]))/f2_func(par)
    C1 = (Fb_G2_func(par)*exp(2*par[3])+Fb_G1_func(par)*exp(par[3]))/Fb_func(par)
    A2 = ((f1_G1_func(par)*exp(par[3]))/f1_func(par))^2
    B2 = ((f2_G1_func(par)*exp(par[3]))/f2_func(par))^2
    C2 = ((Fb_G1_func(par)*exp(par[3]))/Fb_func(par))^2
    log_f1_G2 = A1-A2
    log_f2_G2 = B1-B2
    log_Fb_G2 = C1-C2

    h33 = sum(event1*log_f1_G2)+sum(event2*log_f2_G2)+sum((1-event1-event2)*log_Fb_G2)

    A1 = (f1_AB_func(par)*exp(par[1]+par[2]))/f1_func(par)
    B1 = (f2_AB_func(par)*exp(par[1]+par[2]))/f2_func(par)
    C1 = (Fb_AB_func(par)*exp(par[1]+par[2]))/Fb_func(par)
    A2 = ((f1_A1_func(par)*exp(par[1]))/f1_func(par))*((f1_B1_func(par)*exp(par[2]))/f1_func(par))
    B2 = ((f2_A1_func(par)*exp(par[1]))/f2_func(par))*((f2_B1_func(par)*exp(par[2]))/f2_func(par))
    C2 = ((Fb_A1_func(par)*exp(par[1]))/Fb_func(par))*((Fb_B1_func(par)*exp(par[2]))/Fb_func(par))
    log_f1_AB = A1-A2
    log_f2_AB = B1-B2
    log_Fb_AB = C1-C2

    h12 = h21 = sum(event1*log_f1_AB)+sum(event2*log_f2_AB)+sum((1-event1-event2)*log_Fb_AB)

    A1 = (f1_AG_func(par)*exp(par[1]+par[3]))/f1_func(par)
    B1 = (f2_AG_func(par)*exp(par[1]+par[3]))/f2_func(par)
    C1 = (Fb_AG_func(par)*exp(par[1]+par[3]))/Fb_func(par)
    A2 = ((f1_A1_func(par)*exp(par[1]))/f1_func(par))*((f1_G1_func(par)*exp(par[3]))/f1_func(par))
    B2 = ((f2_A1_func(par)*exp(par[1]))/f2_func(par))*((f2_G1_func(par)*exp(par[3]))/f2_func(par))
    C2 = ((Fb_A1_func(par)*exp(par[1]))/Fb_func(par))*((Fb_G1_func(par)*exp(par[3]))/Fb_func(par))
    log_f1_AG = A1-A2
    log_f2_AG = B1-B2
    log_Fb_AG = C1-C2

    h13 = h31 = sum(event1*log_f1_AG)+sum(event2*log_f2_AG)+sum((1-event1-event2)*log_Fb_AG)

    A1 = (f1_BG_func(par)*exp(par[2]+par[3]))/f1_func(par)
    B1 = (f2_BG_func(par)*exp(par[2]+par[3]))/f2_func(par)
    C1 = (Fb_BG_func(par)*exp(par[2]+par[3]))/Fb_func(par)
    A2 = ((f1_B1_func(par)*exp(par[2]))/f1_func(par))*((f1_G1_func(par)*exp(par[3]))/f1_func(par))
    B2 = ((f2_B1_func(par)*exp(par[2]))/f2_func(par))*((f2_G1_func(par)*exp(par[3]))/f2_func(par))
    C2 = ((Fb_B1_func(par)*exp(par[2]))/Fb_func(par))*((Fb_G1_func(par)*exp(par[3]))/Fb_func(par))
    log_f1_BG = A1-A2
    log_f2_BG = B1-B2
    log_Fb_BG = C1-C2

    h23 = h32 = sum(event1*log_f1_BG)+sum(event2*log_f2_BG)+sum((1-event1-event2)*log_Fb_BG)

    return(matrix(c(h11,h12,h13,h21,h22,h23,h31,h32,h33),3,3))

  }

  t.event = t.event-eta
  xbar = sum(t.event*event1)/sum(event1) # initial guess
  ybar = sum(t.event*event2)/sum(event2) # initial guess

  par_old = c(log(xbar),log(ybar),log(Gamma.0))
  count  = 0
  random = 0

  repeat {

    check1 = f1_func(par_old)
    check2 = f2_func(par_old)
    check3 = Fb_func(par_old)

    if (length(check1[check1 < epsilon.2]) > 0|
        length(check2[check2 < epsilon.2]) > 0|
        length(check3[check3 < epsilon.2]) > 0) {

      par_old = c(log(xbar*exp(runif(1,-r.1,r.1))),
                  log(ybar*exp(runif(1,-r.2,r.2))),
                  log(Gamma.0*exp(runif(1,-r.3,r.3))))

      count  = 0
      random = random+1
      next

    }

    par_new = par_old-solve(Hessian_func(par_old))%*%score_func(par_old)
    count   = count+1

    if (min(exp(par_new)) < epsilon.1) {

      par_old = c(log(xbar*exp(runif(1,-r.1,r.1))),
                  log(ybar*exp(runif(1,-r.2,r.2))),
                  log(Gamma.0*exp(runif(1,-r.3,r.3))))

      count  = 0
      random = random+1
      next

    }

    Dif = max(abs(exp(par_new)-exp(par_old)))

    if (Dif > D) {

      par_old = c(log(xbar*exp(runif(1,-r.1,r.1))),
                  log(ybar*exp(runif(1,-r.2,r.2))),
                  log(Gamma.0*exp(runif(1,-r.3,r.3))))

      count  = 0
      random = random+1
      next

    }
    if (Dif < epsilon.0) {break}

    par_old = par_new

  }

  ### MLE of alpha, beta, gamma ###
  Alpha_hat = exp(par_new[1])
  Beta_hat  = exp(par_new[2])
  Gamma_hat = exp(par_new[3])

  ### standard error of of alpha, beta, gamma ###
  Trans = matrix(c(exp(par_new[1]),0,0,0,exp(par_new[2]),0,0,0,exp(par_new[3])),3,3)
  Info  = t(Trans)%*%solve(-Hessian_func(par_new))%*%Trans
  Alpha_se = sqrt(Info[1,1])
  Beta_se  = sqrt(Info[2,2])
  Gamma_se = sqrt(Info[3,3])

  ### standard error of log-alpha, log-beta, log-gamma and CP ###
  log_Info = solve(-Hessian_func(par_new))
  log_Alpha_se = sqrt(log_Info[1,1])
  log_Beta_se  = sqrt(log_Info[2,2])
  log_Gamma_se = sqrt(log_Info[3,3])
  CI_Alpha = c(Alpha_hat*exp(-qnorm(0.975)*log_Alpha_se),
               Alpha_hat*exp(qnorm(0.975)*log_Alpha_se))
  CI_Beta  = c(Beta_hat*exp(-qnorm(0.975)*log_Beta_se),
               Beta_hat*exp(qnorm(0.975)*log_Beta_se))
  CI_Gamma = c(Gamma_hat*exp(-qnorm(0.975)*log_Gamma_se),
               Gamma_hat*exp(qnorm(0.975)*log_Gamma_se))

  ### results ###
  Alpha.res = c(Estimate = Alpha_hat,SE = Alpha_se,CI.lower = CI_Alpha[1],CI.upper = CI_Alpha[2])
  Beta.res  = c(Estimate = Beta_hat, SE = Beta_se, CI.lower = CI_Beta[1], CI.upper = CI_Beta[2] )
  Gamma.res = c(Estimate = Gamma_hat,SE = Gamma_se,CI.lower = CI_Gamma[1],CI.upper = CI_Gamma[2])

  if (Gamma_hat <= 1) {

    return(list(n = n,Iteration = count,Randomization = random,
                Alpha = Alpha.res,Beta = Beta.res,Gamma = Gamma.res,
                MeanX = "unavailable",MeanY = "unavailable",logL = logL_func(par_new)))

  } else {

    ### standard error of mux and muy ###
    mux_Trans0 = exp(par_new[1])*beta(exp(par_new[1])+1/exp(par_new[3]),1-1/exp(par_new[3]))
    mux_Trans1 = mux_Trans0*(digamma(exp(par_new[1])+1/exp(par_new[3]))-digamma(exp(par_new[1])))
    mux_Trans2 = mux_Trans0*exp(par_new[3])^-2*(digamma(1-1/exp(par_new[3]))-digamma(exp(par_new[1])+1/exp(par_new[3])))
    muy_Trans0 = exp(par_new[2])*beta(exp(par_new[2])+1/exp(par_new[3]),1-1/exp(par_new[3]))
    muy_Trans1 = muy_Trans0*(digamma(exp(par_new[2])+1/exp(par_new[3]))-digamma(exp(par_new[2])))
    muy_Trans2 = muy_Trans0*exp(par_new[3])^-2*(digamma(1-1/exp(par_new[3]))-digamma(exp(par_new[2])+1/exp(par_new[3])))
    mux_Trans  = c(mux_Trans1,0,mux_Trans2)
    muy_Trans  = c(0,muy_Trans1,muy_Trans2)
    mux_Info   = t(mux_Trans)%*%Info%*%mux_Trans
    muy_Info   = t(muy_Trans)%*%Info%*%muy_Trans
    mux_se = sqrt(mux_Info)
    muy_se = sqrt(muy_Info)

    ### standard error of log-mux, log-muy and CP ###
    mux_hat = mux_Trans0
    muy_hat = muy_Trans0

    log_mux_Trans1 = digamma(exp(par_new[1])+1/exp(par_new[3]))-digamma(exp(par_new[1]))
    log_mux_Trans2 = (digamma(1-1/exp(par_new[3]))-digamma(exp(par_new[1])+1/exp(par_new[3])))/exp(par_new[3])^2
    log_muy_Trans1 = digamma(exp(par_new[2])+1/exp(par_new[3]))-digamma(exp(par_new[2]))
    log_muy_Trans2 = (digamma(1-1/exp(par_new[3]))-digamma(exp(par_new[2])+1/exp(par_new[3])))/exp(par_new[3])^2
    log_mux_Trans  = c(log_mux_Trans1,0,log_mux_Trans2)
    log_muy_Trans  = c(0,log_muy_Trans1,log_muy_Trans2)
    log_mux_Info   = t(log_mux_Trans)%*%Info%*%log_mux_Trans
    log_muy_Info   = t(log_muy_Trans)%*%Info%*%log_muy_Trans
    log_mux_se = sqrt(log_mux_Info)
    log_muy_se = sqrt(log_muy_Info)

    CI_mux = c(mux_hat*exp(-qnorm(0.975)*log_mux_se),
               mux_hat*exp(qnorm(0.975)*log_mux_se))
    CI_muy = c(muy_hat*exp(-qnorm(0.975)*log_muy_se),
               muy_hat*exp(qnorm(0.975)*log_muy_se))

    MeanX.res = c(Estimate = mux_hat+eta,SE = mux_se,CI.lower = CI_mux[1]+eta,CI.upper = CI_mux[2]+eta)
    MeanY.res = c(Estimate = muy_hat+eta,SE = muy_se,CI.lower = CI_muy[1]+eta,CI.upper = CI_muy[2]+eta)

    return(list(n = n,Iteration = count,Randomization = random,
                Alpha = Alpha.res,Beta  = Beta.res ,Gamma = Gamma.res,
                MeanX = MeanX.res,MeanY = MeanY.res,logL = logL_func(par_new)))

  }

}

