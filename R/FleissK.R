#' Fleiss' Kappa
#'
#' Function to calculate Fleiss' Kappa (Fleiss, 1971).
#'
#' @param ratings A dataframe or matrix of N x P with N the number of observations and P the number of raters.
#' @param weights The weight matrix that has to be used to calculate the reliability. Default is \code{unweighted},
#' @param conflev Confidence level to calculate the confidence intervals of Fleiss' Kappa.
#' @param N Final population correction factor for the variance of Fleiss' Kappa (Gwet, 2008).
#'
#' @return \item{method}{Method that was used to calculate reliability}
#' @return \item{subjects}{Number of subjects in the dataframe}
#' @return \item{nraters}{Number of raters}
#' @return \item{value}{Value for Fleiss' Kappa}
#' @return \item{StdErr}{The standard error of the estimated Kappa value}
#' @return \item{stat.name}{The corresponding test statistic}
#' @return \item{statistic}{The value of the test statistic}
#' @return \item{p.value}{The p-value for the test}
#' @return \item{alpha}{Alpha-level used to calculate the confidence intervals}
#' @return \item{Po}{Overall proportion of agreement}
#' @return \item{Pe}{Proportion of agreement expected by chance}
#' @return \item{ratings}{The original dataframe with the ratings}
#' @return \item{WeightsMatrix}{The weights matrix used.}
#'
#' @references Fleiss, J. L. (1971). Measuring nominal scale agreement among many raters. \emph{Psychological
#' Bulletin}, Vol.76(5), pp.378-382
#' @references Gwet, K.L. (2008). Variance Estimation of Nominal-Scale Inter-Rater Reliability with Random
#' Selection of Raters \emph{Psychometrika}, Vol.73(3), pp.407-430
#'
#' @examples
#' # Load data
#' data(PsychMorbid)
#'
#' FleissK(PsychMorbid)
FleissK <- function(ratings,weights="unweighted",conflev=0.95,N=Inf){
  ratings.mat <- as.matrix(ratings)
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # final population correction

  # creating a vector containing all categories used by the raters

  categ.init <- unique(as.vector(ratings.mat))
  if (is.numeric(categ.init)){
    categ <- sort(as.vector(na.omit(categ.init)))
  }else{
    categ.init <- trim(categ.init) #trim vector elements to remove leading and trailing blanks
    categ <- na.omit(categ.init[nchar(categ.init)>0])
  }
  q <- length(categ)

  # creating the weights matrix

  if (is.character(weights)){
    if (weights=="quadratic")
      weights.mat<-quadratic.weights(categ)
    else if (weights=="ordinal")
      weights.mat<-ordinal.weights(categ)
    else if (weights=="linear")
      weights.mat<-linear.weights(categ)
    else if (weights=="radical")
      weights.mat<-radical.weights(categ)
    else if (weights=="ratio")
      weights.mat<-ratio.weights(categ)
    else if (weights=="circular")
      weights.mat<-circular.weights(categ)
    else if (weights=="bipolar")
      weights.mat<-bipolar.weights(categ)
    else weights.mat<-identity.weights(categ)
  }else weights.mat= as.matrix(weights)

  # creating the nxq agreement matrix representing the distribution of raters by subjects and category

  agree.mat <- matrix(0,nrow=n,ncol=q)
  if (is.character(ratings.mat)) ratings.mat<-trim(ratings.mat)
  for(k in 1:q){
    k.mis <-(ratings.mat==categ[k])
    in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
    agree.mat[,k] <- in.categ.k%*%rep(1,r)
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))

  # calculating fleiss's generalized kappa coefficient

  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more

  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/(ri.vec%*%t(rep(1,q)))))
  pe <- sum(weights.mat * (pi.vec%*%t(pi.vec)))
  fleiss.kappa <- (pa-pe)/(1-pe)

  # calculating variance, stderr & p-value of gwet's ac1 coefficient

  den.ivec <- ri.vec*(ri.vec-1)
  den.ivec <- den.ivec - (den.ivec==0) # this operation replaces each 0 value with -1 to make the next ratio calculation always possible.
  pa.ivec <- sum.q/den.ivec

  pe.r2 <- pe*(ri.vec>=2)
  kappa.ivec <- (n/n2more)*(pa.ivec-pe.r2)/(1-pe)
  pi.vec.wk. <- weights.mat%*%pi.vec
  pi.vec.w.k <- t(weights.mat)%*%pi.vec
  pi.vec.w <- (pi.vec.wk. + pi.vec.w.k)/2

  pe.ivec <- (agree.mat%*%pi.vec.w)/ri.vec
  kappa.ivec.x <- kappa.ivec - 2*(1-fleiss.kappa) * (pe.ivec-pe)/(1-pe)

  var.fleiss <- ((1-f)/(n*(n-1))) * sum((kappa.ivec.x - fleiss.kappa)^2)
  stderr <- sqrt(var.fleiss)# kappa's standard error
  p.value <- 2*(1-pt(abs(fleiss.kappa/stderr),n-1))

  sum(apply(pi.vec,1,FUN=function(x){x*(1-x)}))^2 -> one
  sum(apply(pi.vec,1,FUN=function(x){x*(1-x)*((1-x)-x)})) -> two

  se <- sqrt((2/(n*r*(r-1)*one))*((one)-two))

  lcb <- fleiss.kappa - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
  ucb <- min(1,fleiss.kappa + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
  Results <- structure(list(method = "Fleiss Kappa", subjects = n,
                         nraters = r, value = fleiss.kappa, StdErr = se,stat.name = "z",
                         statistic = fleiss.kappa/stderr,
                         lcl=lcb, ucl=ucb, p.value = p.value,
                         alpha = 1-conflev,
                         Po=pa,Pe=pe,ratings=ratings,weights=weights, WeightsMatrix = weights.mat),
                          class = "FleissK")
  return(Results)
}
