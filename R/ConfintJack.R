#' Confidence intervals using the jackknife method.
#'
#' This function can be used to calculate the confidence intervals of Cohen's Kappa or Fleiss' Kappa.
#'
#' @param object Object of type CohenK or FleissK, see \code{\link{CohenK}} and \code{\link{FleissK}}.
#' @param alpha Alpha-level for computing the confidence intervals
#'
#' @return Object with the statistic and its confidence intervals
#' @references Fleiss J.L., Davies M. (1982). Jackknifing functions of multinomial frequencies, with an application to a
#' measure of concordance. \emph{Am J Epidemiol};115: 841-845.
#'
#' @seealso \code{\link{CohenK}} and \code{\link{FleissK}}
#'
#' @examples
#' data(PsychMorbid)
#' KappaStat = FleissK(PsychMorbid)
#' ConfintJack(KappaStat)
ConfintJack <- function(object,alpha=0.05){
  if(!any(class(object)%in%c("CohenK", "FleissK"))) stop("This function can only be used for objects of type CohenK or FleissK.")
  MethAGREL <- class(object)
  theta <- function(x,xdata){
    do.call(MethAGREL,args=list(ratings=xdata[x,]))$value
  }
  n <- nrow(object$ratings)
  ResJack <- jackknife(1:n,theta,xdata=object$ratings)
  ValueStat <- object$value
  LCL <- ValueStat - qnorm(1-alpha/2)*ResJack$jack.se
  UCL <- ValueStat + qnorm(1-alpha/2)*ResJack$jack.se
  Results <- c(ValueStat,LCL,UCL)
  names(Results) <- c("Statistic","LCL","UCL")
  return(Results)
}
