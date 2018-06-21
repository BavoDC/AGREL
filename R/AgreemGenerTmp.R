#' Function to calculate the overall proportion agreement and specific agreement.
#'
#' This function can be used in case of multinomial variables and when there are more than two raters. Note that it can also be used for
#'  dichotomous variables and when there are two raters
#'
#' @param df A dataframe or matrix of N x P with N the number of observations and P the number of raters.
#'
#' @return A list containing the overall proportion of agreement, the specific agreement, the number of raters, the number of observations
#' , the original dataframe and the (marginal) prevalences.
#'
#' @note This is a temporary function! A more advanced function with more options and output will be made public as soon as the article is published.
#'
#' @examples
#' data(PsychMorbid)
#'
#' AgreemGeneralizedVector(PsychMorbid)
AgreemGeneralizedVector <- function(df){
  if(is.matrix(df)) df = as.data.frame(df)
  if(!is.data.frame(df)) stop("Must be of type dataframe.")
  UniqValues  = as.character(unique(unlist(df)))
  if(length(UniqValues)==1) stop("Only one unique value.")
  if(any(table(unlist(df))<10)){
    warning("The frequency of some levels of the variable are very low", immediate. = T)
    LowPrev = names(table(unlist(df)))[table(unlist(df))<10]
  }else{
    LowPrev = NULL
  }

  R           = ncol(df)
  N           = nrow(df)
  S           = NULL
  Ps          = NULL
  Prev        = NULL

  for(j in UniqValues){
    Pj      = sum(apply(df,1,function(x) sum(x==j)))/(N*R)
    S[j]    = sum(apply(df,1,function(x) sum(x==j)*(sum(x==j) - 1)))
    Ps[j]   = S[j]/sum(apply(df,1,function(x) sum(x==j)*(R - 1)))
    Prev[j] = sum(apply(df,1,function(x) sum(x==j)*(R - 1))) / (N * R * (R - 1))
  }
  Oposs          = R*(R-1)*N
  Po             = sum(S)/Oposs

  Results        = list(OverallAgreement = Po, Ps = data.frame("Ps" = Ps), nraters = R, N = N, data = df,
                        Prevalence = Prev)
  Results        = Results[!unlist(lapply(Results, is.null))]
  class(Results) = append(class(Results), "AgreemGeneralizedVector")
  return(Results)
}
