#' Calculating agreement for dichotomous variables.
#'
#' Function to calculate specific agreement and overall proportion of agreement for dichotomous variables.
#'
#' @param ratings A dataframe or matrix of N x P with N the number of observations and P the number of raters.
#' @param CI Logical, indicates if confidence intervals have to be calculated
#' @param ConfLevel The confidence level to be used in calculating the confidence intervals. Possible values
#' are \code{"continuity"}, \code{"Fleiss"} and \code{"bootstrap"}.
#' @param correction Method of calculating the confidence intervals (de Vet et al., 2017). The confidence intervals (CI) can be calculated using
#' continuity correction, Fleiss correction or by use of bootstrap samples.
#' @param NrBoot In case of bootstrap methodology to calculate the CI, the number of bootstrap samples.
#' @param Parallel Logical, indicates if parallel computing has to be used to compute the confidence intervals.
#' Implemented only when using bootstrapping to calculate the confidence intervals.
#' @param no_cores Number of cores if parallel computing is used. Default is 1 core less than the number of
#' cores present.
#'
#' @return A list with the following components:
#' @return \item{SumTable}{The summed table as defined in the article of de Vet et al. (2017)}
#' @return \item{ObservedAgreem}{The overall proportion of agreement (with CI).}
#' @return \item{SpecificAgreem}{The specific agreement for each of the categories (with CIs).}
#' @references De Vet HCW, Mokkink LB, Terwee CB, Hoekstra OS, Knol DL. Clinicians are
#' right not to like Cohen’s k. \emph{BMJ} 2013;346:f2125.
#' @references de Vet, H.C.W., Terwee C.B., Knol, D.L., Bouter, L.M. (2006). When to use agreement
#' versus reliability measures. \emph{Journal of Clinical Epidemiology}, Vol.59(10), pp.1033-1039
#' @references de Vet, H.C.W., Dikmans, R.E., Eekhout, I. (2017). Specific agreement on dichotomous outcomes can be
#' calculated for more than two raters. \emph{Journal of Clinical Epidemiology}, Vol.83, pp.85-89
#'
#' @details This function is based on the functions as given in the appendix of the article of de Vet et al. (2017).
#'
#' @examples
#' # Load data
#' data(Agreement_deVetArticle)
#' Df = Agreement_deVetArticle
#'
#' DiagnosticAgreement.deVet(Df)
DiagnosticAgreement.deVet <- function(ratings, CI = T, ConfLevel = 0.95, correction=c("continuity","Fleiss", "bootstrap"), NrBoot = 1e3,
                                      Parallel = F, no_cores = detectCores() - 1){
  if(length(unique(unlist(ratings)))>2) stop("Multinomial variable (number of unique values > 2). Consider using PA.matrix")
  if(is.matrix(ratings)) ratings = as.data.frame(ratings)
  stopifnot(is.numeric(NrBoot))
  if(NrBoot < 200) stop("200 is the minimum number of bootstrap samples.")
  correction = match.arg(correction)
  SumTab     = sumtable(ratings)
  AgreemTab  = Agreement(SumTab)
  AgreemTab  = SpecAgreem(AgreemTab)
  if(CI) AgreemTab  = CIAgreem(AgreemTab, level=ConfLevel, AgreemStat = "all", correction=correction, NrBoot = NrBoot,
                               Parallel=Parallel, no_cores = co_cores)
  return(AgreemTab)
}
