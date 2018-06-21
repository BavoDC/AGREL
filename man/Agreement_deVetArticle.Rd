\name{Agreement_deVetArticle}
\alias{Agreement_deVetArticle}
\docType{data}
\title{
Breast reconstruction dataset as used in de Vet et al. (2017)}
\description{
In this dataset, the ratings of the surgeons are dichotomized (see \code{\link{Agreement_deVet}}). The scores 1-3 are now categorized as dissatisfied and score 4-5 are categorized as satisfied.
}
\usage{data("Agreement_deVetArticle")}
\format{
  A data frame with 50 observations containing the dichotomized version of the ratings of the surgeons.
  \describe{
    \item{Surgeon1...Surgeon4}{These variables contain the dichotomized ratings of each of the surgeons.}
  }
}
\source{
de Vet, H.C.W., Dikmans, R.E., Eekhout, I. (2017). Specific agreement on dichotomous outcomes can be
calculated for more than two raters. \emph{Journal of Clinical Epidemiology}, Vol.83, pp.85-89
}
\examples{
data(Agreement_deVetArticle)
library(AGREL)
DiagnosticAgreement.deVet(Agreement_deVetArticle)
}
\keyword{datasets}
