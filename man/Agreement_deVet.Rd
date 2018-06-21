\name{Agreement_deVet}
\alias{Agreement_deVet}
\docType{data}
\title{
Breast reconstruction dataset
}
\description{
This dataset contains the ratings of 4 surgeons on the breast reconstruction of 50 women. Each of the surgeons were given photographs of the breast reconstruction and the surgeons had to rate the reconstruction quality on a 5-point ordinal scale (very dissatisfied = 1 to very satisfied = 5).
}
\usage{data("Agreement_deVet")}
\format{
  A data frame with 50 observations containing the ratings of each of the 4 surgeons.
  \describe{
    \item{\code{ID}}{Patient ID}
    \item{Surgeon1...Surgeon4}{These variables contain the ratings of each of the surgeons.}
  }
}
\source{
de Vet, H.C.W., Dikmans, R.E., Eekhout, I. (2017). Specific agreement on dichotomous outcomes can be
calculated for more than two raters. \emph{Journal of Clinical Epidemiology}, Vol.83, pp.85-89}
\examples{
data(Agreement_deVet)
library(AGREL)
CohenK(Agreement_deVet[, 2:3], weight = "squared")
}
\keyword{datasets}
