\name{DentalStudy}
\alias{DentalStudy}
\docType{data}
\title{
Dental reliability study}
\description{
For the dental reliability study, 4 different dentists were asked to give 10 patients a DMFS score. This dataset was used to illustrate the modified Bland-Altman plot in Jones et al. (2011).
}
\usage{data("DentalStudy")}
\format{
  A data frame with 40 observations on the following 3 variables.
  \describe{
    \item{\code{patient}}{The patient ID}
    \item{\code{dentist}}{The variable indicating the dentist that rated the patient.}
    \item{\code{DMFS}}{The DMFS score}
  }
}
\source{
Jones M., Dobson A., O’Brian S. (2011). A graphical method for assessing agreement with the
mean between multiple observers using continuous measures. \emph{Int J Epidemiol}. Vol.40: pp. 1308-1313.
}
\examples{
data(DentalStudy)
library(AGREL)
BAplotMultipleR(dentist, patient, DMFS, DentalStudy)
}
\keyword{datasets}
