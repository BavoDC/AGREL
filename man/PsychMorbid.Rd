\name{PsychMorbid}
\alias{PsychMorbid}
\docType{data}
\title{
Psychological morbidity dataset
}
\description{
This dataset contains the ratings of 6 psychiatrists who were asked to diagnose 30 patients (Fleiss, 1971). Possible diagnoses were depression, personality disorder (PD), schizophrenia, neurosis and other. Note that this dataset is different from the original dataset documented in Fleiss (1971) in that it was transformed into a format where the ratings of the psychiatrists are put in separate columns.
}
\usage{data("PsychMorbid")}
\format{
  A data frame with 30 observations containing the ratings of the 6 psychiatrists.
  \describe{
     \item{Psychiatrist1...Psychiatrist6}{These variables contain the ratings of each of the psychiatrists.}
  }
}
\details{
The original dataset of Fleiss is \code{\link{PsychMorbidOrig}}
}
\source{
Fleiss, J. L. (1971). Measuring nominal scale agreement among many raters. \emph{Psychological Bulletin}, Vol.76(5), pp.378-382
}
\examples{
data(PsychMorbid)
library(AGREL)
FleissK(PsychMorbid)
}
\keyword{datasets}
