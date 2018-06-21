\name{PsychMorbidOrig}
\alias{PsychMorbidOrig}
\docType{data}
\title{
Psychological morbidity dataset as used in Fleiss (1971)}
\description{
This is the original dataset as described in Fleiss (1971). 6 psychiatrists were asked to diagnose 30 patients and the possible diagnoses were depression, personality disorder (PD), schizophrenia, neurosis and other. }
\usage{data("PsychMorbidOrig")}
\format{
  A data frame with information on the 30 diagnosed patients. The variables are
  \describe{
    \item{\code{depression}}{For each of the patients, the number of times that the patient was diagnosed with depression.}
    \item{\code{PD}}{For each of the patients, the number of times that the patient was diagnosed with a personality disorder.}
    \item{\code{schizophrenia}}{For each of the patients, the number of times that the patient was diagnosed with schizophrenia.}
    \item{\code{neurosis}}{For each of the patients, the number of times that the patient was diagnosed with neurosis.}
    \item{\code{other}}{Variable when none of the other categories was applicable.}
  }
}
\source{
Fleiss, J. L. (1971). Measuring nominal scale agreement among many raters. \emph{Psychological Bulletin}, Vol.76(5), pp.378-382
}
\examples{
data(PsychMorbidOrig)
head(PsychMorbidOrig)
}
\keyword{datasets}
