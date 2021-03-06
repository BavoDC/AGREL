% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CohenK.R
\name{CohenK}
\alias{CohenK}
\title{Cohen's Kappa}
\usage{
CohenK(ratings, weight = c("unweighted", "equal", "squared"),
  sort.levels = FALSE)
}
\arguments{
\item{ratings}{A dataframe or matrix of N x 2 with N the number of observations. The columns contain the ratings
of the 2 raters.}

\item{weight}{The weight matrix that has to be used to calculate the reliability. Default is \code{unweighted},
\code{'squared'} can be used to calculate Cohen's Weighted Kappa (Cohen, 1968).}

\item{sort.levels}{Sort if the levels are numeric.}
}
\value{
\item{method}{Method that was used to calculate reliability and weights used}

\item{subjects}{Number of subjects in the dataframe}

\item{nraters}{Number of raters}

\item{irr.name}{Type of reliability measure}

\item{value}{Value for Cohen's Kappa}

\item{StdErr}{The standard error of the estimated Kappa value}

\item{stat.name}{The corresponding test statistic}

\item{statistic}{The value of the test statistic}

\item{p.value}{The p-value for the test}

\item{Po}{Overall proportion of agreement}

\item{Pe}{Proportion of agreement expected by chance}

\item{ratings}{The original dataframe with the ratings}
}
\description{
Slightly adjusted function from \code{\link[irr]{kappa2}} to calculate Cohen's Kappa (Cohen, 1960).
}
\examples{
# Load data
data(PsychMorbid)
Df = PsychMorbid[,1:2]

# Unweighted kappa
CohenK(Df)

# Weighted kappa
data(Agreement_deVet)
CohenK(Agreement_deVet[,2:3], weight = "squared")
}
\references{
Cohen, J. (1960). A Coefficient of Agreement for Nominal Scales. \emph{Educational and Psychological
Measurement}, Vol.20(1), pp.37-46

Cohen, J. (1968).Weighted kappa: Nominal scale agreement provision for scaled
disagreement or partial credit. \emph{Psychological Bulletin}, Vol.70(4), pp.213-220
}
