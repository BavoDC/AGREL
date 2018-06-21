\name{AGREL-package}
\alias{AGREL-package}
\alias{AGREL}
\docType{package}
\title{
\packageTitle{AGREL}
}
\description{
\packageDescription{AGREL}
}
%\details{

%The DESCRIPTION file:
%\packageDESCRIPTION{AGREL}
%\packageIndices{AGREL}
%~~ An overview of how to use the package, including the most important functions ~~
%}
\author{
\packageAuthor{AGREL}

Maintainer: \packageMaintainer{AGREL}
}
\references{
Cohen, J. (1960). A Coefficient of Agreement for Nominal Scales. \emph{Educational and Psychological
Measurement}, Vol.20(1), pp.37-46

Cohen, J. (1968).Weighted kappa: Nominal scale agreement provision for scaled
disagreement or partial credit. \emph{Psychological Bulletin}, Vol.70(4), pp.213-220

De Vet HCW, Mokkink LB, Terwee CB, Hoekstra OS, Knol DL. Clinicians are
right not to like Cohen’s k. \emph{BMJ} 2013;346:f2125.

de Vet, H.C.W., Terwee C.B., Knol, D.L., Bouter, L.M. (2006). When to use agreement
versus reliability measures. \emph{Journal of Clinical Epidemiology}, Vol.59(10), pp.1033-1039

de Vet, H.C.W., Dikmans, R.E., Eekhout, I. (2017). Specific agreement on dichotomous outcomes can be
calculated for more than two raters. \emph{Journal of Clinical Epidemiology}, Vol.83, pp.85-89

Fleiss, J. L. (1971). Measuring nominal scale agreement among many raters. \emph{Psychological
Bulletin}, Vol.76(5), pp.378-382

Fleiss J.L., Davies M. (1982). Jackknifing functions of multinomial frequencies, with an application to a
measure of concordance. \emph{Am J Epidemiol};115: 841-845.

Uebersax, J.S. (1982). A design-independent method for measuring the reliability of psychiatric
diagnosis. \emph{Journal of Psychiatric Research}, Vol.17(4), pp.335-342

}
\keyword{ package }
\seealso{
  \code{\link{AgreemGeneralizedVector}}
  \code{\link{CohenK}}
  \code{\link{ConfintJack}}
  \code{\link{DiagnosticAgreement.deVet}}
  \code{\link{FleissK}}
  \code{\link{Kappa.Matrix}}
}
\examples{
# Vignette of package
vignette(package = "AGREL")

# Load data
data(Agreement_deVetArticle)

# Specific and overall agreement, dichotomous variable
DiagnosticAgreement.deVet(Agreement_deVetArticle)

# Specific and overall agreement, multinomial variable
data(PsychMorbid)
AgreemGeneralizedVector(PsychMorbid)

# Cohen's Kappa
Df = Agreement_deVetArticle[,1:2]
Results1 = CohenK(Df)
Results1
ConfintJack(Results1)

# Fleiss' Kappa
FleissK(Agreement_deVetArticle)
}
