#' Cohen's Kappa
#'
#' Slightly adjusted function from \code{\link[irr]{kappa2}} to calculate Cohen's Kappa (Cohen, 1960).
#'
#' @param ratings A dataframe or matrix of N x 2 with N the number of observations. The columns contain the ratings
#' of the 2 raters.
#' @param weight The weight matrix that has to be used to calculate the reliability. Default is \code{unweighted},
#' \code{'squared'} can be used to calculate Cohen's Weighted Kappa (Cohen, 1968).
#' @param sort.levels Sort if the levels are numeric.
#'
#' @return \item{method}{Method that was used to calculate reliability and weights used}
#' @return \item{subjects}{Number of subjects in the dataframe}
#' @return \item{nraters}{Number of raters}
#' @return \item{irr.name}{Type of reliability measure}
#' @return \item{value}{Value for Cohen's Kappa}
#' @return \item{StdErr}{The standard error of the estimated Kappa value}
#' @return \item{stat.name}{The corresponding test statistic}
#' @return \item{statistic}{The value of the test statistic}
#' @return \item{p.value}{The p-value for the test}
#' @return \item{Po}{Overall proportion of agreement}
#' @return \item{Pe}{Proportion of agreement expected by chance}
#' @return \item{ratings}{The original dataframe with the ratings}
#'
#' @references Cohen, J. (1960). A Coefficient of Agreement for Nominal Scales. \emph{Educational and Psychological
#' Measurement}, Vol.20(1), pp.37-46
#' @references Cohen, J. (1968).Weighted kappa: Nominal scale agreement provision for scaled
#' disagreement or partial credit. \emph{Psychological Bulletin}, Vol.70(4), pp.213-220
#'
#' @examples
#' # Load data
#' data(PsychMorbid)
#' Df = PsychMorbid[,1:2]
#'
#' # Unweighted kappa
#' CohenK(Df)
#'
#' # Weighted kappa
#' data(Agreement_deVet)
#' CohenK(Agreement_deVet[,2:3], weight = "squared")

CohenK <- function(ratings, weight = c("unweighted", "equal", "squared"),
                    sort.levels = FALSE)
{
  ratings <- as.matrix(na.omit(ratings))
  if (is.character(weight))
    weight = match.arg(weight)
  ns <- nrow(ratings)
  nr <- ncol(ratings)
  if (nr > 2) {
    stop("Number of raters exceeds 2. Try FleissK")
  }
  r1 <- ratings[, 1]
  r2 <- ratings[, 2]
  if ((is.numeric(r1)) | (is.numeric(r2)))
    sort.levels <- TRUE
  if (!is.factor(r1))
    r1 <- factor(r1)
  if (!is.factor(r2))
    r2 <- factor(r2)
  if (length(levels(r1)) >= length(levels(r2))) {
    lev <- c(levels(r1), levels(r2))
  }
  else {
    lev <- c(levels(r2), levels(r1))
  }
  if (sort.levels)
    lev <- sort(lev)
  lev <- lev[!duplicated(lev)]
  r1 <- factor(ratings[, 1], levels = lev)
  r2 <- factor(ratings[, 2], levels = lev)
  ttab <- table(r1, r2)
  nc <- ncol(ttab)
  if (is.numeric(weight))
    w <- 1 - (weight - min(weight))/(max(weight) - min(weight))
  else if (weight == "equal")
    w <- (nc - 1):0/(nc - 1)
  else if (weight == "squared")
    w <- 1 - (0:(nc - 1))^2/(nc - 1)^2
  else w <- c(1, rep(0, nc - 1))
  wvec <- c(sort(w, decreasing = FALSE), w[2:length(w)])
  nw <- length(w)
  weighttab <- matrix(0, nrow = nw, ncol = nw)
  for (i in 1:nw) {
    weighttab[i, ] <- wvec[(nw - (i - 1)):(2 * nw - i)]
  }
  agreeP <- sum(ttab * weighttab)/ns
  tm1 <- apply(ttab, 1, sum)
  tm2 <- apply(ttab, 2, sum)
  eij <- outer(tm1, tm2)/ns
  chanceP <- sum(eij * weighttab)/ns
  value <- (agreeP - chanceP)/(1 - chanceP)
  w.i <- apply(rep(tm2/ns, nc) * weighttab, 2, sum)
  w.j <- apply(rep(tm1/ns, each = nc) * weighttab, 1, sum)
  var.matrix <- (eij/ns) * (weighttab - outer(w.i, w.j, "+"))^2
  varkappa <- (sum(var.matrix) - chanceP^2)/(ns * (1 - chanceP)^2)
  SEkappa <- sqrt(varkappa)
  u <- value/SEkappa
  p.value <- 2 * (1 - pnorm(abs(u)))
  Results <- structure(list(method = paste("Cohen's Kappa for 2 Raters (Weights: ",
                                        paste(weight, collapse = ","), ")", sep = ""), subjects = ns,
                         nraters = nr, irr.name = "Kappa", value = value, weights = weight, StdErr = SEkappa, stat.name = "z",
                         statistic = u, p.value = p.value,Po=agreeP,Pe=chanceP,ratings=ratings),
                       class = "CohenK")
  return(Results)
}
