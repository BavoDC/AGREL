#' Modified Bland-Altman plot for multiple raters
#'
#' @param rater Variable indicating which rater made the rating.
#' @param subject Variable for the subject.
#' @param variable Variable containing the ratings of the subjects.
#' @param data The dataframe
#' @param LoA The limits of agreement that have to be used. See Details.
#' @param ArgzLoess Arguments for \code{loess.smooth} if \code{LoA = "loess"}.
#' @param symbols.graph Optional, the symbols to be used in the graph for the raters.
#' @param colSymbols Optional, the color of the symbols
#' @param Legend Logical, indicates whether a legend has to be printed. Default is \code{TRUE}.
#' @param ArgzLegend The arguments for the legend, see \code{\link{legend}}.
#' @param xlab The label for the x-axis.
#' @param ylab The label for the y-axis.
#' @param ... Additional arguments to be passed to \code{\link{plot}}.
#'
#' @return Returns a list with the following objects
#' @return \item{AnovaSummary}{The summary of the two-way ANOVA}
#' @return \item{AvgPerSubject}{The average score per subject.}
#' @return \item{LoA}{The limits of agreement calculated according to Jones et al. (2011).}
#' @return \item{ICC}{The intraclass correlation coefficients. See \code{\link{ICC}}}
#'
#' @references Bland, J.M., Altman, D.G.(1999). Measuring agreement in method comparison studies. \emph{Statistical Methods In Medical Research},
#' Vol.8(2), pp. 135-160
#'
#' @references Jones, M., Dobson, A., O'Brian, S. (2011). A graphical method for assessing agreement with the
#' mean between multiple observers using continuous measures. \emph{Int J Epidemiol}. Vol.40: pp. 1308-1313.
#'
#' @references Royston, P., Wright, E.M. (1988). How to construct 'normal ranges' for fetal variables. \emph{Ultrasound Obstet Gynecol}, Vol.11: pp. 30-38.
#'
#' @details The limits of agreement (LoA) have to be interpreted differently than in Bland-Altman plots with 2 raters. In the modified Bland-Altman
#' plots, the LoA indicate how different an individual rater can be when compared with the mean of all the raters (Jones et al., 2011). However, as with the
#' Bland-Altman plots for 2 raters, it may be that the variability of the differences increase as the magnitude of the measurement increases (Bland and
#' Altman, 1999). Alternative LoA can then be plotted using a method based on the method of Royston and Wright (1998). Specifying \code{LoA="loess"}, we get
#' approximate 2.5 and 97.5 percentile curves. In this method, the standard deviation of the difference scores per subject is calculated and the LoA
#' per subject are calculated as -/+ 1.96 times the standard deviation. A loess fit is then used to connect these LoA per subject. Note that we assume
#' that the difference scores are normally distributed (with mean 0 as we are working with centered values).
#'
#' @examples
#' data(DentalStudy)
#' BAplotMultipleR(dentist, patient, DMFS, DentalStudy)
BAplotMultipleR <- function(rater, subject, variable, data, LoA = c("Default", "loess"), ArgzLoess = list(span = 1, degree = 1),
                            symbols.graph=NULL, colSymbols=NULL,
                            Legend = T, ArgzLegend = list(x = "topright"),
                            xlab = "Average", ylab = "Difference",
                            ...){
  if(missing(data)|!is.data.frame(data)) stop(paste(data,"must be of type dataframe"))
  arguments =  as.list(match.call())
  LoA = match.arg(LoA)
  if(Legend && !is.list(ArgzLegend)) stop("ArgzLegend must be of type list.")

  data$rater = eval(arguments$rater,data)
  if(!is.factor(data$rater)) data$rater = tryCatch(as.factor(data$rater), warning = function(w) stop("Conversion of rater to factor failed"),
                                                       error = function(e) stop("Conversions of rater to factor failed"))
  data$subject = eval(arguments$subject,data)
  if(!is.factor(data$subject)) data$subject = tryCatch(as.factor(data$subject), warning = function(w) stop("Conversion of subject to factor failed"),
                                                         error = function(e) stop("Conversions of subject to factor failed"))
  data$var = eval(arguments$variable,data)
  if(!is.numeric(data$var)) data$var = tryCatch(as.numeric(data$var),
                                                              warning = function(w) stop("Conversion of variable to numeric failed"),
                                                          error = function(e) stop("Conversions of variable to numeric failed"))
  if(anyNA(data$var)) stop("Missing values detected.")

  summ = summary(aov(var ~ subject + rater, data = data))

  MSE  = summ[[1]]$`Mean Sq`[3]
  m = length(unique(data$subject))

  DfICC = TransfData(subject, var, rater, data)
  ICCv  = ICC(DfICC)


  # Plot
  mean.per.subj  = with(data, aggregate(var, by = list(subject), mean))
  diff.raters    = unique(data$rater)
  diff.per.rater = list()
  data = data[with(data,order(subject)),]

  for(i in diff.raters){
    temp = data[data$rater==i,]
    diff.per.rater[[i]] = with(temp, var - mean.per.subj[, 2])
  }


  min.diff.raters = lapply(diff.per.rater, min)
  max.diff.raters = lapply(diff.per.rater, max)
  if(is.null(symbols.graph)) symbols = 1:length(diff.raters) else symbols = symbols.graph
  if(length(symbols)!=length(diff.raters)) stop("Length of vector of symbols.graph is not equal to number of raters")


  plot(1 , 1, xlim = range(data$var),ylim= c(min(unlist(min.diff.raters)),
                                                 max(unlist(max.diff.raters))), col="white", xlab = xlab, ylab = ylab, ...)

  for(i in unique(diff.raters)){
    if(!is.null(colSymbols)) points(mean.per.subj[, 2] ,diff.per.rater[[i]], pch = symbols[which(i==diff.raters)],
                                    col = colSymbols[which(unique(diff.raters=i))])
      else points(mean.per.subj[, 2], diff.per.rater[[i]], pch = symbols[which(i==diff.raters)])
  }

  abline(h = 0, lty = 2)

  if(LoA=="Default"){
    abline(h = sqrt(MSE) * 1.96)
    abline(h = sqrt(MSE) * -1.96)
  }else{
  temp2 = matrix(NA, m, 1)
  for(i in seq_len(m)){
    temp = NULL
    for(j in 1:length(diff.raters)){
      temp = c(temp, (diff.per.rater[[j]][i]))
    }
    temp2[i,] = sd(temp)
  }
  ArgzLoessLCL = append(ArgzLoess, list(x = mean.per.subj[order(mean.per.subj[,2]),2], y = (-1.96*(temp2))[order(mean.per.subj[,2])]))
  ArgzLoessUCL = append(ArgzLoess, list(x = mean.per.subj[order(mean.per.subj[,2]),2], y = (1.96*(temp2))[order(mean.per.subj[,2])]))

  LoessLCL = do.call("loess.smooth", ArgzLoessLCL)
  lines(LoessLCL$x, LoessLCL$y, lwd = 2)
  LoessUCL = do.call("loess.smooth", ArgzLoessUCL)
  lines(LoessUCL$x, LoessUCL$y, lwd = 2)
  }

  if(Legend && length(ArgzLegend)==1){
    ArgzLegend = append(ArgzLegend, list(legend = diff.raters, ncol = 2, pch = symbols, bty = "n"))
  }
  colnames(mean.per.subj) = c("Subject", "Average")

  if(Legend) do.call("legend", ArgzLegend)
  Results = list(AnovaSummary=summ, AvgPerSubject = mean.per.subj, LoA = c("LCL" = sqrt(MSE) * -1.96, "UCL" = sqrt(MSE) * 1.96),
                 ICC = ICCv)
  class(Results) = append(class(Results), "BAplotMultipleR")
  return(Results)
}
