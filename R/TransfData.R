#' Transposing data to use the functions in AGREL.
#'
#' Function to transpose dataframes if they are not yet in the desired format for the functions.
#'
#' @param id Variable indicating the ID of the rated observations.
#' @param var Variable indicating the ratings.
#' @param rater Variable indicating the raters.
#' @param data Dataframe
#'
#' @return Returns a transposed dataframe that can be used for the functions in this package.
#'
#' @examples
#' # Sample data
#' Df = cbind.data.frame(ID = sort(rep(1:10, 4)), var = sample(letters[1:2], 40, TRUE), raters = rep(paste("rater",1:4), 10))
#' NewDf = TransfData(ID, var, raters, Df)
#'
#' FleissK(NewDf)
#'
TransfData <- function(id,var,rater,data){
  Argz      = as.list(match.call())[-1]
  df        = data
  df$id     = eval(Argz$id,data)
  df$var    = eval(Argz$var,data)
  df$rater  = eval(Argz$rater,data)
  tmp       = df[,c("id","var","rater")]
  tmp       = reshape(tmp,timevar = "rater",idvar = "id",direction="wide")
  colnames(tmp) = gsub("var.","",colnames(tmp))
  return(tmp[,-1])
}
