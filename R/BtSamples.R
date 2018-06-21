BtSamples <- function(df,fun,...){
  Nboot <- sample(1:nrow(df),nrow(df),replace=T)
  df.tmp <- df[Nboot,]
  Result <- fun(df.tmp,...)
  return(Result)
}

