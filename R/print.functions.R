print.AgreementTable <- function(x,...){
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#\n")
  cat("#          Agreement           #\n")
  cat("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#\n")
  cat(paste("\n- Number of raters =",x$nraters,"\n\n"))
  cat(paste("- Sample size =",x$N,"\n\n"))
  cat("\n- Summed table:\n\n")
  print(x$SumTable)
  if(any(class(x)%in%"ObservedAgreem")){
    cat("\n\n- Observed Agreement:\n\n")
    print(x$ObservedAgreem)
  }
  if(any(class(x)%in%"SpecificAgreem")){
    cat("\n\n- Specific Agreement:\n\n")
    print(x$SpecificAgreem)
  }
}

print.AgreemGeneralizedVector <- function(x,...){
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#\n")
  cat("# Agreement for multinomial variables #\n")
  cat("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#\n")
  cat(paste("\n- Number of raters =",x$nraters,"\n\n"))
  cat(paste("- Sample size =",x$N,"\n\n"))
  cat(paste("- Overall agreement =", round(x$OverallAgreement, 3)))
  cat("\n\n- Specific agreement:\n\n")
  print(x$Ps)
}


print.CohenK <- function(x,...){
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#\n")
  cat("#       Cohen's Kappa          #\n")
  cat("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#\n")
  cat(paste("\n- Number of raters =",x$nraters,"\n\n"))
  cat(paste("- Sample size =",nrow(x$ratings),"\n\n"))
  cat(paste("- Cohen's Kappa =",round(x$value,3),"\n\n"))
  cat('- Overall proportion of agreement:',x$Po,'\n\n')
  cat('- Proportion chance agreement:',x$Pe,'\n\n')
  cat(paste("- Used weights = ", x$weights,"\n\n"))
}

print.FleissK <- function(x,...){
    cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#\n")
    cat("#  Fleiss' Kappa Coefficient  #\n")
    cat("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#\n")
    cat(paste("\n- Number of raters =",x$nraters,"\n\n"))
    cat(paste("- Sample size =",nrow(x$ratings),"\n\n"))
    cat(paste("- Fleiss' Kappa =",round(x$value,3),"\n\n"))
    cat('- Overall proportion of agreement:',x$Po,'\n\n')
    cat('- Proportion chance agreement:',x$Pe,'\n\n')
    if (x$weights!="unweighted") {
      cat('\n')
      if (!is.numeric(x$weights)) {
        cat('Weights: ', x$weights,'\n')
        cat('---------------------------\n')
      }
      else{
        cat('Weights: Custom Weights\n')
        cat('---------------------------\n')
      }
      print(x$WeightsMatrix)
    }

}

print.BAplotMultipleR <- function(x, ...){
  cat("\nAnova summary:\n")
  print(x$AnovaSummary)
  print(x$ICC)
}
