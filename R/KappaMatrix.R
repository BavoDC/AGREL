#' Matrix of Kappa-type coefficients
#'
#' Function to compute a matrix of kappa-type coefficients to assess the reliability of nominal variables
#' (Roberts & McNamee, 1998).
#'
#' @param df A dataframe or matrix of N x P with N the number of observations and P the number of raters.
#' @param FleissKappa Logical, indicates if Fleiss' Kappa has to be calculated.
#' @param CI.est Logical, indicates if confidence intervals using bootstrapping have to be calculated.
#' @param N.bt Number of bootstrap samples.
#' @param parallel Logical, indicates if parallel computing has to be used.
#' @param NrCores Number of cores if parallel computing is used. Default is 1 core less than the number of
#' cores present.
#'
#' @return A list containing the kappa-matrix and Fleiss' Kappa if asked for.
#' @references Roberts C., McNamee R. (1998). A matrix of kappa-type coefficients to assess the reliability of
#' nominal scales. \emph{Statistics in medicine}, Vol.17(4), pp.471-88
#'
#' @examples
#' data(PsychMorbid)
#' Kappa.Matrix(PsychMorbid)
Kappa.Matrix <- function(df,FleissKappa=F,CI.est=F,N.bt=1e3,parallel=F,NrCores=detectCores - 1){
if(is.matrix(df)) df = as.data.frame(df)
if(!is.data.frame(df)) stop("Must be of type dataframe.")
UniqValues <- as.character(unique(unlist(df)))
if(length(UniqValues)==1) stop("Only one unique value.")
if(!parallel){
R <- ncol(df)
N <- nrow(df)
c <- length(UniqValues)
KappaMatrix <- matrix(NA,c,c)
MinLambda = matrix(NA,c,c)
MaxLambda = matrix(NA,c,c)
colnames(KappaMatrix) <- rownames(KappaMatrix) <- UniqValues
colnames(MinLambda) = rownames(MinLambda) = UniqValues
colnames(MaxLambda) = rownames(MaxLambda) = UniqValues

for(j in UniqValues){
  # Kappa coefficients Kj
  Pj <- sum(apply(df,1,function(x) sum(x==j)))/(N*R)
  Kj <- ((sum(apply(df,1,function(x) sum(x==j)*(sum(x==j) - 1)))/(R*(R-1)*N)) - Pj^2)/(Pj*(1-Pj))
    # CI
    if(CI.est){
      Bt.j <- j
      BtResKj <- replicate(N.bt,BtSamples(df,fun=function(z,j=Bt.j){
        Pj <- sum(apply(z,1,function(x) sum(x==j)))/(N*R)
        Kj <- ((sum(apply(z,1,function(x) sum(x==j)*(sum(x==j) - 1)))/(R*(R-1)*N)) - Pj^2)/(Pj*(1-Pj))
        return(Kj)
      }))
      CI <- quantile(BtResKj,c(0.025,0.975))
      Kj <- paste(round(Kj,3)," (",round(CI[1],3)," to ",round(CI[2],3),")",sep="")
    }
  KappaMatrix[which(UniqValues==j),which(UniqValues==j)] <- Kj


  # Off-diagonal elements
  for(k in UniqValues){
    Pk <- sum(apply(df,1,function(x) sum(x==k)))/(N*R)
    if(which(UniqValues==j)<which(UniqValues==k)){
      CorrCoef = F
      Kjk <- (Pj*Pk - (sum(apply(df,1,function(x) sum(x==j)*sum(x==k)))/(R*(R-1)*N)))/(Pj*Pk)
      KappaMatrix[which(UniqValues==j),which(UniqValues==k)] <- Kjk
    }else if(which(UniqValues==j)>which(UniqValues==k)){
      CorrCoef = T
      Lambda_jk <- ((sum(apply(df,1,function(x) sum(x==j)*sum(x==k)))/(R*(R-1)*N)) - Pj*Pk)/(sqrt(Pj*(1-Pj)*Pk*(1-Pk)))
      KappaMatrix[which(UniqValues==j),which(UniqValues==k)] <- Lambda_jk
      MinLambda_jk = - ((Pj * Pk)/((1 - Pj) * (1 - Pk)))^0.5
      MaxLambda_jk = ((Pj * (1 - Pk))/(Pk * (1 - Pj)))^0.5
      MinLambda[which(UniqValues==j),which(UniqValues==k)] = MinLambda_jk
      MaxLambda[which(UniqValues==j),which(UniqValues==k)] = MaxLambda_jk
    }
    # CI
    if(CI.est & j!=k){
     Bt.k <- k
     Pj.o <- Pj
     BtResOffD <- replicate(N.bt,BtSamples(df,fun=function(z,j=Bt.j,k=Bt.k,Lambda=CorrCoef){
       Pj <- sum(apply(z,1,function(x) sum(x==j)))/(N*R)
       Pk <- sum(apply(z,1,function(x) sum(x==k)))/(nrow(z)*ncol(z))
       if(!Lambda) return((Pj*Pk - (sum(apply(z,1,function(x) sum(x==j)*sum(x==k)))/(R*(R-1)*N)))/(Pj*Pk))
       if(Lambda) return(((sum(apply(z,1,function(x) sum(x==j)*sum(x==k)))/(R*(R-1)*N)) - Pj*Pk)/(sqrt(Pj*(1-Pj)*Pk*(1-Pk))))
     }))
     CI2 <- quantile(BtResOffD,c(0.025,0.975),na.rm=T)
     if(!CorrCoef){
       KappaMatrix[which(UniqValues==j),which(UniqValues==k)] <- paste(round(Kjk,3)," (",round(CI2[1],3)," to ",
                                                                                         round(CI2[2],3),")",sep="")
     }else if(CorrCoef){
       KappaMatrix[which(UniqValues==j),which(UniqValues==k)] <- paste(round(Lambda_jk,3)," (",round(CI2[1],3)," to ",
                                                                                        round(CI2[2],3),")",sep="")
     }
    }
  }
}
}else{
KappaMatrix <- KappaMPar(df,CI.est=CI.est,N.bt=N.bt)
}
if(FleissKappa) Fl.k <- FleissK(df)[3]
if(!parallel){
if(FleissKappa) return(list(FleissKappa = Fl.k,KappaMatrix = KappaMatrix, MinLambda = MinLambda, MaxLambda = MaxLambda))
  else return(list(KappaMatrix = KappaMatrix, MinLambda = MinLambda, MaxLambda = MaxLambda))
}else{
  return(KappaMatrix)
}
}







