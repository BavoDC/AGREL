KappaMPar <- function(df,CI.est=F,N.bt=1e3){
UniqValues <- as.character(unique(unlist(df)))
R <- ncol(df)
N <- nrow(df)

NrCores = parallel::detectCores() - 1
doParallel::registerDoParallel(cores=NrCores)
KappaMatrix <- foreach(j=UniqValues,.combine="rbind") %do% {
  # Intraclass Kappa
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
    CI <- quantile(BtResKj,c(0.025,0.975),na.rm=T)
    Kj <- paste(round(Kj,3)," (",round(CI[1],3)," to ",round(CI[2],3),")",sep="")
  }

  # Off-diagonal elements
  InnerLoop <- foreach(k=UniqValues,.combine="c",.export="Bt.j") %dopar% {
    Pk <- sum(apply(df,1,function(x) sum(x==k)))/(N*R)
    if(which(UniqValues==j)<which(UniqValues==k)){
      CorrCoef = F
      Kjk <- (Pj*Pk - (sum(apply(df,1,function(x) sum(x==j)*sum(x==k)))/(R*(R-1)*N)))/(Pj*Pk)
      ResultsSecLoop <- Kjk
    }else if(which(UniqValues==j)>which(UniqValues==k)){
      CorrCoef = T
      Lambda_jk <- ((sum(apply(df,1,function(x) sum(x==j)*sum(x==k)))/(R*(R-1)*N)) - Pj*Pk)/(sqrt(Pj*(1-Pj)*Pk*(1-Pk)))
      ResultsSecLoop <- Lambda_jk
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
        ResultsSecLoop <- paste(round(Kjk,3)," (",round(CI2[1],3)," to ",round(CI2[2],3),")",sep="")
      }else if(CorrCoef){
        ResultsSecLoop <- paste(round(Lambda_jk,3)," (",round(CI2[1],3)," to ",round(CI2[2],3),")",sep="")
      }
    }
    if(j!=k) ResultsSecLoop else NA
  }
  InnerLoop[is.na(InnerLoop)] <- Kj
  InnerLoop
}
rownames(KappaMatrix) <- colnames(KappaMatrix) <- UniqValues
doParallel::stopImplicitCluster()
return(KappaMatrix)
}



