#include <Rcpp.h>
// [[Rcpp::export]]
Rcpp::NumericMatrix KappaDiag(Rcpp::NumericMatrix df, Rcpp::NumericMatrix ResultsK,int NrCat,int Nj, int N, int R){
  Rcpp::NumericMatrix PjMat;
  Rcpp::NumericMatrix Kj(NrCat,1);
  int nij;
  Rcpp::NumericMatrix Varnij;
  int sumnj;
  double Pj;
  double sumVarnij;

  for(int j; j < NrCat; j++){
    for(int i; i < Nj; i++){
      for(int col; col < R; col++){
        nij    += (df(i,col)==j);
      }
      PjMat(i,j)  = nij;
      Varnij(i,j) = nij*(1-nij);
    }
    for(int i; i < Nj; i++){
      sumnj     += PjMat(i,j);
      sumVarnij += Varnij(i,j);
    }
    Pj      = sumnj/(N*R);
    Kj(j,1) = (sumVarnij/(R*(R-1)*N) - pow(Pj,2))/(Pj*(1-Pj));
  }
  return(Kj);
//for(j in UniqValues){
//  # Kappa coefficients Kj
//  Pj <- sum(apply(df,1,function(x) sum(x==j)))/(N*R)
//  Kj <- ((sum(apply(df,1,function(x) sum(x==j)*(sum(x==j) - 1)))/(R*(R-1)*N)) - Pj^2)/(Pj*(1-Pj))
//  # CI
//  if(CI.est){
//    Bt.j <- j
//    BtResKj <- replicate(N.bt,BtSamples(df,fun=function(z,j=Bt.j){
//      Pj <- sum(apply(z,1,function(x) sum(x==j)))/(N*R)
//      Kj <- ((sum(apply(z,1,function(x) sum(x==j)*(sum(x==j) - 1)))/(R*(R-1)*N)) - Pj^2)/(Pj*(1-Pj))
//      return(Kj)
//    }))
//    CI <- quantile(BtResKj,c(0.025,0.975))
//    Kj <- paste(round(Kj,3)," (",round(CI[1],3)," to ",round(CI[2],3),")",sep="")
//  }
//  KappaMatrix[which(UniqValues==j),which(UniqValues==j)] <- Kj//

//  # Off-diagonal elements
//  for(k in UniqValues){
//    Pk <- sum(apply(df,1,function(x) sum(x==k)))/(N*R)
//    if(which(UniqValues==j)<which(UniqValues==k)){
//      CorrCoef = F
//      Kjk <- (Pj*Pk - (sum(apply(df,1,function(x) sum(x==j)*sum(x==k)))/(R*(R-1)*N)))/(Pj*Pk)
//      KappaMatrix[which(UniqValues==j),which(UniqValues==k)] <- Kjk
//    }else if(which(UniqValues==j)>which(UniqValues==k)){
//      CorrCoef = T
//      Lambda_jk <- ((sum(apply(df,1,function(x) sum(x==j)*sum(x==k)))/(R*(R-1)*N)) - Pj*Pk)/(sqrt(Pj*(1-Pj)*Pk*(1-Pk)))
//      KappaMatrix[which(UniqValues==j),which(UniqValues==k)] <- Lambda_jk
//    }
//    # CI
//    if(CI.est & j!=k){
//      Bt.k <- k
//      Pj.o <- Pj
//      BtResOffD <- replicate(N.bt,BtSamples(df,fun=function(z,j=Bt.j,k=Bt.k,Lambda=CorrCoef){
//        Pj <- sum(apply(z,1,function(x) sum(x==j)))/(N*R)
//        Pk <- sum(apply(z,1,function(x) sum(x==k)))/(nrow(z)*ncol(z))
//        if(!Lambda) return((Pj*Pk - (sum(apply(z,1,function(x) sum(x==j)*sum(x==k)))/(R*(R-1)*N)))/(Pj*Pk))
//        if(Lambda) return(((sum(apply(z,1,function(x) sum(x==j)*sum(x==k)))/(R*(R-1)*N)) - Pj*Pk)/(sqrt(Pj*(1-Pj)*Pk*(1-Pk))))
//      }))
//      CI2 <- quantile(BtResOffD,c(0.025,0.975),na.rm=T)
//      if(!CorrCoef){
//        KappaMatrix[which(UniqValues==j),which(UniqValues==k)] <- paste(round(Kjk,3)," (",round(CI2[1],3)," to ",
//                                                                        round(CI2[2],3),")",sep="")
//      }else if(CorrCoef){
//        KappaMatrix[which(UniqValues==j),which(UniqValues==k)] <- paste(round(Lambda_jk,3)," (",round(CI2[1],3)," to ",
//                                                                        round(CI2[2],3),")",sep="")
//      }
//    }
//  }
//}
}
