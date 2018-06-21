sumtable <- function(df, ratings=NULL, levels=NULL, offdiag=NULL){
  stopifnot(is.data.frame(df))
  if(is.null(ratings)){ratings=colnames(df)}
  stopifnot(is.character(ratings))
  stopifnot(all(unlist(lapply(df[ratings], function(x) any(class(x)%in%"factor")))))
  if(is.null(levels)){levels=levels(df[,ratings[1]])}
  if(is.null(offdiag)){if(length(levels==2)){offdiag=FALSE}
    if(length(levels>3)){offdiag=TRUE}}
  df <- df[ratings]
  nraters <- length(ratings)
  crosval <- t(combn(ratings,2))
  for (z in 1:nrow(crosval)) {
    v1 <- df[crosval[z,1]][,1]
    v1.1 <- factor(v1, levels = levels)
    v2 <- df[crosval[z,2]][,1]
    v2.1 <- factor(v2, levels = levels)
    Table <- table(v1.1, v2.1, dnn=NULL)
    crostab <- Table
    if (z == 1) {
      sumtable <- crostab
    }
    if (z > 1) {
      sumtable <- sumtable + crostab
    }
  }
  if(offdiag==TRUE){
    ## off diagonal means in matrix
    mat1 <-sumtable*0
    for (i in levels){ for (j in levels){
      mat1[i,j] <- (sumtable[i,j]+sumtable[j,i])/2
      mat1[j,i] <- (sumtable[i,j]+sumtable[j,i])/2
    }}
    sumtable <- mat1
  }
  Results        = list(SumTable=sumtable, nraters = nraters, N = nrow(df), data = df)
  class(Results) = append(class(Results),"AgreementTable")
  return(Results)
}

Agreement <- function(AgreemTable){
  stopifnot(any(class(AgreemTable)%in%"AgreementTable"))
  Table     = AgreemTable$SumTable
  stopifnot(nrow(Table)==ncol(Table))
  agreement = sum(diag(Table)/sum(Table))
  AgreemTable[["ObservedAgreem"]] = agreement
  Results          = AgreemTable
  class(Results)   = append(class(Results), "ObservedAgreem")
  return(Results)
}

SpecAgreem <- function(AgreemTable, specific=c("both","positive","negative")){
  stopifnot(any(class(AgreemTable)%in%"AgreementTable"))
  Table       = AgreemTable$SumTable
  specific    = match.arg(specific)
  stopifnot(nrow(Table)==2 & ncol(Table)==2)
  SpecAgr.Pos = if(specific%in%c("both","positive")) (2*(Table[1,1]))/((2*(Table[1,1]))+(Table[1,2])+(Table[2,1])) else NA
  SpecAgr.Neg = if(specific%in%c("both","negative")) (2*(Table[2,2]))/((2*(Table[2,2]))+(Table[1,2])+(Table[2,1])) else NA
  SpecAgr     = c(SpecAgr.Pos, SpecAgr.Neg)
  names(SpecAgr) = colnames(Table)
  SpecAgr        = as.data.frame(SpecAgr[!is.na(SpecAgr)])
  colnames(SpecAgr) = "SpecificAgreement"
  class(SpecAgr)    = append(class(SpecAgr),specific)
  AgreemTable[["SpecificAgreem"]] = SpecAgr
  Results            = AgreemTable
  class(Results)     = append(class(Results), "SpecificAgreem")
  return(Results)
}

CIAgreem <- function(AgreemTable, level=0.95, AgreemStat=c("all","observed", "positive", "negative"), correction=c("continuity","Fleiss", "bootstrap"),
                     NrBoot = 1e3, Parallel=F, no_cores = detectCores() - 1){
  stopifnot(any(class(AgreemTable)%in%"AgreementTable"))
  AgreemStat = match.arg(AgreemStat)
  correction = match.arg(correction)
  alpha      = 1 - level
  TypeAgr    = if(identical(AgreemStat,"observed")) "ObservedAgreem" else "SpecificAgreem"
  if(!any(class(AgreemTable)%in%TypeAgr)) stop(paste(TypeAgr," has not yet been calculated."))
  m        = AgreemTable$nraters

  # Observed agreement
  if(any(AgreemStat%in%c("all","observed"))){
    p        = AgreemTable$ObservedAgreem
    n        = AgreemTable$N
    Agrm.CI  = if(!identical(correction, "bootstrap")){
      CI.Agreem(p, n, m, level, correction)
    }else{
        CI.Agreem.bt(AgreemTable, "observed", alpha, NrBoot, Parallel, no_cores)
    }
    names(Agrm.CI) = c("LowerCI","ObservedAgreement","UpperCI")
    # cat("\n\nObserved agreement:\n\n")
    # print(Agrm.CI)
    AgreemTable$ObservedAgreem = Agrm.CI
  }

  ResT = NULL
  # Positive agreement
  if(any(AgreemStat%in%c("all","positive"))){
    if(any(class(AgreemTable$SpecificAgreem)%in%"negative")) stop("Positive agreement has not yet been calculated")
    p        = AgreemTable$SpecificAgreem[1,]
    NewTab   = addmargins(AgreemTable$SumTable)
    n.total  = AgreemTable$N
    n        = NewTab[1,3]/NewTab[3,3] * n.total
    Agrm.CI  = if(!identical(correction, "bootstrap")){
      CI.Agreem(p, n, m, level, correction)
    }else{
      CI.Agreem.bt(AgreemTable, "positive", alpha, NrBoot, Parallel, no_cores)
    }
    TempF = AgreemTable$SpecificAgreem
    ResT  = if(nrow(TempF)!=1) rbind.data.frame(Agrm.CI,c(NA,TempF[2,],NA)) else as.data.frame(t(Agrm.CI))
    colnames(ResT) = c("LowerCI","SpecificAgreement","UpperCI")
    rownames(ResT) = rownames(TempF)
  }

  # Negative agreement
  if(any(AgreemStat%in%c("all","negative"))){
    if(any(class(AgreemTable$SpecificAgreem)%in%"positive")) stop("Negative agreement has not yet been calculated")
    TypeNr   = if(any(class(AgreemTable$SpecificAgreem)%in%"both")) 2 else 1
    p        = AgreemTable$SpecificAgreem[TypeNr,]
    NewTab   = addmargins(AgreemTable$SumTable)
    n.total  = AgreemTable$N
    n        = NewTab[2,3]/NewTab[3,3] * n.total
    Agrm.CI  = if(!identical(correction, "bootstrap")){
      CI.Agreem(p, n, m, level, correction)
    }else{
      CI.Agreem.bt(AgreemTable, "negative", alpha, NrBoot, Parallel, no_cores)
    }
    if(identical(AgreemStat, "all")){
      ResT[2,] = Agrm.CI
    }else{
      TempF = AgreemTable$SpecificAgreem
      ResT  = if(nrow(TempF)!=1) rbind.data.frame(c(NA,TempF[1,],NA),Agrm.CI) else as.data.frame(Agrm.CI)
      rownames(ResT) = rownames(TempF)
    }
    colnames(ResT) = c("LowerCI","SpecificAgreement","UpperCI")
  }
  if(!is.null(ResT)){
    # cat("\n\nSpecific agreement:\n\n")
    # print(ResT)
    AgreemTable$SpecificAgreem = ResT
  }
  invisible(AgreemTable)
}

CI.Agreem <- function(p, n, m, level=0.95, correction=c("continuity","Fleiss")){
  stopifnot(p>=0 | p<=1)
  stopifnot(is.numeric(p))
  stopifnot(is.numeric(n))
  stopifnot(is.numeric(m))
  stopifnot(level>=0|level<=1)
  correction = match.arg(correction)
  a <- qnorm(1-((1-level)/2))
  n <- n*sqrt(m-1)
  if (correction=="continuity"){
    CIagreement<- p+c(-(a*(sqrt(1/n*(p*(1-p))))-1/(2*n)),0, (a*(sqrt(1/n*(p*
                                                                            (1-p))))+1/(2*n)))
    CIagreement <- c(CIlow=CIagreement[1], p=CIagreement[2], CIhigh=CIagreement[3])
  }else{
    FCIlow <- ((2*n*p+(a*a)-1)-a*sqrt((a*a)-(2+(1/n))+4*p*(n*(1-p)+1)))/(2*
                                                                           ((a*a)+n))
    FCIhigh <- ((2*n*p+(a*a)-1)+a*sqrt((a*a)-(2+(1/n))+4*p*(n*(1-p)+1)))/(2
                                                                          *((a*a)+n))
    CIagreement <- c(CIlow=FCIlow, p=p, CIhigh=FCIhigh)
  }
  CIagreement
}


CI.Agreem.bt <- function(AgreemTable, Agreem.stat = c("observed", "positive", "negative"), alpha = 0.05, NrBoot = 1e3,
                         Parallel = F, no_cores = detectCores() - 1){
  df  = AgreemTable$data
  p   = if(identical(Agreem.stat,"observed")) AgreemTable$ObservedAgreem else SpecAgreem(AgreemTable, Agreem.stat)$SpecificAgreem[1,]
  if(!Parallel){
  BtAgrm <- function(df, Agreem.stat){
    df.bt       = df[sample(1:nrow(df),nrow(df),T),]
    TabSum      = sumtable(df.bt)
    Agreem      = if(identical(Agreem.stat,"observed")) Agreement(TabSum)$ObservedAgreem else SpecAgreem(TabSum, Agreem.stat)$SpecificAgreem
    if(is.data.frame(Agreem)) Agreem = Agreem[1,]
    return(Agreem)
  }
  BtSamples   = replicate(NrBoot, BtAgrm(df,Agreem.stat))
  }else{
    cl       = makeCluster(no_cores)
    clusterExport(cl,c("df","Agreem.stat"), envir = environment())
    clusterEvalQ(cl,c("library(AGREL)"))
    BtSamples = parSapply(cl, 1:NrBoot,function(x){
      df.bt       = df[sample(1:nrow(df),nrow(df),T),]
      TabSum      = sumtable(df.bt)
      Agreem      = if(identical(Agreem.stat,"observed")) Agreement(TabSum)$ObservedAgreem else SpecAgreem(TabSum, Agreem.stat)$SpecificAgreem
      if(is.data.frame(Agreem)) Agreem = Agreem[1,]
      return(Agreem)
    })
    stopCluster(cl)
  }
  CIagreement = quantile(as.vector(BtSamples),c(alpha/2,1-alpha/2))
  CIagreement = c(CIagreement[1], p=p, CIagreement[2])
  names(CIagreement)[c(1,3)] = c("CIlow","CIhigh")
  return(CIagreement)
}

