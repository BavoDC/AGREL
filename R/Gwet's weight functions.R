identity.weights<-function(categ){
  weights<-diag(length(categ))
  return (weights)
}

quadratic.weights<-function(categ){
  q<-length(categ)
  weights <- diag(q)
  if (is.numeric(categ)) {
    categ.vec <- sort(categ)
  }
  else {
    categ.vec<-1:length(categ)
  }
  xmin<-min(categ.vec)
  xmax<-max(categ.vec)
  for(k in 1:q){
    for(l in 1:q){
      weights[k,l] <- 1-(categ.vec[k]-categ.vec[l])^2/(xmax-xmin)^2
    }
  }
  return (weights)
}

linear.weights<-function(categ){
  q<-length(categ)
  weights <- diag(q)
  if (is.numeric(categ)) {
    categ.vec <- sort(categ)
  }
  else {
    categ.vec<-1:length(categ)
  }
  xmin<-min(categ.vec)
  xmax<-max(categ.vec)
  for(k in 1:q){
    for(l in 1:q){
      weights[k,l] <- 1-abs(categ.vec[k]-categ.vec[l])/abs(xmax-xmin)
    }
  }
  return (weights)
}
#--------------------------------
radical.weights<-function(categ){
  q<-length(categ)
  weights <- diag(q)
  if (is.numeric(categ)) {
    categ.vec <- sort(categ)
  }
  else {
    categ.vec<-1:length(categ)
  }
  xmin<-min(categ.vec)
  xmax<-max(categ.vec)
  for(k in 1:q){
    for(l in 1:q){
      weights[k,l] <- 1-sqrt(abs(categ.vec[k]-categ.vec[l]))/sqrt(abs(xmax-xmin))
    }
  }
  return (weights)
}

#--------------------------------
ratio.weights<-function(categ){
  q<-length(categ)
  weights <- diag(q)
  if (is.numeric(categ)) {
    categ.vec <- sort(categ)
  }
  else {
    categ.vec<-1:length(categ)
  }
  xmin<-min(categ.vec)
  xmax<-max(categ.vec)
  for(k in 1:q){
    for(l in 1:q){
      weights[k,l] <- 1-((categ.vec[k]-categ.vec[l])/(categ.vec[k]+categ.vec[l]))^2 / ((xmax-xmin)/(xmax+xmin))^2
    }
  }
  return (weights)
}

#--------------------------------
circular.weights<-function(categ){
  q<-length(categ)
  weights <- diag(q)
  if (is.numeric(categ)) {
    categ.vec <- sort(categ)
  }
  else {
    categ.vec<-1:length(categ)
  }
  xmin<-min(categ.vec)
  xmax<-max(categ.vec)
  U = xmax-xmin+1
  for(k in 1:q){
    for(l in 1:q){
      weights[k,l] <- (sin(pi*(categ.vec[k]-categ.vec[l])/U))^2
    }
  }
  weights <- 1-weights/max(weights)
  return (weights)
}

#--------------------------------
bipolar.weights<-function(categ){
  q<-length(categ)
  weights <- diag(q)
  if (is.numeric(categ)) {
    categ.vec <- sort(categ)
  }
  else {
    categ.vec<-1:length(categ)
  }
  xmin<-min(categ.vec)
  xmax<-max(categ.vec)
  for(k in 1:q){
    for(l in 1:q){
      if (k!=l)
        weights[k,l] <- (categ.vec[k]-categ.vec[l])^2 / (((categ.vec[k]+categ.vec[l])-2*xmin)*(2*xmax-(categ.vec[k]+categ.vec[l])))
      else weights[k,l] <- 0
    }
  }
  weights <- 1-weights/max(weights)
  return (weights)
}


#--------------------------------
ordinal.weights<-function(categ){
  q<-length(categ)
  weights <- diag(q)
  categ.vec<-1:length(categ)
  for(k in 1:q){
    for(l in 1:q){
      nkl <- max(k,l)-min(k,l)+1
      weights[k,l] <- nkl * (nkl-1)/2
    }
  }
  weights <- 1-weights/max(weights)
  return (weights)
}
