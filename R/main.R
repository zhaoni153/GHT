GHT = function(X, u = 0, nsim = 999, targetMatrix){  
  centered = F  # whether the data has been centered with mean 0, should be FALSE.
  B = targetMatrix
  getS = function(X, B, centered){    
    if (B == "equal"){
      S = shrinkcovmat.equal(X, centered = centered)
    } 
    if (B == "identity"){
      S = shrinkcovmat.identity(X, centered = centered)
    } 
    if (B == "unequal"){
      S = shrinkcovmat.unequal(X, centered = centered)
    } 
    return(S)
  }
  tX = t(X)
  S = getS(tX, B = B, centered = centered)  
  
  Sigmahat = S$Sigmahat
  lambda = S$lambda
  Sigmasample = S$Sigmasample
  Target = S$Target
  p=ncol(X) 
  N=nrow(X)
  n=N-1
  means =colMeans(X, na.rm=TRUE)
  xx = (means-u)%*% t(means-u)
  T1 = sum(N*xx*solve(Sigmahat))
  #  T1=sum((means-u)^2)*N
  # direct permutation approach. 
  stat0 = rep(NA, nsim)
  for (j in 1:nsim){
    permute = sample(c(-1,1), size = N, replace = T, prob = c(0.5, 0.5)) 
    X0 = X * permute
    tX0 = t(X0)
    S2 = getS(tX0, B = B, centered = centered)
    means =colMeans(X0, na.rm=TRUE)    
    stat0[j]=as.numeric(N*t(means-u)%*%solve(S2$Sigmahat)%*%(means-u))    
    
  }
  
  pval <- (sum(stat0>=T1)+ 1)/nsim 
  pval <- ifelse(pval > 1, 1, pval)  # pvalues within 1/nsim and 1
  return(pval)
} 