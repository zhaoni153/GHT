log_geo_mean <- function(dat) {
  log_data <- log(dat)
  log_gm <- mean(log_data[is.finite(log_data)])
  return(log_gm)
}

#' getCLR
#' Obtain the centered log ratio transformation of compositional data. 
#'
#' @param dat n x p matrix of numerical variables. For microbiome data, n: sample size, p: number of taxa 
#' @return n x p matrix of the log-ratio transformed data.
#' @details A pseudocount is added to all taxa if zero values exist. A pseudocount of 0.01 is added if "dat" is read counts (rarefied or not), 
#' and a pseudocount that is 0.01 times the smallest nonzero relative abundances is added if "dat" is already proportions. 
#' @examples n = 10; p = 100
#' set.seed(1)
#' dat = matrix(rbinom(n*p, size = 1000,  0.005),n,p)
#' X = getCLR(dat)
#' @export
getCLR = function(dat){
  # dat: n X p matrix of relative abundances, n: sample size, p: number of taxa
  # data should be rarefied or proportional 
  n = nrow(dat); p = ncol(dat)
  
  if (any(dat == 0)){
    if (all(rowSums(dat) > 1)){ 
      # read counts. Assuming the total read count for all samples are large. 
      delta = 0.01
      dat = dat + delta
      # t(apply(dat, 1, imputeZero, delta = delta))
      print("A pseudocount of 0.01 added to all reads")   
    }else if(all(rowSums(dat) == 1)){  # already proportions
      delta = min(dat[dat > 0])*0.01
      dat = dat  + delta
      #  t(apply(dat, 1, imputeZero, delta = delta))
      print(sprintf("A pseudocount of %s, which is
                    0.01 times the smallest nonzero values, added to all reads", delta))
    }
    }
  
  
  # print("data is converted to percentages.")
  dat = dat/rowSums(dat)
  log_geoMean = apply(dat, 1, log_geo_mean)
  logclr = log(dat) - log_geoMean
  # logclr = logclr[,-NCOL(logclr)]
  return(logclr)  
  }

#' Generalized Hoteling's test 
#' The fuction testing whether the mean of X is differnt from u, a hypothesized population average. 
#' GHT replaces the sample covariance matrix in classical Hoteling's test with a shrinkage based (positve definite) covariance matrix.
#' Significance is evaluated via permuation. 
#' The method is designed for paired microbiome studies, in which X is the paired differences after log-ratio transformation. 
#' However, the method is equally applicable to other high dimensional settings.  
#' 
#' @param X n x p matrix of numerical variables
#' @param u a vector of numerical variables indicating the true value of the mean
#' @param nsim number of permutations. "equal": diagonal matrix with equal diagonal element;"unequal", diagonal matrix with unequal diagonal element;"identity": identity matrix
#' @param target target matrix for covariance estimate. 
#' @return p value
#' @examples set.seed(1)
#' n=10; p=100
#' dat = matrix(rnorm(n*p),n,p)
#' test1 = GHT(dat)
#' # A test similar to paired microbiome data
#' set.seed(1)
#' dat1 = matrix(rbinom(n*p, size = 1000,  0.005),n,p)
#' dat2 = matrix(rbinom(n*p, size = 1000,  0.005),n,p)
#' X1 = CLR(dat1); X2 = CLR(dat2)
#' X = X1 - X2
#' test2 = GHT(X)
#' @export
GHT = function(X, u = 0, nsim = 1000, target = "equal", centered = F){  
  B = target
  getS = function(X, B, centered){   
    # Touloumis 2015 approach for covariance estimate 
    # centered: whether the data centered to have mean zero. Becase we aim to test the means, centered = F
    # B: target matrix
    if (B == "equal"){
      S = ShrinkCovMat::shrinkcovmat.equal(X, centered = centered)
    } 
    if (B == "identity"){
      S = ShrinkCovMat::shrinkcovmat.identity(X, centered = centered)
    } 
    if (B == "unequal"){
      S = ShrinkCovMat::shrinkcovmat.unequal(X, centered = centered)
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
  pval <- ifelse(pval > 1, 1, pval)  
  return(pval)
} 