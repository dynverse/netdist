## Normal correlation (default pearson)
Adjcor <- function(x,method='pearson',...){

    Adj <- abs(doCall(cor, x=x, method=method, ...))
    diag(Adj) <- 0
    
    return(Adj)
}

## Function for check the variance by features
checkvar <- function(x, tol=1e-5, ...){
  
  ## Compute the variance by columns
  feat.var <- apply(x,2,var)
  
  ## Get the indexes of the feature with low variance
  idx <- which(feat.var<tol)
  if (length(idx) == 0L){
    return (NULL)
  } else {
    return (idx)
  }
}
