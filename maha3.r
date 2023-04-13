### this is a modified version of R mahalanobis distance calculation
maha3<-function (x, center, cov, inverted = FALSE, ...) 
{
 # library(MASS)
  require(MASS)
  x <- if (is.vector(x)) 
    matrix(x, ncol = length(x))
  else as.matrix(x)
  if (!isFALSE(center)) 
    x <- sweep(x, 2L, center)
  if (!inverted) 
    pcov <- ginv(cov)  # here we use ginv instead of inv
  #setNames(rowSums(x %*% pcov %*% t(x)), rownames(x))
    pp<-x %*% pcov %*% t(x)
    result<-diag(pp)
}
