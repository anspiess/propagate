cor2cov <- function(C, var = NULL)
{
  if (is.null(var)) stop("cor2cov: cannot calculate covariance matrix without variances")
  if (ncol(C) != nrow(C)) stop("cor2cov: 'C' is not a square numeric matrix!")
  if (length(var) != ncol(C)) stop("cor2cov: length of 'var' and dimension of 'C' are not equal!")
  if (any(!is.finite(var))) warning("cor2cov: 'var' had 0 or NA entries; result is doubtful!")
  d <- sqrt(var)
  V <- outer(d, d) * C
  return(V) 
}

