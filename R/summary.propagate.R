summary.propagate <- function(object, ...)
{
  ## print error propagation results
  message("Results from error propagation:")
  print(object$prop)
  
  ## print simulation results
  if (is.finite(object$sim[1])) {
    message("Results from Monte Carlo simulation:")
    print(object$sim)
  }
  
  ## print WS degrees-of-freedom
  message("Welch-Satterthwaite degrees of freedom:")
  print(object$ws.df)
  
  ## print coverage factor
  message("Coverage factor (k):")
  print(object$k)
  
  ## print coverage factor
  message("Expanded uncertainty:")
  print(object$u.exp)
  
  ## print symbolic derivatives of gradient
  if (is.matrix(object$gradient)) {
    message("Symbolic gradient matrix:")
    print(as.character(object$gradient))
  }
  
  ## print evaluated derivatives of gradient
  if (is.matrix(object$gradient)) {
    message("Evaluated gradient matrix (sensitivity):")
    print(object$evalGrad)
  }
   
  ## print symbolic derivatives of hessian
  if (is.matrix(object$hessian)) {
    message("Symbolic hessian matrix:")
    NCOL <- ncol(object$hessian)    
    STR <- as.character(object$hessian) 
    SPLIT <- split(STR, as.factor(rep(1:NCOL, NCOL)))
    for (i in 1:length(SPLIT)) print(SPLIT[[i]])
  }
    
  ## print evaluated derivatives of hessian
  if (is.matrix(object$hessian)) {
    message("Evaluated hessian matrix:")
    NCOL <- ncol(object$hessian)
    EVAL <- object$evalHess 
    SPLIT <- split(EVAL, as.factor(rep(1:NCOL, NCOL)))
    for (i in 1:length(SPLIT)) print(SPLIT[[i]]) 
  }
  
  ## print covariance matrix
  message("Covariance matrix:")
  print(object$covMat)
  
  ## print rel. contribution
  message("Relative contribution:")
  print(object$rel.contr)
  
  ## Skewness and excess kurtosis of evaluated MC simulations
  if (length(object$resSIM) > 1) {
    message("Skewness / Excess Kurtosis of MC evaluations:")
    cat(skewness(object$resSIM), "/", kurtosis(object$resSIM), "\n")
  }
    
  ## Shapiro test for normality of MC distribution
  if (length(object$resSIM) > 1) {
    message("Shapiro-Wilk test for normality: ")
    if (length(object$resSIM) > 5000) 
    DAT <- object$resSIM[1:5000] else DAT <- object$resSIM 
    PVAL <- shapiro.test(DAT)$p.value
    cat(PVAL)
    if (PVAL >= 0.05) cat(" => normal\n") else cat(" => non-normal\n")
  }
    
  ## Kolmogorov-Smirnov test for normality of MC distribution
  if (length(object$resSIM) > 1) {
    message("Kolmogorov-Smirnov test for normality: ")
    DAT <- object$resSIM
    simDAT <- rnorm(length(DAT), mean(DAT, na.rm = TRUE), sd(DAT, na.rm = TRUE))
    PVAL <- ks.test(DAT, simDAT)$p.value
    cat(PVAL)
    if (PVAL >= 0.05) cat(" => normal\n") else cat(" => non-normal\n")  
  }
}