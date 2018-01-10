propagate <- function(
expr, 
data, 
second.order = TRUE,
do.sim = TRUE, 
cov = TRUE, 
df = NULL,
nsim = 100000,
alpha = 0.05,
...
)
{            
  op <- options(warn = -1)
  on.exit(options(op))
  
  ## version 1.0-4: convert function to expression
  if (is.function(expr)) {
    ARGS <- as.list(args(expr))
    ARGS <- ARGS[-length(ARGS)]
    VARS <- names(ARGS)    
    expr <- body(expr)
    class(expr) <- "expression"
    isFun <- TRUE
  } else isFun <- FALSE
  
  ## check for correct expression and number of simulations
  if (!is.expression(expr)) stop("propagate: 'expr' must be an expression")
  if (nsim < 10000) stop("propagate: 'nsim' should be >= 10000 !")
  
  ## check for matching variable names
  if (!isFun) VARS <- all.vars(expr)  
  m <- match(VARS, colnames(data))
  if (any(is.na(m))) stop("propagate: variable names of input dataframe and expression do not match!")
  if (length(unique(m)) != length(m)) stop("propagate: some variable names are repetitive!")
  
  DATA <- as.matrix(data)
  EXPR <- expr
  
  ## create variables from input data
  ## version 1.0-5: check for simulated data => large nrow(data)
  if (nrow(data) > 3) {
    meanVAL <- apply(DATA, 2, function(x) mean(x, na.rm = TRUE))
    sdVAL <- apply(DATA, 2, function(x) sd(x, na.rm = TRUE))
    dfVAL <- NULL
    isRaw <- TRUE
  } else {
    meanVAL <- DATA[1, ]
    sdVAL <- DATA[2, ] 
    dfVAL <- if (nrow(DATA) == 3) DATA[3, ] else NULL
    isRaw <- FALSE
  }
       
  ## stat data: if no covariance matrix is supplied, create one with diagonal variances
  if (is.logical(cov) & !isRaw) {
    SIGMA <- diag(sdVAL^2, nrow = length(VARS), ncol = length(VARS))    
    colnames(SIGMA) <- rownames(SIGMA) <- colnames(DATA)
  } 
  
  ## raw data: if no covariance matrix is supplied, create one with off-diagonals or not
  if (cov & isRaw) {
    SIGMA <- cov(data)   
  } 
  if (!cov & isRaw) {
    SIGMA <- cov(data) 
    SIGMA[upper.tri(SIGMA)] <- SIGMA[lower.tri(SIGMA)] <- 0 
  } 
  
  ## if covariance matrix is supplied, check for symmetry and matching names
  if (is.matrix(cov)) {
    if (NCOL(cov) != NROW(cov)) stop("propagate: 'cov' is not a symmetric matrix!")
    m <- match(colnames(cov), colnames(DATA))            
    if (any(is.na(m))) stop("propagate: names of input dataframe and var-cov matrix do not match!")             
    if (length(unique(m)) != length(m)) stop("propagate: some names of the var-cov matrix are repetitive!")             
    SIGMA <- cov
  }
  
  ## version 1.0-5: replace possible NA's in covariance matrix with 0's
  SIGMA[is.na(SIGMA)] <- 0  
  
  ## version 1.0-5: No diagonals with 0, 
  ## otherwise tmvtnorm:::checkSymmetricPositiveDefinite throws an error!
  if (any(diag(SIGMA) == 0)) {
    DIAG <- diag(SIGMA)
    DIAG[DIAG == 0] <- 2E-16
    diag(SIGMA) <- DIAG
  }
  
  ## This will bring the variables in 'data' and 'expr' in the same 
  ## order as in the covariance matrix
  m1 <- match(colnames(SIGMA), colnames(DATA))
  meanVAL <- meanVAL[m1]
  m2 <- match(colnames(SIGMA), VARS)
  
  ############ first- and second-order Taylor expansion-based error propagation ################
  ## first-order mean: eval(EXPR)
  ## version 1.0-4: continue with NA's when differentiation not possible
  MEAN1 <- try(eval(EXPR, envir = as.list(meanVAL)), silent = TRUE)
  if (!is.numeric(MEAN1)) {
    message("propagate: there was an error in calculating the first-order mean")
    MEAN1 <- NA
  }  
  
  ## evaluate gradient vector
  GRAD <- try(makeGrad(EXPR, m2), silent = TRUE)  
  if (!inherits(GRAD, "try-error")) evalGRAD <- try(sapply(GRAD, eval, envir = as.list(meanVAL)), silent = TRUE)
  if (inherits(GRAD, "try-error")) evalGRAD <- try(numGrad(EXPR, as.list(meanVAL)), silent = TRUE)  
  if (!inherits(evalGRAD, "try-error")) evalGRAD <- as.vector(evalGRAD) else evalGRAD <- NA
  
  ## first-order variance: g.S.t(g) 
  VAR1 <- try(as.numeric(t(evalGRAD) %*% SIGMA %*% matrix(evalGRAD)), silent = TRUE)  
  if (inherits(VAR1, "try-error")) {
    message("propagate: there was an error in calculating the first-order variance")
    VAR1 <- NA
  }
 
  ## second-order mean: firstMEAN + 0.5 * tr(H.S) 
  if (second.order) {
    HESS <- try(makeHess(EXPR, m2), silent = TRUE)
    if (!inherits(HESS, "try-error")) evalHESS <- try(sapply(HESS, eval, envir = as.list(meanVAL)), silent = TRUE)
    if (inherits(HESS, "try-error")) evalHESS <- try(numHess(EXPR, as.list(meanVAL)), silent = TRUE)
    if (!inherits(evalHESS, "try-error")) evalHESS <- matrix(evalHESS, ncol = length(meanVAL), byrow = TRUE) else evalHESS <- NA  
    
    valMEAN2 <- try(0.5 * tr(evalHESS %*% SIGMA), silent = TRUE)
    if (!inherits(valMEAN2, "try-error")) {
      MEAN2 <- MEAN1 + valMEAN2
    } else {
      message("propagate: there was an error in calculating the second-order mean")
      MEAN2 <- NA
    }
    
    ## second-order variance: firstVAR + 0.5 * tr(H.S.H.S)
    valVAR2 <- try(0.5 * tr(evalHESS %*% SIGMA %*% evalHESS %*% SIGMA), silent = TRUE)
    if (!inherits(valVAR2, "try-error")) {
      VAR2 <- VAR1 + valVAR2
    } else {
      message("propagate: there was an error in calculating the second-order variance")
      VAR2 <- NA
    }
  } else MEAN2 <- VAR2 <- HESS <- evalHESS <- NA
  
  ## total mean and variance  
  if (second.order) totalVAR <- VAR2 else totalVAR <- VAR1
  if (second.order) totalMEAN <- MEAN2 else totalMEAN <- MEAN1
  errorPROP <- sqrt(totalVAR)  
  
  ## sensitivity index/contribution/relative contribution
  if (is.numeric(evalGRAD)) {
    sensitivity <- evalGRAD
    contribution <- outer(sensitivity, sensitivity, "*") * SIGMA
    rel.contribution <- abs(contribution)/sum(abs(contribution), na.rm = TRUE)
  } else sensitivity <- contribution <- rel.contribution <- NA
  
  ## WS degrees of freedom, coverage factor and expanded uncertainty
  if (!is.null(dfVAL)) dfVAL[is.na(dfVAL)] <- 1E6
  if (is.null(dfVAL)) dfVAL <- rep(1E6, ncol(DATA))
  ws <- WelchSatter(ui = sqrt(diag(SIGMA)), ci = sensitivity, df = dfVAL, dftot = df, uc = errorPROP, alpha = alpha)
  
  ## confidence interval based on either first- or second-order mean
  if (is.na(MEAN2)) confMEAN <- MEAN1 else confMEAN <- MEAN2
  confPROP <- confMEAN + c(-1, 1) * ws$u.exp
  names(confPROP) <- paste(c(alpha/2, 1 - alpha/2) * 100, "%", sep = "")
  
  ################## Monte-Carlo simulation using multivariate t-distribution #####################
  if (do.sim) {  
    if (is.na(ws$ws.df) | is.infinite(ws$ws.df)) DF <- 1E6 else DF <- ws$ws.df
    if (is.numeric(df)) DF <- 1E6
    
    ## if raw data, don't create Monte Carlo data
    if (!isRaw) {
      datSIM <- rtmvt(nsim, mean = meanVAL, sigma = SIGMA, df = floor(DF)) 
      colnames(datSIM) <- colnames(DATA)      
    } else datSIM <- DATA
      
    ## try vectorized evaluation, which is much faster  
    resSIM <- try(eval(EXPR, as.data.frame(datSIM)), silent = TRUE) 
    
    ## use 'row-wise' method if 'vectorized' throws an error
    if (inherits(resSIM, "try-error")) {
      message("propagate: using 'vectorized' evaluation gave an error. Switching to 'row-wise' evaluation...")
      resSIM <- apply(datSIM, 1, function(x) eval(EXPR, envir = as.list(x)))     
    }
    
    ## alpha-based confidence interval of MC simulations
    confSIM <- quantile(resSIM, c(alpha/2, 1 - (alpha/2)), na.rm = TRUE) 
    
    ## warning in case of single evaluated result
    if (length(unique(resSIM)) == 1) message("propagate: Monte Carlo simulation gave unique repetitive values! Are all derivatives constants?")   
  } else resSIM <- datSIM <- confSIM <- allSIM <- NA 
  
  
  outPROP <- c(Mean.1 = MEAN1, Mean.2 = MEAN2, sd.1 = sqrt(VAR1), sd.2 = sqrt(VAR2), 
               confPROP[1], confPROP[2])   
  
  outSIM <- c(Mean = mean(resSIM, na.rm = TRUE), sd = sd(resSIM, na.rm = TRUE), 
              Median = median(resSIM, na.rm = TRUE), MAD = mad(resSIM, na.rm = TRUE),
              confSIM[1], confSIM[2])
  
  OUT <- list(gradient = GRAD, evalGrad = evalGRAD,
              hessian = HESS, evalHess = evalHESS,
              rel.contr  = rel.contribution, covMat = SIGMA, ws.df = floor(ws$ws.df), 
              k = ws$k, u.exp = ws$u.exp, resSIM = resSIM, datSIM = datSIM, 
              prop = outPROP, sim = outSIM, expr = EXPR, data = DATA, alpha = alpha)
  
  class(OUT) <- "propagate"
  return(OUT)                                     
}

