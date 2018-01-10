fitDistr <- function(
  object, 
  nbin = 100,
  weights = FALSE,
  verbose = TRUE, 
  brute = c("fast", "slow"),
  plot = c("hist", "qq"),  
  distsel = NULL,
  ...)
{
  op <- options(warn = -1)
  on.exit(options(op))
  plot <- match.arg(plot)
  brute <- match.arg(brute)
  
  if (is.vector(object)) X <- object
  else if (class(object) == "propagate") X <- object$resSIM
  else stop("fitDistr: object must be either a numeric vector of an object of class 'propagate'!")
  
  MEAN <- mean(X, na.rm = TRUE)
  VAR <- var(X, na.rm = TRUE)
  SD <- sd(X, na.rm = TRUE)
  MIN <- min(X, na.rm = TRUE)
  MAX <- max(X, na.rm = TRUE)
  
  ## version 1.0-5: get histogram density
  DENS <- hist(X, freq = FALSE, breaks = nbin, plot = FALSE, ...)
  DENS$x <- DENS$mids
  DENS$y <- DENS$density    
  
  ## unweighted fitting or weighted fitting 
  wts <- rep(1, length(DENS$x)) 
  if (weights) {
    wts <- 1/DENS$counts
    wts[!is.finite(wts)] <- 1
  }
  if (is.numeric(weights)) {
    if (length(weights) != length(DENS$x)) stop("fitDistr: 'weights' must be a vector of length ", length(DENS$x), "!")
    wts <- weights
  }
  
  ## optimization function, minimum residual sum-of-squares is criterion
  optFun <- function(start, densfun, quantiles, density, eval = FALSE) {
    START <- as.list(start)
    START$x <- quantiles
    
    ## get density values from density function
    EVAL <- try(do.call(densfun, START), silent = TRUE) 
    if (inherits(EVAL, "try-error")) return(NA) 
    EVAL[is.nan(EVAL)] <- 0
    
    ## residual sum-of-squares to density values of object
    RSS <- sqrt(wts) * (density - EVAL)       
    if (eval) return(EVAL) else return(RSS)   
  }
  
  ## AIC/BIC function for optFun output
  ## version 1.0-6: added BIC 
  fitBIC <- function(fitobj) {
    ## taken and modified from stats:::logLik.nls
    RESID <- fitobj$fvec
    N <- length(RESID)
    W <- wts
    ZW <- W == 0    
    VAL <- -N * (log(2 * pi) + 1 - log(N) - sum(log(W + ZW)) + log(sum(W * RESID^2)))/2
    attr(VAL, "nobs") <- sum(!ZW)
    attr(VAL, "df") <- 1L + length(fitobj$par)
    class(VAL) <- "logLik"
    BIC(VAL)
  }
  
  ## define distribution names
  distNAMES <- c("Normal", 
                 "Skewed-normal", 
                 "Generalized normal", 
                 "Log-normal",
                 "Scaled/shifted t-",
                 "Logistic", 
                 "Uniform", 
                 "Triangular", 
                 "Trapezoidal",
                 "Curvilinear Trapezoidal",
                 "Gamma",  
                 "Inverse Gamma",
                 "Cauchy", 
                 "Laplace",
                 "Gumbel", 
                 "Johnson SU",
                 "Johnson SB",
                 "3P Weibull",
                 "2P Beta",
                 "4P Beta",
                 "Arcsine",
                 "von Mises",
                 "Inverse Gaussian",
                 "Generalized Extreme Value",
                 "Rayleigh",
                 "Chi-Square",
                 "Exponential",
                 "F-",
                 "Burr",
                 "Chi",
                 "Inverse Chi-Square",
                 "Cosine"
  )
  
  ## define distribution functions
  funLIST <- list(dnorm, 
                  dsn, 
                  dgnorm, 
                  dlnorm, 
                  dst, 
                  dlogis, 
                  dunif, 
                  dtriang, 
                  dtrap,
                  dctrap, 
                  dgamma,   
                  dinvgamma,
                  dcauchy, 
                  dlaplace,
                  dgumbel, 
                  dJSU, 
                  dJSB,
                  dweibull2, 
                  dbeta,
                  dbeta2,
                  darcsin,
                  dmises,
                  dinvgauss,
                  dgevd,
                  drayleigh,
                  dchisq,
                  dexp,
                  df,
                  dburr,
                  dchi,
                  dinvchisq,
                  dcosine
  )
  
  ## define start parameter list
  parLIST <- list(norm = c(mean = MEAN, sd = SD), 
                  sn = c(location = MEAN, scale = SD, shape = 1),
                  gnorm = c(alpha = 1, xi = MEAN, kappa = -0.1), 
                  lnorm = c(meanlog = 0, sdlog = 1),
                  st = c(mean = MEAN, sd = SD, df = 10),
                  logis = c(location = MEAN, scale = SD), 
                  unif = c(min = MIN, max = MAX), 
                  triang = c(a = MIN, b = (MIN + MAX)/2, c = MAX),
                  trap = c(a = 1.01 * MIN, b = MIN + 0.5 * (MEAN - MIN), 
                           c = MEAN + 0.5 * (MAX - MEAN), d = 0.99 * MAX),
                  ctrap = c(a = 1.01 * MIN, b = 0.99 * MAX, d = 0.01),
                  gamma = c(shape = MEAN^2/VAR, rate = MEAN/VAR),   
                  invgamma = c(shape = 1, scale = 10),
                  cauchy = c(location = MEAN, scale = SD), 
                  laplace = c(mean = MEAN, sigma = SD),                   
                  gumbel = c(location = MEAN - (sqrt(6) * SD/pi) * 0.5772, scale = sqrt(6) * SD/pi),
                  jsu = c(xi = MEAN, lambda = 1,  gamma = -MEAN, delta = 1),
                  jsb = c(xi = MIN, lambda = MAX - MIN,  gamma = 0, delta = 1),
                  weib = c(location = min(X, na.rm = TRUE), shape = 3, scale = 1),
                  beta = c(shape1 = 10, shape2 = 10),
                  beta2 = c(alpha1 = 10, alpha2 = 10, a = 0.9 * MIN, b = 1.1 * MAX),
                  arcsin = c(a = MIN, b = MAX),
                  mises = c(mu = MEAN, kappa = 1),
                  invgauss = c(mean = MEAN, dispersion = 0.1),
                  gevd = c(loc = MEAN, scale = SD, shape = 0),
                  rayleigh = c(mu = MEAN, sigma = SD),
                  chisq = c(df = 5),
                  exp = c(rate = 1),
                  fdist = c(df1 = 10, df2 = 10),
                  burr = c(k = 1),
                  chi = c(nu = 5),
                  invchisq = c(nu = 5),
                  cosine = c(mu = MEAN, sigma = SD)
  )
                  
  ## version 1-0-6: check for function selection
  if (!is.null(distsel) & !all(distsel %in% 1:length(distNAMES))) 
    stop(paste("Selection of distributions must be in 1 to", length(distNAMES), "!"))
  if (!is.null(distsel)) iter <- distsel else iter <- 1:length(distNAMES)
  
  ## preallocate list and vectors
  fitLIST <- vector("list", length = length(iter))
  BICS <- RSS <- MSE <- rep(NA, length(iter))
  
  ## fit all distributions 
  for (i in 1:length(iter)) {
    sel <- iter[i]
    if (verbose) message(sel, " of ", length(distNAMES), ": Fitting ", distNAMES[sel]," distribution...", sep = "")
    
    ## version 1.0-6: use gridded 'optFun' for all distributions
    ## create grid of starting parameters or single parameter
    if (brute == "slow") SEQ <- sapply(parLIST[[sel]], function(x) x * c(0.01, 0.1, 1, 10, 100))
    else SEQ <- sapply(parLIST[[sel]], function(x) x * c(0.1, 1, 10))
    GRID <- do.call(expand.grid, split(SEQ, 1:ncol(SEQ)))  
    colnames(GRID) <- names(parLIST[[sel]])
    
    ## preallocate empty vector for RSS
    rssVEC <- rep(NA, nrow(GRID))
    
    ## collect RSS for all grid values by calling 'optFun'
    for (j in 1:nrow(GRID)) {
      if (verbose) counter(j)
      PARS <- GRID[j, ]
      FIT <- try(nls.lm(par = PARS, fn = optFun, densfun = funLIST[[sel]], quantiles = DENS$x,
                        density = DENS$y, control = nls.lm.control(maxiter = 10000, maxfev = 10000)), silent = TRUE)
      if (inherits(FIT, "try-error")) rssVEC[j] <- NA else rssVEC[j] <- FIT$deviance     
    }
    
    ## select parameter combination with lowest RSS and re-fit
    WHICH <- which.min(rssVEC) 
    bestPAR <- GRID[WHICH, ]
    FIT <- try(nls.lm(par = bestPAR, fn = optFun, densfun = funLIST[[sel]], quantiles = DENS$x,
                      density = DENS$y, control = nls.lm.control(maxiter = 10000, maxfev = 10000)), silent = TRUE)      
    ## calculate GOF measures
    if (inherits(FIT, "try-error")) {
      FIT <- NA
      if (verbose) message("fitDistr: error in calculating BIC.\n")
    } else {
      fitLIST[[i]] <- FIT
      BICS[i] <- tryCatch(fitBIC(FIT), error = function(e) NA)
      RSS[i] <- FIT$fvec^2
      MSE[i] <- mean(FIT$fvec^2, na.rm = TRUE)
    }    
    cat("\n\n")
  } 
  
  ## aggregate and sort ascending by AIC
  ORDER <- order(BICS)
  statDAT <- data.frame("Distribution" = distNAMES[iter], "BIC" = BICS, "RSS" = RSS, "MSE" = MSE)
  statDAT <- statDAT[ORDER, ]
  
  ## select best fit
  SEL <- ORDER[1]
  bestFIT <- fitLIST[[SEL]]
  evalLIST <- as.list(bestFIT$par)
  evalLIST$x <- DENS$x
  evalY <- do.call(funLIST[[iter[SEL]]], evalLIST)  
  
  ## version 1.0-6: sort fitLIST by AIC
  fitLIST <- fitLIST[ORDER]
  names(fitLIST) <- distNAMES[ORDER]
  
  ## version 1.0-6: create parLIST and s.e. list
  parLIST <- lapply(fitLIST, function(x) coef(x))
  seLIST <- lapply(fitLIST, function(x) tryCatch(sqrt(diag(vcov(x))), error = function(e) NA))
                   
  ## version 1.0-6: create sorted funLIST
  funLIST <- funLIST[ORDER]
  names(funLIST) <- distNAMES[iter[ORDER]]
  
  ## plot best fit as histogram
  if (plot == "hist") {     
    TEXT <- paste(distNAMES[iter[SEL]], "distribution,\n BIC =", round(BICS[SEL], 3))
    hist(X, freq = FALSE, breaks = nbin, cex.axis = 1.5, las = 0, col = "dodgerblue3",
         cex.lab = 1.5, xlab = "Bin", main = TEXT, ...)
    lines(DENS$x, evalY, col = "red3", lty = 1, lwd = 3, ...)    
  }  
  
  ## version 1.0-6: plot best fit as qqplot
  if (plot == "qq") {
    TEXT <- paste(distNAMES[iter[SEL]], "distribution\n BIC =", round(BICS[SEL], 3))
    qqplot(DENS$y, evalY, asp = 1, pch = 16, col = "dodgerblue3", xlab = "Density", ylab = "Fitted", main = TEXT, ...)
    grid()
    segments(-0.1, -0.1, 1.2 * max(evalY, na.rm = TRUE), 1.2 * max(evalY, na.rm = TRUE), col = "red3", lty = 1, lwd = 2, ...)
  }
  
  OUT <- list(stat = statDAT, fit = fitLIST, par = parLIST, se = seLIST, dens = funLIST, 
              bestfit = bestFIT, bestpar = parLIST[[1]], bestse = seLIST[[1]], 
              fitted = evalY, residuals = DENS$y - evalY)
  
  class(OUT) <- "fitDistr"
  return(OUT)
}

