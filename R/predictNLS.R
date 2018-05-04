predictNLS <- function(
  model, 
  newdata,
  newerror,
  interval = c("confidence", "prediction", "none"),
  alpha = 0.05, 
  ...
)
{
  interval <- match.arg(interval)
  
  ## get right-hand side of formula
  RHS <- as.list(eval(model$call$formula))[[3]]
  EXPR <- as.expression(RHS)
  
  ## all variables in model
  VARS <- all.vars(EXPR)
  
  ## coefficients
  COEF <- coef(model)
  
  ## extract predictor variable
  predVAR <- setdiff(VARS, names(COEF))  
  
  ## version 1.0-6: take original predictor values, if 'newdata' is missing,
  ## as in classical 'predict.nls'
  if (missing(newdata)) {
    newdata <- try(eval(model$data)[, predVAR, drop = FALSE], silent = TRUE)
    if (inherits(newdata, "try-error")) {
      newdata <- sapply(predVAR, function(i) get(i))
      names(newdata) <- predVAR
    }
  }
  
  ## check that 'newdata' and 'newerror' has same name as predVAR
  if (length(setdiff(colnames(newdata), predVAR)) != 0) stop("predictNLS: 'newdata' should have column name(s): ", predVAR, "!\n")
  if (!missing(newerror)) {
    if (length(setdiff(colnames(newerror), predVAR)) != 0) stop("predictNLS: 'newerror' should have column name(s): ", predVAR, "!\n")
  }
  
  ## get variance-covariance matrix
  VCOV <- vcov(model)
  
  ## create and check 'newerror' for same dimension as 'newdata'
  if (missing(newerror)) {
    newerror <- newdata
    newerror[] <- 0
  } else {
    if (!identical(dim(newdata), dim (newerror))) stop("'newdata' and 'newerror' should have the same dimensions!")
    if (!identical(colnames(newdata), colnames(newerror))) stop("'newdata' and 'newerror' should have the same column names!")
  }
  
  ## iterate over all entries in 'newdata' and 'newerror' as in usual 'predict.' functions
  NR <- NROW(newdata)
  outMAT <- matrix(nrow = NR, ncol = 12)
  propLIST <- vector("list", length = NR)
  
  ## version 1.06: if interval = "prediction" add "eps" variable for residual standard error
  if (interval == "prediction") {
    form <- formula(model)
    form <- as.formula(paste(form[2], form[1], paste(form[3], " + rv", sep = ""), sep = " "))
    EXPR <- as.expression(form[[3]])
  }
  
  for (i in 1:NR) {
    message("predictNLS: ", paste0("Propagating predictor value #", i, "..."))
    tempDATA <- newdata[i, , drop = FALSE]
    errorDATA <- newerror[i, , drop = FALSE]
    
    ## create dataframe of variables for 'propagate'
    ## version 1.0-6: if interval = "prediction", add zero mean for residual variance
    DAT <- as.numeric(c(COEF, tempDATA)) 
    names(DAT) <- c(names(COEF), predVAR)
    DAT <- rbind(DAT, 0)  
    row.names(DAT) <- NULL 
    if (interval == "prediction") DAT <- cbind(DAT, rv = c(0, 0))
    
    ## create (augmented) covariance matrix
    if (interval == "confidence") COV <- mixCov(VCOV, errorDATA^2)
    else {
      r <- residuals(model)
      w <- weights(model)
      rss <- sum(if (is.null(w)) r^2 else r^2 * w)
      n <- length(residuals(model))
      p <- length(coef(model))
      rv <- rss/(n - p)
      COV <- mixCov(VCOV, errorDATA^2, rv)
    }
    
    ## version 1.06: include degrees of freedom
    DF <- df.residual(model)
    
    ## call 'propagate'   
    PROP <- propagate(expr = EXPR, data = DAT, cov = COV, alpha = alpha, df = DF, ...)
    propLIST[[i]] <- PROP
    
    ## populate outMAT
    outPROP <- PROP$prop
    outSIM <- PROP$sim
    OUT <- c(outPROP, outSIM)
    outMAT[i, ] <- OUT
  }
  
  outMAT <- as.data.frame(outMAT)
  colnames(outMAT) <- c(paste("Prop.", names(outPROP), sep = ""), 
                        paste("Sim.", names(outSIM), sep = ""))
  return(list(summary = outMAT, prop = propLIST))
}