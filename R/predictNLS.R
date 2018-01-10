predictNLS <- function(
model, 
newdata,
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
  
  ## check that 'newdata' has same name as predVAR
  if (!identical(names(newdata)[1:length(predVAR)], predVAR)) stop("predictNLS: newdata should have name(s) ", predVAR, "!\n")
    
  ## get variance-covariance matrix
  VCOV <- vcov(model)
  
  ## iterate over all entries in 'newdata' as in usual 'predict.' functions
  NR <- NROW(newdata)
  outMAT <- matrix(nrow = NR, ncol = 12)
  propLIST <- vector("list", length = NR)
  
  for (i in 1:NR) {
    #cat("Propagating predictor value #", i, " ...\n", sep = "")
    message("predictNLS: ", paste0("Propagating predictor value #", i, "..."))
    tempDATA <- newdata[i, , drop = FALSE]
    colnames(tempDATA) <- colnames(newdata) 
    
    ## create dataframe of variables for 'propagate'
    SEL <- which(names(tempDATA) == predVAR) 
    predVEC <- c(COEF, as.numeric(tempDATA[SEL]))  
    names(predVEC) <- c(names(COEF), predVAR)
    DF <- rbind(predVEC, 0)  
    row.names(DF) <- NULL 
    
    ## create covariance matrix for 'propagate'       
    if (NCOL(tempDATA) == length(predVAR)) forCOV <- matrix(rep(0, length(predVAR)), nrow = 1) 
    else forCOV <- tempDATA[, -(1:length(predVAR)), drop = FALSE] 
    colnames(forCOV) <- predVAR
   
    ## create covariance matrix
    COV <- mixCov(VCOV, forCOV^2)
    colnames(COV) <- rownames(COV) <- colnames(DF)
    
    ## call 'propagate'   
    PROP <- propagate(expr = EXPR, data = DF, cov = COV, alpha = alpha, ...)
    propLIST[[i]] <- PROP
    
    ## populate outMAT and override confidence/prediction values from 'propagate'
    outPROP <- PROP$prop
      
    if (is.na(outPROP[4])) {
      SD <- outPROP[3] 
      MEAN <- outPROP[1]
    } else {        
      SD <-  outPROP[4] 
      MEAN <- outPROP[2]
    }   
        
    if (interval != "none") {
      TQUAN <- qt(1 - alpha/2, df.residual(model))
      
      if (interval == "confidence") {
        outPROP[5] <- MEAN - TQUAN * SD 
        outPROP[6] <- MEAN + TQUAN * SD         
      }      
      if (interval == "prediction") {
        rss <- sum(residuals(model)^2)
        n <- length(residuals(model))
        p <- length(coef(model))
        rv <- rss/(n - p)
        outPROP[5] <- MEAN - TQUAN * sqrt(SD^2 + rv) 
        outPROP[6] <- MEAN + TQUAN * sqrt(SD^2 + rv)          
      }        
    } else {
      outPROP[5] <- NA
      outPROP[6] <- NA
    }     
            
    outSIM <- PROP$sim
    OUT <- c(outPROP, outSIM)
    outMAT[i, ] <- OUT
  }
  
  outMAT <- as.data.frame(outMAT)
  colnames(outMAT) <- c(paste("Prop.", names(outPROP), sep = ""), 
                        paste("Sim.", names(outSIM), sep = ""))
    
  return(list(summary = outMAT, prop = propLIST))
}