stochContr <- function(prop, plot = TRUE)
{
  ## check for correct input
  if (class(prop) != "propagate") stop("This is no 'propagate' object!")
  
  ## get data from 'prop'
  DATA <- prop$data
  EXPR <- prop$expr
  datSIM<- prop$datSIM
  resSIM <- prop$resSIM
  
  ## get means either from original data or from simulated input
  MEAN <- if (nrow(DATA) > 3) colMeans(datSIM, na.rm = TRUE) else DATA[1, ]
 
  ## preallocate result list
  resList <- vector("list", length(MEAN))
  
  ## iterate over original result data
  for (i in 1:length(MEAN)) {
    tempDATA <- datSIM
    
    ## substitute simulated data by mean without variance
    tempDATA[, i] <- rep(MEAN[i], nrow(tempDATA))
    
    ## try vectorized evaluation, which is much faster  
    tempRES <- try(eval(EXPR, as.data.frame(tempDATA)), silent = TRUE) 
    
    ## use 'row-wise' method if 'vectorized' throws an error
    if (inherits(tempRES, "try-error")) {
      print("Using 'vectorized' evaluation gave an error. Switching to 'row-wise' evaluation...")
      tempRES <- apply(tempDATA, 1, function(x) eval(EXPR, envir = as.list(x)))     
    }
    
    ## store in result list
    resList[[i]] <- tempRES
  }
  
  ## plot distributions
  allDAT <- cbind(all = resSIM, as.data.frame(resList))
  colnames(allDAT) <- c("all", colnames(datSIM))
  COL <- c("cornflowerblue", rep("firebrick", length(resList)))
  if (plot) boxplot(allDAT, outline = FALSE, col = COL, las = 1, pars = list(bty = "c"))
  
  ## calculate contribution as differences in variance
  VARS <- colVarsC(allDAT)
  DIFF <- VARS[1] - VARS[-1]
  CONTRIB <- DIFF/sum(DIFF, na.rm = TRUE)
  
  return(CONTRIB)
}








