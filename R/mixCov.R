mixCov <- function(...)
{
  ## convert arguments to named list
  covLIST <- list(...)
  
  ## extract function call names
  allVars <- all.vars(sys.call(0))
  
  ## get dimensions
  dimLIST <- lapply(covLIST, function(x) NCOL(x))
  
  ## get columns names
  nameVec <- NULL
  
  ## create final covariance matrix
  DIMS <- do.call("sum", dimLIST)
  covMAT <- matrix(NA, nrow = sum(DIMS), ncol = sum(DIMS))
  
  ## populate matrix and augment name vector
  counter <- 1
  for (i in 1:length(covLIST)) {
    if (!is.null(colnames(covLIST[[i]]))) nameVec <- c(nameVec, colnames(covLIST[[i]]))
    else nameVec <- c(nameVec, allVars[i])
    POS <- counter:(counter + dimLIST[[i]] - 1)
    INS <- covLIST[[i]]
    if (NCOL(INS) > NROW(INS)) INS <- diag(as.numeric(INS)) ## check if vector or covariance matrix
    if (NCOL(INS) == 1) INS <- as.numeric(INS)
    covMAT[POS, POS] <- INS
    counter <- counter + dimLIST[[i]]
  }
  
  ## fill with 0's
  covMAT[is.na(covMAT)] <- 0
  
  ## create colnames/rownames
  rownames(covMAT) <- colnames(covMAT) <- nameVec
  
  return(covMAT)  
}


