colVarsC <- function(x)
{
  if (!is.matrix(x)) x <- as.matrix(x)
  CN <- colnames(x)
  if (NCOL(x) == 1) stop("colVarsC: 'x' should have more than one column!")
  OUT <- .Call("C_colVarsC", x, PACKAGE = "propagate")
  names(OUT) <- CN
  OUT
}

rowVarsC <- function(x)
{
  if (!is.matrix(x)) x <- as.matrix(x)
  RN <- rownames(x)
  if (NROW(t(x)) == 1) stop("rowVarsC: 'x' should have more than one row!")
  OUT <- .Call("C_rowVarsC", x, PACKAGE = "propagate")
  names(OUT) <- RN
  OUT
}
