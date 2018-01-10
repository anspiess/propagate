WelchSatter <- function(ui, ci = NULL, df = NULL, dftot = NULL, uc = NULL, alpha = 0.05)
{
  if (is.null(ci)) ci <- rep(1, length(ui))
  if (is.null(df)) stop("WelchSatter: Please supply 'df's!")
  if (length(ui) != length(df)) stop("WelchSatter: Different number of values in 'ui' and 'df'!")
  if (!is.null(dftot) & length(dftot) != 1) stop("WelchSatter: Total degrees of freedom must be a single number!")
  if (is.null(uc)) uc <- sqrt(sum((ci * ui)^2))
  ws.df <- (uc^4)/sum(((ci * ui)^4)/df)
  if (is.numeric(dftot)) ws.df <- dftot
  k <- qt(1 - alpha/2, floor(ws.df))
  u.exp <- k * uc
  return(list(ws.df = ws.df, k = k, u.exp = u.exp))   
}