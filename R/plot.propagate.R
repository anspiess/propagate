plot.propagate <- function(x, logx = FALSE, ...)
{
  object <- x
  
  ## plot setup
  par(mar = c(4, 4, 3, 1))   
  
  ## cut off histogram extreme values
  ## and plot histogram
  resSIM <- object$resSIM
  FILTER <- quantile(resSIM, c(0.001, 0.999), na.rm = TRUE) 
  plotDATA <- resSIM[resSIM > FILTER[1] & resSIM < FILTER[2]]
  plotDATA <- plotDATA[!is.na(plotDATA)]
  
  if (logx) plotDATA <- suppressWarnings(log(plotDATA)) 
    
  ## histogram with fitted density curve and confidence intervals
  HIST <- hist(plotDATA, col = "dodgerblue3", breaks = 100, 
               main = NULL, cex.main = 1, freq = FALSE,  xlab = "Bin", ylab = "Density", ...)
  DENS <- density(plotDATA)
  lines(DENS, col = "red3", lwd = 3)
  title(main = paste("Histogram of Monte Carlo simulation results with density curve (red),\n", 
                     (1 - object$alpha) * 100, "% confidence interval (black), median (orange) and mean (green)", sep = ""), cex.main = 1)
  abline(v = object$sim[c(5, 6)], col = "black", lwd = 3)
  abline(v = object$sim[3], col = "darkorange", lwd = 3, lty = 1)
  abline(v = object$sim[1], col = "green3", lwd = 3, lty = 1)
}
  
  