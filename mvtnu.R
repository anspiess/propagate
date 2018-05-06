mvtvnu <- function (
  n = 1000000, 
  mu = c(0, 0), 
  sigma = diag(length(mu)),
  nu = rep(10, length(mu)),
  cut = TRUE
)
{
  ## scaled/shifted t-distribution
  qt.scaled <- function (p, df, mean = 0, sd = 1, ncp, lower.tail = F, log.p = FALSE) 
  {
    mean + sd * stats::qt(p, df, ncp = ncp, log.p = log.p)
  }
  
  ## input checks
  if (!all.equal(length(mu), length(nu), length(diag(sigma))))
   stop("'mu', 'nu' and 'sigma' need to be of same dimension!")
  len <- length(mu)
  n0 <- n
  if (cut) n <- n + len * 10
  
  ## Step 1: map Pearson to Spearman correlations
  sigma2 <- sigma
  sigma2 <- 2 * sin(pi * sigma2 / 6)
  diag(sigma2) <- 1
  
  ## STEP 2: sample from multivariate centered/scaled normal distribution
  X <- MASS:::mvrnorm(n, mu = rep(0, len), Sigma = sigma2, empirical = TRUE)
  print(cor(X), method = "spearman")
  
  ## STEP 3: convert to uniform distribution
  P <- pnorm(X, lower.tail = TRUE)
  print(cor(P), method = "spearman")
  
  ## transform the marginals to t-distributions with individual df's
  Z <- sapply(1:len, function(i) qt.scaled(P[, i], df = nu[i], mean = mu[i], 
                                 sd = sqrt(diag(sigma)[i]), lower.tail = FALSE, log.p = FALSE))
  print(cor(Z), method = "spearman")
  #plot(Z[, 1], Z[, 2], pch = 16, col = "#33555533")
  
  ## optionally remove extreme points
  if (cut) {
    CUT <- apply(Z, 2, function(x) quantile(x, c(5/n, 1-5/n), na.rm = TRUE))
    for (i in 1:ncol(Z)) {
      SEL <- which(Z[, i] < CUT[1, i] | Z[, i] > CUT[2, i])
      Z[SEL, i] <- NA
    }
    Z <- na.omit(Z)
    if (nrow(Z) > n0) Z <- Z[1:n0, ]
  }
  
  return(Z)
}

####################################################
SIGMA <- matrix(c(1, 0.2, 0.5, 0.2, 2, 0, 0.5, 0, 4), nrow = 3)
MU <- c(10, 2, 20)
NU <- c(10, 5, 3)
set.seed(123)
res <- mvtvnu(n = 1000000, mu = MU, sigma = SIGMA, nu = NU, cut = T)
cor(res)
#plot(res[, 1:2], pch = 16, cex = 0.5, asp = 1)
#for (i in 1:ncol(res)) print(MASS:::fitdistr(res[, i], "t"))
#library(rgl)
#plot3d(res[, 1:3])