########### Generating random deviates from distributions #############
## Random samples from the skew-normal distribution, taken from package 'VGAM'
rsn <- function(n, location = 0, scale = 1, shape = 0) 
{
  rho <- shape/sqrt(1 + shape^2)
  u0 <- rnorm(n)
  v <- rnorm(n)
  u1 <- rho * u0 + sqrt(1 - rho^2) * v
  res <- location + scale * ifelse(u0 >= 0, u1, -u1)
  res[scale <= 0] <- NA
  res
}

## Random samples from the generalized normal distribution, taken from 
## inverseCDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
rgnorm <- function(n, alpha = 1, xi = 1, kappa = -0.1) 
{
  ## 'invcerf' taken from erf.R in package 'pracma'
  invcerf <- function(x) -qnorm(x/2)/sqrt(2) 
  
  res <- xi + (alpha/kappa) * (1 - exp((sqrt(2) * kappa) * invcerf(2 * runif(n))))
  res
}

## Random samples from the scales and shifted t-distribution, 
## taken from the 'metRology' package
rst <- function(n, mean = 0, sd = 1, df = 2) 
{
  mean + sd * rt(n, df = df)
}

## Random samples from the Gumbel distribution, taken from 
## inverseCDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
rgumbel <- function(n, location = 0, scale = 1) 
{
  res <- location - scale * log(-log(runif(n)))
  res[scale <= 0] <- NaN
  res
}

## Random samples from the Johnson SU distribution, taken from 
## inverseCDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
rJSU <- function(n, xi = 0, lambda = 1, gamma = 1, delta = 1) 
{
  erf.inv <- function(x) qnorm((x + 1)/2)/sqrt(2)
  res <- xi + lambda * sinh((-gamma + sqrt(2) * erf.inv(-1 + 2 * runif(n)))/delta)
  res
}

## Random samples from the Johnson SB distribution, taken from 
## inverseCDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
rJSB <- function(n, xi = 0, lambda = 1, gamma = 1, delta = 1) 
{
  erf.inv <- function(x) qnorm((x + 1)/2)/sqrt(2)
  res <- xi + lambda/(1 + exp((gamma - sqrt(2) * erf.inv(-1 + 2 * runif(n)))/delta))
  res
}

## Random samples from the three-parameter weibull distribution, taken from 
## inverseCDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
rweibull2 <- function(n, location = 0, shape = 1, scale = 1) 
{
  res <- location + scale * (-log(1 - runif(n))^(1/shape))
  res
}

## Random samples from the four-parameter beta distribution, taken from 
## inverseCDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
rbeta2 <- function(n, alpha1 = 1, alpha2 = 1, a = 0, b = 0)  
{
  res <- a + qbeta(runif(n), shape1 = alpha1, shape2 = alpha2, lower.tail = TRUE, log.p = FALSE)
  res <- rescale(res, a, b)
  res
}

## Random samples from the triangular distribution, taken from 
## inverseCDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
rtriang <- function(n, a = 0, b = 0.5, c = 1)  
{
  res <- numeric(length(n))
  x <- runif(n)
  CRIT <- (b - a)/(c - a)
  res[x >= 0 & x <= CRIT] <- a + sqrt(((b - a) * (c - a)) * x[x >= 0 & x <= CRIT])
  res[CRIT < x & x <= 1] <- c - sqrt(((c - b) * (c - a)) * (1 - x[CRIT < x & x <= 1]))
  res
}

## Random samples from the trapezoidal distribution, taken from 
## inverseCDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
rtrap <- function(n, a = 0, b = 1, c = 2, d = 3)  
{
  res <- numeric(length(n))
  x <- runif(n)
  CRIT1 <- (a - b)/(a + b - c - d)
  CRIT2 <- (a + b - 2 * c)/(a + b - c - d)
  res[x > 0 & x < CRIT1] <- a + sqrt(((b - a) * (-a - b + c + d)) * x[x > 0 & x < CRIT1])
  res[CRIT1 <= x & x < CRIT2] <- 0.5 * (a + b + (- a - b + c + d) * x[CRIT1 <= x & x < CRIT2])
  res[CRIT2 <= x & x < 1] <- d - sqrt(((d - c) * (-a - b + c + d)) * (1 - x[CRIT2 <= x & x < 1]))
  res
}

## Random samples from the Laplacian distribution, taken from 
## inverseCDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
rlaplace <- function(n, mean = 0, sigma = 1)   
{
  RUNIF <- runif(n)
  res <- mean - (sigma * log(1 - (-1 + 2 * RUNIF) * sign(-1 + 2 * RUNIF))) * sign(-1 + 2 * RUNIF) 
  res
}

## Random samples from the Arcsine distribution, taken from 
## inverseCDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
rarcsin <- function(n, a = 2, b = 1)   
{
  res <- a + (b - a) * (sin(pi * runif(n)/2)^2)
  res 
}

## Random samples from the von Mises distribution,  
## using "Rejection Sampling" and runif as the enveloping distribution
rmises <- function(n, mu = 1, kappa = 3)  
{
  OUT <- numeric(n)
  LOWER <- mu - 2 * pi
  UPPER <- mu + 2 * pi  
  POS1 <- 1
  nSAMPLE <- 0
  N <- n
  maxDUNIF <- 1/(UPPER - LOWER) 
  
  while(nSAMPLE < n) {
    ## (1) sample from rectangular distribution 
    u <- runif(N)   
    
    ## (2) Sample from enveloping rectangular distribution from mu - pi/mu + pi
    x <- runif(N, LOWER, UPPER)  
    
    ## (3) calculate scaling factor
    DMISES <- dmises(x, mu = mu, kappa = kappa)    
    if (nSAMPLE == 0) {      
      maxDMISES <- max(DMISES, na.rm = TRUE)
      A <- 1.1 * maxDMISES/maxDUNIF       
    }       
        
    ## (4) Accept candidate value, if probability u is smaller or equal than the density g(x)
    ## at x divided by the density of f(x) * A (scale) 
    VAL <- 1/A * DMISES/maxDUNIF    
    SEL <- x[u <= VAL]
    
    ## (5) fill result vector and increase start position,
    ## calculate ratio of true events
    LEN <- length(SEL)
    if (nSAMPLE == 0) RATIO <- N/LEN/10
    nSAMPLE <- nSAMPLE + LEN
    POS2 <- POS1 + LEN - 1
    OUT[POS1:POS2] <- SEL
    POS1 <- POS2 + 1    
              
    ## if n(events) < n, restart with ratio as factor
    if (nSAMPLE < n) N <- RATIO * n    
  }  
  
  return(OUT[1:n])
}

## Random samples from the curvilinear trapezoidal distribution,  
## using "Rejection Sampling" and runif as the enveloping distribution
rctrap <- function(n, a = 0, b = 1, d = 0.1)  
{
  OUT <- numeric(n)
  LOWER <- a - 2 * d
  UPPER <- b + 2 * d  
  POS1 <- 1
  nSAMPLE <- 0
  N <- n
  maxDUNIF <- 1/(UPPER - LOWER) 
  
  while(nSAMPLE < n) {
    ## (1) sample from rectangular distribution 
    u <- runif(N)   
    
    ## (2) Sample from enveloping rectangular distribution from mu - pi/mu + pi
    x <- runif(N, LOWER, UPPER)  
    
    ## (3) calculate scaling factor
    DCTRAP <- dctrap(x, a = a, b = b, d = d)    
    if (nSAMPLE == 0) {      
      maxDCTRAP <- max(DCTRAP, na.rm = TRUE)
      A <- 1.1 * maxDCTRAP/maxDUNIF       
    }       
    
    ## (4) Accept candidate value, if probability u is smaller or equal than the density g(x)
    ## at x divided by the density of f(x) * A (scale) 
    VAL <- 1/A * DCTRAP/maxDUNIF    
    SEL <- x[u <= VAL]
    
    ## (5) fill result vector and increase start position,
    ## calculate ratio of true events
    LEN <- length(SEL)
    if (nSAMPLE == 0) RATIO <- N/LEN/10
    nSAMPLE <- nSAMPLE + LEN
    POS2 <- POS1 + LEN - 1
    OUT[POS1:POS2] <- SEL
    POS1 <- POS2 + 1    
    
    ## if n(events) < n, restart with ratio as factor
    if (nSAMPLE < n) N <- RATIO * n    
  }  
  
  return(OUT[1:n])
}

## Random samples from the generalized trapezoidal distribution,  
## using "Rejection Sampling" and runif as the enveloping distribution
rgtrap <- function(n, min = 0, mode1 = 1/3, mode2 = 2/3, max = 1, n1 = 2, 
                   n3 = 2, alpha = 1) 
{
  OUT <- numeric(n)
  LOWER <- min
  UPPER <- max  
  POS1 <- 1
  nSAMPLE <- 0
  N <- n
  maxDUNIF <- 1/(UPPER - LOWER) 
    
  while(nSAMPLE < n) {
    ## (1) sample from rectangular distribution 
    u <- runif(N)   
      
    ## (2) Sample from enveloping rectangular distribution from mu - pi/mu + pi
    x <- runif(N, LOWER, UPPER)  
      
    ## (3) calculate scaling factor
    DGTRAP <- dgtrap(x, min = min, mode1 = mode1, mode2 = mode2, max = max, n1 = n1, 
                         n3 = n3, alpha = 1, log = FALSE)    
    if (nSAMPLE == 0) {      
      maxDGTRAP <- max(DGTRAP, na.rm = TRUE)
      A <- 1.1 * maxDGTRAP/maxDUNIF       
    }       
      
    ## (4) Accept candidate value, if probability u is smaller or equal than the density g(x)
    ## at x divided by the density of f(x) * A (scale) 
    VAL <- 1/A * DGTRAP/maxDUNIF    
    SEL <- x[u <= VAL]
      
    ## (5) fill result vector and increase start position,
    ## calculate ratio of true events
    LEN <- length(SEL)
    if (nSAMPLE == 0) RATIO <- N/LEN/10
    nSAMPLE <- nSAMPLE + LEN
    POS2 <- POS1 + LEN - 1
    OUT[POS1:POS2] <- SEL
    POS1 <- POS2 + 1    
      
    ## if n(events) < n, restart with ratio as factor
    if (nSAMPLE < n) N <- RATIO * n    
  }  
    
  return(OUT[1:n])
}

## Random samples from the Inverse Gaussian distribution, taken from 
## package 'statmod'
rinvgauss <- function(n, mean = 1, dispersion = 1) 
{
  mu <- rep_len(mean, n)
  phi <- rep_len(dispersion, n)
  r <- rep_len(0, n)
  i <- (mu > 0 & phi > 0)
  if (!all(i)) {
    r[!i] <- NA
    n <- sum(i)
  }
  phi[i] <- phi[i] * mu[i]
  Y <- rchisq(n, df = 1)
  X1 <- 1 + phi[i]/2 * (Y - sqrt(4 * Y/phi[i] + Y^2))
  firstroot <- as.logical(rbinom(n, size = 1L, prob = 1/(1 + X1)))
  r[i][firstroot] <- X1[firstroot]
  r[i][!firstroot] <- 1/X1[!firstroot]
  mu * r
}

## Random samples from the Generalized Extreme Value distribution, taken from 
## package 'extRemes'
rgevd <- function(n, loc = 0, scale = 1, shape = 0) 
{
  z <- rexp(n)
  out <- numeric(n) + NA
  loc <- rep(loc, n)
  sc <- rep(scale, n)
  sh <- rep(shape, n)
  id <- sh == 0
  if (any(id)) out[id] <- loc[id] - sc[id] * log(z[id])
  if (any(!id)) out[!id] <- loc[!id] + sc[!id] * (z[!id]^(-sh[!id]) - 1)/sh[!id]
  return(out)
}

## Random samples from the Inverse Gamma distribution, taken from 
## package 'MCMCpack' 
rinvgamma <- function(n, shape = 1, scale = 5) 
{
  return(1/rgamma(n = n, shape = shape, rate = scale))
}

## Random samples from the Rayleigh distribution, taken from 
## inverseCDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
rrayleigh <- function(n, mu = 1, sigma = 1) 
{
  RUNIF <- runif(n)
  mu + sigma * sqrt(log(1/(1 - RUNIF)^2))
}

## Random samples from the Burr VIII distribution, taken from 
## inverseCDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
rburr <- function(n, k = 1) 
{
  RUNIF <- runif(n)
  log(tan(0.5 * pi * RUNIF^(1/k)))
}

## Random samples from the Chi distribution, by ANS
rchi <- function(n, nu = 5) sqrt(rchisq(n, nu))

## Random samples from the Inverse Chi-Square distribution, taken from 
## inverseCDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
rinvchisq <- function(n, nu = 5) 1/rchisq(n, nu)

## Random samples from the Cosine distribution,  
## using "Rejection Sampling" and runif as the enveloping distribution
rcosine <- function(n, mu = 5, sigma = 1) 
{
  OUT <- numeric(n)
  LOWER <- mu - sigma 
  UPPER <- mu + sigma    
  POS1 <- 1
  nSAMPLE <- 0
  N <- n
  maxDUNIF <- 1/(UPPER - LOWER) 
  
  while(nSAMPLE < n) {
    ## (1) sample from rectangular distribution 
    u <- runif(N)   
    
    ## (2) Sample from enveloping rectangular distribution from mu - pi/mu + pi
    x <- runif(N, LOWER, UPPER)  
    
    ## (3) calculate scaling factor
    DCOSINE <- (1 + cos((pi * (x - mu))/sigma))/(2 * sigma)    
    
    if (nSAMPLE == 0) {      
      maxDCOSINE <- max(DCOSINE, na.rm = TRUE)
      A <- 1.1 * maxDCOSINE/maxDUNIF       
    }       
    
    ## (4) Accept candidate value, if probability u is smaller or equal than the density g(x)
    ## at x divided by the density of f(x) * A (scale) 
    VAL <- 1/A * DCOSINE/maxDUNIF    
    SEL <- x[u <= VAL]
    
    ## (5) fill result vector and increase start position,
    ## calculate ratio of true events
    LEN <- length(SEL)
    if (nSAMPLE == 0) RATIO <- N/LEN/10
    nSAMPLE <- nSAMPLE + LEN
    POS2 <- POS1 + LEN - 1
    OUT[POS1:POS2] <- SEL
    POS1 <- POS2 + 1    
    
    ## if n(events) < n, restart with ratio as factor
    if (nSAMPLE < n) N <- RATIO * n    
  }  
  
  return(OUT[1:n])
}