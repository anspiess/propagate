## data for Harris model
DAT <- data.frame(x = c(3.4, 7.1, 16.1, 20, 23.1, 34.4, 40, 44.7, 65.9, 78.9, 96.8, 115.4, 120),
                  y = c(9.59, 5.29, 3.63, 3.42, 3.46, 3.06, 3.25, 3.31, 3.50, 3.86, 4.24, 4.62, 4.67))

## unweighted fit
LM <- nlsLM(y ~ A * x + B/x + C, data = DAT, start = list(A = 1, B = 1, C = 1))
summary(LM)

## weighted fit (5% of y)
wts <- 1/c(0.48, 0.26, 0.18, 0.17, 0.17, 0.15, 0.16, 0.17, 0.18, 0.19, 0.21, 0.23, 0.23)^2
LM <- nlsLM(y ~ A * x + B/x + C, data = DAT, start = list(A = 1, B = 1, C = 1), weights = wts)
summary(LM)

## EIV2
wtsx <- 1/(0.03 * DAT[, 1])^2
wtsy <- 1/(0.02 * DAT[, 2])^2
LM <- nlsLM(y ~ A * x + B/x + C, data = DAT, start = list(A = 1, B = 1, C = 1), weights = wtsy, weights.x = wtsx)
summary(LM)

#######################################################################################

## Original Harris paper
xx <- c(3.4, 7.1, 16.1, 20, 23.1, 34.4, 40, 44.7, 65.9, 78.9, 96.8, 115.4, 120)
yy <- c(9.59, 5.29, 3.63, 3.42, 3.46, 3.06, 3.25, 3.31, 3.50, 3.86, 4.24, 4.62, 4.67)
sigy <- sqrt(1/c(4.34, 14.79, 30.86, 34.6, 34.6, 44.44, 39.06, 34.60, 30.86, 27.7, 22.68, 18.9, 18.9))
fct <- function(par) {
  calc <- par[1] * xx + par[2]/xx + par[3]
  sigt <- sqrt(sigy^2)
  fitf <- (calc - yy)/sigt
  fitf
} 
NLS <- nls.lm(par = c(1, 1, 1), fn = fct, control = nls.lm.control(maxiter = 1024))
summary(NLS)

## Harris Data with Joel's EIV2
xx <- c(3.4, 7.1, 16.1, 20, 23.1, 34.4, 40, 44.7, 65.9, 78.9, 96.8, 115.4, 120)
yy <- c(9.59, 5.29, 3.63, 3.42, 3.46, 3.06, 3.25, 3.31, 3.50, 3.86, 4.24, 4.62, 4.67)
sigx <- 0.03 * xx
fct <- function(par) {
  calc <- par[1] * xx + par[2]/xx + par[3]
  sigy <- 0.02 * calc
  sigt <- sqrt(sigy^2 + (par[1] - par[2]/xx^2)^2 * sigx^2)
  fitf <- (calc - yy)/sigt
  fitf
} 
NLS <- nls.lm(par = c(1, 1, 1), fn = fct, control = nls.lm.control(maxiter = 1024))
summary(NLS)


## York model with TV
x <- c(0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)
wx <- c(1000, 1000, 500, 800, 200, 80, 60, 20, 1.8, 1)
wy <- c(1, 1.8, 4, 8, 20, 20, 70, 70, 100, 500)

sigx <- sqrt(1/wx)
sigy <- sqrt(1/wy) 

fct <- function(par) {
  calc1 <- par[1] + par[2] * x
  calc2 <- (calc1 - par[1])/par[2]
  fitf1 <- (calc1 - y)/sigy
  fitf2 <- (calc2 - x)/sigx
  fitf <- fitf1 + fitf2
  fitf
} 
NLS <- nls.lm(par = c(6, -1), fn = fct, control = nls.lm.control(maxiter = 1024))
summary(NLS)


## York model with Joel's EIV2
x <- c(0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)
wx <- c(1000, 1000, 500, 800, 200, 80, 60, 20, 1.8, 1)
wy <- c(1, 1.8, 4, 8, 20, 20, 70, 70, 100, 500)

sigx <- sqrt(1/wx)
sigy <- sqrt(1/wy) 

fct <- function(par) {
  calc <- par[1] + par[2] * x
  sigt <- sqrt(sigy^2 + par[2]^2 * sigx^2)
  fitf <- (calc - y)/sigt
  fitf
} 
NLS <- nls.lm(par = c(6, -1), fn = fct, control = nls.lm.control(maxiter = 1024))
summary(NLS)

deming(y ~ x, xstd = 1/sqrt(wx), yst = 1/sqrt(wy))

