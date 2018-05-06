\name{propagate}
\alias{propagate}
\encoding{latin1}

\title{Propagation of uncertainty using higher-order Taylor expansion and Monte Carlo simulation}

\description{
A general function for the calculation of uncertainty propagation by first-/second-order Taylor expansion and Monte Carlo simulation including covariances. Input data can be any symbolic/numeric differentiable expression and data based on summaries (mean & s.d.) or sampled from distributions. Uncertainty propagation is based completely on matrix calculus accounting for full covariance structure. Monte Carlo simulation is conducted using a multivariate t-distribution with covariance structure. Propagation confidence intervals are calculated from the expanded uncertainties by means of the degrees of freedom obtained from \code{\link{WelchSatter}}, or from the [\eqn{\frac{\alpha}{2}, 1-\frac{\alpha}{2}}] quantiles of the MC evaluations.
}

\usage{
propagate(expr, data, second.order = TRUE, do.sim = TRUE, cov = TRUE, 
          df = NULL, nsim = 1000000, alpha = 0.05, ...)  
}

\arguments{
  \item{expr}{an expression, such as \code{expression(x/y)}.}
  \item{data}{a dataframe or matrix containing either a) the means \eqn{\mu_i}, standard deviations \eqn{\sigma_i} and degrees of freedom \eqn{\nu_i} (optionally) in the first, second and third (optionally) row, or b) sampled data generated from any of \R's \code{\link{distributions}} or those implemented in this package (\code{\link{rDistr}}). If \code{nrow(data)} > 3, sampled data is assumed. The column names must match the variable names.}
  \item{second.order}{logical. If \code{TRUE}, error propagation will be calculated with first- and second-order Taylor expansion. See 'Details'.}  
  \item{do.sim}{logical. Should Monte Carlo simulation be applied?}
  \item{cov}{logical or variance-covariance matrix with the same column names as \code{data}. See 'Details'.}
  \item{df}{an optional scalar with the total degrees of freedom \eqn{\nu_{\mathrm{tot}}} of the system.}
  \item{nsim}{the number of Monte Carlo simulations to be performed, minimum is 10000.}  
  \item{alpha}{the 1 - confidence level.}
  \item{...}{other parameters to be supplied to future methods.}
}

\details{
The implemented methods are:\cr\cr
1) \bold{Monte Carlo simulation:}\cr
For each variable \eqn{m} in \code{data}, simulated data \eqn{X = [x_1, x_2, \ldots, x_n]} with \eqn{n} = \code{nsim} samples is generated from a multivariate t-distribution \eqn{X_{m, n} \sim t(\mu, \Sigma, \nu)} using means \eqn{\mu_i} and covariance matrix \eqn{\boldsymbol{\Sigma}} constructed from the standard deviations \eqn{\sigma_i} of each variable. All data is coerced into a new dataframe that has the same covariance structure as the initial \code{data}: \eqn{\boldsymbol{\Sigma}(\mathtt{data}) = \boldsymbol{\Sigma}(X_{m, n})}. Each row \eqn{i = 1, \ldots, n} of the simulated dataset \eqn{X_{m, n}} is evaluated with \code{expr}, \eqn{y_i = f(x_{m, i})}, and summary statistics (mean, sd, median, mad, quantile-based confidence interval based on [\eqn{\frac{\alpha}{2}, 1-\frac{\alpha}{2}}]) are calculated on \eqn{y}. 

2) \bold{Error propagation:}\cr
The propagated error is calculated by first-/second-order Taylor expansion accounting for full covariance structure using matrix algebra.\cr
The following transformations based on two variables \eqn{x_1, x_2} illustrate the equivalence of the matrix-based approach with well-known classical notations:\cr
\bold{First-order mean:} \eqn{\rm{E[y]} = f(\bar{x}_i)}\cr
\bold{First-order variance:} \eqn{\sigma_y^2 = {\color{red} \nabla \mathbf{\Sigma} \nabla^T}}:\cr
\deqn{{ \color{red}[\rm{j_1}\; \rm{j_2}] \left[ \begin{array}{cc} \sigma_1^2 & \sigma_1\sigma_2 \\ \sigma_2\sigma_1 & \sigma_2^2 \end{array} \right] \left[ \begin{array}{c} \rm{j_1} \\ \rm{j_2} \end{array} \right]} = \rm{j_1}^2 \sigma_1^2 + \rm{2 j_1 j_2} \sigma_1 \sigma_2 + \rm{j_2}^2 \sigma_2^2}
\deqn{= \underbrace{\sum_{i=1}^2 \rm{j_i}^2 \sigma_i^2 + 2\sum_{i=1\atop i \neq k}^2\sum_{k=1\atop k \neq i}^2 \rm{j_i j_k} \sigma_{ik}}_{\rm{classical\;notation}} = \frac{1}{1!} \left(\sum_{i=1}^2 \frac{\partial f}{\partial x_i} \sigma_i \right)^2}
\cr
\bold{Second-order mean:} \eqn{\rm{E}[y] = f(\bar{x}_i) + {\color{blue} \frac{1}{2}\rm{tr}(\mathbf{H\Sigma)}}}:\cr
\deqn{{ \color{blue} \frac{1}{2} \rm{tr} \left[ \begin{array}{cc} \rm{h_1} & \rm{h_2} \\ \rm{h_3} & \rm{h_4} \end{array} \right] \left[ \begin{array}{cc} \sigma_1^2 & \sigma_1\sigma_2 \\ \sigma_2\sigma_1 & \sigma_2^2 \end{array} \right]} = \frac{1}{2} \rm{tr} \left[ \begin{array}{cc} \rm{h_1} \sigma_1^2 + \rm{h_2}\sigma_1\sigma_2 & \rm{h_1}\sigma_1\sigma_2 + \rm{h_2}\sigma_2^2 \\ \rm{h_3} \sigma_1^2 + \rm{h_4} \sigma_1\sigma_2 & \rm{h_3} \sigma_1\sigma_2 + \rm{h_4} \sigma_2^2 \end{array} \right]}
\deqn{ = \frac{1}{2}(\rm{h_1}\sigma_1^2 + \rm{h_2}\sigma_1\sigma_2 + \rm{h_3}\sigma_1\sigma_2 + \rm{h_4}\sigma_2^2) = \frac{1}{2!} \left(\sum_{i=1}^2 \frac{\partial}{\partial x_i} \sigma_i \right)^2 \it f} 
\cr
\bold{Second-order variance:} \eqn{\sigma_y^2 = {\color{red} \nabla\mathbf{\Sigma}\nabla^T} + {\color{blue} \frac{1}{2}\rm{tr}(\mathbf{H\Sigma H\Sigma)}}}:\cr
\deqn{{\color{blue}\frac{1}{2} \rm{tr} \left[ \begin{array}{cc} \rm{h_1} & \rm{h_2} \\ \rm{h_3} & \rm{h_4} \end{array} \right] \left[ \begin{array}{cc} \rm{\sigma_1^2} & \rm{\sigma_1\sigma_2} \\ \rm{\sigma_2\sigma_1} & \rm{\sigma_2^2} \end{array} \right] \left[ \begin{array}{cc} \rm{h_1} & \rm{h_2} \\ \rm{h_3} & \rm{h_4} \end{array} \right] \left[ \begin{array}{cc} \rm{\sigma_1^2} & \rm{\sigma_1\sigma_2} \\ \rm{\sigma_2\sigma_1} & \rm{\sigma_2^2} \end{array} \right]} = \ldots}
\deqn{= \frac{1}{2} (\rm{h_1}^2\sigma_1^4 + \rm{2h_1h_2}\sigma_1^3\sigma_2 + \rm{2h_1h_3}\sigma_1^3\sigma_2 + \rm{h_2}^2\sigma_1^2\sigma_2^2 + \rm{2h_2h_3}\sigma_1^2\sigma_2^2 + \rm{h_3}^2\sigma_1^2\sigma_2^2 + \rm{2h_1h_4}\sigma_1^2\sigma_2^2}
\deqn{+ \rm{2h_2h_4}\sigma_1\sigma_2^3 + \rm{2h_3h_4}\sigma_1\sigma_2^3 + \rm{h_4}^2\sigma_2^4 = \frac{1}{2} (\rm{h_1}\sigma_1^2 + \rm{h_2}\sigma_1\sigma_2 + \rm{h_3}\sigma_1\sigma_2 + \rm{h_4}\sigma_2^2)^2}
\deqn{= \frac{1}{2!} \left( \left(\sum_{i=1}^2 \frac{\partial}{\partial x_i} \sigma_i \right)^2 \it f \right)^2}
\cr
with \eqn{\mathrm{E}(y)} = expectation of \eqn{y}, \eqn{\mathbf{\sigma_y^2}} = variance of \eqn{y}, \eqn{{\color{red} \nabla}} = the p x n gradient matrix with all partial first derivatives \eqn{{\color{red} \rm{j_i}}}, \eqn{\mathbf{\Sigma}} = the p x p covariance matrix, \eqn{{\color{blue}\mathbf{H}}} the Hessian matrix with all partial second derivatives \eqn{{\color{blue} \rm{h_i}}}, \eqn{\sigma_i} = the uncertainties and \eqn{\rm{tr}(\cdot)} = the trace (sum of diagonal) of a matrix. Note that because the Hessian matrices are symmetric, \eqn{{\color{blue} \rm{h_2}} = {\color{blue} \rm{h_3}}}. For a detailed derivation, see 'References'.\cr
The second-order Taylor expansion corrects for bias in nonlinear expressions as the first-order Taylor expansion assumes linearity around \eqn{\bar{x}_i}. There is also a Python library available for second-order error propagation ('soerp', \url{https://pypi.python.org/pypi/soerp}). The 'propagate' package gives \bold{exactly} the same results, see last example under "Examples".\cr
Depending on the input expression, the uncertainty propagation may result in an error that is not normally distributed. The Monte Carlo simulation, starting with a symmetric t-distributions of the variables, can clarify this. For instance, a high tendency from deviation of normality is encountered in formulas in which the error of the denominator is relatively large or in exponential models with a large error in the exponent.\cr 

For setups in which there is no symbolic derivation possible (i.e. \code{e <- expression(abs(x))} => "Function 'abs' is not in the derivatives table") the function automatically switches from symbolic (using \code{\link{makeGrad}} or \code{\link{makeHess}}) to numeric (\code{\link{numGrad}} or \code{\link{numHess}}) differentiation.\cr

The function will try to evaluate the expression in an environment using \code{\link{eval}} which results in a significant speed enhancement (~ 10-fold). If that fails, evaluation is done over the rows of the simulated data using \code{\link{apply}}.

\code{cov} is used in the following ways:\cr
1) If \eqn{\mu_i, \sigma_i} are supplied, a covariance matrix is built with diagonals \eqn{\sigma_i^2}, independent of \code{cov = TRUE, FALSE}.\cr
2) When simulated data is supplied, a covariance matrix is constructed that either has (\code{cov = TRUE}) or has not (\code{cov = FALSE}) off-diagonal covariances.\cr
3) The user can supply an own covariance matrix \eqn{\Sigma}, with the same column/row names as in \code{data}.

The expanded uncertainty used for constructing the confidence interval is calculated from the Welch-Satterthwaite degrees of freedom \eqn{\nu_{\mathrm{WS}}} of the \code{\link{WelchSatter}} function.
}

\value{
A list with the following components:   
\item{gradient}{the symbolic gradient vector \eqn{\nabla} of partial first-order derivatives.}
\item{evalGrad}{the evaluated gradient vector \eqn{\nabla} of partial first-order derivatives, also known as the "sensitivity". See \code{\link{summary.propagate}}.}
\item{hessian}{the symbolic Hessian matrix \eqn{\mathbf{H}} of partial second-order derivatives.}
\item{evalHess}{the evaluated Hessian matrix \eqn{\mathbf{H}} of partial second-order derivatives.}
\item{rel.contr}{the relative contribution matrix, see \code{\link{summary.propagate}}.}
\item{covMat}{the covariance matrix \eqn{\mathbf{\Sigma}} used for Monte Carlo simulation and uncertainty propagation.} 
\item{ws.df}{the Welch-Satterthwaite degrees of freedom \eqn{\nu_{\mathrm{ws}}}, as obtained from \code{\link{WelchSatter}}.}
\item{k}{the coverage factor \eqn{k}, as calculated by \eqn{t(1-(\alpha/2), \nu_{\mathrm{ws}})}.}
\item{u.exp}{the expanded uncertainty, \eqn{k\sigma(y)}, where \eqn{\sigma(y)} is derived either from the second-order uncertainty, if successfully calculated, or first-order otherwise.}
\item{resSIM}{a vector containing the \code{nsim} values obtained from the row-wise expression evaluations \eqn{f(x_{m, i})} of the simulated data in \code{datSIM}.}
\item{datSIM}{a vector containing the \code{nsim} simulated multivariate values for each variable in column format.}  
\item{prop}{a summary vector containing first-/second-order expectations and uncertainties as well as the confidence interval based on \code{alpha}.} 
\item{sim}{a summary vector containing the mean, standard deviation, median, MAD as well as the confidence interval based on \code{alpha}.}
\item{expr}{the original expression \code{expr}.}
\item{data}{the original data \code{data}.}
\item{alpha}{the otiginal \code{alpha}.}
}

\author{
Andrej-Nikolai Spiess
}   

\references{
\bold{Error propagation (in general):}\cr
An Introduction to error analysis.\cr
Taylor JR.\cr
University Science Books (1996), New York.

Evaluation of measurement data - Guide to the expression of uncertainty in measurement.\cr
JCGM 100:2008 (GUM 1995 with minor corrections).\cr
\url{http://www.bipm.org/utils/common/documents/jcgm/JCGM_100_2008_E.pdf}.

Evaluation of measurement data - Supplement 1 to the Guide to the expression of uncertainty in measurement - Propagation of distributions using a Monte Carlo Method.\cr
JCGM 101:2008.\cr
\url{http://www.bipm.org/utils/common/documents/jcgm/JCGM_101_2008_E.pdf}.

\bold{Higher-order Taylor expansion:}\cr
On higher-order corrections for propagating uncertainties.\cr
Wang CM & Iyer HK.\cr
\emph{Metrologia} (2005), \bold{42}: 406-410.

Propagation of uncertainty: Expressions of second and third order uncertainty with third and fourth moments.\cr
Mekid S & Vaja D.\cr
\emph{Measurement} (2008), \bold{41}: 600-609.

\bold{Matrix algebra for error propagation:}\cr
An Introduction to Error Propagation: Derivation, Meaning and Examples of Equation Cy = FxCxFx^t.\cr
\url{www.nada.kth.se/~kai-a/papers/arrasTR-9801-R3.pdf}.

Second order nonlinear uncertainty modeling in strapdown integration using MEMS IMUs.\cr
Zhang M, Hol JD, Slot L, Luinge H.\cr
2011 Proceedings of the 14th International Conference on Information Fusion (FUSION) (2011).

Uncertainty propagation in non-linear measurement equations.\cr
Mana G & Pennecchi F.\cr
\emph{Metrologia} (2007), \bold{44}: 246-251.

A compact tensor algebra expression of the law of propagation of uncertainty.\cr
Bouchot C, Quilantan JLC, Ochoa JCS.\cr
\emph{Metrologia} (2011), \bold{48}: L22-L28.

Nonlinear error propagation law.\cr
Kubacek L.\cr
\emph{Appl Math} (1996), \bold{41}: 329-345.

\bold{Monte Carlo simulation (normal- and t-distribution):}\cr
MUSE: computational aspects of a GUM supplement 1 implementation.\cr
Mueller M, Wolf M, Roesslein M.\cr
\emph{Metrologia} (2008), \bold{45}: 586-594.

Copulas for uncertainty analysis.\cr
Possolo A.\cr
\emph{Metrologia} (2010), \bold{47}: 262-271.

\bold{Multivariate normal distribution:}\cr
Stochastic Simulation.\cr
Ripley BD.\cr
Stochastic Simulation (1987). Wiley. Page 98.

\bold{Testing for normal distribution:}\cr
Testing for  Normality.\cr
Thode Jr. HC.\cr
Marcel Dekker (2002), New York.

Approximating the Shapiro-Wilk W-test for non-normality.\cr
Royston P.\cr 
\emph{Stat Comp} (1992), \bold{2}: 117-119.
}

\examples{
## In these examples, 'nsim = 100000' to save
## Rcmd check time (CRAN). It is advocated
## to use at least 'nsim = 1000000' though...

## Example without given degrees-of-freedom.
EXPR1 <- expression(x/y)
x <- c(5, 0.01)
y <- c(1, 0.01)
DF1 <- cbind(x, y)
RES1 <- propagate(expr = EXPR1, data = DF1, type = "stat", 
                  do.sim = TRUE, verbose = TRUE, 
                  nsim = 100000)
RES1

## Same example with given degrees-of-freedom
## => third row in input 'data'.
EXPR2 <- expression(x/y)
x <- c(5, 0.01, 12)
y <- c(1, 0.01, 5)
DF2 <- cbind(x, y)
RES2 <- propagate(expr = EXPR2, data = DF2, type = "stat", 
                  do.sim = TRUE, verbose = TRUE,
                  nsim = 100000)
RES2

## With the 'summary' function, we can get the
## Welch-Satterthwaite DF's, coverage, expanded uncertainty,
## Gradient and Hessian matrix etc.
summary(RES2)

## Example using a recursive function:
## no Taylor expansion possible, only Monte-Carlo.
a <- c(5, 0.1)
b <- c(100, 2)
DAT <- cbind(a, b)

f <- function(a, b) {
  N <- 0
  for (i in 1:100) {
    N <- N + i * log(a) + b^(1/i)
  }
  return(N)
}

propagate(f, DAT, nsim = 100000)

\dontrun{
################# GUM 2008 (1) ########################
## Example in Annex H.1 from the GUM 2008 manual
## (see 'References'), an end gauge calibration
## study. We use only first-order error propagation,
## with total df = 16 and alpha = 0.01, 
## as detailed in GUM H.1.6.
EXPR3 <- expression(ls + d - ls * (da * the + as * dt))
ls <- c(50000623, 25)
d <- c(215, 9.7)
da <- c(0, 0.58E-6)
the <- c(-0.1, 0.41)
as <- c(11.5E-6, 1.2E-6)
dt <- c(0, 0.029)
DF3 <- cbind(ls, d, da, the, as, dt)
RES3 <- propagate(expr = EXPR3, data = DF3, second.order = FALSE,
                  df = 16, alpha = 0.01)
RES3
## propagate: sd.1 = 31.71 
## GUM H.1.4/H.6c: u = 32  

## Expanded uncertainty, from summary function.
summary(RES3)
## propagate: 92.62
## GUM H.1.6: 93

## Proof that covariance of Monte-Carlo
## simulated dataset is "fairly"" the same 
## as from initial data.
RES3$covMat
cov(RES3$datSIM)
all.equal(RES3$covMat, cov(RES3$datSIM))

## Now using second-order Taylor expansion.
RES4 <- propagate(expr = EXPR3, data = DF3)
RES4
## propagate: sd.2 = 33.91115
## GUM H.1.7: u = 34.
## Also similar to the non-matrix-based approach
## in Wang et al. (2005, page 408): u1 = 33.91115.
## NOTE: After second-order correction ("sd.2"), 
## uncertainty is more similar to the uncertainty
## obtained from Monte Carlo simulation!

#################### GUM 2008 (2) #################
## Example in Annex H.2 from the GUM 2008 manual
## (see 'References'), simultaneous resistance
## and reactance measurement.
data(H.2)

## This gives exactly the means, uncertainties and
## correlations as given in Table H.2:
colMeans(H.2)
sqrt(colVarsC(H.2))/sqrt(5)
cor(H.2)

## H.2.3 Approach 1 using mean values and
## standard uncertainties:
EXPR6a <- expression((V/I) *  cos(phi)) ## R
EXPR6b <- expression((V/I) *  sin(phi)) ## X
EXPR6c <- expression(V/I) ## Z
MEAN6 <- colMeans(H.2)
SD6 <- sqrt(colVarsC(H.2))
DF6 <- rbind(MEAN6, SD6)
COV6ab <- cov(H.2) ## covariance matrix of V, I, phi
COV6c <- cov(H.2[, 1:2])  ## covariance matrix of V, I

RES6a <- propagate(expr = EXPR6a, data = DF6, cov = COV6ab)
RES6b <- propagate(expr = EXPR6b, data = DF6, cov = COV6ab)
RES6c <- propagate(expr = EXPR6c, data = DF6[, 1:2], 
                   cov = COV6c)

## This gives exactly the same values of mean and sd/sqrt(5)
## as given in Table H.4.
RES6a$prop # 0.15892/sqrt(5) = 0.071
RES6b$prop # 0.66094/sqrt(5) = 0.296
RES6c$prop # 0.52846/sqrt(5) = 0.236

######### GUM 2008 Supplement 1 (1) #######################
## Example from 9.2.2 of the GUM 2008 Supplement 1
## (see 'References'), normally distributed input
## quantities. Assign values as in 9.2.2.1.
EXPR7 <- expression(X1 + X2 + X3 + X4)
X1 <- c(0, 1)
X2 <- c(0, 1)
X3 <- c(0, 1)
X4 <- c(0, 1)
DF7 <- cbind(X1, X2, X3, X4)
RES7 <- propagate(expr = EXPR7, data = DF7, nsim = 1E5)
## This will give exactly the same results as in 
## 9.2.2.6, Table 2.
RES7

######### GUM 2008 Supplement 1 (2) #######################
## Example from 9.3 of the GUM 2008 Supplement 1
## (see 'References'), mass calibration.
## Formula 24 in 9.3.1.3 and values as in 9.3.1.4, Table 5.
EXPR8 <- expression((Mrc + dMrc) * (1 + (Pa - Pa0) * ((1/Pw) - (1/Pr))) - Mnom)
Mrc <- rnorm(1E5, 100000, 0.050)
dMrc <- rnorm(1E5, 1.234, 0.020)
Pa <- runif(1E5, 1.10, 1.30)  ## E(Pa) = 1.2, (b-a)/2 = 0.1 
Pw <- runif(1E5, 7000, 9000)  ## E(Pw) = 8000, (b-a)/2 = 1000
Pr <- runif(1E5, 7950, 8050) ## E(Pr) = 8000, (b-a)/2 = 50
Pa0 <- 1.2 
Mnom <- 100000
DF8 <- cbind(Mrc, dMrc, Pa, Pw, Pr, Pa0, Mnom)
RES8 <- propagate(expr = EXPR8, data = DF8, nsim = 1E5)
## This will give exactly the same results as in 
## 9.3.2.3, Table 6
RES8
RES8
 
######### GUM 2008 Supplement 1 (3) #######################
## Example from 9.4 of the GUM 2008 Supplement 1
## (see 'References'), comparioson loss in microwave
## power meter calibration, zero covariance.
## Formula 28 in 9.4.1.5 and values as in 9.4.1.7.
EXPR9 <- expression(X1^2 - X2^2)
X1 <- c(0.050, 0.005)
X2 <- c(0, 0.005)
DF9 <- cbind(X1, X2)
RES9a <- propagate(expr = EXPR9, data = DF9, nsim = 1E5)
## This will give exactly the same results as in 
## 9.4.2.2.7, Table 8, x1 = 0.050.
RES9a

## Using covariance matrix with r(x1, x2) = 0.9
## We convert to covariances using cor2cov.
COR9 <- matrix(c(1, 0.9, 0.9, 1), nrow = 2)
COV9 <- cor2cov(COR9, c(0.005^2, 0.005^2))
colnames(COV9) <- c("X1", "X2")
rownames(COV9) <- c("X1", "X2")
RES9b <- propagate(expr = EXPR9, data = DF9, cov = COV9)
## This will give exactly the same results as in 
## 9.4.3.2.1, Table 9, x1 = 0.050.
RES9b

######### GUM 2008 Supplement 1 (4) #######################
## Example from 9.5 of the GUM 2008 Supplement 1
## (see 'References'), gauge block calibration.
## Assignment of PDF's as in Table 10 of 9.5.2.1.
EXPR10 <- expression(Ls + D + d1 + d2 - Ls *(da *(t0 + Delta) + as * dt) - Lnom)
Lnom <- 50000000
Ls <- propagate:::rst(1000000, mean = 50000623, sd  = 25, df = 18)
D <- propagate:::rst(1000000, mean = 215, sd = 6, df = 25)
d1 <- propagate:::rst(1000000, mean = 0, sd = 4, df = 5)
d2 <- propagate:::rst(1000000, mean = 0, sd = 7, df = 8)
as <- runif(1000000, 9.5E-6, 13.5E-6)
t0 <- rnorm(1000000, -0.1, 0.2)
Delta <- propagate:::rarcsin(1000000, -0.5, 0.5)
da <- propagate:::rctrap(1000000, -1E-6, 1E-6, 0.1E-6)
dt <- propagate:::rctrap(1000000, -0.050, 0.050, 0.025)
DF10 <- cbind(Ls, D, d1, d2, as, t0, Delta, da, dt, Lnom)
RES10 <- propagate(expr = EXPR10, data = DF10, cov = FALSE, alpha = 0.01)
RES10
## This gives the same results as in 9.5.4.2, Table 11.
## However: results are exacter than in the GUM 2008
## manual, especially when comparing sd(Monte Carlo) with sd.2!
## GUM 2008 gives 32 and 36, respectively.
RES10

########## Comparison to Pythons 'soerp' ###################
## Exactly the same results as under 
## https://pypi.python.org/pypi/soerp ! 
EXPR11 <- expression(C * sqrt((520 * H * P)/(M *(t + 460))))
H <- c(64, 0.5)
M <- c(16, 0.1)
P <- c(361, 2)
t <- c(165, 0.5)
C <- c(38.4, 0) 
DAT11 <- makeDat(EXPR11)
RES11 <- propagate(expr = EXPR11, data = DAT11) 
RES11
}   
}

\keyword{array}
\keyword{algebra}
\keyword{multivariate}

