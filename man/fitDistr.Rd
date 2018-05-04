\name{fitDistr}
\alias{fitDistr}

\title{Fitting distributions to observations/Monte Carlo simulations}

\description{
This function fits 32 different continuous distributions by (weighted) NLS to the histogram of Monte Carlo simulation results as obtained by \code{\link{propagate}} or any other vector containing large-scale observations. Finally, the fits are sorted by ascending \code{\link{BIC}}. 
}

\usage{
fitDistr(object, nbin = 100, weights = FALSE, verbose = TRUE, 
         brute = c("fast", "slow"), plot = c("hist", "qq"),  distsel = NULL, ...)
}

\arguments{
  \item{object}{an object of class 'propagate' or a vector containing observations.}
  \item{nbin}{the number of bins in the histogram.}
  \item{weights}{numeric or logical. Either a numeric vector of weights, or if \code{TRUE}, the distributions are fitted with weights = 1/(counts per bin).}
  \item{verbose}{logical. If \code{TRUE}, steps of the analysis are printed to the console.}
  \item{brute}{complexity of the brute force approach. See 'Details'.}
  \item{plot}{if \code{"hist"}, a plot with the "best" distribution (in terms of lowest BIC) on top of the histogram is displayed. If \code{"qq"}, a QQ-Plot will display the difference between the observed and fitted quantiles.}
  \item{distsel}{a vector of distribution numbers to select from the complete cohort as listed below, e.g. \code{c(1:10, 15)}.}
  \item{...}{other parameters to be passed to the plots.}  
}

\details{
Fits the following 32 distributions using (weighted) residual sum-of-squares as the minimization criterion for \code{\link{nls.lm}}:\cr
1) Normal distribution (\code{dnorm}) => \url{https://en.wikipedia.org/wiki/Normal_distribution}\cr
2) Skewed-normal distribution (\code{propagate:::dsn}) => \url{https://en.wikipedia.org/wiki/Skew_normal_distribution}\cr
3) Generalized normal distribution (\code{propagate:::dgnorm}) => \url{https://en.wikipedia.org/wiki/Generalized_normal_distribution}\cr 
4) Log-normal distribution (\code{dlnorm}) => \url{https://en.wikipedia.org/wiki/Log-normal_distribution}\cr 
5) Scaled and shifted t-distribution (\code{propagate:::dst}) => GUM 2008, Chapter 6.4.9.2.\cr
6) Logistic distribution (\code{dlogis}) => \url{https://en.wikipedia.org/wiki/Logistic_distribution}\cr 
7) Uniform distribution (\code{dunif}) => \url{https://en.wikipedia.org/wiki/Uniform_distribution_(continuous)}\cr 
8) Triangular distribution (\code{propagate:::dtriang}) => \url{https://en.wikipedia.org/wiki/Triangular_distribution}\cr 
9) Trapezoidal distribution (\code{propagate:::dtrap}) => \url{https://en.wikipedia.org/wiki/Trapezoidal_distribution}\cr 
10) Curvilinear Trapezoidal distribution (\code{propagate:::dctrap}) => GUM 2008, Chapter 6.4.3.1\cr
11) Gamma distribution (\code{dgamma}) => \url{https://en.wikipedia.org/wiki/Gamma_distribution}\cr 
12) Inverse Gamma distribution (\code{propagate:::dinvgamma}) => \url{https://en.wikipedia.org/wiki/Inverse-gamma_distribution}\cr 
13) Cauchy distribution (\code{dcauchy}) => \url{https://en.wikipedia.org/wiki/Cauchy_distribution}\cr 
14) Laplace distribution (\code{propagate:::dlaplace}) => \url{https://en.wikipedia.org/wiki/Laplace_distribution}\cr 
15) Gumbel distribution (\code{propagate:::dgumbel}) => \url{https://en.wikipedia.org/wiki/Gumbel_distribution}\cr 
16) Johnson SU distribution (\code{propagate:::dJSU}) => \url{https://en.wikipedia.org/wiki/Johnson_SU_distribution}\cr 
17) Johnson SB distribution (\code{propagate:::dJSB}) => \url{https://www.mathwave.com/articles/johnson_sb_distribution.html}\cr 
18) Three-parameter Weibull distribution (\code{propagate:::dweibull2}) => \url{https://en.wikipedia.org/wiki/Weibull_distribution}\cr 
19) Two-parameter beta distribution (\code{dbeta2}) => \url{https://en.wikipedia.org/wiki/Beta_distribution#Two_parameters_2}\cr
20) Four-parameter beta distribution (\code{propagate:::dbeta2}) => \url{https://en.wikipedia.org/wiki/Beta_distribution#Four_parameters_2}\cr 
21) Arcsine distribution (\code{propagate:::darcsin}) => \url{https://en.wikipedia.org/wiki/Arcsine_distribution}\cr 
22) Von Mises distribution (\code{propagate:::dmises}) => \url{https://en.wikipedia.org/wiki/Von_Mises_distribution}\cr 
23) Inverse Gaussian distribution (\code{propagate:::dinvgauss}) => \url{https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution}\cr 
24) Generalized Extreme Value distribution (\code{propagate:::dgevd}) => \url{https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution}\cr 
25) Rayleigh distribution (\code{propagate:::drayleigh}) => \url{https://en.wikipedia.org/wiki/Rayleigh_distribution}\cr
26) Chi-square distribution (\code{dchisq}) => \url{https://en.wikipedia.org/wiki/Chi-squared_distribution}\cr
27) Exponential distribution (\code{dexp}) => \url{https://en.wikipedia.org/wiki/Exponential_distribution}\cr
28) F-distribution (\code{df}) => \url{https://en.wikipedia.org/wiki/F-distribution}\cr
29) Burr distribution (\code{dburr}) => \url{https://en.wikipedia.org/wiki/Burr_distribution}\cr
30) Chi distribution (\code{dchi}) => \url{https://en.wikipedia.org/wiki/Chi_distribution}\cr
31) Inverse Chi-square distribution (\code{dinvchisq}) => \url{https://en.wikipedia.org/wiki/Inverse-chi-squared_distribution}\cr
32) Cosine distribution (\code{dcosine}) => \url{https://en.wikipedia.org/wiki/Raised_cosine_distribution}\cr

All distributions are fitted with a brute force approach, in which the parameter space is extended over three orders of magnitude \eqn{(0.1, 1, 10)\times \beta_i} when \code{brute = "fast"}, or five orders \eqn{(0.01, 0.1, 1, 10, 100)\times \beta_i} when \code{brute = "slow"}. Approx. 20-90s are needed to fit for the fast version, depending mainly on the number of bins. 

The goodness-of-fit (GOF) is calculated with \code{\link{BIC}} from the (weighted) log-likelihood of the fit:
\deqn{\rm{ln}(L) = 0.5 \cdot \left(-N \cdot \left(\rm{ln}(2\pi) + 1 + \rm{ln}(N) - \sum_{i=1}^n log(w_i) + \rm{ln}\left(\sum_{i=1}^n w_i \cdot x_i^2\right) \right) \right)} 
\deqn{\rm{BIC} = - 2\rm{ln}(L) + (N - k)ln(N)}
with \eqn{x_i} = the residuals from the NLS fit, \eqn{N} = the length of the residual vector, \eqn{k} = the number of parameters of the fitted model and \eqn{w_i} = the weights.

In contrast to some other distribution fitting softwares (i.e. Easyfit, Mathwave) that use residual sum-of-squares/Anderson-Darling/Kolmogorov-Smirnov statistics as GOF measures, the application of BIC accounts for increasing number of parameters in the distribution fit and therefore compensates for overfitting. Hence, this approach is more similar to ModelRisk (Vose Software) and as employed in \code{fitdistr} of the 'MASS' package.
Another application is to identify a possible distribution for the raw data prior to using Monte Carlo simulations from this distribution. However, a decent number of observations should be at hand in order to obtain a realistic estimate of the proper distribution. See 'Examples'.\cr
The code for the density functions can be found in file "distr-densities.R".

IMPORTANT: It can be feasible to set \code{weights = TRUE} in order to give more weight to bins with low counts. See 'Examples'.
ALSO: Distribution fitting is highly sensitive to the number of defined histogram bins, so it is advisable to change this parameter and inspect if the order of fitted distributions remains stable.
}

\value{
A list with the following items:

\code{stat}: the by BIC value ascendingly sorted distribution names, including RSS and MSE.\cr
\code{fit}: a list of the results from \code{\link{nls.lm}} for each distribution model, also sorted ascendingly by BIC values.\cr
\code{par}: a list of the estimated parameters of the models in \code{fit}.\cr
\code{se}: a list of the parameters' standard errors, calculated from the square root of the covariance matrices diagonals.\cr
\code{dens}: a list with all density function used for fitting, sorted as in \code{fit}.\cr
\code{bestfit}: the best model in terms of lowest BIC.\cr
\code{bestpar}: the parameters of \code{bestfit}.\cr
\code{bestse}: the parameters' standard errors of \code{bestfit}.\cr
\code{fitted}: the fitted values of \code{bestfit}.\cr
\code{residuals}: the residuals of \code{bestfit}.\cr
}
     
\author{
Andrej-Nikolai Spiess
}  

\references{
Continuous univariate distributions, Volume 1.\cr
Johnson NL, Kotz S and Balakrishnan N.\cr
\emph{Wiley Series in Probability and Statistics, 2.ed} (2004).

A guide on probability distributions.\cr
R-forge distributions core team.\cr
\url{http://dutangc.free.fr/pub/prob/probdistr-main.pdf}.

Univariate distribution relationships.\cr
Leemis LM and McQueston JT.\cr
\emph{The American Statistician} (2008), \bold{62}: 45-53.
}

\examples{
\dontrun{
## Linear example, small error
## => Normal distribution.
EXPR1 <- expression(x + 2 * y)
x <- c(5, 0.01)
y <- c(1, 0.01)
DF1 <- cbind(x, y)
RES1 <- propagate(expr = EXPR1, data = DF1, type = "stat", 
                  do.sim = TRUE, verbose = TRUE)
fitDistr(RES1)

## Ratio example, larger error
## => Gamma distribution.
EXPR2 <- expression(x/2 * y)
x <- c(5, 0.1)
y <- c(1, 0.02)
DF2 <- cbind(x, y)
RES2 <- propagate(expr = EXPR2, data = DF2, type = "stat", 
                  do.sim = TRUE, verbose = TRUE)
fitDistr(RES2)

## Exponential example, large error
## => Log-Normal distribution.
EXPR3 <- expression(x^(2 * y))
x <- c(5, 0.1)
y <- c(1, 0.1)
DF3 <- cbind(x, y)
RES3 <- propagate(expr = EXPR3, data = DF3, type = "stat", 
                  do.sim = TRUE, verbose = TRUE)
fitDistr(RES3)

## Rectangular input distributions result
## in Curvilinear Trapezoidal output distribution.
A <- runif(100000, 20, 25)
B <- runif(100000, 3, 3.5)
DF4 <- cbind(A, B)
EXPR4 <- expression(A + B)
RES4 <- propagate(EXPR4, data = DF4, type = "sim", 
                 use.cov = FALSE, do.sim = TRUE)
fitDistr(RES4)        

## Fitting with 1/counts as weights.
EXPR5 <- expression(x + 2 * y)
x <- c(5, 0.05)
y <- c(1, 0.05)
DF5 <- cbind(x, y)
RES5 <- propagate(expr = EXPR5, data = DF5, type = "stat", 
                  do.sim = TRUE, verbose = TRUE, weights = TRUE)
fitDistr(RES5)

## Using only selected distributions.
EXPR6 <- expression(x + sin(y))
x <- c(5, 0.1)
y <- c(1, 0.2)
DF6 <- cbind(x, y)
RES6 <- propagate(expr = EXPR6, data = DF6, type = "stat", 
                  do.sim = TRUE)
fitDistr(RES6, distsel = c(1:10, 15, 28))
}
}   

\keyword{algebra}
\keyword{univariate}