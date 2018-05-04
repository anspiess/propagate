\name{predictNLS}
\alias{predictNLS}
\encoding{latin1}

\title{Confidence/prediction intervals for (weighted) nonlinear models based on uncertainty propagation}

\description{A function for calculating confidence/prediction intervals of (weighted) nonlinear models for the supplied or new predictor values, by using first-/second-order Taylor expansion and Monte Carlo simulation. This approach can be used to construct more realistic error estimates and confidence/prediction intervals for nonlinear models than what is possible with only a simple linearization (first-order Taylor expansion) approach. Another application is when there is an "error in x" setup with uncertainties in the predictor variable (See 'Examples'). This function will also work in the presence of multiple predictors with/without errors.
}

\usage{
predictNLS(model, newdata, newerror, interval = c("confidence", "prediction", "none"),
           alpha = 0.05, ...)
}

\arguments{
  \item{model}{a model obtained from \code{\link{nls}} or \code{nlsLM} (package 'minpack.lm').}
  \item{newdata}{a data frame with new predictor values, having the same column names as in \code{model}. See \code{\link{predict.nls}} and 'Examples'. If omitted, the model's predictor values are employed.}
  \item{newerror}{a data frame with optional error values, having the same column names as in \code{model} and in the same order as in \code{newdata}. See 'Examples'.}
  \item{interval}{A character string indicating if confidence/prediction intervals are to be calculated or not.}
  \item{alpha}{the \eqn{\alpha} level.}
  \item{...}{other parameters to be supplied to \code{\link{propagate}}.}
}

\details{
Calculation of the propagated uncertainty \eqn{\sigma_y} using \eqn{\nabla \Sigma \nabla^T} is called the "Delta Method" and is widely applied in NLS fitting. However, this method is based on first-order Taylor expansion and thus assummes linearity around \eqn{f(x)}. The second-order approach as implemented in the \code{\link{propagate}} function can partially correct for this restriction by using a second-order polynomial around \eqn{f(x)}.\cr
Confidence and prediction intervals are calculated in a usual way using \eqn{t(1 - \frac{\alpha}{2}, \nu) \cdot \sigma_y} (1) or \eqn{t(1 - \frac{\alpha}{2}, \nu) \cdot \sqrt{\sigma_y^2 + \textcolor{red}{\sigma_r^2}}} (2), respectively, where the residual variance \eqn{\textcolor{red}{\sigma_r^2} = \frac{\sum_{i=1}^n (y_i - \hat{y}_i)^2}{n - \nu}} (3).
The inclusion of \eqn{\textcolor{red}{\sigma_r^2}} in the prediction interval is implemented as an extended gradient and "augmented" covariance matrix. So instead of using \eqn{y = f(x, \beta)} (4) we take \eqn{y = f(x, \beta) + \textcolor{red}{\sigma_r^2}} (5) as the expression and augment the \eqn{n \times n} covariance matrix \eqn{C} to an \eqn{n+1 \times n+1} covariance matrix, where \eqn{C_{n+1, n+1} = \textcolor{red}{\sigma_r^2}}. Partial differentiation and matrix multiplication will then yield, for example with two coefficients \eqn{\beta_1} and \eqn{\beta_2} and their corresponding covariance matrix \eqn{\Sigma}:     
\deqn{\left[\frac{\partial f}{\partial \beta_1}\; \frac{\partial f}{\partial \beta_2}\; \textcolor{red}{1}\right] \left[ \begin{array}{ccc} \sigma_1^2 & \sigma_1\sigma_2 & 0 \\ \sigma_2\sigma_1 & \sigma_2^2 & 0 \\ 0 & 0 & \textcolor{red}{\sigma_r^2} \end{array} \right] \left[ \begin{array}{c} \frac{\partial f}{\partial \beta_1} \\ \frac{\partial f}{\partial \beta_2} \\ \textcolor{red}{1} \end{array} \right] = \left(\frac{\partial f}{\partial \beta_1}\right)^2\sigma_1^2 + 2 \frac{\partial f}{\partial \beta_1} \frac{\partial f}{\partial \beta_2} \sigma_1 \sigma_2 + \left(\frac{\partial f}{\partial \beta_2}\right)^2\sigma_2^2 + \textcolor{red}{\sigma_r^2}}
\eqn{\equiv \sigma_y^2 + \textcolor{red}{\sigma_r^2}}, where \eqn{\sigma_y^2 + \textcolor{red}{\sigma_r^2}} then goes into (2).\cr
The advantage of the augmented covariance matrix is that it can be exploited for creating Monte Carlo simulation-based prediction intervals. This is new from version 1.0-6 and is based on the paradigm that we simply add another dimension with \eqn{\mu = 0} and \eqn{\sigma^2 = \textcolor{red}{\sigma_r^2}} to the multivariate t-distribution random number generator (in our case \code{\link{rtmvt}}). All \eqn{n} simulations are then evaluated with (5) and the usual \eqn{[1 - \frac{\alpha}{2}, \frac{\alpha}{2}]} quantiles calculated.\cr 
If errors are supplied to the predictor values in \code{newerror}, they need to have the same column names and order than the new predictor values. 
}

\value{
A list with the following items:\cr
\code{summary}: The mean/error estimates obtained from first-/second-order Taylor expansion and Monte Carlo simulation, together with calculated confidence/prediction intervals based on asymptotic normality.\cr
\code{prop}: the complete output from \code{\link{propagate}} for each value in \code{newdata}.
}

\author{
Andrej-Nikolai Spiess
}   

\references{
Nonlinear Regression.\cr
Seber GAF & Wild CJ.\cr
John Wiley & Sons; 1ed, 2003.

Nonlinear Regression Analysis and its Applications.\cr
Bates DM & Watts DG.\cr
Wiley-Interscience; 1ed, 2007.

Statistical Error Propagation.\cr
Tellinghuisen J.\cr
\emph{J. Phys. Chem. A} (2001), \bold{47}: 3917-3921.

Least-squares analysis of data with uncertainty in x and y: A Monte Carlo
methods comparison.\cr
Tellinghuisen J.\cr
\emph{Chemometr Intell Lab} (2010), \bold{47}: 160-169.

From the author's blog:\cr
http://rmazing.wordpress.com/2013/08/14/predictnls-part-1-monte-carlo-simulation-confidence-intervals-for-nls-models/\cr
http://rmazing.wordpress.com/2013/08/26/predictnls-part-2-taylor-approximation-confidence-intervals-for-nls-models/
}

\examples{
## Example from ?nls.
DNase1 <- subset(DNase, Run == 1)
fm3DNase1 <- nls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
                 data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1))

## Using a single predictor value without error.
PROP1 <- predictNLS(fm3DNase1, newdata = data.frame(conc = 2))
PRED1 <- predict(fm3DNase1, newdata = data.frame(conc = 2))
PROP1$summary
PRED1
## => Prop.Mean.1 equal to PRED1

## Using a single predictor value with error.
PROP2 <- predictNLS(fm3DNase1, newdata = data.frame(conc = 2), newerror = data.frame(conc = 0.5))
PROP2$summary

\dontrun{
## Using a sequence of predictor values without error.
CONC <- seq(1, 12, by = 1)
PROP3 <- predictNLS(fm3DNase1, newdata = data.frame(conc = CONC))
PRED3 <- predict(fm3DNase1, newdata = data.frame(conc = CONC))
PROP3$summary
PRED3
## => Prop.Mean.1 equal to PRED3

## Plot mean and confidence values from first-/second-order 
## Taylor expansion and Monte Carlo simulation.
plot(DNase1$conc, DNase1$density)
lines(DNase1$conc, fitted(fm3DNase1), lwd = 2, col = 1)
points(CONC, PROP3$summary[, 1], col = 2, pch = 16)
lines(CONC, PROP3$summary[, 5], col = 2)
lines(CONC, PROP3$summary[, 6], col = 2)
lines(CONC, PROP3$summary[, 11], col = 4)
lines(CONC, PROP3$summary[, 12], col = 4)

## Using a sequence of predictor values with error.
PROP4 <- predictNLS(fm3DNase1, newdata = data.frame(conc = 1:5), 
                    newerror = data.frame(conc = (1:5)/10))
PROP4$summary

## Using multiple predictor values.
## 1: Setup of response values with gaussian error of 10\%.
x <- seq(1, 10, by = 0.01)
y <- seq(10, 1, by = -0.01)
a <- 2
b <- 5
c <- 10
z <- a * exp(b * x)^sin(y/c)
z <- z + sapply(z, function(x) rnorm(1, x, 0.10 * x))

## 2: Fit 'nls' model.
MOD <- nls(z ~ a * exp(b * x)^sin(y/c), 
           start = list(a = 2, b = 5, c = 10))

## 3: single newdata without errors.
DAT1 <- data.frame(x = 4, y = 3)
PROP5 <- predictNLS(MOD, newdata = DAT1)
PROP5$summary

## 4: single newdata with errors.
DAT2 <- data.frame(x = 4, y = 3)
ERR2 <- data.frame(x = 0.2, y = 0.1)
PROP6 <- predictNLS(MOD, newdata = DAT2, newerror = ERR2)
PROP6$summary

## 5: multiple newdata with errors.
DAT3 <- data.frame(x = 1:4, y = 3)
ERR3 <- data.frame(x = rep(0.2, 4), y = seq(1:4)/10)
PROP7 <- predictNLS(MOD, newdata = DAT3, newerror = ERR3)
PROP7$summary

## 6: Linear model to compare conf/pred intervals.
set.seed(123)
X <- 1:20
Y <- 3 + 2 * X + rnorm(20, 0, 2)
plot(X, Y)
LM <- lm(Y ~ X)
NLS <- nlsLM(Y ~ a + b * X, start = list(a = 3, b = 2)) 
predict(LM, newdata = data.frame(X = 14.5), interval = "conf") 
predictNLS(NLS, newdata = data.frame(X = 14.5), interval = "conf")$summary
predict(LM, newdata = data.frame(X = 14.5), interval = "pred")
predictNLS(NLS, newdata = data.frame(X = 14.5), interval = "pred")$summary
} 
}

\keyword{array}
\keyword{algebra}
\keyword{multivariate}

