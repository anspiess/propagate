\name{predictNLS}
\alias{predictNLS}
\encoding{latin1}

\title{Confidence intervals for nonlinear models based on uncertainty propagation}

\description{A function for calculating confidence intervals for the fitted values of nonlinear models by using first-/second-order Taylor expansion and Monte Carlo simulation. This approach can be used to construct more realistic error estimates and confidence/prediction intervals for nonlinear models than what is possible with only a simple linearization (first-order Taylor expansion) approach. Another application is when there is an "error in x" setup with uncertainties in the predictor variable (See 'Examples'). This function will also work in the presence of multiple predictors with/without errors.
}

\usage{
predictNLS(model, newdata, interval = c("confidence", "prediction", "none"),
           alpha = 0.05, ...)
}

\arguments{
  \item{model}{a model obtained from \code{\link{nls}} or \code{nlsLM} (package 'minpack.lm').}
  \item{newdata}{A named list or data frame in which to look for variables with which to predict. The first \eqn{n} columns must contain the predictor values, the following \eqn{n} columns can contain errors. See \code{\link{predict.nls}} and 'Examples'.}
  \item{interval}{A character string indicating if confidence/prediction intervals on the mean responses are to be calculated or not.}
  \item{alpha}{the 1 - confidence level.}
  \item{...}{other parameters to be supplied to \code{\link{propagate}}.}
}

\details{
Calculation of the propagated uncertainty \eqn{\sigma_y^2} using \eqn{\nabla_x C_x \nabla_x^T} is called the "Delta Method" and is widely applied in NLS fitting. However, this method is based on first-order Taylor expansion and thus assummes linearity around \eqn{f(x)}. The second-order approach as implemented in the \code{\link{propagate}} function can partially correct for this restriction by using a second-order polynomial around \eqn{f(x)}.\cr
Confidence/prediction intervals are calculated in a usual way using \eqn{t(1 - \frac{\alpha}{2}, \nu) \cdot \sigma_{prop}} or \eqn{t(1 - \frac{\alpha}{2}, \nu) \cdot \sqrt{\sigma_{prop}^2 + rv}}, respectively, where \eqn{rv} = the residual variance \eqn{\frac{\sum_{i=1}^n (y_i - \hat{y}_i)^2}{n - \nu}}.
If errors are supplied to the predictor values in \code{newdata}, for \eqn{n} predictor values they have to be in the next \eqn{n} columns, i.e. predictors in \code{[, 1:3]}, and errors in \code{[, 4:6]}, see 'Examples'. 
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
PROP2 <- predictNLS(fm3DNase1, newdata = data.frame(conc = 2, error = 0.2))
PRED2 <- predict(fm3DNase1, newdata = data.frame(conc = 2))
PROP2$summary
PRED2

\dontrun{
## Using a sequence of predictor values without error.
CONC <- seq(1, 12, by = 1)
PROP3 <- predictNLS(fm3DNase1, newdata = data.frame(conc = CONC))
PRED3 <- predict(fm3DNase1, newdata = data.frame(conc = CONC))
PROP3$summary
PRED3
## => Prop.Mean.1 equal to PRED2

## Plot predicted and confidence values from 
## first-/second-order Taylor expansion 
## and Monte Carlo simulation.
plot(DNase1$conc, DNase1$density)
lines(DNase1$conc, fitted(fm3DNase1), lwd = 2, col = 1)
points(CONC, PROP3$summary[, 1], col = 2, pch = 16)
lines(CONC, PROP3$summary[, 5], col = 2)
lines(CONC, PROP3$summary[, 6], col = 2)

## Using a sequence of predictor values with error.
DAT <- data.frame(conc = CONC, error = rnorm(12, 0, 0.1))
PROP4 <- predictNLS(fm3DNase1, newdata = DAT)
PRED4 <- predict(fm3DNase1, newdata = DAT)
PROP4$summary
PRED4
## => Prop.Mean.1 equal to PRED4

## Using multiple predictor values
## 1: Setup of response values
## with gaussian error of 10\%.
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
## 3: newdata without errors and prediction.
DAT1 <- data.frame(x = 4, y = 3)
PROP5 <- predictNLS(MOD, newdata = DAT1)
PROP5$summary
## 4: newdata with errors and prediction.
DAT2 <- data.frame(x = 4, y = 3, error.x = 0.2, error.y = 0.1)
PROP6 <- predictNLS(MOD, newdata = DAT2)
PROP6$summary
} 
}

\keyword{array}
\keyword{algebra}
\keyword{multivariate}

