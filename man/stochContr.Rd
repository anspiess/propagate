\name{stochContr}
\alias{stochContr}
\encoding{latin1}

\title{Stochastic contribution analysis of Monte Carlo simulation-derived propagated uncertainty}

\description{
Conducts a "stochastic contribution analysis" by calculating the change in propagated uncertainty when each of the simulated variables is kept constant at its mean, i.e. the uncertainty is removed.
}

\usage{
stochContr(prop, plot = TRUE)
}

\arguments{
  \item{prop}{a \code{propagate} object.}
  \item{plot}{logical. If \code{TRUE}, a boxplot with the original and mean-value propagated distribution.}  
}

\details{
This function takes the Monte Carlo simulated data \eqn{X_n} from a \code{\link{propagate}} object (\code{...$datSIM}), sequentially substitutes each variable \eqn{\beta_i} by its mean \eqn{\bar{\beta_i}} and then re-evaluates the output distribution \eqn{Y_n = f(\beta, X_n)}. Optional boxplots are displayed that compare the original \eqn{Y_n\mathrm{(orig)}} to those obtained from removing \eqn{\sigma} from each \eqn{\beta_i}. Finally, the relative contribution \eqn{C_i} for all \eqn{\beta_i} is calculated by \eqn{C_i = \sigma(Y_n\mathrm{(orig)})-\sigma(Y_n)}, and divided by its sum so that \eqn{\sum_{i=1}^n C_i = 1}.
}

\value{
The relative contribution \eqn{C_i} for all variables.
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
a <- c(15, 1)
b <- c(100, 5)
c <- c(0.5, 0.02)
DAT <- cbind(a, b, c)
EXPR <- expression(a * b^sin(c))
RES <- propagate(EXPR, DAT)
stochContr(RES)
}   

\keyword{matrix}
\keyword{multivariate}
\keyword{algebra}

