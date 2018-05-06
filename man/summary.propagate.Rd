\name{summary.propagate}
\alias{summary.propagate}

\title{Summary function for 'propagate' objects}

\description{
Provides a printed summary of the results obtained by \code{\link{propagate}}, such as statistics of the first/second-order uncertainty propagation, Monte Carlo simulation, the covariance matrix, symbolic as well as evaluated versions of the Gradient ("sensitivity") and Hessian matrices, relative contributions, the coverage factor and the Welch-Satterthwaite degrees of freedom. If \code{do.sim = TRUE} was set in \code{propagate}, skewness/kurtosis and Shapiro-Wilks/Kolmogorov-Smirnov tests for normality are calculated on the Monte-Carlo evaluations. 
}

\usage{
\method{summary}{propagate}(object, ...)
}

\arguments{
  \item{object}{an object returned from \code{\link{propagate}}.}     
  \item{...}{other parameters for future methods.}       
}

\value{
A printed output with the items listed in 'Description'.
}

\details{
Calculates the "sensitivity"" \eqn{S_i} of each variable \eqn{x_i} to the propagated uncertainty, as defined in the \emph{Expression of the Uncertainty of Measurement in Calibration, Eqn 4.2, page 9} (see 'References'):\cr
\deqn{S_i = \mathrm{eval}\left(\frac{\partial f}{\partial x_i}\right)}\cr 
The "contribution" matrix is then \eqn{\mathbf{C} = \mathbf{SS}^T\mathbf{\Sigma}}, where \eqn{\mathbf{\Sigma}} is the covariance matrix. In the implementation here, the "relative contribution" matrix \eqn{\mathbf{C}_{\mathrm{rel}}} is rescaled to sum up to 1.
}

\references{
Expression of the Uncertainty of Measurement in Calibration.\cr
European Cooperation for Accreditation (EA-4/02), 1999.
}  
  
\author{
Andrej-Nikolai Spiess
}
        
\examples{
EXPR1 <- expression(x^2 * sin(y))
x <- c(5, 0.01)
y <- c(1, 0.01)
DF1 <- cbind(x, y)
RES1 <- propagate(expr = EXPR1, data = DF1, type = "stat", 
                  do.sim = TRUE, verbose = TRUE, nsim = 100000)
summary(RES1)
}

\keyword{models}
\keyword{nonlinear}
