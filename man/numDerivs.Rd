\name{numDerivs}
\alias{numGrad}
\alias{numHess}
\encoding{latin1}

\title{Functions for creating Gradient and Hessian matrices by numerical differentiation (Richardson's method) of the partial derivatives}

\description{
These two functions create Gradient and Hessian matrices by Richardson's central finite difference method of the partial derivatives for any expression. 
}

\usage{
numGrad(expr, envir = .GlobalEnv)
numHess(expr, envir = .GlobalEnv)
}

\arguments{
  \item{expr}{an expression, such as \code{expression(x/y)}.}
  \item{envir}{the \code{\link{environment}} to evaluate in.}  
}

\details{
Calculates first- and second-order numerical approximation using Richardson's \bold{central difference formula}:\cr
\deqn{f'_i(x) \approx \frac{f(x_1, \ldots, x_i + d, \ldots, x_n) - f(x_1, \ldots, x_i - d, \ldots, x_n)}{2d}}
\deqn{f''_i(x) \approx \frac{f(x_1, \ldots, x_i + d, \ldots, x_n) - 2f(x_1, \ldots, x_n) + f(x_1, \ldots, x_i - d, \ldots, x_n)}{d^2}}
}

\value{The numeric Gradient/Hessian matrices.}

\note{The two functions are modified versions of the \code{genD} function in the 'numDeriv' package, but a bit more easy to handle because they use expressions and the function's \code{x} value must not be defined as splitted scalar values \code{x[1], x[2], ... x[n]} in the body of the function.}

\author{
Andrej-Nikolai Spiess
}   

\examples{
## Check for equality of symbolic  
## and numerical derivatives.
EXPR <- expression(2^x + sin(2 * y) - cos(z))
x <- 5
y <- 10
z <- 20

symGRAD <- evalDerivs(makeGrad(EXPR))
numGRAD <- numGrad(EXPR)
all.equal(symGRAD, numGRAD)

symHESS <- evalDerivs(makeHess(EXPR))
numHESS <- numHess(EXPR)
all.equal(symHESS, numHESS)
}   

\keyword{array}
\keyword{algebra}
\keyword{multivariate}

