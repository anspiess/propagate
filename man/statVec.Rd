\name{statVec}
\alias{statVec}
\encoding{latin1}

\title{Transform an input vector into one with defined mean and standard deviation}

\description{
Transforms an input vector into one with defined \eqn{\mu} and \eqn{\sigma} by using a scaled-and-shifted Z-transformation.
}

\usage{
statVec(x, mean, sd)
}

\arguments{
  \item{x}{the input vector to be transformed.}
  \item{mean}{the desired mean of the created vector.}   
  \item{sd}{the desired standard deviation of the created vector.}   
}

\details{
Calculates vector \eqn{V} using a Z-transformation of the input vector \eqn{X} and subsequent scaling by \code{sd} and shifting by \code{mean}:
\deqn{V = \frac{X - \mu_X}{\sigma_X} \cdot \rm{sd} + \rm{mean}}
}

\value{
A vector with defined \eqn{\mu} and \eqn{\sigma}.
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
## Create a 10-sized vector with mean = 10 and s.d. = 1.
x <- rnorm(10, 5, 2)
mean(x) ## => mean is not 5!
sd(x) ## => s.d. is not 2!

z <- statVec(x, 5, 2)
mean(z) ## => mean is 5!
sd(z) ## => s.d. is 2!
}   

\keyword{matrix}
\keyword{multivariate}
\keyword{algebra}

