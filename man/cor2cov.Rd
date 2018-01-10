\name{cor2cov}
\alias{cor2cov}
\encoding{latin1}

\title{Converting a correlation matrix into a covariance matrix}

\description{
Converts a correlation matrix into a covariance matrix using variance information. It is therefore the opposite of \code{\link{cov2cor}}.
}

\usage{
cor2cov(C, var)
}

\arguments{
  \item{C}{a symmetric numeric correlation matrix \eqn{\mathbf{C}}.}
  \item{var}{a vector of variances \eqn{\sigma_n^2}.}   
}

\details{
Calculates the covariance matrix \eqn{\mathbf{\Sigma}} using a correlation matrix \eqn{\mathbf{C}} and outer products of the standard deviations \eqn{\sigma_n}:
\deqn{\mathbf{\Sigma} = \mathbf{C} \cdot \sigma_n \otimes \sigma_n}
}

\value{
The corresponding covariance matrix.
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
## Example in Annex H.2 from the GUM 2008 manual
## (see 'References'), simultaneous resistance
## and reactance measurement.
data(H.2)
attach(H.2)

## Original covariance matrix.
COV <- cov(H.2)
## extract variances
VAR <- diag(COV)

## cor2cov covariance matrix.
COV2 <- cor2cov(cor(H.2), VAR) 

## Equal to original covariance matrix.
all.equal(COV2, COV)
}   

\keyword{matrix}
\keyword{multivariate}
\keyword{algebra}

