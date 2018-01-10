\name{moments}
\alias{skewness}
\alias{kurtosis}
\encoding{latin1}

\title{Skewness and (excess) Kurtosis of a vector of values}

\description{
These functions calculate skewness and excess kurtosis of a vector of values. They were taken from the package 'moments'.
}

\usage{
skewness(x, na.rm = FALSE) 
kurtosis(x, na.rm = FALSE) 
}

\arguments{
  \item{x}{a numeric vector, matrix or data frame.}
  \item{na.rm}{logical. Should missing values be removed?}  
}

\details{
Skewness:
\deqn{\frac{\frac{1}{n} \sum_{i=1}^n (x_i-\overline{x})^3}{\left(\frac{1}{n} \sum_{i=1}^n (x_i-\overline{x})^2\right)^{3/2}}}

(excess) Kurtosis:
\deqn{\frac{\frac{1}{n} \sum_{i=1}^n (x_i - \overline{x})^4}{\left(\frac{1}{n} \sum_{i=1}^n (x_i - \overline{x})^2\right)^2} - 3}
}

\value{The skewness/kurtosis values.}

\author{
Andrej-Nikolai Spiess
}   

\examples{
X <- rnorm(100, 20, 2)
skewness(X)
kurtosis(X)
}   

\keyword{array}
\keyword{algebra}
\keyword{multivariate}

