\name{makeDat}
\alias{makeDat}

\title{Create a dataframe from the variables defined in an expression}

\description{
Creates a dataframe from the variables defined in an expression by \code{\link{cbind}}ing the corresponding data found in the workspace. This is a convenience function for creating a dataframe to be passed to \code{\link{propagate}}, when starting with data which was simulated from distributions, i.e. when \code{type = "sim"}. Will throw an error if a variable is defined in the expression but is not available from the workspace.
}

\usage{
makeDat(expr)
}

\arguments{
  \item{expr}{an expression to be use for \code{\link{propagate}}.}  
}

\value{
A dataframe containing the data defined in \code{expr} in columns.
}
     
\author{
Andrej-Nikolai Spiess
}  

\examples{
## Simulating from uniform
## and normal distribution,
## run 'propagate'.
EXPR1 <- expression(a + b^c)
a <- rnorm(100000, 12, 1)
b <- rnorm(100000, 5, 0.1)
c <- runif(100000, 6, 7)

DAT1 <- makeDat(EXPR1)
propagate(EXPR1, DAT1, type = "sim", cov = FALSE)
}   

\keyword{algebra}
\keyword{univariate}

