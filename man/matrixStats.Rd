\name{matrixStats}
\alias{colVarsC}
\alias{rowVarsC}
\encoding{latin1}

\title{Fast column- and row-wise versions of variance coded in C++}

\description{
These two functions are fast C++ versions for column- and row-wise \code{\link{var}}iance calculation on matrices/data.frames and are meant to substitute the classical \code{apply(mat, 1, var)} approach. 
}

\usage{
  colVarsC(x)
  rowVarsC(x)  
}

\arguments{
  \item{x}{a matrix or data.frame} 
}

\details{
They are coded in a way that they automatically remove \code{NA} values, so they behave like \code{na.rm = TRUE}.
}

\value{
A vector with the variance values.
}

\author{
  Andrej-Nikolai Spiess
}   

\examples{
## Speed comparison on large matrix.
## ~ 110x speed increase!
\dontrun{
MAT <- matrix(rnorm(10 * 500000), ncol = 10)
system.time(RES1 <- apply(MAT, 1, var))
system.time(RES2 <- rowVarsC(MAT))
all.equal(RES1, RES2)
}
}  
   
\keyword{univar}

