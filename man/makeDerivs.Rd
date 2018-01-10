\name{makeDerivs}
\alias{makeGrad}
\alias{makeHess}
\alias{evalDerivs}
\encoding{latin1}

\title{Utility functions for creating Gradient- and Hessian-like matrices with symbolic derivatives and evaluating them in an environment}

\description{
These are three different utility functions that create matrices containing the symbolic partial derivatives of first (\code{makeGrad}) and second (\code{makeHess}) order and a function for evaluating these matrices in an environment. The evaluations of the symbolic derivatives are used within the \code{\link{propagate}} function to calculate the propagated uncertainty, but these functions may also be useful for other applications. 
}

\usage{
makeGrad(expr, order = NULL)
makeHess(expr, order = NULL)
evalDerivs(deriv, envir)
}

\arguments{
  \item{expr}{an expression, such as \code{expression(x/y)}.}
  \item{order}{order of creating partial derivatives, i.e. \code{c(2, 1)}. See 'Examples'.} 
  \item{deriv}{a matrix returned from \code{makeGrad} or \code{makeHess}.}
  \item{envir}{an environment to evaluate in. By default the workspace.}
}

\details{
Given a function \eqn{f(x_1, x_2, \ldots, x_n)}, the following matrices containing symbolic derivatives of \eqn{f} are returned:\cr

\bold{makeGrad:}\cr
\deqn{\nabla(f) = \left[\frac{\partial f}{\partial x_1}, \ldots, \frac{\partial f}{\partial x_n}\right]} 
\bold{makeHess:}\cr

\deqn{H(f) = \left[ \begin{array}{cccc} \frac{\partial^2 f}{\partial x_1^2} & \frac{\partial^2 f}{\partial x_1\,\partial x_2} & \cdots & \frac{\partial^2 f}{\partial x_1\,\partial x_n} \\ \frac{\partial^2 f}{\partial x_2\,\partial x_1} & \frac{\partial^2 f}{\partial x_2^2} & \cdots & \frac{\partial^2 f}{\partial x_2\,\partial x_n} \\ \vdots & \vdots & \ddots & \vdots \\ \frac{\partial^2 f}{\partial x_n\,\partial x_1} & \frac{\partial^2 f}{\partial x_n\,\partial x_2} & \cdots & \frac{\partial^2 f}{\partial x_n^2} \end{array} \right]
}
}

\value{The symbolic or evaluated Gradient/Hessian matrices.}


\author{
Andrej-Nikolai Spiess
}   

\references{
The Matrix Cookbook (Version November 2012).\cr
Petersen KB & Pedersen MS.\cr
\url{http://ais.informatik.uni-freiburg.de/teaching/ws12/mapping/pdf/imm3274.pdf}
}

\examples{
EXPR <- expression(a^b + sin(c))
ENVIR <- list(a = 2, b = 3, c = 4)

## First-order partial derivatives: Gradient.
GRAD <- makeGrad(EXPR) 

## This will evaluate the Gradient.
evalDerivs(GRAD, ENVIR)

## Second-order partial derivatives: Hessian.
HESS <- makeHess(EXPR) 

## This will evaluate the Hessian.
evalDerivs(HESS, ENVIR)

## Change derivatives order.
GRAD <- makeGrad(EXPR, order = c(2,1,3)) 
evalDerivs(GRAD, ENVIR)
}   

\keyword{array}
\keyword{algebra}
\keyword{multivariate}

