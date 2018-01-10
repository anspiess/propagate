\name{interval}
\alias{interval}
\encoding{latin1}

\title{Uncertainty propagation based on interval arithmetics}

\description{
Calculates the uncertainty of a model by using interval arithmetics based on a "combinatorial sequence grid evaluation" approach, thereby avoiding the classical dependency problem that inflates the result interval.
}

\usage{
interval(df, expr, seq = 10, plot = TRUE)
}

\arguments{
  \item{df}{a 2-row dataframe/matrix with lower border values \eqn{A_i} in the first row and upper border values \eqn{B_i} in the second row. Column names must correspond to the variable names in \code{expr}.}
  \item{expr}{an expression, such as \code{expression(x/y)}.}   
  \item{seq}{the sequence length from \eqn{A_i} to \eqn{B_i} in \eqn{[A_i, B_i]}.}   
  \item{plot}{logical. If \code{TRUE}, plots the evaluations and min/max values as blue border lines.}
}

\details{
For two variables \eqn{{\color{red}x}, {\color{blue}y}} with intervals \eqn{[{\color{red}x_1}, {\color{red}x_2}]} and \eqn{[{\color{blue}y_1}, {\color{blue}y_2}]}, the four basic arithmetic operations \eqn{\langle \mathrm{op} \rangle \in \{+, -, \cdot, /\}} are
\deqn{[{\color{red}x_1}, {\color{red}x_2}] \,\langle\!\mathrm{op}\!\rangle\, [{\color{blue}y_1}, {\color{blue}y_2}] = }
\deqn{\left[ \min({\color{red}x_1} {\langle\!\mathrm{op}\!\rangle} {\color{blue}y_1}, {\color{red}x_1} \langle\!\mathrm{op}\!\rangle {\color{blue}y_2}, {\color{red}x_2} \langle\!\mathrm{op}\!\rangle {\color{blue}y_1}, {\color{red}x_2} \langle\!\mathrm{op}\!\rangle {\color{blue}y_2}),  \max({\color{red}x_1} {\langle\!\mathrm{op}\!\rangle} {\color{blue}y_1}, {\color{red}x_1} {\langle\!\mathrm{op}\!\rangle} {\color{blue}y_2}, {\color{red}x_2} {\langle\!\mathrm{op}\!\rangle} {\color{blue}y_1}, {\color{red}x_2} {\langle\!\mathrm{op}\!\rangle} {\color{blue}y_2})\right]
}
So for a function \eqn{f([{\color{red}x_1}, {\color{red}x_2}], [{\color{blue}y_1}, {\color{blue}y_2}], [{\color{green}z_1}, {\color{green}z_2}], ...)} with \eqn{k} variables, we have to create all combinations \eqn{C_i = {\{\{{\color{red}x_1}, {\color{red}x_2}\}, \{{\color{blue}y_1}, {\color{blue}y2}\}, \{{\color{green}z_1}, {\color{green}z_2}\}, ...\} \choose k}}, evaluate their function values \eqn{R_i = f(C_i)} and select \eqn{R = [\min R_i, \max R_i]}.\cr 
The so-called \emph{dependency problem} is a major obstacle to the application of interval arithmetic and arises when the same variable exists in several terms of a complicated and often nonlinear function. In these cases, over-estimation can cover a range that is significantly larger, i.e. \eqn{\min R_i \ll \min f(x, y, z, ...) , \max R_i \gg \max f(x, y, z, ...)}. For an example, see \url{http://en.wikipedia.org/w/index.php?title=Interval_arithmetic} under "Dependency problem". A partial solution to this problem is to refine \eqn{R_i} by dividing \eqn{[{\color{red}x_1}, {\color{red}x_2}]} into \eqn{i} smaller subranges to obtain sequence \eqn{({\color{red}x_1}, x_{1.1}, x_{1.2}, x_{1.i}, {\color{red}x_2})}. Again, all combinations are evaluated as described above, resulting in a larger number of \eqn{R_i} in which \eqn{\min R_i} and \eqn{\max R_i} may be closer to \eqn{\min f(x, y, z, ...)} and \eqn{\max f(x, y, z, ...)}, respectively. This is the "combinatorial sequence grid evaluation" approach which works quite well in scenarios where monotonicity changes direction (see 'Examples'), obviating the need to create multivariate derivatives (Hessians) or use some multivariate minimization algorithm.\cr
If intervals are of type \eqn{[{\color{red}x_1} < 0, {\color{red}x_2} > 0]}, a zero is included into the middle of the sequence to avoid wrong results in case of even powers, i.e. \eqn{[-1, 1]^2 = [-1, 1][-1, 1] = [-1, 1]} when actually the right interval is \eqn{[0, 1]}, see \code{curve(x^2, -1, 1)}.
}

\value{
A 2-element vector with the resulting interval and an (optional) plot of all evaluations.
}

\author{
Andrej-Nikolai Spiess
}   

\references{
\bold{Wikipedia entry is quite good, especially the section
on the 'dependency problem':}\cr
\url{http://en.wikipedia.org/wiki/Interval_arithmetic}

\bold{Comparison to Monte Carlo and error propagation:}\cr
Interval Arithmetic in Power Flow Analysis.\cr
Wang Z & Alvarado FL.\cr
Power Industry Computer Application Conference (1991): 156-162.

\bold{Computer implementation}\cr
Interval arithmetic: From principles to implementation.\cr
Hickey T, Ju Q, Van Emden MH.\cr
\emph{JACM} (2001), \bold{48}: 1038-1068.

Complete Interval Arithmetic and its Implementation on the Computer.\cr
Kulisch UW.\cr
In: Numerical Validation in Current Hardware Architectures.
Lecture Notes in Computer Science \bold{5492} (2009): 7-26. 
}

\examples{
## Example 1: even squaring of negative interval.
EXPR1 <- expression(x^2)
DAT1 <- data.frame(x = c(-1, 1))
interval(DAT1, EXPR1)

## Example 2: A complicated nonlinear model.
## Reduce sequence length to 2 => original interval
## for quicker evaluation.
EXPR2 <- expression(C * sqrt((520 * H * P)/(M *(t + 460))))
H <- c(64, 65)
M <- c(16, 16.2)
P <- c(361, 365)
t <- c(165, 170)
C <- c(38.4, 38.5)
DAT2 <- makeDat(EXPR2)
interval(DAT2, EXPR2, seq = 2)

## Example 3: Body Mass Index taken from
## http://en.wikipedia.org/w/index.php?title=Interval_arithmetic
EXPR3 <- expression(m/h^2)
m <- c(79.5, 80.5)
h <- c(1.795, 1.805)
DAT3 <- makeDat(EXPR3)
interval(DAT3, EXPR3)

## Example 4: Linear model.
EXPR4 <- expression(a * x + b)
a <- c(1, 2)
b <- c(5, 7)
x <- c(2, 3)
DAT4 <- makeDat(EXPR4)
interval(DAT4, EXPR4)

## Example 5: Overestimation from dependency problem.
# Original interval with seq = 2 => [1, 7]
EXPR5 <- expression(x^2 - x + 1)
x <- c(-2, 1)
DAT5 <- makeDat(EXPR5)
interval(DAT5, EXPR5, seq = 2)

# Refine with large sequence => [0.75, 7]
interval(DAT5, EXPR5, seq = 100)
# Tallies with curve function.
curve(x^2 - x + 1, -2, 1)

## Example 6: Underestimation from dependency problem.
# Original interval with seq = 2 => [0, 0]
EXPR6 <- expression(x - x^2)
x <- c(0, 1)
DAT6 <- makeDat(EXPR6)
interval(DAT6, EXPR6, seq = 2)

# Refine with large sequence => [0, 0.25]
interval(DAT6, EXPR6, seq = 100)
# Tallies with curve function.
curve(x - x^2, 0, 1)
}   

\keyword{matrix}
\keyword{multivariate}
\keyword{algebra}

