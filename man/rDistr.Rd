\name{rDistr}
\alias{rDistr}

\title{Creating random samples from a variety of useful distributions}

\description{
These are random sample generators for 22 different continuous distributions which are not readily available as \code{\link{Distributions}} in \R. Some of them are implemented in other specialized packages (i.e. \code{rsn} in package 'sn' or \code{rtrapezoid} in package 'trapezoid'), but here they are collated in a way that makes them easily accessible for Monte Carlo-based uncertainty propagation.
}

\details{
Random samples can be drawn from the following distributions:\cr
1) Skewed-normal distribution: \code{propagate:::rsn(n, location = 0, scale = 1, shape = 0)}\cr
2) Generalized normal distribution: \code{propagate:::rgnorm(n, alpha = 1, xi = 1, kappa = -0.1)}\cr 
3) Scaled and shifted t-distribution: \code{propagate:::rst(n, mean = 0, sd = 1, df = 2)}\cr 
4) Gumbel distribution: \code{propagate:::rgumbel(n, location = 0, scale = 1)}\cr
5) Johnson SU distribution: \code{propagate:::rJSU(n, xi = 0, lambda = 1, gamma = 1, delta = 1)}\cr
6) Johnson SB distribution: \code{propagate:::rJSB(n, xi = 0, lambda = 1, gamma = 1, delta = 1)}\cr
7) 3P Weibull distribution: \code{propagate:::rweibull2(n, location = 0, shape = 1, scale = 1)}\cr
8) 4P Beta distribution: \code{propagate:::rbeta2(n, alpha1 = 1, alpha2 = 1, a = 0, b = 0)}\cr
9) Triangular distribution: \code{propagate:::rtriang(n, a = 0, b = 1, c = 0.5)}\cr
10) Trapezoidal distribution: \code{propagate:::rtrap(n, a = 0, b = 1, c = 2, d = 3)}\cr
11) Laplacian distribution: \code{propagate:::rlaplace(n, mean = 0, sigma = 1)}\cr
12) Arcsine distribution: \code{propagate:::rarcsin(n, a = 2, b = 1)}\cr
13) von Mises distribution: \code{propagate:::rmises(n, mu = 1, kappa = 3)}\cr
14) Curvilinear Trapezoidal distribution: \code{propagate:::rctrap(n, a = 0, b = 1, d = 0.1)}\cr
15) Generalized trapezoidal distribution:\cr
\code{propagate:::rgtrap(n, min = 0, mode1 = 1/3, mode2 = 2/3, max = 1, n1 = 2, n3 = 2, alpha = 1)}\cr
16) Inverse Gaussian distribution: \code{propagate:::rinvgauss(n, mean = 1, dispersion = 1)}
17) Generalized Extreme Value distribution: \code{propagate:::rgevd(n, loc = 0, scale = 1, shape = 0)}
with \code{n} = number of samples.\cr
18) Inverse Gamma distribution: \code{propagate:::rinvgamma(n, shape = 1, scale = 5)}\cr
19) Rayleigh distribution: \code{propagate:::rrayleigh(n, mu = 1, sigma = 1)}\cr
20) Burr distribution: \code{propagate:::rburr(n, k = 1)}\cr
21) Chi distribution: \code{propagate:::rchi(n, nu = 5)}\cr
22) Inverse Chi-Square distribution: \code{propagate:::rinvchisq(n, nu = 5)}\cr
23) Cosine distribution:  \code{propagate:::rcosine(n, mu = 5, sigma = 1)}\cr

1) - 12), 17) - 22) use the inverse cumulative distribution function as mapping functions for \code{\link{runif}} (\bold{Inverse Transform Method}):\cr
(1) \eqn{U \sim \mathcal{U}(0, 1)}\cr
(2) \eqn{Y = F^{-1}(U, \beta)}

16) uses binomial selection from a \eqn{\chi^2}-distribution.

13) - 15), 23) employ "Rejection Sampling" using a uniform envelope distribution (\bold{Acceptance Rejection Method}):\cr
(1) Find \eqn{F_{max} = \max(F([x_{min}, x_{max}], \beta)}\cr
(2) \eqn{U_{max} = 1/(x_{max} - x_{min})}\cr
(3) \eqn{A = F_{max}/U_{max}}\cr
(4) \eqn{U \sim \mathcal{U}(0, 1)}\cr
(5) \eqn{X \sim \mathcal{U}(x_{min}, x_{max})}\cr
(6) \eqn{Y \iff U \le A \cdot \mathcal{U}(X, x_{min}, x_{max})/F(X, \beta)}\cr

These four distributions are coded in a vectorized approach and are hence not much slower than implementations in C/C++ (0.2 - 0.5 sec for 100000 samples; 3 GHz Quadcore processor, 4 GByte RAM). The code for the random generators is in file "distr-samplers.R".
}

\value{
A vector with \code{n} samples from the corresponding distribution.
}
     
\author{
Andrej-Nikolai Spiess
}  

\seealso{
See also \code{\link{propagate}}, in which GUM 2008 Supplement 1 examples use these distributions.
}

\references{
\bold{Inverse CDFs were taken from:}\cr
"The Ultimate Univariate Probability Distribution Explorer"\cr
\url{http://blog.wolfram.com/data/uploads/2013/02/ProbabilityDistributionExplorer.zip}.

\bold{Rejection Sampling in R:}\cr
Rejection Sampling.\cr
\url{https://www.r-bloggers.com/rejection-sampling/}.

An example of rejection sampling.\cr
\url{http://www.mas.ncl.ac.uk/~ndjw1/teaching/sim/reject/circ.html}.

\bold{Rejection Sampling in general:}\cr
Non-uniform random variate generation.\cr
Devroye L.\cr
Springer-Verlag, New York (1986).

\bold{Distributions:}\cr
Continuous univariate distributions, Volume 1.\cr
Johnson NL, Kotz S and Balakrishnan N.\cr
\emph{Wiley Series in Probability and Statistics, 2.ed} (2004).
}

\examples{
\dontrun{
## First we create random samples from the
## von Mises distribution.
X <- propagate:::rmises(1000000, mu = 1, kappa = 2)

## then we fit all available distributions
## with 'fitDistr'.
fitDistr(X, nbin = 200)
## => von Mises wins! (lowest BIC)
}
}

\keyword{algebra}
\keyword{univariate}

