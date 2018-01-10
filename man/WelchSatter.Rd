\name{WelchSatter}
\alias{WelchSatter}
\encoding{latin1}

\title{Welch-Satterthwaite approximation to the 'effective degrees of freedom'}

\description{
Calculates the Welch-Satterthwaite approximation to the 'effective degrees of freedom' by using the samples' uncertainties and degrees of freedoms, as described in Welch (1947) and Satterthwaite (1946). External sensitivity coefficients can be supplied optionally.
}

\usage{
WelchSatter(ui, ci = NULL, df = NULL, dftot = NULL, uc = NULL, alpha = 0.05)
}

\arguments{
  \item{ui}{the uncertainties \eqn{u_i} for each variable \eqn{x_i}.}
  \item{ci}{the sensitivity coefficients \eqn{c_i = \partial y/\partial x_i}.}   
  \item{df}{the degrees of freedom for the samples, \eqn{\nu_i}.}
  \item{dftot}{an optional known total degrees of freedom for the system, \eqn{\nu_{\mathrm{tot}}}. Overrides the internal calculation of \eqn{\nu_{\mathrm{ws}}}.} 
  \item{uc}{the combined uncertainty, u(y).}
  \item{alpha}{the significance level for the t-statistic. See 'Details'.}
}

\details{
\deqn{\nu_{\rm{eff}} \approx \frac{u(y)^4}{\sum_{i = 1}^n \frac{(c_iu_i)^4}{\nu_i}},  \quad k = t(1 - (\alpha/2), \nu_{\rm{eff}}), \quad u_{\rm{exp}} = ku(y)}
}

\value{
A list with the following items:\cr
\item{ws.df}{the 'effective degrees of freedom'.}
\item{k}{the coverage factor for calculating the expanded uncertainty.}
\item{u.exp}{the expanded uncertainty \eqn{u_{\rm{exp}}}.}
}

\author{
Andrej-Nikolai Spiess
}   

\references{
An Approximate Distribution of Estimates of Variance Components.\cr
Satterthwaite FE.\cr
\emph{Biometrics Bulletin} (1946), \bold{2}: 110-114.

The generalization of "Student's" problem when several different population variances are involved.\cr
Welch BL.\cr
\emph{Biometrika} (1947), \bold{34}: 28-35.
}

\examples{
## Taken from GUM H.1.6, 4).
WelchSatter(ui = c(25, 9.7, 2.9, 16.6), df = c(18, 25.6, 50, 2), uc = 32, alpha = 0.01)
}   

\keyword{matrix}
\keyword{multivariate}
\keyword{algebra}

