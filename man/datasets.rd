\name{datasets}
\alias{H.2}
\alias{H.3}
\alias{H.4}
\encoding{latin1}

\title{Datasets from the GUM "Guide to the expression of uncertainties in measurement" (2008)}

\description{
Several datasets found in "Annex H" of the GUM that are used in illustrating the different approaches to quantifying measurement uncertainty.
}

\usage{
H.2
H.3
H.4
}

\details{
\bold{H.2: Simultaneous resistance and reactance measurement, Table H.2}\cr
This example demonstrates the treatment of multiple measurands or output quantities determined simultaneously in the same measurement and the correlation of their estimates. It considers only the random variations of the observations; in actual practice, the uncertainties of corrections for systematic effects would also contribute to the uncertainty of the measurement results. The data are analysed in two different ways with each yielding essentially the same numerical values.\cr 
H.2.1  The measurement problem:\cr 
The resistance \emph{R} and the reactance \emph{X} of a circuit element are determined by measuring the amplitude \emph{V} of a sinusoidally-alternating  potential  difference  across  its  terminals,  the  amplitude  \emph{I}  of  the  alternating  current passing through it, and the phase-shift angle \eqn{\phi} of the alternating potential difference relative to the alternating current. Thus the three input quantities are \emph{V}, \emph{I}, and \eqn{\phi} and the three output quantities -the measurands- are the three impedance components \emph{R}, \emph{X}, and \emph{Z}. Since \eqn{Z^2 = R^2 + X^2}, there are only two independent output 
quantities.\cr
H.2.2  Mathematical model and data:\cr
The measurands are related to the input quantities by Ohm's law: 
\deqn{R = \frac{V}{I}\cos\phi;\quad X = \frac{V}{I}\sin\phi;\quad Z = \frac{V}{I} \qquad (\mathrm{H.7})}

\bold{H.3: Calibration of a thermometer, Table H.6}\cr
This example illustrates the use of the method of least squares to obtain a linear calibration curve and how the parameters of the fit, the intercept and slope, and their estimated variances and covariance, are used to 
obtain from the curve the value and standard uncertainty of a predicted correction.\cr 
H.3.1 The measurement problem:\cr
A thermometer is calibrated by comparing \emph{n} = 11 temperature readings \eqn{t_k} of the thermometer, each having negligible uncertainty, with corresponding known reference temperatures \eqn{t_{R,k}} in the temperature range 21°C to 27°C to obtain the corrections \eqn{b_k = t_{R,k} - t_k} to the readings. The \emph{measured} corrections \eqn{b_k} and \emph{measured} temperatures \eqn{t_k} are the input quantities of the evaluation. A linear calibration curve \deqn{b(t) = y_1 + y_2(t-t_0) \qquad (\mathrm{H.12})} is fitted to the measured corrections and temperatures by the method of least squares. The parameters \eqn{y_1} and \eqn{y_2}, which are respectively the intercept and slope of the calibration curve, are the two measurands or output 
quantities to be determined. The temperature \eqn{t_0} is a conveniently chosen exact reference temperature; it is not an independent parameter to be determined by the least-squares fit. Once \eqn{y_1} and \eqn{y_2} are found, along with their  estimated  variances  and  covariance, Equation (H.12)  can  be  used  to  predict  the  value  and  standard uncertainty of the correction to be applied to the thermometer for any value \eqn{t} of the temperature. 

\bold{H.4: Measurement of activity, Table H.7}\cr
This example is similar to example H.2, the simultaneous measurement of resistance and reactance, in that 
the data can be analysed in two different ways but each yields essentially the same numerical result. The first approach  illustrates  once  again  the  need  to  take  the  observed  correlations  between  input  quantities  into account.\cr 
H.4.1  The measurement problem:\cr 
The  unknown  radon (\eqn{{}^{222}\mathrm{Rn}}) activity concentration in a water sample is determined by liquid-scintillation counting against a radon-in-water standard sample having a known activity concentration. The unknown activity concentration is obtained by measuring three counting sources consisting of approximately 5g of water and 12g of organic emulsion scintillator in vials of volume 22ml:\cr
Source (a) a  \emph{standard}  consisting  of  a  mass \eqn{m_S} of  the  standard  solution  with  a  known  activity concentration;\cr
Source (b) a  matched  \emph{blank}  water  sample  containing  no  radioactive  material,  used  to  obtain  the background counting rate;\cr
Source (c) the \emph{sample} consisting of an aliquot of mass \eqn{m_x} with unknown activity concentration.\cr
Six cycles of measurement of the three counting sources are made in the order standard - blank - sample; and  each  dead-time-corrected counting  interval \eqn{T_0} for each source  during  all six cycles is 60 minutes. Although the background counting rate cannot be assumed to be constant over the entire counting interval (65 hours), it is assumed that the number of counts obtained for each blank may be used as representative of the background counting rate during the measurements of the standard and sample in the same cycle. The 
data are given in Table H.7, where\cr
\eqn{t_S, t_B, t_x} are  the  times  from  the  reference  time  \eqn{t} = 0  to  the  midpoint  of  the  dead-time-corrected counting intervals \eqn{T_0} = 60 min for the standard, blank, and sample vials, respectively; although \eqn{t_B} is given for completeness, it is not needed in the analysis;\cr
\eqn{C_S, C_B, C_x} are the number of counts recorded in the dead-time-corrected counting intervals \eqn{T_0} = 60 min for the standard, blank, and sample vials, respectively.\cr
The observed counts may be expressed as\cr
\deqn{C_S = C_B + \varepsilon A_S T_0 m_S e^{-\lambda t_S} \qquad (\mathrm{H.18a})}
\deqn{C_x = C_B + \varepsilon A_x T_0 m_x e^{-\lambda t_x} \qquad (\mathrm{H.18b})}

where\cr
\eqn{\varepsilon} is the liquid scintillation detection efficiency for \eqn{{}^{222}\mathrm{Rn}} for a given source composition, assumed to be independent of the activity level;\cr
\eqn{A_S} is the activity concentration of the standard at the reference time \eqn{t} = 0;\cr
\eqn{A_x} is the measurand and is defined as the unknown activity concentration of the sample at the reference time \eqn{t} = 0;\cr
\eqn{m_S} is the mass of the standard solution;\cr
\eqn{m_x} is the mass of the sample aliquot;\cr
\eqn{\lambda} is the decay constant for \eqn{{}^{222}\mathrm{Rn}}: \eqn{\lambda = (ln 2)/T_{1/2} = 1.25894 \cdot 10^{-4}\ \mathrm{min}^{-1} (T_{1/2} = 5505.8\ \mathrm{min})}.\cr
(...) This suggests combining Equations (H.18a) and (H.18b) to obtain the following expression for the unknown concentration in terms of the known quantities:\cr
\deqn{... = A_S \frac{m_S}{m_x}\frac{C_x - C_B}{C_S - C_B}e^{\lambda(t_x - t_S)} \qquad (\mathrm{H.19})}
where \eqn{(C_x - C_B)e^{\lambda t_x}} and \eqn{(C_S - C_B)e^{\lambda t_S}} are, respectively, the background-corrected counts of the sample and the standard  at the reference time \eqn{t} = 0 and for the time interval \eqn{T_0} = 60 min.
}

\author{
Andrej-Nikolai Spiess, taken mainly from the GUM 2008 manual.
}

\references{
Evaluation of measurement data - Guide to the expression of uncertainty in measurement.\cr
JCGM 100:2008 (GUM 1995 with minor corrections).\cr
\url{http://www.bipm.org/utils/common/documents/jcgm/JCGM_100_2008_E.pdf}.

Evaluation of measurement data - Supplement 1 to the Guide to the expression of uncertainty in measurement - Propagation of distributions using a Monte Carlo Method.\cr
JCGM 101:2008.\cr
\url{http://www.bipm.org/utils/common/documents/jcgm/JCGM_101_2008_E.pdf}.
}

\examples{
## See "Examples" in 'propagate'.
}

\keyword{datasets}

