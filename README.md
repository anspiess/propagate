## propagate: Propagation of uncertainty using higher-order Taylor expansion and Monte Carlo simulation (R-package)

*propagate* is an R package that can conduct uncertainty propagation based on first- and second-order Taylor expansion as well as Monte Carlo simulation. It also houses functionality to estimate confidence and prediction intervals for nonlinear models, create large correlation matrices and automatic distribution fitting.

#### Installation
You can install the latest development version of the code using the [devtools](https://cran.r-project.org/package=devtools) R package.

```R
# Install devtools, if you haven't already.
install.packages("devtools")
library(devtools)
install_github("anspiess/propagate")
source("https://install-github.me/anspiess/propagate")
```