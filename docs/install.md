# Download and Installation

Download
-------------------------------------------------------------------------------

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

`enchantr` is currently not available from CRAN and must be installed from the
GitHub repository directly by first cloning it from [https://github.com/immcantation/enchantr](https://github.com/immcantation/enchantr)

Then building it using the following R commands from the package root:

```R
install.packages(c("devtools", "roxygen2"))
library(devtools)
install_deps(dependencies=TRUE)
document()
install()
```

Alternatively, you can install the package directly from the GitHub repository, but this
will not build the documentation:

```R
library(devtools)
install_github("immcantation/enchantr@master")
```
