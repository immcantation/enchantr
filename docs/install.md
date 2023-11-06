Download
-------------------------------------------------------------------------------

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

`enchantr` is currently not available from CRAN and must be installed from the
bitbucket repo directly by first cloning the bitbucket repository [https://bitbucket.org/kleinstein/enchantr](https://bitbucket.org/kleinstein/enchantr)

Then build using the following R commands from the package root:

```R
install.packages(c("devtools", "roxygen2"))
library(devtools)
install_deps(dependencies=TRUE)
document()
install()
```

Alternatively, you can install directly form the bitbucket repository, but this
will not build the documentation:

```R
library(devtools)
install_bitbucket("kleinstein/enchantr@master")
```
