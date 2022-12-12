Download
-------------------------------------------------------------------------------

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

`enchantr` is current not available from CRAN and must be installed from the
bitbucket repo directly by first cloning the bitbucket repository:

`https://bitbucket.org/kleinstein/enchantr <https://bitbucket.org/kleinstein/enchantr>`_

Then build using the following R commands from the package root::

    install.packages(c("devtools", "roxygen2"))
    library(devtools)
    install_deps(dependencies=T)
    document()
    install()

Alternatively, you can install directly form the bitbucket repository, but this
will not build the documentation::

    library(devtools)
    install_bitbucket("kleinstein/enchantr@master")