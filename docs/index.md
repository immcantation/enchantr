[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)


enchantR
-------------------------------------------------------------------------------

enchantR is part of the [Immcantation](http://immcantation.readthedocs.io) 
analysis framework for Adaptive Immune Receptor Repertoire sequencing 
(AIRR-seq). Immcantation provides a set of tools to investigate lymphocyte 
receptor clonal lineages, diversity, gene usage, and other repertoire level 
properties, with a focus on high-throughput immunoglobulin (Ig) sequencing.

enchantR:

1. Provides additional quality control and reporting functions for other R
   packages in the Immcantation framework. 
2. Provides parametrized `.Rmd` reports for standard analysis.
3. Is part of the magic behind [nf-co.re/airrflow](https://nf-co.re/airrflow),
   a start-to-finish Nextflow pipeline to analyze AIRR-seq data.


Contact
-------------------------------------------------------------------------------

For help and questions, please contact the [Immcantation Group](mailto:immcantation@googlegroups.com)
or use the [issue tracker](https://bitbucket.org/kleinstein/enchantr/issues?status=new&status=open).


# Dependencies

**Depends:** FALSE  
**Imports:** airr, alakazam, bookdown, ComplexHeatmap, data.table, doParallel, dowser, dplyr, DT, foreach, ggraph, ggplot2, gridExtra, igraph, knitr, plotly, RColorBrewer, reshape2, rmarkdown, scales, scoper, shazam, stringi, stringr, tidyr  
**Suggests:** optparse, testthat


# Authors

[Susanna Marquez](mailto:susanna.marquez@yale.edu) (aut, cre)  
[Steven Kleinstein](mailto:steven.kleinstein@yale.edu) (aut, cph)


# License

AGPL-3
