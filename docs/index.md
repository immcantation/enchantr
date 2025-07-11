# 
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

**IMPORTANT!** 
enchantR has moved to https://github.com/immcantation/enchantr

To update Git configuration settings use:

```
   git config user.email "your-gh-user@email.com"
   git config user.name "your-gh-user-name"
   git remote set-url origin git@github.com:immcantation/enchantr.git
```

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

If you need help or have any questions, please contact the [Immcantation Group](mailto:immcantation@googlegroups.com).

If you have discovered a bug or have a feature request, you can open an issue using the [issue tracker](https://github.com/immcantation/enchantr/issues).

To receive alerts about Immcantation releases, news, events, and tutorials, join the [Immcantation News](https://groups.google.com/g/immcantation-news) Google Group. [Membership settings](https://groups.google.com/g/immcantation-news/membership) can be adjusted to change the frequency of email updates.


## Dependencies

**Depends:** FALSE  
**Imports:** airr, alakazam, bookdown, ComplexHeatmap, data.table, doParallel, dowser, dplyr, DT, foreach, ggraph, ggplot2, gridExtra, igraph, iterators, knitr, plotly, RColorBrewer, R.utils, reshape2, rlang, rmarkdown, scales, scoper, shazam, stringi, stringr, tidyr, xfun  
**Suggests:** optparse, testthat


## Authors

[Susanna Marquez](mailto:susanna.marquez@yale.edu) (aut, cre)  
[Edel Aron](mailto:edel.aron@yale.edu) (aut)  
[Gisela Gabernet](mailto:gisela.gabernet@yale.edu) (aut)  
[Kenneth Hoehn](mailto:kenneth.hoehn@yale.edu) (aut)  
[Cole Jensen](mailto:cole.jensen@yale.edu) (aut)  
[Edward Lee](mailto:edward.lee@yale.edu) (aut)  
[Noah Yann Lee](mailto:noah.yann.lee@yale.edu) (aut)  
[Hailong Meng](mailto:hailong.meng@yale.edu) (aut)  
[Steven Kleinstein](mailto:steven.kleinstein@yale.edu) (aut, cph)


## License

AGPL-3
