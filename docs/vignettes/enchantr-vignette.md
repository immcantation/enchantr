# enchantR reports

The 'enchantr' package provides template parametrized reports to analize AIRR-seq
data with the Immcantation framework. The different reports can be used from the 
command line, to be incorporated into complex workflows, or can be used 
from RStudio to create project templates.


## Quick start

### Command line

```
Rscript -e "enchantr:::enchantr_report('report_name', \
  report_params=list('input'='report_input', 'outdir'=getwd()))"
```

### RStudio

Got to `File > New Project... > New Directory...` and select an Immcantation project.




