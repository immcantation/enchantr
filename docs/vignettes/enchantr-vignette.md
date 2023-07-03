# enchantR reports

The 'enchantr' package provides template parametrized reports to analyze AIRR-seq
data with the [Immcantation](https://immcantation.readthedocs.io/en/stable/) framework. The different reports can be used from the command line, from RStudio to create project templates or incorporated into complex workflows.

## Quick start

### Command line

```
Rscript -e "enchantr:::enchantr_report('report_name', \
  report_params=list('input'='report_input', 'outdir'=getwd()))"
```

### RStudio

Go to `File > New Project... > New Directory...` and select an Immcantation project.
