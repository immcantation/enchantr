**readInput** - *Load input repertertoires*

Description
--------------------

Load input repertertoires


Usage
--------------------
```
readInput(input, pattern = "pass.tsv", col_select = NULL)
```

Arguments
-------------------

input
:   Path to repertoire(s). Can be a directory or comma separated paths.
The target files can be repertoire files or file-of-files (
tabulated text files with one or more paths to repertoires in the first column).

pattern
:   If input is a directory, the pattern to select the input
repertoire files to be loaded.

col_select
:   Columns to select from the repertoire. Will be passed to airr::read_rearrangement











