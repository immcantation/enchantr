**countSequencesPerCell** - *countSequencesPerCell*

Description
--------------------

`countSequencesPerCell` counts the number of sequences of each isotype in each sample's cell


Usage
--------------------
```
countSequencesPerCell(
db,
sample_id = "sample_id",
cell_id = "cell_id",
locus = "locus",
c_call = "c_call"
)
```

Arguments
-------------------

db
:   data.frame with AIRR-format style columns.

sample_id
:   column in `db` containing sample identifiers

cell_id
:   column in `db` containing cell identifiers

locus
:   column in `db` containing locus assignments

c_call
:   column in `db` containing constant region assignments




Value
-------------------

A data.frame with cell counts



Examples
-------------------

```R
db <- Example10x

```

**Error in eval(expr, envir, enclos)**: object 'Example10x' not found
```R
db[['sample_id']] <- 'example_sample'

```

**Error in db[["sample_id"]] <- "example_sample"**: object 'db' not found
```R
countSequencesPerCell(db[1:10,])
```

**Error in group_by(., sample_id, cell_id, locus)**: object 'db' not found

See also
-------------------

See also [plotSequencesPerCell](plotSequencesPerCell.md).






