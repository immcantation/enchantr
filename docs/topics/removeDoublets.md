**removeDoublets** - *TODO: example
removeDoublets*

Description
--------------------

`removeDoublets` removes cells with multiple heavy chain sequences.


Usage
--------------------
```
removeDoublets(
db,
cell_id = "cell_id",
locus = "locus",
sequence_id = "sequence_id",
fields = NULL
)
```

Arguments
-------------------

db
:   data.frame with AIRR-format style columns.

cell_id
:   column in `db` containing cell identifiers

locus
:   column in `db` containing locus assignments

sequence_id
:   column in `db` containing locus assignments

fields
:   Columns in `db`, in addition to `sample_id`,
that should be used to group sequences to be
analyzed independently.




Value
-------------------

The input data.frame (`db`) with doublets removed.









