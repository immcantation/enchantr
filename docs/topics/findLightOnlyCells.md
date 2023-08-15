**findLightOnlyCells** - *TODO: example
findLightOnlyCells*

Description
--------------------

`findLightOnlyCells` identifies cells with only light chains.


Usage
--------------------
```
findLightOnlyCells(
db,
sample_id = "sample_id",
cell_id = "cell_id",
locus = "locus",
fields = NULL
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

fields
:   Columns in `db`, in addition to `sample_id`,
that should be used to group sequences to be
analyzed independently.




Value
-------------------

The input `db` with an additional column named `light_only_cell`
with values TRUE/FALSE.









