**findSingleCellDuplicates** - *findSingleCellDuplicates*

Description
--------------------

Finds sequences that share the same identifier and sequence between groups.


Usage
--------------------
```
findSingleCellDuplicates(
db,
fields,
cell_id = "cell_id",
seq = "sequence_alignment",
sequence_id = "sequence_id"
)
```

Arguments
-------------------

db
:   data.frame with AIRR-format style columns.

fields
:   Columns in `db`, in addition to `sample_id`,
that should be used to group sequences to be 
analyzed independently.

cell_id
:   column in `db` containing cell identifiers

seq
:   column in `db` containing sequences to be compared

sequence_id
:   column in `db` containing sequence identifiers




Value
-------------------

A list with fields: 

+ `dups`:    a data.frame with the column `sc_duplicate`
with values TRUE/FALSE to indicate whether the 
the row corresponds to a duplicated entry.
+ `fields`:  a data.frame showing the input fields used
+ `cell_id`: column in `db` containing cell identifiers
+ `seq`:     column in `db` containing sequence data
+ `sequence_id` column in `db` containin sequence identifiers










