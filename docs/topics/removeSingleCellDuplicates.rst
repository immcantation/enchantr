removeSingleCellDuplicates
--------------------------

**removeSingleCellDuplicates**

Description
~~~~~~~~~~~

``removeSingleCellDuplicates`` removes cells sharing the same identifier
and sequence between groups.

Usage
~~~~~

::

   removeSingleCellDuplicates(
   db,
   fields,
   cell_id = "cell_id",
   seq = "sequence_alignment",
   sequence_id = "sequence_id",
   mode = c("sequences", "cells")
   )

Arguments
~~~~~~~~~

db
   data.frame with AIRR-format style columns.
fields
   Columns in ``db``, in addition to ``sample_id``, that should be used
   to group sequences to be analyzed independently.
cell_id
   column in ``db`` containing cell identifiers
seq
   column in ``db`` containing sequences to be compared
sequence_id
   column in ``db`` containing sequence identifiers
mode
   Use ``sequences`` to remove duplicated sequences and ``cells`` to
   remove cells with duplicated sequences.

Value
~~~~~

A data.frame: a modified ``db`` without the duplicated sequences or
cells.
