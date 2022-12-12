isHeavyChain
------------

**isHeavyChain**

Description
~~~~~~~~~~~

``isHeavyChain`` checks if a locus is a heavy or light chain locus.

Usage
~~~~~

::

   isHeavyChain(loci)

Arguments
~~~~~~~~~

loci
   a vector of valid loci
cell_id
   column in ``db`` containing cell identifiers
locus
   column in ``db`` containing locus assignments
sequence_id
   column in ``db`` containing locus assignments
fields
   Columns in ``db``, in addition to ``sample_id``, that should be used
   to group sequences to be analyzed independently.

Value
~~~~~

The input data.frame (``db``) with doublets removed.

Examples
~~~~~~~~

.. code:: r

   isHeavyChain(c("IGH","igh","TRA"))

::

   [1]  TRUE  TRUE FALSE
