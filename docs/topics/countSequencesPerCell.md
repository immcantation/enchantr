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
data(Example10x, package="alakazam")
db <- Example10x
db[["sample_id"]] <- "example_sample"
countSequencesPerCell(db[1:10,])
```


```
# A tibble: 10 × 6
   sample_id    cell_id locus cell_num_sequences cell_num_isotypes cell_isotypes
   <chr>        <chr>   <chr>              <int>             <int> <chr>        
 1 example_sam… AAAGCA… IGH                    1                 1 IGHD         
 2 example_sam… ATTGGT… IGH                    1                 1 IGHM         
 3 example_sam… CATGGC… IGH                    1                 1 IGHM         
 4 example_sam… GAACAT… IGH                    1                 1 IGHD         
 5 example_sam… GACGTG… IGH                    1                 1 IGHM         
 6 example_sam… GCTGGG… IGH                    1                 1 IGHM         
 7 example_sam… GGTATT… IGH                    1                 1 IGHM         
 8 example_sam… TCTCTA… IGH                    1                 1 IGHM         
 9 example_sam… TTAGTT… IGH                    1                 1 IGHA1        
10 example_sam… TTTCCT… IGH                    1                 1 IGHM         

```



See also
-------------------

See also [plotSequencesPerCell](plotSequencesPerCell.md).






