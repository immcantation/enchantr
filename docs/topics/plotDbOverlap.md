**plotDbOverlap** - *Plot heatmap of one or two db features overlap*

Description
--------------------

Plot a matrix to visualize the the number of objects (e.g. clones and sequences) shared 
between groups.


Usage
--------------------
```
plotDbOverlap(
db,
group = "sample",
features = c("clone_id", "sequence_alignment"),
heatmap_colors = c("white", "orange", "grey80"),
print_zero = FALSE,
long_x_angle = 90,
title = NULL,
xlab = NULL,
ylab = NULL,
plot_order = NULL,
silent = F,
similarity = c("min", "jaccard"),
na.rm = FALSE,
identity = c("exact", "ambiguous", "ham_nt", "ham_aa"),
threshold = 0,
geom_text_size = 3
)
```

Arguments
-------------------

db
:   Changeo db data.frame

group
:   Vector with column names for grouping. Overlap will
be calculated across the groups. e.g `SAMPLE`

features
:   Column name of the feature column(s) shared across group. 
e.g. `CLONE`, `SEQUENCE_INPUT`

heatmap_colors
:   Vector of colors representing low and high values on 
the heatmap and the diagonal. Default is 
c("white","orange", "grey80")

print_zero
:   Show labels on zero overlap cells or not.

long_x_angle
:   Angle to rotate x axis labels when any of the labels is longer
than 6 characters

title
:   A string that will be used for the heatmap title. If `NULL`,
`group` and `featureS` will be used.

xlab
:   Text to be used as x axis title

ylab
:   Text to be used as y axis title

plot_order
:   A vector to reorder the grouped columns. Can contain either 
the specific names of the grouped columns (e.g. 
`c("Donor1","Donor2")` ); the names 
of columns in `db` (e.g. `"DONOR"`, which
has factor levels "Donor1" and "Donor2" or `NULL` 
if ordering is not relevant.

silent
:   If T, the plot will not be printed, it will be found in the
returned list.

similarity
:   vector of the same length as `features`. For each
`feature`, method used to quantify the overlap. 
"min" will use the number of shared features over the 
number of features in the smaller set. "jaccard" will 
use the jaccard index, expressed as a percent, 
and defined as the intersection of features over the 
union. Can also be one of the possible values in the last grouping
column. For example to make the percent always relative to
the memory samples, use `group="SAMPLE", "SORT")` and 
`similarity=c("memory")`

na.rm
:   logical. If TRUE, NA values will be removed and not 
considered

identity
:   Vector of the same length as `features` spcifying 
how to establish identity. For each `features`, 
compare the exact value of the features (`identity='exact'`),
allow ambiguous characters in DNA sequences (`identity='ambiguous'`) 
using the function `seqEqual`, or use hamming distance and
a threshold (`identity='ham'`). See `details`.

threshold
:   Identity threshold to be used when `identity='ham'`.

geom_text_size
:   Plot text size




Value
-------------------

A list with the plot object and a data.frame with the values.


Details
-------------------

Can be used to visualize the number of clones and/or sequences shared between 
compartments, which are potentially antigen-specific

**Note**: When using `exact=FALSE` and `distance="min"` the overlap 
will be calculated using the number of shared features from the smallest set. In case 
of ties, the smallest number of shared features.



Examples
-------------------

```R
### Not run:
data("ExampleDb", package="alakazam")
# db <- ExampleDb
# 
# ## Plot the number of sequences that overlap across samples
# overlap <- plotDbOverlap(db, group="sample_id", features="clone_id", identity="exact", similarity="jaccard")
# overlap <- plotDbOverlap(db, group="sample_id", features=c("clone_id","junction"))
# 
# ## The returned plot can be modified
# ## To edit the axis labels. the title and change the color scale and change
# ## the theme
# overlap$p + xlab("Sequence") + ylab("Clone") + ggtitle("New title") +
#  scale_fill_gradient(low="white", high="orange", na.value="black") + theme_enchantr()
```








