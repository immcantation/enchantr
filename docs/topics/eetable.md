**eetable** - *eeplot: a convenience function to save tables*

Description
--------------------

This function takes a `data.frame`, and saves it as a .tsv file in
the `outdir` directory.


Usage
--------------------
```
eetable(df, caption = NULL, outdir = NULL, file = NULL, show_max = NULL)
```

Arguments
-------------------

df
:   data.frame

caption
:   caption for the table

outdir
:   directory where the output table be saved.

file
:   filename. If null, `p` will be used.

show_max
:   Number of lines to show. All data will be saved in the .tsv file.




Value
-------------------

It returns a list with a `DT::datatable` and the caption, 
updated with additional text with a link to the destination file.




See also
-------------------

See print_table_caption, the use of which is required for
correct placement of the caption.






