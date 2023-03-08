Version 0.1.1: March 8, 2023
-------------------------------------------------------------------------------

Bug fixes:

+ Added the field `cell_id` to the set of selected columns used to subset the input repertoire in the find threshold report.

General:

+ Use plotly in clone size plots.

+ Resize abundande plots.
 

Version 0.1.0: February 28, 2023
-------------------------------------------------------------------------------

Bug fixes:

+ Fixed a bug in collapse duplicates report, where it was trying to select columns
  that sometimes did not exist.

Breaks:

+ Renamed param `skip_overlap` to `skip_convergence` in the Define Clones report.

General:

+ Added dependencies: ComplexHeatmap, plotly