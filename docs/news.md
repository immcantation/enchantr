Version 0.1.3: June 22, 2023
-------------------------------------------------------------------------------

General:

+ Find threshold report. Parametrized distance and threshold subsampling.
+ Define clones project. Added example data, available from url. Made `singlecell`
  optional, with default value `NULL`. `species` can be set to `auto` to auto
  detect from the field `species`.

Bug fixes:

+ The convergence UpSet was not showing in the report because it was missing
  a grid.draw statement.
  
+ Updated the logic in the find threshold report to show the distance to 
  nearest plot only if it is not `NULL`. 
  
+ Define clones project. Run chucks that require the `tissue` field only if
  it exists.
  
+ Updated figures widths and heights.

Version 0.1.2: June 2, 2023
-------------------------------------------------------------------------------

General:

+ Lineage report. Added parameters: `tips`, `traits`. All trees are saved to a .pdf
  and exported in graphml format.
  
+ Find threshold report. Added parameters: `findthreshold_method`, 
  `findthreshold_model`, `findthreshold_edge`, `findthreshold_cutoff`, 
   `findthreshold_spc`.
  
+ Define clones report. Limit the size of UpSet plots to 31 sets. 31 is the
  largest size allowed by `ComplexHeatmap::make_comb_mat`.
  
+ Create single page, self contained reports.

Bug fixes:

+ Fixed bug in collapse duplicates that caused the function to fail when all 
  values in d_call were NA.
  
+ In the define clones report. Added conditional chunk evaluation to chunks of 
  code that don't need to be executed if only light chain data is available. 
  
+ In the find threshold report. Skip distance plots and findThreshold if all
  distance values are NA and show a message.
  

Version 0.1.1: March 8, 2023
-------------------------------------------------------------------------------

Bug fixes:

+ Added the field `cell_id` to the set of selected columns used to subset the 
  input repertoire in the find threshold report.

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
