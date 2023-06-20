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
