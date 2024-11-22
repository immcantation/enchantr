# Release Notes

Version 0.1.18: November 25, 2024
-------------------------------------------------------------------------------

General:

+ Fixed CollapseDuplicates to annotate isotype information in c_region or c_primer fields in the c_call field when c_call is NA.


Version 0.1.17: October 14, 2024
-------------------------------------------------------------------------------

General:

+ Fixed bug in single-cell QC report when database was empty.
+ Fixed bug in collapseDuplicates to provide as numeric field "duplicate_count" and "consensus_count".


Version 0.1.16: May 29, 2024
-------------------------------------------------------------------------------

General:

+ Fixed bug plotting mutation frequency.


Version 0.1.14: April 21, 2024
-------------------------------------------------------------------------------

General:

+ Fixed bug to run Dowser lineages with raxml
+ Added saving RDS object for Dowser formatted clones and trees

Version 0.1.13: March 12, 2024
-------------------------------------------------------------------------------

General:

+ Improved defineClones report and added mutation frequency calculation
+ Fixed bug that eliminated light chains from bulk or mixed sequencing data from
final repertoires.


Version 0.1.10: January 11, 2024
-------------------------------------------------------------------------------

General:

+ Added a copy of `fuse.js@6.4.6` to be able to render the reports when there is 
no internet access.


Version 0.1.5: September 15, 2023
-------------------------------------------------------------------------------

General:

+ Bumped versions to `dowser >= 1.2.0`, `ggplot2 >= 3.4.0` and `shazam >= 1.1.2`.
+ Updated some reports' text to be more concise and clear.
+ In the define clones report, added the option to focus the convergence analysis on 
  a particular V gene by setting `convergence_vgene`. Added the option to specify
  a threshold for the convergence analysis (`convergence_threshold`). Added a step
  to remove light chain only cells after `createGermlines`.
+ Updated the updated single cell qc report to remove cells with duplicated 
  sequences, not just the sequences.
  
+ Added additional parameters to the dowser lineages report: `num_fields`, `chain`, `cell`, 
  `heavy`, `collapse`, and `columns`. In the same report, added additional 
  validation steps to check that the fields passed to dowser functions have data.
  Updated the example input data and the default values of some parameters (
  `traits`: from `c_call` to `day`, `tips`: from `c_call` to `day`). Fixed the default value
  of `outname`, from `define-clones` to `dowser`.

Version 0.1.4: July 31, 2023
-------------------------------------------------------------------------------

General:

+ Aesthetic adjusments to figures and tables.

+ Added conditional evaluation logic to some chucks that don't need to be
  evaluated always.
  
Bug fixes:

+ In find threshold report, updated output path used to save the final tables
  to use the parameter `outname`. Changed the output file suffix `threshold_summary`
  to `threshold-summary`.  


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
