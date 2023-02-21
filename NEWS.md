Version 0.0.6.999: February 21, 2023
-------------------------------------------------------------------------------

Bug fixes:

+ Fixed a bug in collapse duplicates report, where it was trying to select columns
  that sometimes did not exist.

Breaks:

+ Renamed param `skip_overlap` to `skip_convergence` in the Define Clones report.

General:

+ Added dependencies: ComplexHeatmap, plotly