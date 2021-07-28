#!/usr/bin/env Rscript
# 
# Validate --input. Check: 
# - filenames are .fasta, .fa or .tsv
# - column names are MiAIRR compliant
#
# Arguments:
#   --input   Tabulated data
#   --collapseby Names of the metadata columns to be used to collapse
#              duplicated sequences. Default: sample_id
#   --cloneby Names of the metadata columns that need to be used to group
#             files to identify clonal groups. Default: subject_id
#   --output  Output name. Default:validate_input
#   -h  Display help.
# Example: ./validate_input.R --input ../../test-datasets/metadata.tsv --collapseby sample_id --cloneby subject_id

# TODO: should validate metadata using AIRR schema, and properly tested in the
# framework, so that we incorporate any updates in the standard
# https://github.com/airr-community/airr-standards/blob/v1.3.1/specs/airr-schema.yaml
# https://github.com/airr-community/airr-standards/blob/v1.3.1/NCBI_implementation/mapping_MiAIRR_BioSample.tsv
# metadata data.frame
# miairr_mapping. Local version of https://github.com/airr-community/airr-standards/blob/v1.3.1/NCBI_implementation/mapping_MiAIRR_BioSample.tsv

# Imports
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("enchantr"))
suppressPackageStartupMessages(library("DT"))
suppressPackageStartupMessages(library("dplyr"))

# Define commmandline arguments
opt_list <- list(make_option(c("--input"), dest="input", default=NULL,
                             help="Input file."),
                 make_option(c("--collapseby"), dest="collapseby", default='sample_id',
                             help="Grouping fields to collapse duplicated sequences."),                    
                 make_option(c("--cloneby"), dest="cloneby", default='subject_id',
                             help="Grouping fields to identify clonally related sequences."),                 
                 make_option(c("--output"), dest="output", default="validate_input",
                             help="Output name."),
                 make_option(c("--miairr"), dest="miairr", default=NULL,
                             help="Output name."))
# Parse arguments
opt <- parse_args(OptionParser(option_list=opt_list))

tryCatch(
   report <- validate_input(opt, outdir=opt$outdir),
   error = function(e) {
      stop(safeError(e))
   }
)

report

