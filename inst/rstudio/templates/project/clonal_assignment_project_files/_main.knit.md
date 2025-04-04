--- 
title: "Immcantation - enchantR"
subtitle: "Clonal Assignment"
author: ""
date: "Updated: Fri Apr  4 12:00:41 2025"
knit: enchantr::render_book
site: bookdown::bookdown_site
documentclass: book
bibliography: "references.bib"
biblio-style: apalike
link-citations: yes
description: "enchantr define clones"
output: enchantr::immcantation
params:
   input: 
      label: "`input`: Path to repertoires files"
      input: file
      value: "https://yale.box.com/shared/static/dajikxtvbc3c5yw1o2v69fjtx5bh58r2.tsv"
   imgt_db:
      label: "`imgt_db`: Path to IMGT reference  database."
      input: file
      value: "https://raw.githubusercontent.com/nf-core/test-datasets/airrflow/database-cache/imgtdb_base.zip"
   species:
      label: "`species`"
      input: text
      value: "human"        
   force:
      label: "`force`: if clone_id already exists, overwrite (force=TRUE) or skip clone assignment and render the report (force=FALSE)."
      input: logical        
      value: !r FALSE       
   cloneby:
      label: "`cloneby`"
      input: text
      value: "subject_id"    
   singlecell:
      label: "`singlecell`"
      input: text
      value: NULL          
   threshold:
      label: "`threshold`"
      input: numeric
      value: 1
   min_n:
      label: "`min_n` minimum number of observations to sample to run alakazam::estimateAbundance."
      input: numeric
      value: 30
   outputby:
      label: "`outputby` name of the column in `input` that contains sample identifiers that will be used to split the output db." 
      input: text
      value: "sample_id"
   model:
      label: "`model`"
      input: text
      value: "hierarchical" 
   method:
       label: "`method`"
       input: text
       value: "nt"
   linkage:
       label: "`linkage`"
       input: text
       value: "single"
   nboot:
       label: "`nboot`: number of bootstrap realizations to generate in `estimateAbundance`"
       input: numeric
       value: 200               
   outname:
      label: "`outname`"
      input: text
      value: "define-clones"
   nproc:
      label: "`nproc`"
      input: numeric
      value: 1      
   log:
      label: "`log`"
      input: text
      value: "command_log"
   outdir:
      label: '`outdir`: Output directory'
      input: text
      value: !r file.path(getwd(),'enchantr')
   date: 
      label: '`date`: Run date'
      input: date
      value: !r format(Sys.time(), "%Y-%m-%d")
   logo:
      label: "`logo`: Path to report logo"
      input: file
      value: !r file.path("assets", "logo.png")
   logolink:
      label: "`logolink`: URL to be added to the logo"
      input: text
      value: "immcantation.org"
   echo: 
      label: '`echo`: Show code in the report.'
      input: checkbox
      value: TRUE
   cache:
      label: '`cache`: Use cached results'
      input: checkbox
      value: FALSE
editor_options: 
  chunk_output_type: console
---

<div style="margin-top: -10px; font-size: 16px; color: gray;">
Report of clonal assignment, clone size distribution, clonal abundance, diversity and mutation frequency created by immcantation enchantR package. Immcantation package SCOPer is used to define clones. Package Alakazam is used for clonal abundance and diversity analysis. Package SHazaM is used for advanced analysis of somatic hypermutation (SHM). 
</div>





# Input parameters


``` r
printParams(params)
```

```{=html}
<div class="datatables html-widget html-fill-item" id="input20250404120042" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="input20250404120042">{"x":{"filter":"top","vertical":false,"filterHTML":"<tr>\n  <td data-type=\"factor\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"width: 100%; display: none;\">\n      <select multiple=\"multiple\" style=\"width: 100%;\" data-options=\"[&quot;input&quot;,&quot;imgt_db&quot;,&quot;species&quot;,&quot;force&quot;,&quot;cloneby&quot;,&quot;threshold&quot;,&quot;min_n&quot;,&quot;outputby&quot;,&quot;model&quot;,&quot;method&quot;,&quot;linkage&quot;,&quot;nboot&quot;,&quot;outname&quot;,&quot;nproc&quot;,&quot;log&quot;,&quot;outdir&quot;,&quot;date&quot;,&quot;logo&quot;,&quot;logolink&quot;,&quot;echo&quot;,&quot;cache&quot;]\"><\/select>\n    <\/div>\n  <\/td>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n<\/tr>","fillContainer":false,"data":[["input","imgt_db","species","force","cloneby","threshold","min_n","outputby","model","method","linkage","nboot","outname","nproc","log","outdir","date","logo","logolink","echo","cache"],["https://yale.box.com/shared/static/dajikxtvbc3c5yw1o2v69fjtx5bh58r2.tsv","https://raw.githubusercontent.com/nf-core/test-datasets/airrflow/database-cache/imgtdb_base.zip","human","FALSE","subject_id","1","30","sample_id","hierarchical","nt","single","200","define-clones","1","command_log","D:/Yale_Keinstein_lab/Airrflow_report_Improvement/enchantr/inst/rstudio/templates/project/define_clones_project_files/enchantr","2025-04-04","assets/logo.png","immcantation.org","TRUE","FALSE"]],"container":"<table class=\"stripe hover order-column row-border compact\">\n  <thead>\n    <tr>\n      <th>parameter<\/th>\n      <th>value<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"scrollX":true,"pageLength":5,"initComplete":"function(settings, json) {\n$(this.api().table().header()).css({'font-size': '50% !important'});\n}","columnDefs":[{"name":"parameter","targets":0},{"name":"value","targets":1}],"order":[],"autoWidth":false,"orderClasses":false,"orderCellsTop":true,"lengthMenu":[5,10,25,50,100]}},"evals":["options.initComplete"],"jsHooks":[]}</script>
```
<table>
<caption>(#tab:input) Input parameters.</caption>
</table>

# Read repertoires
<div style="margin-top: -10px; font-size: 15px; color: gray;">
Tables summarizing repertoire distributions for each sample.
</div>


``` r
# Read repertoire
db <- readInput(params[['input']], col_select = NULL)

if (params$species =="auto") {
    if ("species" %in% colnames(db)) {
        species <- unique(db[['species']])
        if (length(species)>1) {
            stop("Multiple species detected. Only one allowed.")
        }
    } else  {
        stop("Can't detect species. The column `species` does not exist in `db`.")
    }
} else {
    species <- params$species
}

# Check for single cell label
if (!is.null(singlecell)) {
    if (singlecell %in% colnames(db)) {
         db[[ singlecell ]] <- as.logical(db[[ singlecell ]] )   
    } else {
        stop("`",singlecell, "` is not a valid field in `db`.")
    }
} else {
   singlecell <- "single_cell"
   db[[singlecell]] <- F
   if ("cell_id" %in% colnames(db)) {
       message("Setting `singlecell` using `cell_id`.") 
       db[[singlecell]][!db[['cell_id']] %in% c(NA, '')] <- T
   }

}

if (!"locus" %in% colnames(db)) {
    db[['locus']] <- getLocus(db[['v_call']])
}
heavy_chains <- isHeavyChain(db[['locus']])

if ("clone_id" %in% colnames(db)) {
    if (params$force) {
        # Reset if force
        warning("Overwritting clone_id.")
        db$clone_id <- NULL
    }   
    # Reset always if clone_id exists
    if ("clone_size_count" %in% colnames(db)) {    
        warning("Overwritting clone_size_count.")    
        db$clone_size_count <- NULL
    }
    if ("clone_size_freq" %in% colnames(db)) {    
        warning("Overwritting clone_size_freq.")    
        db$clone_size_freq <- NULL
    }
}
```


``` r
input_size <- nrow(db)
input_sizes <- db %>%
    group_by(!!!rlang::syms(unique(c("input_file")))) %>%
    summarize(input_size=n())

# Track input size of the expected output groups, for the
# end of report summary
input_sizes_byoutput <- db %>%
    group_by(!!!rlang::syms(unique(params$outputby))) %>%
    mutate(
        num_input_files = length(unique(input_file)),
        input_files = paste(unique(input_file), collapse=",")
    ) %>%
    ungroup() %>%
    group_by(!!!rlang::syms(unique(c(params$outputby, "input_file", "num_input_files", "input_files")))) %>%
    summarize(
        input_size = n()
    )
```


Number of sequences loaded: 3000. Number of heavy chain sequences loaded: 3000.


``` r
input_samples_summary <- db %>%
   group_by(sample_id, subject_id, tissue) %>%
   summarize(size = n()) %>%
   arrange(subject_id)

eetable(input_samples_summary)$table
```



## Sequences per locus


``` r
input_locus_summary <- db %>%
   group_by(!!!rlang::syms(unique(c("sample_id",params$cloneby, "locus")))) %>%
   summarize(n=n(), .groups="drop") %>%
   pivot_wider(names_from=locus, values_from=n) %>%
   rowwise() %>%
   mutate(Total = sum(!!!rlang::syms(unique(db[['locus']]))))

total <- data.frame(list("sample_id"="Total", 
                    t(input_locus_summary %>%
   select(!!!rlang::syms(c(unique(db[['locus']]), "Total"))) %>%
   colSums(na.rm = T))))

input_locus_summary <- bind_rows(input_locus_summary, total)

tab_caption <- paste0("Input data. Number of sequences in each ", 
                   paste("sample", params$cloneby, sep=", "),
                   " and locus."
      )
eetable(input_locus_summary)$table
```

```{=html}
<div class="datatables html-widget html-fill-item" id="input_locus_summary20250404120049" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="input_locus_summary20250404120049">{"x":{"filter":"top","vertical":false,"filterHTML":"<tr>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"180\" data-max=\"3000\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"180\" data-max=\"3000\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n<\/tr>","fillContainer":false,"data":[["P05_FNA_0_Y1","P05_FNA_12_Y1","P05_FNA_2_0_Y1","P05_FNA_2_12_Y1","P05_FNA_2_28_Y1","P05_FNA_2_5_Y1","P05_FNA_2_60_Y1","P05_FNA_3_12_Y1","P05_FNA_3_28_Y1","P05_FNA_5_Y1","P05_FNA_60_Y1","Total"],["P05","P05","P05","P05","P05","P05","P05","P05","P05","P05","P05",null],[218,208,293,180,221,245,465,245,189,283,453,3000],[218,208,293,180,221,245,465,245,189,283,453,3000]],"container":"<table class=\"stripe hover order-column row-border compact\">\n  <thead>\n    <tr>\n      <th>sample_id<\/th>\n      <th>subject_id<\/th>\n      <th>IGH<\/th>\n      <th>Total<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"scrollX":true,"pageLength":5,"initComplete":"function(settings, json) {\n$(this.api().table().header()).css({'font-size': '50% !important'});\n}","columnDefs":[{"className":"dt-right","targets":[2,3]},{"name":"sample_id","targets":0},{"name":"subject_id","targets":1},{"name":"IGH","targets":2},{"name":"Total","targets":3}],"order":[],"autoWidth":false,"orderClasses":false,"orderCellsTop":true,"lengthMenu":[5,10,25,50,100]}},"evals":["options.initComplete"],"jsHooks":[]}</script>
```
<table>
<caption>(#tab:input-locus-summary) Input data. Number of sequences in each sample, subject_id and locus.</caption>
</table>


``` r
cat("## Sequences per constant region\n")
```


``` r
# Create c_gene column only for plotting (will be removed before storing the df again)
db <- db %>% filter(locus %in% c("IGH", "IGK", "IGL")) %>%
    mutate(c_gene = alakazam::getGene(c_call, first=TRUE))

input_cgene_summary <- db %>%
   group_by(!!!rlang::syms(unique(c("sample_id",params$cloneby, "c_gene")))) %>%
   summarize(n=n(), .groups="drop") %>%
       pivot_wider(names_from=c_gene, values_from=n)

tab_caption <- paste0("Input data. Number of sequences in each ", 
                   paste("sample", params$cloneby, sep=", "),
                   " and c_gene"
      )
eetable(input_cgene_summary)$table
```



# Clonal assignment
<div style="margin-top: -10px; font-size: 15px; color: gray;">
Clonal assignment performed with scoper::hierarchicalClones.This approach clusters B or T cell receptor sequences based on junction region sequence similarity within partitions that share the same V gene, J gene, and junction length, allowing for ambiguous V or J gene annotations.By default, 'subject_id'is used for grouping. Sequences with different subject_id will be classified as separate clones.
</div>





Clonal assignment performed with `scoper::hierarchicalClones`, version 1.3.0 within subject_id.

To know more details about the method, visit the documentation website [https://scoper.readthedocs.io/en/](https://scoper.readthedocs.io/en/1.3.0/topics/hierarchicalClones/)


``` r
if (sum(heavy_chains)>0) {
    if (params$model == "hierarchical") {
        if (all(db[[ singlecell ]] == T )) {
            cell_id <- 'cell_id'
        } else {
            db_light <- db %>% filter(!isHeavyChain(locus))
            cell_id <- NULL
            if (all(c(T,F) %in% db[[ singlecell ]])) {
                warning("Mix of single and bulk data. Setting cell_id=`NULL`.")
            }
        }
        db <- hierarchicalClones(db, 
                                 params$threshold, 
                                 method=params$method,
                                 linkage=params$linkage, 
                                 normalize="len",
                                 junction="junction", 
                                 v_call="v_call", j_call="j_call", 
                                 clone="clone_id", 
                                 fields=params$cloneby,
                                 cell_id=cell_id, 
                                 locus="locus", 
                                 only_heavy=TRUE, 
                                 split_light=TRUE,
                                 first=FALSE, 
                                 cdr3=FALSE, mod3=FALSE, 
                                 max_n=0, nproc=params$nproc,
                                 verbose=FALSE, log=NULL,
                                 summarize_clones=FALSE) 
    } else {
        stop("Unsuported model requested. Supported models: hierarchical")
    }
} else {
    warning("No heavy chain sequences found.")
    db$clone_id <- NA
}
```

```
## Warning in defineClonesScoper(db = db, threshold = threshold, model =
## "hierarchical", : Single cell mode requested, but `db` doesn't contain light
## chain data. Skipping.
```

3000 sequences passed the clonal assignment step 
and 0 were removed. 0 sequences
have `clone_id==NA`.

## Create germlines


``` r
dowser_v <- packageVersion("dowser")
```


Identification of the V(D)J germline sequences from which each of the observed 
sequences is derived is performed with `dowser::createGermlines`, version 2.3.
These reconstructed germlines will be used in downstream analysis to infer somatic 
mutations and reconstruct lineages. 

`dowser::createGermlines` takes the alignment information 
in the rearrangement file as well as the reference database used by the 
alignment software and generates a germline sequence for each individual observed sequence.
Because clonal relations have already been inferred, the function assigns the 
same germline to all sequences belonging to the same clone. 

Two types of germlines
are created, `germline_alignment` and `germline_alignment_d_mask`. The last one
has the D region masked, meaning that all nucleotides in the N/P and D-segments
are replaced with N's. This is often done because the germline base calls from 
this region are unreliable for B cell receptor alignments.

Documentation for `dowser::createGermlines` is available here: [https://dowser.readthedocs.io/en/latest/topics/createGermlines](https://dowser.readthedocs.io/en/latest/topics/createGermlines/).


``` r
pre_germ_size <- nrow(db)
db[['tmp_nrow']] <- 1:nrow(db)

# download and unzip if needed
imgt_db <- prepareIMGT(params$imgt_db)

references <- dowser::readIMGT(file.path(imgt_db,
                                species, "vdj"),
                               quiet=TRUE)
```

```
## [1] "Read in 1199 from 17 fasta files"
```

``` r
# tmp fix (add unique sequence id ) for
# Error in dowser::createGermlines(db_sp, references, locus = "locus",  :
#   Sequence IDs are not unique!
dup_ids <- db %>%
    ungroup() %>%
    group_by(!!!rlang::syms(c(params$cloneby))) %>%
    mutate(dup=duplicated(sequence_id)) %>%
    pull(dup) %>% any()

if (dup_ids) {
    db <- db %>%
        mutate(sequence_id = paste(sample_id, sequence_id, sep = '_'))
}

is_na_clone <- is.na(db[["clone_id"]])
not_na_clone_db <- data.frame()
na_clone_db <- data.frame()
if (any(is_na_clone)) {
    na_clone_db <- dowser::createGermlines(
        db[is_na_clone,,drop=FALSE],
        references,
        locus = "locus",
        nproc = params$nproc,
        seq = "sequence_alignment",
        v_call = "v_call",
        d_call = "d_call",
        j_call = "j_call",
        amino_acid = FALSE,
        id = "sequence_id",
        clone = "tmp_nrow",
        v_germ_start = "v_germline_start",
        v_germ_end = "v_germline_end",
        v_germ_length = "v_germline_length",
        d_germ_start = "d_germline_start",
        d_germ_end = "d_germline_end",
        d_germ_length = "d_germline_length",
        j_germ_start = "j_germline_start",
        j_germ_end = "j_germline_end",
        j_germ_length = "j_germline_length",
        np1_length = "np1_length",
        np2_length = "np2_length",
        na.rm = TRUE,
        fields = params$cloneby
    )    
}

if (any(!is_na_clone)) {
    not_na_clone_db <- dowser::createGermlines(
        db[!is_na_clone,,drop=FALSE],
        references,
        locus = "locus",
        nproc = params$nproc,
        seq = "sequence_alignment",
        v_call = "v_call",
        d_call = "d_call",
        j_call = "j_call",
        amino_acid = FALSE,
        id = "sequence_id",
        clone = "clone_id",
        v_germ_start = "v_germline_start",
        v_germ_end = "v_germline_end",
        v_germ_length = "v_germline_length",
        d_germ_start = "d_germline_start",
        d_germ_end = "d_germline_end",
        d_germ_length = "d_germline_length",
        j_germ_start = "j_germline_start",
        j_germ_end = "j_germline_end",
        j_germ_length = "j_germline_length",
        np1_length = "np1_length",
        np2_length = "np2_length",
        na.rm = TRUE,
        fields = params$cloneby
    )
}

db <- bind_rows(not_na_clone_db, na_clone_db) %>%
    arrange(tmp_nrow) %>%
    select(-tmp_nrow)

any_germ_fail <- pre_germ_size > nrow(db)
```

2998 sequences passed the germline reconstruction step
and 2 failed.


``` r
# Remove cells that after createGermlines
# have only light chains
db <- findLightOnlyCells(db, 
                         sample_id = "sample_id", cell_id = cell_id,
                         locus = "locus", fields = NULL)
num_light_only <- sum(db[["light_only_cell"]], na.rm = T)
if (num_light_only > 0) {
    # Remove the cells
    db <- db %>% dplyr::filter(!light_only_cell)
    cat(paste0("Removed ",num_light_only," sequences in light chain only cells."))
}
db[["light_only_cell"]] <- NULL
```


``` r
# Add again light chains before observed mutations when mixture of sc + bulk or only bulk
if (nrow(db_light) > 0){
    db <- bind_rows(db, db_light)
}
```

# Clone size distribution

Find your [clone sizes table here](tables/clone_sizes_table.tsv). Most real datasets, will have most clones of size 1 (one sequence). Straight sequence count as a measure of the size of the clones is not the best measure to compare clone size between samples due to possible disproportionate sampling. See the [Clonal abundance](#clonal_abundance) section.

Description of terms:

* `clone_size_count`: Clone size as sequence counts. In a sample (`sample_id`), the number of heavy chain 
sequences with the same `clone_id`.

* `clone_size_freq`: Clone size as percent of the repertoire. `clone_size_count` divided by the number of heavy chain sequences in the sample (`sample_id`).


## Number of clones 
<div style="margin-top: -10px; font-size: 15px; color: gray;">
Summary of the number of clones base on heavy chain sequences, and clone size per sample. Singletons which are clone of single sequence are included in the table. 
</div>


``` r
# Add clone_size 
clone_sizes <- countClones(
            db %>% filter(isHeavyChain(locus)), # Keep heavy chains only
            groups=unique(c("sample_id", params$cloneby)))
```



``` r
tmp <- eetable(clone_sizes, 
               caption=NULL, 
               outdir=params$outdir, file="clone_sizes_table")

db <- db %>%
   left_join(clone_sizes) %>%
   rename(
      clone_size_count = seq_count,
      clone_size_freq = seq_freq
   ) %>% 
   mutate_at(vars(starts_with("clone_size")), round,2)
```



``` r
num_clones_table <- db %>%
   filter(isHeavyChain(locus)) %>% # Keep heavy chains only
   group_by(sample_id) %>%
   mutate(sequences=n()) %>%
   group_by(sample_id) %>%
   mutate(
      number_of_clones=length(unique(clone_id)),
   ) %>%
   group_by(!!!rlang::syms(unique(c("sample_id","sequences", params$cloneby, "number_of_clones")))) %>%
   summarize_at(vars(starts_with("clone_size")), list("min"=min, "median"=median, "max"=max)) %>%
   mutate_at(vars(starts_with("clone_size")), round, 2)

tab_caption <- "Summary of the number of clones, and clone size, per sample. Includes singletons (clone_size == 1)."
tab <- eetable(num_clones_table, caption=tab_caption, outdir=params$outdir, file="num_clones_table")
tab$table
```

```{=html}
<div class="datatables html-widget html-fill-item" id="num_clones_table20250404120113" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="num_clones_table20250404120113">{"x":{"filter":"top","vertical":false,"filterHTML":"<tr>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"integer\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"180\" data-max=\"465\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\" disabled=\"\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"integer\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"133\" data-max=\"269\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\" disabled=\"\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0\" data-max=\"1\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0\" data-max=\"0.01\" data-scale=\"2\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"1\" data-max=\"2\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0\" data-max=\"0.01\" data-scale=\"2\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"4\" data-max=\"15\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0.02\" data-max=\"0.03\" data-scale=\"2\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n<\/tr>","fillContainer":false,"data":[["P05_FNA_0_Y1","P05_FNA_12_Y1","P05_FNA_2_0_Y1","P05_FNA_2_12_Y1","P05_FNA_2_28_Y1","P05_FNA_2_5_Y1","P05_FNA_2_60_Y1","P05_FNA_3_12_Y1","P05_FNA_3_28_Y1","P05_FNA_5_Y1","P05_FNA_60_Y1"],[217,208,293,180,221,245,465,244,189,283,453],["P05","P05","P05","P05","P05","P05","P05","P05","P05","P05","P05"],[167,164,201,133,169,177,269,171,146,198,252],[1,1,1,1,1,1,1,1,1,1,1],[0,0,0,0.01,0,0,0,0,0.01,0,0],[1,1,1,1,1,1,2,2,1,1,2],[0,0,0,0.01,0,0,0,0.01,0.01,0,0],[4,4,7,5,5,6,15,7,5,6,9],[0.02,0.02,0.02,0.03,0.02,0.02,0.03,0.03,0.03,0.02,0.02]],"container":"<table class=\"stripe hover order-column row-border compact\">\n  <thead>\n    <tr>\n      <th>sample_id<\/th>\n      <th>sequences<\/th>\n      <th>subject_id<\/th>\n      <th>number_of_clones<\/th>\n      <th>clone_size_count_min<\/th>\n      <th>clone_size_freq_min<\/th>\n      <th>clone_size_count_median<\/th>\n      <th>clone_size_freq_median<\/th>\n      <th>clone_size_count_max<\/th>\n      <th>clone_size_freq_max<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"scrollX":true,"pageLength":5,"initComplete":"function(settings, json) {\n$(this.api().table().header()).css({'font-size': '50% !important'});\n}","columnDefs":[{"className":"dt-right","targets":[1,3,4,5,6,7,8,9]},{"name":"sample_id","targets":0},{"name":"sequences","targets":1},{"name":"subject_id","targets":2},{"name":"number_of_clones","targets":3},{"name":"clone_size_count_min","targets":4},{"name":"clone_size_freq_min","targets":5},{"name":"clone_size_count_median","targets":6},{"name":"clone_size_freq_median","targets":7},{"name":"clone_size_count_max","targets":8},{"name":"clone_size_freq_max","targets":9}],"order":[],"autoWidth":false,"orderClasses":false,"orderCellsTop":true,"lengthMenu":[5,10,25,50,100]}},"evals":["options.initComplete"],"jsHooks":[]}</script>
```
<table>
<caption>(#tab:num-clones-table)  Summary of the number of clones, and clone size, per sample. Includes singletons (clone_size == 1). File can be found here: <a href='tables/num_clones_table.tsv'> num_clones_table.tsv</a></caption>
</table>


## Clone size distribution without singletons

<div style="margin-top: -10px; font-size: 15px; color: gray;">
Plot of clone size distribution. In most cases, the majority of clones are singletons (clones consisting of a single sequence). To improve readability for non-singleton clones, singletons are excluded from the plot.
</div>


``` r
caption <- "Clone size distribution, excluding singletons (subset to clone size > 1). Size is measured as number of heavy chain sequences belonging to the same clone."
clone_size_atleast2 <- ggplot(clone_sizes %>% filter(seq_count>1), aes(x=seq_count, color=sample_id, fill=sample_id))+
    geom_bar() + theme_enchantr() +
    facet_wrap(~sample_id, scales = "free_y", ncol=3) +
    xlab("Clone size (Sequences per clone)")
clone_size_atleast2 <- eeplot(clone_size_atleast2, 
       outdir=params$outdir, 
       file=knitr::opts_current$get('clone-size_atleast2'),
       caption=caption
       )
ggplotly(clone_size_atleast2 + theme(panel.spacing=unit(2, 'lines'), legend.position="right"))
```

<div class="figure">

```{=html}
<div class="plotly html-widget html-fill-item" id="htmlwidget-663ee8a3f8235a4018ac" style="width:768px;height:739.2px;"></div>
<script type="application/json" data-for="htmlwidget-663ee8a3f8235a4018ac">{"x":{"data":[{"orientation":"v","width":[0.90000000000000013,0.90000000000000036,0.90000000000000036],"base":[0,0,0],"x":[2,3,4],"y":[30,7,2],"text":["count: 30<br />seq_count:  2<br />sample_id: P05_FNA_0_Y1<br />sample_id: P05_FNA_0_Y1","count:  7<br />seq_count:  3<br />sample_id: P05_FNA_0_Y1<br />sample_id: P05_FNA_0_Y1","count:  2<br />seq_count:  4<br />sample_id: P05_FNA_0_Y1<br />sample_id: P05_FNA_0_Y1"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(248,118,109,1)","line":{"width":1.8897637795275593,"color":"rgba(248,118,109,1)"}},"name":"P05_FNA_0_Y1","legendgroup":"P05_FNA_0_Y1","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.90000000000000036,0.90000000000000036],"base":[0,0,0],"x":[2,3,4],"y":[30,4,2],"text":["count: 30<br />seq_count:  2<br />sample_id: P05_FNA_12_Y1<br />sample_id: P05_FNA_12_Y1","count:  4<br />seq_count:  3<br />sample_id: P05_FNA_12_Y1<br />sample_id: P05_FNA_12_Y1","count:  2<br />seq_count:  4<br />sample_id: P05_FNA_12_Y1<br />sample_id: P05_FNA_12_Y1"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(219,142,0,1)","line":{"width":1.8897637795275593,"color":"rgba(219,142,0,1)"}},"name":"P05_FNA_12_Y1","legendgroup":"P05_FNA_12_Y1","showlegend":true,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.90000000000000036,0.90000000000000036,0.90000000000000036,0.90000000000000036],"base":[0,0,0,0,0],"x":[2,3,4,5,7],"y":[28,15,6,1,2],"text":["count: 28<br />seq_count:  2<br />sample_id: P05_FNA_2_0_Y1<br />sample_id: P05_FNA_2_0_Y1","count: 15<br />seq_count:  3<br />sample_id: P05_FNA_2_0_Y1<br />sample_id: P05_FNA_2_0_Y1","count:  6<br />seq_count:  4<br />sample_id: P05_FNA_2_0_Y1<br />sample_id: P05_FNA_2_0_Y1","count:  1<br />seq_count:  5<br />sample_id: P05_FNA_2_0_Y1<br />sample_id: P05_FNA_2_0_Y1","count:  2<br />seq_count:  7<br />sample_id: P05_FNA_2_0_Y1<br />sample_id: P05_FNA_2_0_Y1"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(174,162,0,1)","line":{"width":1.8897637795275593,"color":"rgba(174,162,0,1)"}},"name":"P05_FNA_2_0_Y1","legendgroup":"P05_FNA_2_0_Y1","showlegend":true,"xaxis":"x3","yaxis":"y3","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.90000000000000036,0.90000000000000036],"base":[0,0,0],"x":[2,3,5],"y":[27,4,3],"text":["count: 27<br />seq_count:  2<br />sample_id: P05_FNA_2_12_Y1<br />sample_id: P05_FNA_2_12_Y1","count:  4<br />seq_count:  3<br />sample_id: P05_FNA_2_12_Y1<br />sample_id: P05_FNA_2_12_Y1","count:  3<br />seq_count:  5<br />sample_id: P05_FNA_2_12_Y1<br />sample_id: P05_FNA_2_12_Y1"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(100,178,0,1)","line":{"width":1.8897637795275593,"color":"rgba(100,178,0,1)"}},"name":"P05_FNA_2_12_Y1","legendgroup":"P05_FNA_2_12_Y1","showlegend":true,"xaxis":"x","yaxis":"y4","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.90000000000000036,0.90000000000000036],"base":[0,0,0],"x":[2,3,5],"y":[32,8,1],"text":["count: 32<br />seq_count:  2<br />sample_id: P05_FNA_2_28_Y1<br />sample_id: P05_FNA_2_28_Y1","count:  8<br />seq_count:  3<br />sample_id: P05_FNA_2_28_Y1<br />sample_id: P05_FNA_2_28_Y1","count:  1<br />seq_count:  5<br />sample_id: P05_FNA_2_28_Y1<br />sample_id: P05_FNA_2_28_Y1"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,189,92,1)","line":{"width":1.8897637795275593,"color":"rgba(0,189,92,1)"}},"name":"P05_FNA_2_28_Y1","legendgroup":"P05_FNA_2_28_Y1","showlegend":true,"xaxis":"x2","yaxis":"y5","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.90000000000000036,0.90000000000000036,0.90000000000000036],"base":[0,0,0,0],"x":[2,3,4,6],"y":[23,14,4,1],"text":["count: 23<br />seq_count:  2<br />sample_id: P05_FNA_2_5_Y1<br />sample_id: P05_FNA_2_5_Y1","count: 14<br />seq_count:  3<br />sample_id: P05_FNA_2_5_Y1<br />sample_id: P05_FNA_2_5_Y1","count:  4<br />seq_count:  4<br />sample_id: P05_FNA_2_5_Y1<br />sample_id: P05_FNA_2_5_Y1","count:  1<br />seq_count:  6<br />sample_id: P05_FNA_2_5_Y1<br />sample_id: P05_FNA_2_5_Y1"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,193,167,1)","line":{"width":1.8897637795275593,"color":"rgba(0,193,167,1)"}},"name":"P05_FNA_2_5_Y1","legendgroup":"P05_FNA_2_5_Y1","showlegend":true,"xaxis":"x3","yaxis":"y6","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.90000000000000036,0.90000000000000036,0.90000000000000036,0.90000000000000036,0.90000000000000036,0.89999999999999858,0.89999999999999858,0.89999999999999858],"base":[0,0,0,0,0,0,0,0,0],"x":[2,3,4,5,6,7,9,12,15],"y":[54,20,6,10,1,1,1,1,1],"text":["count: 54<br />seq_count:  2<br />sample_id: P05_FNA_2_60_Y1<br />sample_id: P05_FNA_2_60_Y1","count: 20<br />seq_count:  3<br />sample_id: P05_FNA_2_60_Y1<br />sample_id: P05_FNA_2_60_Y1","count:  6<br />seq_count:  4<br />sample_id: P05_FNA_2_60_Y1<br />sample_id: P05_FNA_2_60_Y1","count: 10<br />seq_count:  5<br />sample_id: P05_FNA_2_60_Y1<br />sample_id: P05_FNA_2_60_Y1","count:  1<br />seq_count:  6<br />sample_id: P05_FNA_2_60_Y1<br />sample_id: P05_FNA_2_60_Y1","count:  1<br />seq_count:  7<br />sample_id: P05_FNA_2_60_Y1<br />sample_id: P05_FNA_2_60_Y1","count:  1<br />seq_count:  9<br />sample_id: P05_FNA_2_60_Y1<br />sample_id: P05_FNA_2_60_Y1","count:  1<br />seq_count: 12<br />sample_id: P05_FNA_2_60_Y1<br />sample_id: P05_FNA_2_60_Y1","count:  1<br />seq_count: 15<br />sample_id: P05_FNA_2_60_Y1<br />sample_id: P05_FNA_2_60_Y1"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,186,222,1)","line":{"width":1.8897637795275593,"color":"rgba(0,186,222,1)"}},"name":"P05_FNA_2_60_Y1","legendgroup":"P05_FNA_2_60_Y1","showlegend":true,"xaxis":"x","yaxis":"y7","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.90000000000000036,0.90000000000000036,0.90000000000000036,0.90000000000000036],"base":[0,0,0,0,0],"x":[2,3,4,5,7],"y":[35,11,2,1,1],"text":["count: 35<br />seq_count:  2<br />sample_id: P05_FNA_3_12_Y1<br />sample_id: P05_FNA_3_12_Y1","count: 11<br />seq_count:  3<br />sample_id: P05_FNA_3_12_Y1<br />sample_id: P05_FNA_3_12_Y1","count:  2<br />seq_count:  4<br />sample_id: P05_FNA_3_12_Y1<br />sample_id: P05_FNA_3_12_Y1","count:  1<br />seq_count:  5<br />sample_id: P05_FNA_3_12_Y1<br />sample_id: P05_FNA_3_12_Y1","count:  1<br />seq_count:  7<br />sample_id: P05_FNA_3_12_Y1<br />sample_id: P05_FNA_3_12_Y1"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,166,255,1)","line":{"width":1.8897637795275593,"color":"rgba(0,166,255,1)"}},"name":"P05_FNA_3_12_Y1","legendgroup":"P05_FNA_3_12_Y1","showlegend":true,"xaxis":"x2","yaxis":"y8","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.90000000000000036,0.90000000000000036],"base":[0,0,0],"x":[2,3,5],"y":[21,9,1],"text":["count: 21<br />seq_count:  2<br />sample_id: P05_FNA_3_28_Y1<br />sample_id: P05_FNA_3_28_Y1","count:  9<br />seq_count:  3<br />sample_id: P05_FNA_3_28_Y1<br />sample_id: P05_FNA_3_28_Y1","count:  1<br />seq_count:  5<br />sample_id: P05_FNA_3_28_Y1<br />sample_id: P05_FNA_3_28_Y1"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(179,133,255,1)","line":{"width":1.8897637795275593,"color":"rgba(179,133,255,1)"}},"name":"P05_FNA_3_28_Y1","legendgroup":"P05_FNA_3_28_Y1","showlegend":true,"xaxis":"x3","yaxis":"y9","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.90000000000000036,0.90000000000000036,0.90000000000000036,0.90000000000000036],"base":[0,0,0,0,0],"x":[2,3,4,5,6],"y":[34,14,3,1,2],"text":["count: 34<br />seq_count:  2<br />sample_id: P05_FNA_5_Y1<br />sample_id: P05_FNA_5_Y1","count: 14<br />seq_count:  3<br />sample_id: P05_FNA_5_Y1<br />sample_id: P05_FNA_5_Y1","count:  3<br />seq_count:  4<br />sample_id: P05_FNA_5_Y1<br />sample_id: P05_FNA_5_Y1","count:  1<br />seq_count:  5<br />sample_id: P05_FNA_5_Y1<br />sample_id: P05_FNA_5_Y1","count:  2<br />seq_count:  6<br />sample_id: P05_FNA_5_Y1<br />sample_id: P05_FNA_5_Y1"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(239,103,235,1)","line":{"width":1.8897637795275593,"color":"rgba(239,103,235,1)"}},"name":"P05_FNA_5_Y1","legendgroup":"P05_FNA_5_Y1","showlegend":true,"xaxis":"x","yaxis":"y10","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.90000000000000036,0.90000000000000036,0.90000000000000036,0.90000000000000036,0.90000000000000036,0.89999999999999947,0.89999999999999858],"base":[0,0,0,0,0,0,0,0],"x":[2,3,4,5,6,7,8,9],"y":[58,22,11,5,2,1,2,2],"text":["count: 58<br />seq_count:  2<br />sample_id: P05_FNA_60_Y1<br />sample_id: P05_FNA_60_Y1","count: 22<br />seq_count:  3<br />sample_id: P05_FNA_60_Y1<br />sample_id: P05_FNA_60_Y1","count: 11<br />seq_count:  4<br />sample_id: P05_FNA_60_Y1<br />sample_id: P05_FNA_60_Y1","count:  5<br />seq_count:  5<br />sample_id: P05_FNA_60_Y1<br />sample_id: P05_FNA_60_Y1","count:  2<br />seq_count:  6<br />sample_id: P05_FNA_60_Y1<br />sample_id: P05_FNA_60_Y1","count:  1<br />seq_count:  7<br />sample_id: P05_FNA_60_Y1<br />sample_id: P05_FNA_60_Y1","count:  2<br />seq_count:  8<br />sample_id: P05_FNA_60_Y1<br />sample_id: P05_FNA_60_Y1","count:  2<br />seq_count:  9<br />sample_id: P05_FNA_60_Y1<br />sample_id: P05_FNA_60_Y1"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(255,99,182,1)","line":{"width":1.8897637795275593,"color":"rgba(255,99,182,1)"}},"name":"P05_FNA_60_Y1","legendgroup":"P05_FNA_60_Y1","showlegend":true,"xaxis":"x2","yaxis":"y11","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":35.282134934905258,"r":7.3059360730593621,"b":37.546975117553657,"l":37.260273972602747},"plot_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"ArialMT","size":14.611872146118724},"xaxis":{"domain":[0,0.30855126608551264],"automargin":true,"type":"linear","autorange":false,"range":[0.85500000000000009,16.145],"tickmode":"array","ticktext":["4","8","12","16"],"tickvals":[4,8,12,16],"categoryorder":"array","categoryarray":["4","8","12","16"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"ArialMT","size":11.689497716894984},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y10","title":"","hoverformat":".2f"},"annotations":[{"text":"Clone size (Sequences per clone)","x":0.5,"y":0,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"ArialMT","size":14.611872146118724},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"top","annotationType":"axis","yshift":-21.917808219178088},{"text":"count","x":0,"y":0.5,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"ArialMT","size":14.611872146118724},"xref":"paper","yref":"paper","textangle":-90,"xanchor":"right","yanchor":"center","annotationType":"axis","xshift":-21.917808219178088},{"text":"P05_FNA_0_Y1","x":0.15427563304275632,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(26,26,26,1)","family":"ArialMT","size":11.689497716894984},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"},{"text":"P05_FNA_12_Y1","x":0.5,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(26,26,26,1)","family":"ArialMT","size":11.689497716894984},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"},{"text":"P05_FNA_2_0_Y1","x":0.84572436695724362,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(26,26,26,1)","family":"ArialMT","size":11.689497716894984},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"},{"text":"P05_FNA_2_12_Y1","x":0.15427563304275632,"y":0.72458246056314535,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(26,26,26,1)","family":"ArialMT","size":11.689497716894984},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"},{"text":"P05_FNA_2_28_Y1","x":0.5,"y":0.72458246056314535,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(26,26,26,1)","family":"ArialMT","size":11.689497716894984},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"},{"text":"P05_FNA_2_5_Y1","x":0.84572436695724362,"y":0.72458246056314535,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(26,26,26,1)","family":"ArialMT","size":11.689497716894984},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"},{"text":"P05_FNA_2_60_Y1","x":0.15427563304275632,"y":0.47458246056314529,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(26,26,26,1)","family":"ArialMT","size":11.689497716894984},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"},{"text":"P05_FNA_3_12_Y1","x":0.5,"y":0.47458246056314529,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(26,26,26,1)","family":"ArialMT","size":11.689497716894984},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"},{"text":"P05_FNA_3_28_Y1","x":0.84572436695724362,"y":0.47458246056314529,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(26,26,26,1)","family":"ArialMT","size":11.689497716894984},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"},{"text":"P05_FNA_5_Y1","x":0.15427563304275632,"y":0.22458246056314529,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(26,26,26,1)","family":"ArialMT","size":11.689497716894984},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"},{"text":"P05_FNA_60_Y1","x":0.5,"y":0.22458246056314529,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(26,26,26,1)","family":"ArialMT","size":11.689497716894984},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"}],"yaxis":{"domain":[0.77541753943685465,1],"automargin":true,"type":"linear","autorange":false,"range":[-1.5,31.5],"tickmode":"array","ticktext":["0","10","20","30"],"tickvals":[0,10,20,30],"categoryorder":"array","categoryarray":["0","10","20","30"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"ArialMT","size":11.689497716894984},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":"","hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(51,51,51,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","x0":0,"x1":0.30855126608551264,"y0":0.77541753943685465,"y1":1},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":0.30855126608551264,"y0":0,"y1":23.37899543378996,"yanchor":1,"ysizemode":"pixel"},{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(51,51,51,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","x0":0.35811540058115399,"x1":0.64188459941884601,"y0":0.77541753943685465,"y1":1},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0.35811540058115399,"x1":0.64188459941884601,"y0":0,"y1":23.37899543378996,"yanchor":1,"ysizemode":"pixel"},{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(51,51,51,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","x0":0.69144873391448725,"x1":1,"y0":0.77541753943685465,"y1":1},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0.69144873391448725,"x1":1,"y0":0,"y1":23.37899543378996,"yanchor":1,"ysizemode":"pixel"},{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(51,51,51,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","x0":0,"x1":0.30855126608551264,"y0":0.52541753943685465,"y1":0.72458246056314535},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":0.30855126608551264,"y0":0,"y1":23.37899543378996,"yanchor":0.72458246056314535,"ysizemode":"pixel"},{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(51,51,51,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","x0":0.35811540058115399,"x1":0.64188459941884601,"y0":0.52541753943685465,"y1":0.72458246056314535},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0.35811540058115399,"x1":0.64188459941884601,"y0":0,"y1":23.37899543378996,"yanchor":0.72458246056314535,"ysizemode":"pixel"},{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(51,51,51,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","x0":0.69144873391448725,"x1":1,"y0":0.52541753943685465,"y1":0.72458246056314535},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0.69144873391448725,"x1":1,"y0":0,"y1":23.37899543378996,"yanchor":0.72458246056314535,"ysizemode":"pixel"},{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(51,51,51,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","x0":0,"x1":0.30855126608551264,"y0":0.27541753943685471,"y1":0.47458246056314529},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":0.30855126608551264,"y0":0,"y1":23.37899543378996,"yanchor":0.47458246056314529,"ysizemode":"pixel"},{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(51,51,51,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","x0":0.35811540058115399,"x1":0.64188459941884601,"y0":0.27541753943685471,"y1":0.47458246056314529},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0.35811540058115399,"x1":0.64188459941884601,"y0":0,"y1":23.37899543378996,"yanchor":0.47458246056314529,"ysizemode":"pixel"},{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(51,51,51,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","x0":0.69144873391448725,"x1":1,"y0":0.27541753943685471,"y1":0.47458246056314529},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0.69144873391448725,"x1":1,"y0":0,"y1":23.37899543378996,"yanchor":0.47458246056314529,"ysizemode":"pixel"},{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(51,51,51,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","x0":0,"x1":0.30855126608551264,"y0":0,"y1":0.22458246056314529},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":0.30855126608551264,"y0":0,"y1":23.37899543378996,"yanchor":0.22458246056314529,"ysizemode":"pixel"},{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(51,51,51,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","x0":0.35811540058115399,"x1":0.64188459941884601,"y0":0,"y1":0.22458246056314529},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0.35811540058115399,"x1":0.64188459941884601,"y0":0,"y1":23.37899543378996,"yanchor":0.22458246056314529,"ysizemode":"pixel"}],"xaxis2":{"type":"linear","autorange":false,"range":[0.85500000000000009,16.145],"tickmode":"array","ticktext":["4","8","12","16"],"tickvals":[4,8,12,16],"categoryorder":"array","categoryarray":["4","8","12","16"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"ArialMT","size":11.689497716894984},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"domain":[0.35811540058115399,0.64188459941884601],"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y11","title":"","hoverformat":".2f"},"yaxis2":{"type":"linear","autorange":false,"range":[-1.5,31.5],"tickmode":"array","ticktext":["0","10","20","30"],"tickvals":[0,10,20,30],"categoryorder":"array","categoryarray":["0","10","20","30"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"ArialMT","size":11.689497716894984},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"domain":[0.77541753943685465,1],"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x2","title":"","hoverformat":".2f"},"xaxis3":{"type":"linear","autorange":false,"range":[0.85500000000000009,16.145],"tickmode":"array","ticktext":["4","8","12","16"],"tickvals":[4,8,12,16],"categoryorder":"array","categoryarray":["4","8","12","16"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"ArialMT","size":11.689497716894984},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"domain":[0.69144873391448725,1],"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y9","title":"","hoverformat":".2f"},"yaxis3":{"type":"linear","autorange":false,"range":[-1.4000000000000001,29.399999999999999],"tickmode":"array","ticktext":["0","10","20"],"tickvals":[0,10,20],"categoryorder":"array","categoryarray":["0","10","20"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"ArialMT","size":11.689497716894984},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"domain":[0.77541753943685465,1],"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x3","title":"","hoverformat":".2f"},"yaxis4":{"type":"linear","autorange":false,"range":[-1.3500000000000001,28.350000000000001],"tickmode":"array","ticktext":["0","10","20"],"tickvals":[0,10,20],"categoryorder":"array","categoryarray":["0","10","20"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"ArialMT","size":11.689497716894984},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"domain":[0.52541753943685465,0.72458246056314535],"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":"","hoverformat":".2f"},"yaxis5":{"type":"linear","autorange":false,"range":[-1.6000000000000001,33.600000000000001],"tickmode":"array","ticktext":["0","10","20","30"],"tickvals":[0,10,20,30],"categoryorder":"array","categoryarray":["0","10","20","30"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"ArialMT","size":11.689497716894984},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"domain":[0.52541753943685465,0.72458246056314535],"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x2","title":"","hoverformat":".2f"},"yaxis6":{"type":"linear","autorange":false,"range":[-1.1500000000000001,24.149999999999999],"tickmode":"array","ticktext":["0","5","10","15","20"],"tickvals":[0,5,10,14.999999999999998,20],"categoryorder":"array","categoryarray":["0","5","10","15","20"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"ArialMT","size":11.689497716894984},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"domain":[0.52541753943685465,0.72458246056314535],"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x3","title":"","hoverformat":".2f"},"yaxis7":{"type":"linear","autorange":false,"range":[-2.7000000000000002,56.700000000000003],"tickmode":"array","ticktext":["0","20","40"],"tickvals":[0,20,40],"categoryorder":"array","categoryarray":["0","20","40"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"ArialMT","size":11.689497716894984},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"domain":[0.27541753943685471,0.47458246056314529],"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":"","hoverformat":".2f"},"yaxis8":{"type":"linear","autorange":false,"range":[-1.75,36.75],"tickmode":"array","ticktext":["0","10","20","30"],"tickvals":[0,10,20,30],"categoryorder":"array","categoryarray":["0","10","20","30"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"ArialMT","size":11.689497716894984},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"domain":[0.27541753943685471,0.47458246056314529],"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x2","title":"","hoverformat":".2f"},"yaxis9":{"type":"linear","autorange":false,"range":[-1.05,22.050000000000001],"tickmode":"array","ticktext":["0","5","10","15","20"],"tickvals":[0,5,10,15,20],"categoryorder":"array","categoryarray":["0","5","10","15","20"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"ArialMT","size":11.689497716894984},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"domain":[0.27541753943685471,0.47458246056314529],"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x3","title":"","hoverformat":".2f"},"yaxis10":{"type":"linear","autorange":false,"range":[-1.7000000000000002,35.700000000000003],"tickmode":"array","ticktext":["0","10","20","30"],"tickvals":[2.2204460492503131e-16,10,20,30],"categoryorder":"array","categoryarray":["0","10","20","30"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"ArialMT","size":11.689497716894984},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"domain":[0,0.22458246056314529],"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":"","hoverformat":".2f"},"yaxis11":{"type":"linear","autorange":false,"range":[-2.9000000000000004,60.899999999999999],"tickmode":"array","ticktext":["0","20","40","60"],"tickvals":[0,20,40,60],"categoryorder":"array","categoryarray":["0","20","40","60"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"ArialMT","size":11.689497716894984},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"domain":[0,0.22458246056314529],"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x2","title":"","hoverformat":".2f"},"showlegend":true,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"font":{"color":"rgba(0,0,0,1)","family":"ArialMT","size":11.689497716894984},"title":{"text":"sample_id","font":{"color":"rgba(0,0,0,1)","family":"ArialMT","size":14.611872146118724}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"2e744b0b34d":{"x":{},"colour":{},"fill":{},"type":"bar"}},"cur_data":"2e744b0b34d","visdat":{"2e744b0b34d":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

<p class="caption">(\#fig:clone-size-atleast2)Clone size distribution, excluding singletons (subset to clone size > 1). Size is measured as number of heavy chain sequences belonging to the same clone. <a href='ggplots/clone_size_atleast2.RData'>ggplot file: clone_size_atleast2.RData</a></p>
</div>

``` r
caption <- clone_size_atleast2$enchantr$html_caption
```


``` r
cat("Only singletons detected: there aren't clones of size>1.\n\n")
```


# Clonal abundance {#clonal_abundance}

<div style="margin-top: -10px; font-size: 15px; color: gray;">
Clonal abundance is the size of each clone (as a fraction of the entire repertoire). To correct for the different number of sequences in each of the samples, `estimateAbundance` estimates the clonal abundance distribution along with confidence intervals on these clone sizes using bootstrapping. 200 random bootstrap samples were taken, with size the number of sequences in the sample with less sequences (N). The y-axis shows the clone abundance (i.e., the size as a percent of the repertoire) and the x-axis is a rank of each clone, where the rank is sorted by clone size from larger (rank 1, left) to smaller (right). The shaded areas are confidence intervals.
</div>


``` r
# set empty default abundance
a <- new("AbundanceCurve")
```


``` r
# calculate the rank-abundance curve
a <- estimateAbundance(db %>% filter(isHeavyChain(locus)),
                       group = "sample_id", min_n = params$min_n)
if (nrow(a@abundance)==0) {
    cat("\nAll groups failed to pass the threshold min_n=",params$min_n,". Skipping clonal abundance report.\n\n")
}
```




``` r
# annotate
a@abundance <- a@abundance %>%
   left_join( db %>% 
                 select(any_of(c(unique(c("sample_id", 
                                                "subject_id",
                                                params$cloneby)))
                                       )
                  ) %>%
                 distinct(),
              
              )
p <- plotAbundanceCurve(a, annotate="depth", silent = T)
```


``` r
tab <- eetable(a@abundance, file = "clonal_abundance", 
               outdir=params$outdir, show_max=10,
               caption="Example 10 lines of the clonal abundance file.")
tab$table %>%
        DT::formatRound(columns=c("p","p_sd","lower","upper"),digits=3)
```

```{=html}
<div class="datatables html-widget html-fill-item" id="a@abundance20250404120115" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="a@abundance20250404120115">{"x":{"filter":"top","vertical":false,"filterHTML":"<tr>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\" disabled=\"\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0.007333333333333\" data-max=\"0.017972222222223\" data-scale=\"15\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0.006441382433574\" data-max=\"0.01028096514477\" data-scale=\"15\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0\" data-max=\"0.000102688881199\" data-scale=\"15\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0.019958210913788\" data-max=\"0.037094765854504\" data-scale=\"15\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"integer\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"1\" data-max=\"10\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\" disabled=\"\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n<\/tr>","fillContainer":false,"data":[["P05_FNA_0_Y1","P05_FNA_0_Y1","P05_FNA_0_Y1","P05_FNA_0_Y1","P05_FNA_0_Y1","P05_FNA_0_Y1","P05_FNA_0_Y1","P05_FNA_0_Y1","P05_FNA_0_Y1","P05_FNA_0_Y1"],["201","166","100","13","475","134","655","472","571","578"],[0.01797222222222222,0.01694444444444445,0.013,0.01275,0.01261111111111111,0.0125,0.01236111111111111,0.01186111111111111,0.01172222222222222,0.007333333333333333],[0.009117276379554231,0.01028096514476922,0.008136804109443327,0.008890241033395158,0.008441233336983116,0.007871533507947158,0.00800221563400056,0.008083976318537791,0.008181475082433132,0.00644138243357447],[0.0001026888811981978,0,0,0,0,0,0,0,0,0],[0.03584175556324624,0.03709476585450373,0.02894784300376643,0.03017455223933466,0.02915562443669687,0.02792792217867666,0.02804516555027556,0.02770541354731988,0.02775761872420303,0.01995821091378826],[1,2,3,4,5,6,7,8,9,10],["P05","P05","P05","P05","P05","P05","P05","P05","P05","P05"]],"container":"<table class=\"stripe hover order-column row-border compact\">\n  <thead>\n    <tr>\n      <th>sample_id<\/th>\n      <th>clone_id<\/th>\n      <th>p<\/th>\n      <th>p_sd<\/th>\n      <th>lower<\/th>\n      <th>upper<\/th>\n      <th>rank<\/th>\n      <th>subject_id<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"scrollX":true,"pageLength":5,"initComplete":"function(settings, json) {\n$(this.api().table().header()).css({'font-size': '50% !important'});\n}","columnDefs":[{"targets":2,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 3, 3, \",\", \".\", null);\n  }"},{"targets":3,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 3, 3, \",\", \".\", null);\n  }"},{"targets":4,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 3, 3, \",\", \".\", null);\n  }"},{"targets":5,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 3, 3, \",\", \".\", null);\n  }"},{"className":"dt-right","targets":[2,3,4,5,6]},{"name":"sample_id","targets":0},{"name":"clone_id","targets":1},{"name":"p","targets":2},{"name":"p_sd","targets":3},{"name":"lower","targets":4},{"name":"upper","targets":5},{"name":"rank","targets":6},{"name":"subject_id","targets":7}],"order":[],"autoWidth":false,"orderClasses":false,"orderCellsTop":true,"lengthMenu":[5,10,25,50,100]}},"evals":["options.initComplete","options.columnDefs.0.render","options.columnDefs.1.render","options.columnDefs.2.render","options.columnDefs.3.render"],"jsHooks":[]}</script>
```
<table>
<caption>(#tab:clonal-abundance)  Example 10 lines of the clonal abundance file. File can be found here: <a href='tables/clonal_abundance.tsv'> clonal_abundance.tsv</a></caption>
</table>

## Abundance plot by sample





``` r
abundanceSample <- p + facet_wrap(~ subject_id + sample_id, ncol=3)
abundanceSample <- eeplot(abundanceSample, 
                     outdir=params$outdir, 
                     file=knitr::opts_current$get('abundanceSample'),
                     caption=paste0("Clonal abundance plot per `sample_id`. Clonal abundance was calculated with a sample of N=",
                      a@n,
                      " sequences with ",
                      a@nboot,
                      " boostrapping repetitions."),
                     a)
abundanceSample +
    theme(legend.position = "none")
```

<div class="figure">
<img src="figures/abundanceSample-1.png" alt="Clonal abundance plot per `sample_id`. Clonal abundance was calculated with a sample of N=180 sequences with 200 boostrapping repetitions. &lt;a href='ggplots/abundanceSample.RData'&gt;ggplot file: abundanceSample.RData&lt;/a&gt;" width="576" />
<p class="caption">(\#fig:abundanceSample)Clonal abundance plot per `sample_id`. Clonal abundance was calculated with a sample of N=180 sequences with 200 boostrapping repetitions. <a href='ggplots/abundanceSample.RData'>ggplot file: abundanceSample.RData</a></p>
</div>

``` r
caption <- abundanceSample$enchantr$html_caption  
```

# Diversity

The clonal abundance distribution can be characterized using diversity statistics. Diversity scores (D) are calculated using the generalized diversity index (Hill numbers), which covers many different measures of diversity in a single function with a single varying parameter, the diversity order q.

The function alphaDiversity resamples the sequences (200 random bootstrapping events, with the number of sequences in the sample with less sequences (N)) and calculates diversity scores (D) over a interval of diversity orders (q). The diversity (D) is shown on the y-axis and the x-axis is the parameter q. - q = 0 corresponds to Species Richness - q = 1 corresponds to Shannon Entropy - q = 2 corresponds to Simpson Index

Inspection of this figure is useful to determine whether any difference in diversity between two repertoires depends on the statistic used or if it is a universal property. 


The clonal diversity $D$ of the repertoire was calculated according to the general formula of Hill Diversity
numbers:

$$
\begin{aligned}
    ^{q}D = \left( \sum_{i=1}^Rp_i^q \right)^{1/(1-q)}
\end{aligned}
$$

where:

* $p_i$ is the proportion of unique sequences belonging to clone $i$.
* $q$ are the values of the different diversity numbers.
* $R$ is the Richness, the number of different clones in the sample.

At $q=1$ the function is undefined and the limit to zero equals the exponential of the Shannon Entropy:

$$
\begin{aligned}
    ^{1}D = exp \left(  \sum_{i=1}^Rp_i ln(p_i)  \right)
\end{aligned}
$$

The intuition about the different Hill Diversity values is the following:

* At $q=0$ the diversity index equals the number of clones in the sample.
* At $q=1$ the diversity index is the geometric mean of the clones in the sample,
weighted by their proportion in the sample.
* At $q>1$ more weight is given to the clones with higher proportions in the sample.


## Diversity curves

The following table shows the summary of the diversity calculations per sample.


``` r
# generate the Hill diversity curve
d <- alphaDiversity(db %>% filter(isHeavyChain(locus)), 
                    group = "sample_id", 
                    min_n=params$min_n)

d@diversity <- d@diversity %>%
   left_join( db %>% 
                 select(any_of(c(unique(c("sample_id", 
                                                "subject_id", 
                                                params$cloneby)))
                                       )
                  ) %>%
                 distinct(),
              
              )

diversitySample <- plotDiversityCurve(d, silent = T, annotate="depth")

tab <- eetable(d@diversity, file = "clonal_diversity", 
               outdir=params$outdir, show_max=10,
               caption="Example 10 lines of the clonal diversity file.")
tab$table %>%
        DT::formatRound(columns=c("q","d","d_sd","d_lower", "d_upper", "e", "e_lower", "e_upper"),digits=3)
```

```{=html}
<div class="datatables html-widget html-fill-item" id="d@diversity20250404120142" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="d@diversity20250404120142">{"x":{"filter":"top","vertical":false,"filterHTML":"<tr>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\" disabled=\"\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0\" data-max=\"0.9\" data-scale=\"15\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"119.748833321211\" data-max=\"132.37\" data-scale=\"14\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"4.76725640891506\" data-max=\"5.71911493294869\" data-scale=\"15\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"108.539574029186\" data-max=\"123.026349133459\" data-scale=\"14\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"130.958092613235\" data-max=\"141.713650866541\" data-scale=\"13\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0.904652363233443\" data-max=\"1\" data-scale=\"15\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0.819971096390316\" data-max=\"0.929412624714503\" data-scale=\"15\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0.989333630076569\" data-max=\"1.0705873752855\" data-scale=\"15\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\" disabled=\"\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n<\/tr>","fillContainer":false,"data":[["P05_FNA_0_Y1","P05_FNA_0_Y1","P05_FNA_0_Y1","P05_FNA_0_Y1","P05_FNA_0_Y1","P05_FNA_0_Y1","P05_FNA_0_Y1","P05_FNA_0_Y1","P05_FNA_0_Y1","P05_FNA_0_Y1"],[0,0.1,0.2,0.3,0.4,0.5,0.6000000000000001,0.7000000000000001,0.8,0.9],[132.37,131.1026559487649,129.7998061436617,128.4618531806912,127.0894161486685,125.6833481435267,124.2447521142568,122.7749943325391,121.275714729315,119.7488333212109],[4.76725640891506,4.863570069189956,4.96206504206852,5.062816732137501,5.165904719357658,5.27141052352003,5.379414544388517,5.48999214103982,5.603208866046778,5.719114932948694],[123.0263491334587,121.5702337768657,120.0743373722622,118.538914725375,116.964428951162,115.351573369702,113.7012933493444,112.0148074606931,110.2936271540078,108.5395740291862],[141.7136508665413,140.6350781206642,139.5252749150612,138.3847916360075,137.214403346175,136.0151229173514,134.7882108791692,133.535181204385,132.2578023046222,130.9580926132355],[1,0.9904257456278986,0.9805832601319158,0.9704755849564949,0.960107397058763,0.9494851412217775,0.9386171497639706,0.9275137442965857,0.9161873138121552,0.9046523632334431],[0.9294126247145028,0.9184122820644078,0.9071114102308845,0.8955119341646518,0.8836173525055674,0.8714329029969178,0.8589657275012795,0.8462250318100254,0.8332222342978605,0.8199710963903167],[1.070587375285497,1.06243920919139,1.054055110032947,1.045439235748338,1.036597441611959,1.027537379446637,1.018268572026662,1.008802456783146,0.9991523933264501,0.9893336300765696],["P05","P05","P05","P05","P05","P05","P05","P05","P05","P05"]],"container":"<table class=\"stripe hover order-column row-border compact\">\n  <thead>\n    <tr>\n      <th>sample_id<\/th>\n      <th>q<\/th>\n      <th>d<\/th>\n      <th>d_sd<\/th>\n      <th>d_lower<\/th>\n      <th>d_upper<\/th>\n      <th>e<\/th>\n      <th>e_lower<\/th>\n      <th>e_upper<\/th>\n      <th>subject_id<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"scrollX":true,"pageLength":5,"initComplete":"function(settings, json) {\n$(this.api().table().header()).css({'font-size': '50% !important'});\n}","columnDefs":[{"targets":1,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 3, 3, \",\", \".\", null);\n  }"},{"targets":2,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 3, 3, \",\", \".\", null);\n  }"},{"targets":3,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 3, 3, \",\", \".\", null);\n  }"},{"targets":4,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 3, 3, \",\", \".\", null);\n  }"},{"targets":5,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 3, 3, \",\", \".\", null);\n  }"},{"targets":6,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 3, 3, \",\", \".\", null);\n  }"},{"targets":7,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 3, 3, \",\", \".\", null);\n  }"},{"targets":8,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 3, 3, \",\", \".\", null);\n  }"},{"className":"dt-right","targets":[1,2,3,4,5,6,7,8]},{"name":"sample_id","targets":0},{"name":"q","targets":1},{"name":"d","targets":2},{"name":"d_sd","targets":3},{"name":"d_lower","targets":4},{"name":"d_upper","targets":5},{"name":"e","targets":6},{"name":"e_lower","targets":7},{"name":"e_upper","targets":8},{"name":"subject_id","targets":9}],"order":[],"autoWidth":false,"orderClasses":false,"orderCellsTop":true,"lengthMenu":[5,10,25,50,100]}},"evals":["options.initComplete","options.columnDefs.0.render","options.columnDefs.1.render","options.columnDefs.2.render","options.columnDefs.3.render","options.columnDefs.4.render","options.columnDefs.5.render","options.columnDefs.6.render","options.columnDefs.7.render"],"jsHooks":[]}</script>
```
<table>
<caption>(#tab:clonal-diversity)  Example 10 lines of the clonal diversity file. File can be found here: <a href='tables/clonal_diversity.tsv'> clonal_diversity.tsv</a></caption>
</table>

## Diversity plot by sample


``` r
# plot duplicated cells
diversitySample <- diversitySample + 
    geom_vline(xintercept = c(0,1,2), color = "grey50", linetype = "dashed") +
    facet_wrap(~sample_id + subject_id, ncol=3, scales = "free_x")
diversitySample <- eeplot(diversitySample, 
                     outdir=params$outdir, 
                     file=knitr::opts_current$get('label'),
                     caption=paste0("Clonal diversity per `sample_id`.Clonal diversity was calculated with a sample of N=",
                      d@n,
                      " sequences with ",
                      a@nboot,
                      " boostrapping repetitions."),
                     d)
diversitySample +
    theme(legend.position = "none")
```

<div class="figure">
<img src="figures/diversitySample-1.png" alt="Clonal diversity per `sample_id`.Clonal diversity was calculated with a sample of N=180 sequences with 200 boostrapping repetitions. &lt;a href='ggplots/diversitySample.RData'&gt;ggplot file: diversitySample.RData&lt;/a&gt;" width="576" />
<p class="caption">(\#fig:diversitySample)Clonal diversity per `sample_id`.Clonal diversity was calculated with a sample of N=180 sequences with 200 boostrapping repetitions. <a href='ggplots/diversitySample.RData'>ggplot file: diversitySample.RData</a></p>
</div>

``` r
caption <- diversitySample$enchantr$html_caption
```

## Diversity at different `q`

- `q=0` corresponds to Richness (number of different clones in the sample).
- `q=1` corresponds to Shannon diversity.
- `q=2` corresponds to Simpson diversity.


``` r
# plot duplicated cells
div_data <- d@diversity %>% dplyr::filter(q==0 | q==1 | q==2)

div_at_q <- ggplot(div_data, aes(x=q, y=d, color=sample_id)) +
    geom_point() + theme_enchantr() +
    scale_x_continuous(breaks= c(0,1,2))
    
div_at_q <- eeplot(div_at_q, 
                     outdir=params$outdir, 
                     file=knitr::opts_current$get('label'),
                     caption="Diversity at different q values.",
                     div_at_q)
ggplotly(div_at_q)
```

<div class="figure">

```{=html}
<div class="plotly html-widget html-fill-item" id="htmlwidget-dde37cda1b8a48ed2a5d" style="width:480px;height:288px;"></div>
<script type="application/json" data-for="htmlwidget-dde37cda1b8a48ed2a5d">{"x":{"data":[{"x":[0,1,2],"y":[132.37,118.1981161323904,101.93680679604904],"text":["q: 0<br />d: 132.37000<br />sample_id: P05_FNA_0_Y1","q: 1<br />d: 118.19812<br />sample_id: P05_FNA_0_Y1","q: 2<br />d: 101.93681<br />sample_id: P05_FNA_0_Y1"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(248,118,109,1)","opacity":1,"size":5.6692913385826778,"symbol":"circle","line":{"width":1.8897637795275593,"color":"rgba(248,118,109,1)"}},"hoveron":"points","name":"P05_FNA_0_Y1","legendgroup":"P05_FNA_0_Y1","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0,1,2],"y":[134.46000000000001,121.04874463012031,105.42658204084533],"text":["q: 0<br />d: 134.46000<br />sample_id: P05_FNA_12_Y1","q: 1<br />d: 121.04874<br />sample_id: P05_FNA_12_Y1","q: 2<br />d: 105.42658<br />sample_id: P05_FNA_12_Y1"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(219,142,0,1)","opacity":1,"size":5.6692913385826778,"symbol":"circle","line":{"width":1.8897637795275593,"color":"rgba(219,142,0,1)"}},"hoveron":"points","name":"P05_FNA_12_Y1","legendgroup":"P05_FNA_12_Y1","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0,1,2],"y":[131.90000000000001,115.52337873822303,96.149961948617374],"text":["q: 0<br />d: 131.90000<br />sample_id: P05_FNA_2_0_Y1","q: 1<br />d: 115.52338<br />sample_id: P05_FNA_2_0_Y1","q: 2<br />d:  96.14996<br />sample_id: P05_FNA_2_0_Y1"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(174,162,0,1)","opacity":1,"size":5.6692913385826778,"symbol":"circle","line":{"width":1.8897637795275593,"color":"rgba(174,162,0,1)"}},"hoveron":"points","name":"P05_FNA_2_0_Y1","legendgroup":"P05_FNA_2_0_Y1","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0,1,2],"y":[119.54000000000001,102.24404251050048,83.558264289478231],"text":["q: 0<br />d: 119.54000<br />sample_id: P05_FNA_2_12_Y1","q: 1<br />d: 102.24404<br />sample_id: P05_FNA_2_12_Y1","q: 2<br />d:  83.55826<br />sample_id: P05_FNA_2_12_Y1"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(100,178,0,1)","opacity":1,"size":5.6692913385826778,"symbol":"circle","line":{"width":1.8897637795275593,"color":"rgba(100,178,0,1)"}},"hoveron":"points","name":"P05_FNA_2_12_Y1","legendgroup":"P05_FNA_2_12_Y1","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0,1,2],"y":[132.16,118.04966443659761,101.80616917238585],"text":["q: 0<br />d: 132.16000<br />sample_id: P05_FNA_2_28_Y1","q: 1<br />d: 118.04966<br />sample_id: P05_FNA_2_28_Y1","q: 2<br />d: 101.80617<br />sample_id: P05_FNA_2_28_Y1"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(0,189,92,1)","opacity":1,"size":5.6692913385826778,"symbol":"circle","line":{"width":1.8897637795275593,"color":"rgba(0,189,92,1)"}},"hoveron":"points","name":"P05_FNA_2_28_Y1","legendgroup":"P05_FNA_2_28_Y1","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0,1,2],"y":[132.75999999999999,116.70123666238962,97.773458636909083],"text":["q: 0<br />d: 132.76000<br />sample_id: P05_FNA_2_5_Y1","q: 1<br />d: 116.70124<br />sample_id: P05_FNA_2_5_Y1","q: 2<br />d:  97.77346<br />sample_id: P05_FNA_2_5_Y1"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(0,193,167,1)","opacity":1,"size":5.6692913385826778,"symbol":"circle","line":{"width":1.8897637795275593,"color":"rgba(0,193,167,1)"}},"hoveron":"points","name":"P05_FNA_2_5_Y1","legendgroup":"P05_FNA_2_5_Y1","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0,1,2],"y":[128.21000000000001,110.71585105454079,90.180410693808568],"text":["q: 0<br />d: 128.21000<br />sample_id: P05_FNA_2_60_Y1","q: 1<br />d: 110.71585<br />sample_id: P05_FNA_2_60_Y1","q: 2<br />d:  90.18041<br />sample_id: P05_FNA_2_60_Y1"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(0,186,222,1)","opacity":1,"size":5.6692913385826778,"symbol":"circle","line":{"width":1.8897637795275593,"color":"rgba(0,186,222,1)"}},"hoveron":"points","name":"P05_FNA_2_60_Y1","legendgroup":"P05_FNA_2_60_Y1","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0,1,2],"y":[125.72499999999999,109.60190544095155,91.466058949740997],"text":["q: 0<br />d: 125.72500<br />sample_id: P05_FNA_3_12_Y1","q: 1<br />d: 109.60191<br />sample_id: P05_FNA_3_12_Y1","q: 2<br />d:  91.46606<br />sample_id: P05_FNA_3_12_Y1"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(0,166,255,1)","opacity":1,"size":5.6692913385826778,"symbol":"circle","line":{"width":1.8897637795275593,"color":"rgba(0,166,255,1)"}},"hoveron":"points","name":"P05_FNA_3_12_Y1","legendgroup":"P05_FNA_3_12_Y1","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0,1,2],"y":[131.565,115.87246044209478,97.332965369503924],"text":["q: 0<br />d: 131.56500<br />sample_id: P05_FNA_3_28_Y1","q: 1<br />d: 115.87246<br />sample_id: P05_FNA_3_28_Y1","q: 2<br />d:  97.33297<br />sample_id: P05_FNA_3_28_Y1"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(179,133,255,1)","opacity":1,"size":5.6692913385826778,"symbol":"circle","line":{"width":1.8897637795275593,"color":"rgba(179,133,255,1)"}},"hoveron":"points","name":"P05_FNA_3_28_Y1","legendgroup":"P05_FNA_3_28_Y1","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0,1,2],"y":[131.84999999999999,116.33920124578316,98.123978205851969],"text":["q: 0<br />d: 131.85000<br />sample_id: P05_FNA_5_Y1","q: 1<br />d: 116.33920<br />sample_id: P05_FNA_5_Y1","q: 2<br />d:  98.12398<br />sample_id: P05_FNA_5_Y1"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(239,103,235,1)","opacity":1,"size":5.6692913385826778,"symbol":"circle","line":{"width":1.8897637795275593,"color":"rgba(239,103,235,1)"}},"hoveron":"points","name":"P05_FNA_5_Y1","legendgroup":"P05_FNA_5_Y1","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0,1,2],"y":[126.17,110.19665240010714,92.501953522757276],"text":["q: 0<br />d: 126.17000<br />sample_id: P05_FNA_60_Y1","q: 1<br />d: 110.19665<br />sample_id: P05_FNA_60_Y1","q: 2<br />d:  92.50195<br />sample_id: P05_FNA_60_Y1"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(255,99,182,1)","opacity":1,"size":5.6692913385826778,"symbol":"circle","line":{"width":1.8897637795275593,"color":"rgba(255,99,182,1)"}},"hoveron":"points","name":"P05_FNA_60_Y1","legendgroup":"P05_FNA_60_Y1","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":28.176560121765604,"r":7.3059360730593621,"b":42.130898021308994,"l":43.105022831050235},"plot_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"ArialMT","size":14.611872146118724},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.10000000000000001,2.1000000000000001],"tickmode":"array","ticktext":["0","1","2"],"tickvals":[0,1,2],"categoryorder":"array","categoryarray":["0","1","2"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"ArialMT","size":11.68949771689498},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":{"text":"q","font":{"color":"rgba(0,0,0,1)","family":"ArialMT","size":14.611872146118724}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[81.013177503952136,137.00508678552609],"tickmode":"array","ticktext":["90","100","110","120","130"],"tickvals":[90,100,110,120,130],"categoryorder":"array","categoryarray":["90","100","110","120","130"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"ArialMT","size":11.68949771689498},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":{"text":"d","font":{"color":"rgba(0,0,0,1)","family":"ArialMT","size":14.611872146118724}},"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(51,51,51,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":true,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"font":{"color":"rgba(0,0,0,1)","family":"ArialMT","size":11.68949771689498},"title":{"text":"sample_id","font":{"color":"rgba(0,0,0,1)","family":"ArialMT","size":14.611872146118724}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"2e74df930d0":{"x":{},"y":{},"colour":{},"type":"scatter"}},"cur_data":"2e74df930d0","visdat":{"2e74df930d0":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

<p class="caption">(\#fig:div_at_q)Diversity at different q values. <a href='ggplots/div_at_q.RData'>ggplot file: div_at_q.RData</a></p>
</div>

``` r
caption <- div_at_q$enchantr$html_caption
```


``` r
cat("# Mutation frequency\n")
```

# Mutation frequency


``` r
db <- shazam::observedMutations(db,
                                sequenceColumn = "sequence_alignment",
                                germlineColumn = "germline_alignment_d_mask",
                                frequency = T,
                                combine = T,
                                nproc = params$nproc)
```




``` r
cat("Showing mutation frequency per c_gene only if the `c_call` column is present.\n")
```

Showing mutation frequency per c_gene only if the `c_call` column is present.

``` r
cat("The mutation frequency per sequence is stored in the final dataframes in the `mu_freq` column.\n")
```

The mutation frequency per sequence is stored in the final dataframes in the `mu_freq` column.



``` r
db_cgene <- db %>% filter(!is.na(c_gene), !(is.na(mu_freq)))
db_cgene$mu_freq <- as.numeric(db_cgene$mu_freq)

if (nrow(db_cgene) > 0) {
    mufreqSample <- ggplot(db_cgene, aes(x=c_gene, y=mu_freq, fill=c_gene)) + 
                    geom_boxplot() +
                    facet_wrap(~ subject_id + sample_id, ncol=3, scales = "free") +
                    theme_enchantr() +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))

    mufreqSample <- eeplot(mufreqSample, 
                        outdir=params$outdir, 
                        file='mutation_frequency_sample',
                        caption=paste0("Mutation frequency per C gene per `sample_id`."),
                        mufreqSample)
    print(mufreqSample +
        theme(legend.position = "none"))
    caption <- mufreqSample$enchantr$html_caption  
}
```


``` r
db_locus <- db %>% filter(!is.na(locus), !is.na(mu_freq))
db_locus$mu_freq <- as.numeric(db_locus$mu_freq)

if (nrow(db_locus) > 0) {
    mufreqSampleLocus <- ggplot(db_locus, aes(x=locus, y=mu_freq, fill=locus)) + 
                    geom_boxplot() +
                    facet_wrap(~ subject_id + sample_id, ncol=3, scales = "free") +
                    theme_enchantr() +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))

    mufreqSampleLocus <- eeplot(mufreqSampleLocus, 
                        outdir=params$outdir, 
                        file='mutation_frequency_sample_locus',
                        caption=paste0("Mutation frequency per locus per `sample_id`."),
                        mufreqSampleLocus)
    print(mufreqSampleLocus +
        theme(legend.position = "none"))
    caption <- mufreqSampleLocus$enchantr$html_caption 
}
```

# Final repertoires and tables

Summary tables:

- [Summary of number of clones](tables/num_clones_table.tsv)
- [Clone sizes table](tables/clone_sizes_table.tsv)
- [Clonal abundance](tables/clonal_abundance.tsv)
- [Clonal diversity](tables/clonal_diversity.tsv)

Find your final repertoires here:



- [P05_FNA_0_Y1](repertoires/P05_FNA_0_Y1_define-clones_clone-pass.tsv)
- [P05_FNA_12_Y1](repertoires/P05_FNA_12_Y1_define-clones_clone-pass.tsv)
- [P05_FNA_2_0_Y1](repertoires/P05_FNA_2_0_Y1_define-clones_clone-pass.tsv)
- [P05_FNA_2_12_Y1](repertoires/P05_FNA_2_12_Y1_define-clones_clone-pass.tsv)
- [P05_FNA_2_28_Y1](repertoires/P05_FNA_2_28_Y1_define-clones_clone-pass.tsv)
- [P05_FNA_2_5_Y1](repertoires/P05_FNA_2_5_Y1_define-clones_clone-pass.tsv)
- [P05_FNA_2_60_Y1](repertoires/P05_FNA_2_60_Y1_define-clones_clone-pass.tsv)
- [P05_FNA_3_12_Y1](repertoires/P05_FNA_3_12_Y1_define-clones_clone-pass.tsv)
- [P05_FNA_3_28_Y1](repertoires/P05_FNA_3_28_Y1_define-clones_clone-pass.tsv)
- [P05_FNA_5_Y1](repertoires/P05_FNA_5_Y1_define-clones_clone-pass.tsv)
- [P05_FNA_60_Y1](repertoires/P05_FNA_60_Y1_define-clones_clone-pass.tsv)




# Software versions


``` r
sessionInfo()
```

```
## R version 4.4.3 (2025-02-28 ucrt)
## Platform: x86_64-w64-mingw32/x64
## Running under: Windows 11 x64 (build 26100)
## 
## Matrix products: default
## 
## 
## locale:
## [1] LC_COLLATE=English_United States.utf8 
## [2] LC_CTYPE=English_United States.utf8   
## [3] LC_MONETARY=English_United States.utf8
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.utf8    
## 
## time zone: America/New_York
## tzcode source: internal
## 
## attached base packages:
## [1] grid      stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] plotly_4.10.4         ComplexHeatmap_2.18.0 enchantr_0.1.20      
##  [4] dowser_2.3            scoper_1.3.0          shazam_1.2.0         
##  [7] alakazam_1.3.0        ggplot2_3.5.1         airr_1.5.0           
## [10] tidyr_1.3.1           dplyr_1.1.4           DT_0.33              
## 
## loaded via a namespace (and not attached):
##   [1] RColorBrewer_1.1-3          rstudioapi_0.17.1          
##   [3] jsonlite_2.0.0              shape_1.4.6.1              
##   [5] magrittr_2.0.3              farver_2.1.2               
##   [7] rmarkdown_2.29              GlobalOptions_0.1.2        
##   [9] fs_1.6.5                    zlibbioc_1.52.0            
##  [11] vctrs_0.6.5                 memoise_2.0.1              
##  [13] Rsamtools_2.22.0            ggtree_3.14.0              
##  [15] htmltools_0.5.8.1           S4Arrays_1.6.0             
##  [17] progress_1.2.3              curl_6.2.2                 
##  [19] SparseArray_1.6.2           gridGraphics_0.5-1         
##  [21] sass_0.4.9                  KernSmooth_2.23-26         
##  [23] bslib_0.9.0                 htmlwidgets_1.6.4          
##  [25] cachem_1.1.0                GenomicAlignments_1.42.0   
##  [27] igraph_2.1.4                lifecycle_1.0.4            
##  [29] iterators_1.0.14            pkgconfig_2.0.3            
##  [31] Matrix_1.7-2                R6_2.6.1                   
##  [33] fastmap_1.2.0               GenomeInfoDbData_1.2.13    
##  [35] MatrixGenerics_1.18.1       clue_0.3-66                
##  [37] digest_0.6.37               aplot_0.2.5                
##  [39] colorspace_2.1-1            patchwork_1.3.0            
##  [41] S4Vectors_0.44.0            crosstalk_1.2.1            
##  [43] GenomicRanges_1.58.0        labeling_0.4.3             
##  [45] phylotate_1.3               httr_1.4.7                 
##  [47] polyclip_1.10-7             abind_1.4-8                
##  [49] compiler_4.4.3              bit64_4.6.0-1              
##  [51] withr_3.0.2                 doParallel_1.0.17          
##  [53] BiocParallel_1.40.0         viridis_0.6.5              
##  [55] ggforce_0.4.2               MASS_7.3-64                
##  [57] DelayedArray_0.32.0         rjson_0.2.23               
##  [59] tools_4.4.3                 ape_5.8-1                  
##  [61] quadprog_1.5-8              glue_1.8.0                 
##  [63] nlme_3.1-167                cluster_2.1.8              
##  [65] ade4_1.7-23                 generics_0.1.3             
##  [67] seqinr_4.2-36               gtable_0.3.6               
##  [69] tzdb_0.5.0                  data.table_1.17.0          
##  [71] hms_1.1.3                   tidygraph_1.3.1            
##  [73] XVector_0.46.0              BiocGenerics_0.52.0        
##  [75] ggrepel_0.9.6               foreach_1.5.2              
##  [77] pillar_1.10.1               markdown_2.0               
##  [79] stringr_1.5.1               vroom_1.6.5                
##  [81] yulab.utils_0.2.0           circlize_0.4.16            
##  [83] tweenr_2.0.3                treeio_1.30.0              
##  [85] lattice_0.22-6              bit_4.6.0                  
##  [87] tidyselect_1.2.1            Biostrings_2.74.1          
##  [89] knitr_1.50                  gridExtra_2.3              
##  [91] bookdown_0.42               IRanges_2.40.1             
##  [93] SummarizedExperiment_1.36.0 stats4_4.4.3               
##  [95] xfun_0.51                   graphlayouts_1.2.2         
##  [97] Biobase_2.66.0              diptest_0.77-1             
##  [99] matrixStats_1.5.0           stringi_1.8.7              
## [101] UCSC.utils_1.2.0            lazyeval_0.2.2             
## [103] ggfun_0.1.8                 yaml_2.3.10                
## [105] evaluate_1.0.3              codetools_0.2-20           
## [107] ggraph_2.2.1                tibble_3.2.1               
## [109] ggplotify_0.1.2             cli_3.6.4                  
## [111] munsell_0.5.1               jquerylib_0.1.4            
## [113] Rcpp_1.0.14                 GenomeInfoDb_1.42.3        
## [115] png_0.1-8                   parallel_4.4.3             
## [117] readr_2.1.5                 prettyunits_1.2.0          
## [119] bitops_1.0-9                phangorn_2.12.1            
## [121] viridisLite_0.4.2           tidytree_0.4.6             
## [123] scales_1.3.0                purrr_1.0.4                
## [125] crayon_1.5.3                GetoptLong_1.0.5           
## [127] rlang_1.1.5                 fastmatch_1.1-6
```


<!--chapter:end:index.Rmd-->

