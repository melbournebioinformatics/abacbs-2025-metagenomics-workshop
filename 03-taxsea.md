---
title: "Taxon-set enrichment analysis with TaxSEA"
teaching: 20
exercises: 0
---

:::::::::::::::::::::::::::::::::::::: questions 

- What is taxon-set enrichment analysis?
- How do we get started with the TaxSEA package?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Placeholder
- Placeholder

::::::::::::::::::::::::::::::::::::::::::::::::

## Loading the package






``` r
library(TaxSEA)

# Retrieve taxon sets containing Bifidobacterium longum.
blong.sets <- get_taxon_sets(taxon="Bifidobacterium_longum")

# Run TaxSEA with test data provided
data(TaxSEA_test_data)
taxsea_results <- TaxSEA(taxon_ranks=TaxSEA_test_data)
```

``` warning
Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
approximate in the presence of ties
Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
approximate in the presence of ties
Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
approximate in the presence of ties
Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
approximate in the presence of ties
Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
approximate in the presence of ties
Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
approximate in the presence of ties
Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
approximate in the presence of ties
Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
approximate in the presence of ties
Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
approximate in the presence of ties
Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
approximate in the presence of ties
```

``` r
# Enrichments among metabolite producers from gutMgene and MiMeDB
metabolites.df <- taxsea_results$Metabolite_producers

# Enrichments among health and disease signatures from GMRepoV2 and mBodyMap
disease.df <- taxsea_results$Health_associations

# Enrichments amongh published associations from BugSigDB
bsdb.df <- taxsea_results$BugSigdB
```





``` r
#### Applying TaxSEA functionality to taxonomic ranks  
# This script applies TaxSEA to identify taxonomic enrichment at different taxonomic levels.
# Specifically, we analyze enrichment at the family level using metagenomic data

# Load required libraries
library(TaxSEA)
library(curatedMetagenomicData)
```

``` output
Loading required package: SummarizedExperiment
```

``` output
Loading required package: MatrixGenerics
```

``` output
Loading required package: matrixStats
```

``` output

Attaching package: 'MatrixGenerics'
```

``` output
The following objects are masked from 'package:matrixStats':

    colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    colWeightedMeans, colWeightedMedians, colWeightedSds,
    colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    rowWeightedSds, rowWeightedVars
```

``` output
Loading required package: GenomicRanges
```

``` output
Loading required package: stats4
```

``` output
Loading required package: BiocGenerics
```

``` output
Loading required package: generics
```

``` output

Attaching package: 'generics'
```

``` output
The following objects are masked from 'package:base':

    as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
    setequal, union
```

``` output

Attaching package: 'BiocGenerics'
```

``` output
The following objects are masked from 'package:stats':

    IQR, mad, sd, var, xtabs
```

``` output
The following objects are masked from 'package:base':

    anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
    order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
    unsplit, which.max, which.min
```

``` output
Loading required package: S4Vectors
```

``` output

Attaching package: 'S4Vectors'
```

``` output
The following object is masked from 'package:utils':

    findMatches
```

``` output
The following objects are masked from 'package:base':

    expand.grid, I, unname
```

``` output
Loading required package: IRanges
```

``` output
Loading required package: Seqinfo
```

``` output
Loading required package: Biobase
```

``` output
Welcome to Bioconductor

    Vignettes contain introductory material; view with
    'browseVignettes()'. To cite Bioconductor, see
    'citation("Biobase")', and for packages 'citation("pkgname")'.
```

``` output

Attaching package: 'Biobase'
```

``` output
The following object is masked from 'package:MatrixGenerics':

    rowMedians
```

``` output
The following objects are masked from 'package:matrixStats':

    anyMissing, rowMedians
```

``` output
The following object is masked from 'package:httr':

    content
```

``` output
Loading required package: TreeSummarizedExperiment
```

``` output
Loading required package: SingleCellExperiment
```

``` output
Loading required package: Biostrings
```

``` output
Loading required package: XVector
```

``` output

Attaching package: 'Biostrings'
```

``` output
The following object is masked from 'package:base':

    strsplit
```

``` output
Registered S3 methods overwritten by 'readr':
  method                    from 
  as.data.frame.spec_tbl_df vroom
  as_tibble.spec_tbl_df     vroom
  format.col_spec           vroom
  print.col_spec            vroom
  print.collector           vroom
  print.date_names          vroom
  print.locale              vroom
  str.col_spec              vroom
```

``` r
library(tidyverse)
```

``` output
── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
✔ dplyr     1.1.4     ✔ readr     2.1.5
✔ forcats   1.0.1     ✔ stringr   1.5.2
✔ ggplot2   4.0.0     ✔ tibble    3.3.0
✔ lubridate 1.9.4     ✔ tidyr     1.3.1
✔ purrr     1.1.0     
```

``` output
── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
✖ lubridate::%within%() masks IRanges::%within%()
✖ dplyr::collapse()     masks Biostrings::collapse(), IRanges::collapse()
✖ dplyr::combine()      masks Biobase::combine(), BiocGenerics::combine()
✖ purrr::compact()      masks XVector::compact()
✖ Biobase::content()    masks httr::content()
✖ dplyr::count()        masks matrixStats::count()
✖ dplyr::desc()         masks IRanges::desc()
✖ tidyr::expand()       masks S4Vectors::expand()
✖ dplyr::filter()       masks stats::filter()
✖ dplyr::first()        masks S4Vectors::first()
✖ dplyr::lag()          masks stats::lag()
✖ ggplot2::Position()   masks BiocGenerics::Position(), base::Position()
✖ purrr::reduce()       masks GenomicRanges::reduce(), IRanges::reduce()
✖ dplyr::rename()       masks S4Vectors::rename()
✖ lubridate::second()   masks S4Vectors::second()
✖ lubridate::second<-() masks S4Vectors::second<-()
✖ dplyr::slice()        masks XVector::slice(), IRanges::slice()
ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```

``` r
library(phyloseq)
```

``` output

Attaching package: 'phyloseq'

The following object is masked from 'package:SummarizedExperiment':

    distance

The following object is masked from 'package:Biobase':

    sampleNames

The following object is masked from 'package:GenomicRanges':

    distance

The following object is masked from 'package:IRanges':

    distance
```

``` r
library(MicrobiomeStat)
```

``` output
Registered S3 method overwritten by 'rmutil':
  method         from
  print.response httr
```

``` r
library(dplyr)

# Load sample metadata
metadata_all <- sampleMetadata

# Filter metadata for the specific study (ChngKR_2016)
metadata <- metadata_all %>% 
  filter(study_name == "ChngKR_2016") %>% 
  column_to_rownames('sample_id')

# Extract count data using curatedMetagenomicData
cmd_data <- curatedMetagenomicData(
  pattern = "ChngKR_2016.relative_abundance",
  counts = TRUE,
  dryrun = FALSE
)
```

``` output

$`2021-03-31.ChngKR_2016.relative_abundance`
dropping rows without rowTree matches:
  k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Carnobacteriaceae|g__Granulicatella|s__Granulicatella_elegans
  k__Bacteria|p__Proteobacteria|c__Alphaproteobacteria|o__Caulobacterales|f__Caulobacteraceae|g__Brevundimonas|s__Brevundimonas_abyssalis
```

``` r
# Convert the extracted data to a count matrix
counts_data <- assay(cmd_data[[1]])
counts_data <- counts_data[, rownames(metadata)]  # Subset to relevant samples

# Filter taxa with at least one sample having counts > 100
counts_data <- counts_data[apply(counts_data > 100, 1, sum) > 0, ]

# Extract species names from taxonomic strings
species_names <- gsub("s__", "", sapply(rownames(counts_data), function(y) strsplit(y, "\\|")[[1]][7]))
rownames(counts_data) <- species_names

# Create a taxonomic lineage dataframe
# Remove taxonomic prefixes (k__, p__, c__, etc.) and separate into taxonomic ranks

# make data frame of taxon lineages
taxon_lineages <- data.frame(Name = species_names,
                             Lineage = names(species_names)) %>%
  mutate(Lineage = str_remove_all(Lineage, '[kpcofgs]__')) %>%
  separate(col = Lineage, into = c('kingdom', 'phylum', 'class', 
                                   'order', 'family', 'genus', 'species'), 
           sep = '\\|') %>%
  mutate(name = Name) %>%
  remove_rownames() %>%
  column_to_rownames('name')

# Perform differential abundance testing using LinDA
metadata$study_condition <- factor(metadata$study_condition, levels = c("control", "AD"))

linda_results <- linda(
  feature.dat = counts_data,
  meta.dat = metadata,
  formula = '~study_condition',
  feature.dat.type = 'count',
  prev.filter = 0.05
)
```

``` output
216  features are filtered!
The filtered data has  78  samples and  216  features will be tested!
Pseudo-count approach is used.
Fit linear models ...
Completed.
```

``` r
# Extract log2 fold change values for differential taxa
linda_results <- linda_results$output$study_conditionAD
log2_fold_changes <- linda_results$log2FoldChange
names(log2_fold_changes) <- rownames(linda_results)

# Define the taxonomic rank for enrichment analysis
selected_taxon_level <- 'genus'  # Modify as needed (e.g., genus, phylum)

# Create a named list of species grouped by taxonomic rank
custom_taxon_sets <- taxon_lineages %>%
  group_by(.data[[selected_taxon_level]]) %>% 
  summarise(species = list(species), .groups = "drop") %>%
  deframe()

# Perform enrichment analysis using TaxSEA
custom_taxsea_results <- TaxSEA(taxon_ranks = log2_fold_changes, custom_db = custom_taxon_sets)
custom_taxsea_results <- custom_taxsea_results$custom_sets
```

::::::::::::::::::::::::::::::::::::: callout

Callout sections can highlight information.

They are sometimes used to emphasise particularly important points
but are also used in some lessons to present "asides": 
content that is not central to the narrative of the lesson,
e.g. by providing the answer to a commonly-asked question.

::::::::::::::::::::::::::::::::::::::::::::::::

## Figures

You can also include figures generated from R Markdown:


``` r
pie(
  c(Sky = 78, "Sunny side of pyramid" = 17, "Shady side of pyramid" = 5), 
  init.angle = 315, 
  col = c("deepskyblue", "yellow", "yellow3"), 
  border = FALSE
)
```

<div class="figure" style="text-align: center">
<img src="fig/03-taxsea-rendered-pyramid-1.png" alt="pie chart illusion of a pyramid"  />
<p class="caption">Sun arise each and every morning</p>
</div>



## Math

One of our episodes contains $\LaTeX$ equations when describing how to create
dynamic reports with {knitr}, so we now use mathjax to describe this:

`$\alpha = \dfrac{1}{(1 - \beta)^2}$` becomes: $\alpha = \dfrac{1}{(1 - \beta)^2}$

Cool, right?

::::::::::::::::::::::::::::::::::::: keypoints 

- Use `.md` files for episodes when you want static content
- Use `.Rmd` files for episodes when you need to generate output
- Run `sandpaper::check_lesson()` to identify any issues with your lesson
- Run `sandpaper::build_lesson()` to preview your lesson locally

::::::::::::::::::::::::::::::::::::::::::::::::

[r-markdown]: https://rmarkdown.rstudio.com/
