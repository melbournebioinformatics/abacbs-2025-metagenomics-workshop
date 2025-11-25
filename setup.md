---
title: Setup
---

Please follow the steps below and install the required software **before** the scheduled workshop.

## Datasets

Download the [Franzosa_IDB_2019 RDS data file](../episodes/data/Franzosa_IBD_2019.RDS) for the Multi-Omics block.

You'll load it into R later.

## RStudio Setup

We use RStudio for coding in R.

[Click here and follow the instructions](https://posit.co/download/rstudio-desktop/) to install RStudio Desktop in your system.

### R packages

This workshop requires two key R packages:

* **TaxSEA** (for taxon-set enrichment analysis)
* **mixOmics** (for multi-omics integration and DIABLO)

Please install them **before the workshop** following the instructions below.

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("TaxSEA", "mixOmics"))
```

If the installation completes without error, you can load the packages with:

```r
library(TaxSEA)
library(mixOmics)
```

