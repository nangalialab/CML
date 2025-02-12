# Exports tables with estimates of timings of tree coalescences
source("export_paper_analysis_part_one.R")

## Produces figures
source("export_paper_analysis_part_two.R")

## Produces additional signature plots
# for comparison it requires previously generated growth rates with different prior range (cache/phylores.1e5.RDS))
source("export_paper_analysis_part_three.R")

## Run mixed models
setwd("../mixedmodels")
source("mixedmodel_markdown_generation.R")
