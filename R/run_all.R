#!/usr/bin/env Rscript
# Rebuild every figure, its source data and the QC report:
#   Rscript R/run_all.R            (from the repository root)

ROOT <- normalizePath(".")
if (!file.exists(file.path(ROOT, "R", "helpers.R")))
  stop("Run from the repository root: Rscript R/run_all.R")
if (file.exists(file.path(ROOT, "renv", "activate.R"))) source(file.path(ROOT, "renv", "activate.R"))

source(file.path(ROOT, "R", "helpers.R"))
.bw$fig_sizes <- list()

for (f in c("fig1", "fig3", "fig4", "fig5", "fig6")) {
  message("== ", f)
  source(file.path(ROOT, "R", paste0(f, ".R")))
}

message("== QC report")
source(file.path(ROOT, "R", "qc_report.R"))
message("done: figures/, source_data/, results/figure_qc.md")
