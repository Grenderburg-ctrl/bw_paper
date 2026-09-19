#!/usr/bin/env Rscript
# Verification of blast-vs-control differential expression with reference DESeq2.
# Protocol as specified: round -> drop mito/rRNA -> rowSums>=10 -> ~condition -> Wald, no shrinkage.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
counts_file <- if (length(args) >= 1) args[1] else "salmon.merged.gene_counts_blast2_with_ensembl_ids.csv"
outdir      <- if (length(args) >= 2) args[2] else "results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

log_lines <- c()
say <- function(...) {
  line <- paste0(...)
  cat(line, "\n", sep = "")
  log_lines <<- c(log_lines, line)
}

say("=== Session ===")
say("R:       ", R.version.string)
say("Bioc:    ", as.character(BiocManager::version()))
say("DESeq2:  ", as.character(packageVersion("DESeq2")))
say("")

## ---------------------------------------------------------------- 1. read ---
raw <- read.csv(counts_file, check.names = FALSE, fileEncoding = "UTF-8-BOM",
                stringsAsFactors = FALSE)
say("=== Input ===")
say("file:    ", counts_file)
say("rows:    ", nrow(raw))
say("columns: ", paste(colnames(raw), collapse = ", "))

samples <- c("Bw1", "Bw2", "Bw3", "Bw4", "R004", "R005", "R006", "R007")
stopifnot(all(samples %in% colnames(raw)))
stopifnot("gene_id" %in% colnames(raw))

symbols <- as.character(raw$gene_id)
mat_num <- as.matrix(raw[, samples])
storage.mode(mat_num) <- "double"

# Protocol step 1: round the salmon estimates to integers.
mat <- round(mat_num)
storage.mode(mat) <- "integer"
rownames(mat) <- make.unique(symbols)

dup_n <- sum(duplicated(symbols))
say("duplicate gene symbols: ", dup_n, if (dup_n > 0) "  (made unique with make.unique)" else "")
say("")

## ------------------------------------------------- 2. drop mito / rRNA genes ---
mito <- c("ND1","ND2","ND3","ND4","ND4L","ND5","ND6",
          "COX1","COX2","COX3","CYTB","ATP6","ATP8")
is_mito <- toupper(symbols) %in% toupper(mito)
is_rn   <- grepl("^Rn", symbols) & nchar(symbols) < 8
drop    <- is_mito | is_rn

# The literal "^Rn & nchar<8" rule also catches protein-coding families (Rnf*, Rnase*,
# Rnd*, ...). Flag the true rRNA/snRNA subset so run 3 can spare the collateral genes.
is_true_rrna <- is_rn & grepl("^Rn[0-9]", symbols)          # Rn18s, Rn28s, Rn45s, Rn5s, Rn5-8s, Rn7sl1
collateral   <- is_rn & !is_true_rrna
drop_strict  <- is_mito | is_true_rrna

say("=== Gene exclusion (step 2) ===")
say("mitochondrial removed: ", sum(is_mito), "  -> ", paste(sort(symbols[is_mito]), collapse = ", "))
say("Rn*/<8chr removed:     ", sum(is_rn))
say("  true rRNA/snRNA:     ", sum(is_true_rrna), "  -> ", paste(sort(symbols[is_true_rrna]), collapse = ", "))
say("  COLLATERAL (protein-coding caught by the rule): ", sum(collateral))
say("    ", paste(sort(symbols[collateral]), collapse = ", "))
say("total removed:         ", sum(drop))
say("")

write.csv(data.frame(gene_id = symbols[drop],
                     rule = ifelse(is_mito[drop], "mitochondrial", "Rn_prefix_lt8"),
                     mat[drop, , drop = FALSE], check.names = FALSE),
          file.path(outdir, "excluded_genes.csv"), row.names = FALSE)

## ------------------------------------------------------------- DESeq2 run ---
coldata <- data.frame(
  row.names = samples,
  condition = factor(c(rep("blast", 4), rep("control", 4)),
                     levels = c("control", "blast"))   # control is the reference
)

run_deseq <- function(count_mat, tag) {
  keep <- rowSums(count_mat) >= 10
  say("[", tag, "] genes before filter: ", nrow(count_mat))
  say("[", tag, "] genes after rowSums>=10: ", sum(keep))

  dds <- DESeqDataSetFromMatrix(countData = count_mat[keep, , drop = FALSE],
                                colData = coldata,
                                design = ~ condition)
  dds <- DESeq(dds)

  sf <- sizeFactors(dds)
  say("[", tag, "] size factors:")
  for (s in samples) say("    ", s, "  ", sprintf("%.4f", sf[s]))

  res <- results(dds, contrast = c("condition", "blast", "control"))
  say("[", tag, "] results() description:")
  for (l in mcols(res)$description) say("    ", l)

  sig <- !is.na(res$padj) & res$padj < 0.05
  say("[", tag, "] DEG padj<0.05: ", sum(sig),
      "  (down ", sum(sig & res$log2FoldChange < 0),
      ", up ", sum(sig & res$log2FoldChange > 0), ")")
  say("[", tag, "] padj NA (independent filtering / outliers): ", sum(is.na(res$padj)))
  say("")

  list(dds = dds, res = res, sf = sf, n_kept = sum(keep),
       n_deg = sum(sig),
       n_down = sum(sig & res$log2FoldChange < 0),
       n_up = sum(sig & res$log2FoldChange > 0))
}

say("=== RUN 1: with mito/rRNA exclusion (primary protocol) ===")
r1 <- run_deseq(mat[!drop, , drop = FALSE], "run1")

say("=== RUN 2: WITHOUT mito/rRNA exclusion (sensitivity) ===")
r2 <- run_deseq(mat, "run2")

say("=== RUN 3: mito + TRUE rRNA only (Rnf*/Rnase*/Rnd* spared) ===")
r3 <- run_deseq(mat[!drop_strict, , drop = FALSE], "run3")

## ------------------------------------------------------------- 1. results ---
res1 <- r1$res
out1 <- data.frame(gene_id = rownames(res1), as.data.frame(res1), check.names = FALSE)
out1 <- out1[order(out1$padj, na.last = TRUE), ]
write.csv(out1, file.path(outdir, "deseq2_R_results.csv"), row.names = FALSE)

res2 <- r2$res
out2 <- data.frame(gene_id = rownames(res2), as.data.frame(res2), check.names = FALSE)
out2 <- out2[order(out2$padj, na.last = TRUE), ]
write.csv(out2, file.path(outdir, "deseq2_R_results_no_exclusion.csv"), row.names = FALSE)

## --------------------------------------------------------- 2. panel compare ---
panel <- data.frame(
  gene = c("Tmc6","Mcoln2","Kcnk2","Piezo2","Ano8","Tmem63a","Pkd1","Ano6",
           "Piezo1","Tmem63b","Tmc7","Trpm4","Kcnk3",
           "Myh6","Myh7","Col3a1"),
  ref_lfc  = c(-1.618, 5.544, -1.188, -2.859, -0.630, -1.134, -1.234, -1.098,
               -1.473, -0.580, -0.968, -0.502, 0.207,
               -0.962, 1.701, -2.819),
  ref_padj = c(2.9e-20, 3.6e-08, 9.3e-08, 4.7e-07, 6.6e-06, 1.0e-05, 3.4e-05,
               4.8e-05, 7.7e-05, 1.0e-03, 3.2e-03, 4.0e-03, 0.53,
               4.8e-06, 6.3e-06, 9.5e-07),
  stringsAsFactors = FALSE
)

pull <- function(res, g) {
  if (g %in% rownames(res)) {
    list(lfc = res[g, "log2FoldChange"], padj = res[g, "padj"],
         base = res[g, "baseMean"], pval = res[g, "pvalue"])
  } else list(lfc = NA_real_, padj = NA_real_, base = NA_real_, pval = NA_real_)
}

res3 <- r3$res
cmp <- do.call(rbind, lapply(seq_len(nrow(panel)), function(i) {
  g <- panel$gene[i]
  a <- pull(res1, g); b <- pull(res2, g); cc <- pull(res3, g)
  ref_sig <- panel$ref_padj[i] < 0.05
  my_sig  <- !is.na(a$padj) && a$padj < 0.05
  data.frame(
    gene           = g,
    ref_log2FC     = panel$ref_lfc[i],
    R_log2FC       = a$lfc,
    d_log2FC       = a$lfc - panel$ref_lfc[i],
    lfc_ok         = !is.na(a$lfc) && abs(a$lfc - panel$ref_lfc[i]) < 0.15,
    ref_padj       = panel$ref_padj[i],
    R_padj         = a$padj,
    ref_sig        = ref_sig,
    R_sig          = my_sig,
    call_agrees    = ref_sig == my_sig,
    R_baseMean     = a$base,
    noexcl_log2FC  = b$lfc,
    noexcl_padj    = b$padj,
    noexcl_sig     = !is.na(b$padj) && b$padj < 0.05,
    strict_log2FC  = cc$lfc,
    strict_padj    = cc$padj,
    strict_sig     = !is.na(cc$padj) && cc$padj < 0.05,
    stringsAsFactors = FALSE
  )
}))
write.csv(cmp, file.path(outdir, "panel_comparison.csv"), row.names = FALSE)

say("=== Panel comparison (run 1 vs reference) ===")
print(cmp[, c("gene","ref_log2FC","R_log2FC","d_log2FC","lfc_ok",
              "ref_padj","R_padj","ref_sig","R_sig","call_agrees")])
say(paste(capture.output(print(cmp[, c("gene","ref_log2FC","R_log2FC","d_log2FC","lfc_ok",
              "ref_padj","R_padj","ref_sig","R_sig","call_agrees")])), collapse = "\n"))
say("")
say("=== Panel: run2 (no mito/rRNA exclusion) vs run3 (mito + true rRNA only) ===")
say(paste(capture.output(print(cmp[, c("gene","R_log2FC","R_padj",
              "noexcl_log2FC","noexcl_padj","noexcl_sig",
              "strict_log2FC","strict_padj","strict_sig")])), collapse = "\n"))
say("")
say("max |log2FC| shift run1->run2: ",
    sprintf("%.4f", max(abs(cmp$noexcl_log2FC - cmp$R_log2FC), na.rm = TRUE)))
say("max |log2FC| shift run1->run3: ",
    sprintf("%.4f", max(abs(cmp$strict_log2FC - cmp$R_log2FC), na.rm = TRUE)))
say("significance calls changed run1->run2: ", sum(cmp$R_sig != cmp$noexcl_sig))
say("significance calls changed run1->run3: ", sum(cmp$R_sig != cmp$strict_sig))
say("")
write.csv(data.frame(sample = samples,
                     run1 = as.numeric(r1$sf[samples]),
                     run2 = as.numeric(r2$sf[samples]),
                     run3 = as.numeric(r3$sf[samples])),
          file.path(outdir, "size_factors_all_runs.csv"), row.names = FALSE)

## -------------------------------------------- 1(add). contrast direction ---
say("=== Contrast direction check (raw counts) ===")
nf <- sweep(mat, 2, r1$sf, "/")   # normalised by run1 size factors
for (g in c("Mcoln2", "Col3a1")) {
  if (g %in% rownames(mat)) {
    rb <- mean(mat[g, samples[1:4]]);  rc <- mean(mat[g, samples[5:8]])
    nb <- mean(nf[g, samples[1:4]]);   nc <- mean(nf[g, samples[5:8]])
    say(sprintf("%-8s raw blast=%9.1f  raw control=%9.1f | norm blast=%9.1f  norm control=%9.1f | log2FC=%+.3f  -> %s",
                g, rb, rc, nb, nc, res1[g, "log2FoldChange"],
                ifelse(res1[g, "log2FoldChange"] > 0, "higher in blast", "higher in control")))
    say(sprintf("         per-sample raw: %s", paste(sprintf("%s=%.0f", samples, mat[g, samples]), collapse = "  ")))
  }
}
say("")

## ------------------------------------------------ size factor verification ---
ref_sf <- c(Bw1=1.277, Bw2=0.628, Bw3=1.098, Bw4=1.264,
            R004=1.403, R005=0.826, R006=0.975, R007=0.833)
sf_cmp <- data.frame(sample = samples,
                     ref = ref_sf[samples],
                     R = as.numeric(r1$sf[samples]),
                     pct_dev = 100 * (as.numeric(r1$sf[samples]) - ref_sf[samples]) / ref_sf[samples])
sf_cmp$ok <- abs(sf_cmp$pct_dev) < 2
write.csv(sf_cmp, file.path(outdir, "size_factor_comparison.csv"), row.names = FALSE)
say("=== Size factors vs reference ===")
say(paste(capture.output(print(sf_cmp)), collapse = "\n"))
say("")

## ------------------------------------------------------- global DEG counts ---
say("=== Global counts vs reference ===")
say(sprintf("genes after filter:  ref 18002  | R %d  (dev %+.2f%%)",
            r1$n_kept, 100*(r1$n_kept-18002)/18002))
say(sprintf("DEG padj<0.05:       ref 6371   | R %d  (dev %+.2f%%)",
            r1$n_deg, 100*(r1$n_deg-6371)/6371))
say(sprintf("  down:              ref 3131   | R %d", r1$n_down))
say(sprintf("  up:                ref 3240   | R %d", r1$n_up))
say(sprintf("run2 (no exclusion):        genes %d, DEG %d (down %d, up %d)",
            r2$n_kept, r2$n_deg, r2$n_down, r2$n_up))
say(sprintf("run3 (mito+true rRNA only): genes %d, DEG %d (down %d, up %d)",
            r3$n_kept, r3$n_deg, r3$n_down, r3$n_up))
say("")

res3o <- data.frame(gene_id = rownames(res3), as.data.frame(res3), check.names = FALSE)
write.csv(res3o[order(res3o$padj, na.last = TRUE), ],
          file.path(outdir, "deseq2_R_results_strict_exclusion.csv"), row.names = FALSE)

## ----------------------------------------------------------- 3. diagnostics ---
png(file.path(outdir, "dispersion_estimates.png"), width = 1400, height = 1100, res = 150)
plotDispEsts(r1$dds, main = "DESeq2 dispersion estimates (run 1)")
dev.off()

png(file.path(outdir, "MA_plot.png"), width = 1400, height = 1100, res = 150)
DESeq2::plotMA(res1, ylim = c(-8, 8), main = "MA-plot: blast vs control (MLE log2FC, padj<0.05)")
dev.off()

png(file.path(outdir, "pvalue_histogram.png"), width = 1400, height = 1100, res = 150)
hist(res1$pvalue[!is.na(res1$pvalue)], breaks = 50, col = "grey70", border = "white",
     main = "p-value distribution (run 1)", xlab = "raw p-value")
dev.off()

vsd <- vst(r1$dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, colour = condition, label = rownames(pca))) +
  geom_point(size = 4) + geom_text(vjust = -1.1, size = 3.5, show.legend = FALSE) +
  xlab(paste0("PC1: ", pv[1], "% variance")) +
  ylab(paste0("PC2: ", pv[2], "% variance")) +
  ggtitle("PCA on VST-transformed counts (run 1)") +
  coord_fixed() + theme_bw()
ggsave(file.path(outdir, "PCA_vst.png"), p, width = 8, height = 6.5, dpi = 150)

png(file.path(outdir, "dispersion_estimates_run2.png"), width = 1400, height = 1100, res = 150)
plotDispEsts(r2$dds, main = "DESeq2 dispersion estimates (run 2, no exclusion)")
dev.off()

writeLines(log_lines, file.path(outdir, "analysis_log.txt"))
say("Done. Outputs in: ", normalizePath(outdir))
