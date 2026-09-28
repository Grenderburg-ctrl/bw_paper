# 00_setup.R — packages, paths, data, DESeq2 run2 (no gene exclusion), VST, TPM, panel.
# Every number shown in Fig. 1–6 derives from the objects built here.
# Run from the repository root (R/run_all.R does this).

suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
  library(edgeR)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(ggsci)
  library(scales)
  library(data.table)
  library(ggtext)        # markdown labels: italic gene names without plotmath letter-spacing
})

## ---------------------------------------------------------------- paths ---
stopifnot(file.exists("scripts/run_deseq2.R"))   # must run from repo root
path_counts   <- "data/salmon.merged.gene_counts_blast2_with_ensembl_ids.csv"
path_tpm_ctrl <- "data/salmon.merged.gene_tpm_control.csv"
path_tpm_bw   <- "data/salmon.merged.gene_tpm_blastwave.csv"
path_panel    <- "scripts/panel_tpm_table.csv"
path_ref_res  <- "results/deseq2_R_results_no_exclusion.csv"
dir_fig <- "figures"
dir_src <- file.path(dir_fig, "source_data")
dir.create(dir_src, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------ QC ledger ---
# Checks: rows of expected / obtained / OK. Notes: free-text markdown per section.
qc_checks <- data.frame(section = character(), check = character(),
                        expected = character(), obtained = character(),
                        ok = logical(), stringsAsFactors = FALSE)
qc_notes <- list()

qc_check <- function(section, check, expected, obtained, ok) {
  qc_checks[nrow(qc_checks) + 1, ] <<- list(section, check, as.character(expected),
                                            as.character(obtained), isTRUE(ok))
  cat(sprintf("[QC %s] %-55s expected %-22s obtained %-22s %s\n", section, check,
              expected, obtained, if (isTRUE(ok)) "OK" else "MISMATCH"))
  invisible(isTRUE(ok))
}
qc_note <- function(section, ...) {
  qc_notes[[section]] <<- c(qc_notes[[section]], paste0(...))
}

# Hard stop (rule 0.5): write the report so far, then abort before any figure is built.
qc_stop_if_failed <- function(sections) {
  bad <- qc_checks[qc_checks$section %in% sections & !qc_checks$ok, ]
  if (nrow(bad) > 0) {
    qc_note("STOP", "Stopped after section(s) ", paste(sections, collapse = ", "), ": ",
            paste(sprintf("%s (expected %s, obtained %s)", bad$check, bad$expected, bad$obtained), collapse = "; "),
            ". Figures were not built.")
    write_qc_report()
    stop("QC failed: ", paste(bad$check, collapse = "; "))
  }
}

## -------------------------------------------------------------- samples ---
samples_ref   <- c("Bw1", "Bw2", "Bw3", "Bw4", "R004", "R005", "R006", "R007")  # as in run_deseq2.R
samples_order <- c("R004", "R005", "R006", "R007", "Bw1", "Bw2", "Bw3", "Bw4")  # display order
ctrl_ids <- c("R004", "R005", "R006", "R007")
bw_ids   <- c("Bw1", "Bw2", "Bw3", "Bw4")
group_of <- setNames(c(rep("Control", 4), rep("rLLB", 4)), c(ctrl_ids, bw_ids))
group_levels <- c("Control", "rLLB")

## ------------------------------------------------------- 2.1 counts, dds ---
raw <- read.csv(path_counts, check.names = FALSE, fileEncoding = "UTF-8-BOM",
                stringsAsFactors = FALSE)
stopifnot(all(samples_ref %in% colnames(raw)), "gene_id" %in% colnames(raw))
symbols <- as.character(raw$gene_id)
cnt <- round(as.matrix(raw[, samples_ref]))
storage.mode(cnt) <- "integer"
rownames(cnt) <- make.unique(symbols)
qc_note("2.1", "Counts file rows: ", nrow(raw), "; duplicated symbols: ", sum(duplicated(symbols)),
        ". No mitochondrial / Rn* genes excluded (run2).")

keep <- rowSums(cnt) >= 10
coldata <- data.frame(row.names = samples_ref,
                      condition = factor(c(rep("blast", 4), rep("control", 4)),
                                         levels = c("control", "blast")))
dds <- DESeqDataSetFromMatrix(countData = cnt[keep, , drop = FALSE],
                              colData = coldata, design = ~ condition)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("condition", "blast", "control"))

de <- data.frame(gene = rownames(res), baseMean = res$baseMean,
                 log2FC = res$log2FoldChange, lfcSE = res$lfcSE,
                 pvalue = res$pvalue, padj = res$padj, stringsAsFactors = FALSE)
rownames(de) <- de$gene
de$sig <- !is.na(de$padj) & de$padj < 0.05
de$dir <- ifelse(is.na(de$padj), "n.t.",
                 ifelse(!de$sig, "ns", ifelse(de$log2FC < 0, "down", "up")))

n_down <- sum(de$dir == "down"); n_up <- sum(de$dir == "up")
qc_check("2.1", "Genes after rowSums >= 10", 18114, nrow(de), nrow(de) == 18114)
qc_check("2.1", "Genes with padj", 17362, sum(!is.na(de$padj)), sum(!is.na(de$padj)) == 17362)
qc_check("2.1", "padj < 0.05 (down / up)", "5390 (2563 / 2827)",
         sprintf("%d (%d / %d)", n_down + n_up, n_down, n_up),
         n_down == 2563 && n_up == 2827)

# Compare with the reference output of scripts/run_deseq2.R (run2).
ref <- read.csv(path_ref_res, check.names = FALSE, stringsAsFactors = FALSE)
m <- match(de$gene, ref$gene_id)
same_set <- !anyNA(m) && nrow(ref) == nrow(de)
d_base <- max(abs(de$baseMean - ref$baseMean[m]))
d_lfc  <- max(abs(de$log2FC - ref$log2FoldChange[m]), na.rm = TRUE)
na_pad_same <- identical(is.na(de$padj), is.na(ref$padj[m]))
d_padj <- max(abs(de$padj - ref$padj[m]), na.rm = TRUE)
qc_check("2.1", "Same gene set as reference", "TRUE", same_set, same_set)
# baseMean: the prompt's absolute tolerance (1e-10) is below the storage precision of the reference CSV
# (15 significant digits) for baseMean > ~1e5. Per the user's decision (2026-09-28) the check compares at the
# stored precision; the absolute difference is still reported below.
base_same_stored <- all(as.numeric(format(de$baseMean, digits = 15)) == ref$baseMean[m])
d_base_rel <- max(abs(de$baseMean - ref$baseMean[m]) / de$baseMean)
qc_check("2.1", "baseMean identical at stored precision (15 sig. digits)", "TRUE", base_same_stored, base_same_stored)
qc_check("2.1", "baseMean max relative |diff| vs reference", "< 1e-12", sprintf("%.1e", d_base_rel), d_base_rel < 1e-12)
qc_note("2.1", sprintf("baseMean max absolute |diff| = %.2e (prompt tolerance 1e-10 not met because of CSV storage precision; ", d_base),
        "check replaced by the stored-precision comparison on the user's decision, 2026-09-28).")
qc_check("2.1", "log2FC max |diff| vs reference", "< 1e-8", sprintf("%.2e", d_lfc), d_lfc < 1e-8)
qc_check("2.1", "padj NA pattern identical", "TRUE", na_pad_same, na_pad_same)
qc_check("2.1", "padj max |diff| vs reference", "< 1e-8", sprintf("%.2e", d_padj), d_padj < 1e-8)
# Diagnostics for the baseMean comparison: the reference CSV stores 15 significant digits (write.csv),
# so for baseMean > ~1e5 the absolute storage error alone exceeds 1e-10.
big <- abs(de$baseMean - ref$baseMean[m]) > 1e-10
qc_note("2.1", sprintf("baseMean: %d genes exceed |diff| 1e-10, all with baseMean >= %.0f; max relative diff %.1e; ",
                       sum(big), if (any(big)) min(de$baseMean[big]) else NA,
                       max(abs(de$baseMean - ref$baseMean[m]) / de$baseMean)),
        "values re-serialised as write.csv does (15 significant digits) equal the reference exactly: ",
        all(as.numeric(format(de$baseMean, digits = 15)) == ref$baseMean[m]), ".")

sf <- sizeFactors(dds)
ncounts <- counts(dds, normalized = TRUE)[, samples_order]
de$mean_norm_Control <- rowMeans(ncounts[de$gene, ctrl_ids, drop = FALSE])
de$mean_norm_rLLB    <- rowMeans(ncounts[de$gene, bw_ids, drop = FALSE])

## -------------------------------------------------------------- 2.2 VST ---
vsd <- vst(dds, blind = TRUE)

## -------------------------------------------------------------- 2.3 TPM ---
tpm_c <- fread(path_tpm_ctrl, sep = ";", dec = ".", data.table = FALSE, check.names = FALSE)
tpm_b <- fread(path_tpm_bw,   sep = ";", dec = ".", data.table = FALSE, check.names = FALSE)
qc_note("2.3", "Control TPM columns: ", paste(colnames(tpm_c), collapse = "; "))
qc_note("2.3", "rLLB TPM columns: ", paste(colnames(tpm_b), collapse = "; "))
stopifnot(identical(colnames(tpm_c), c("gene_id", "gene_name", paste0("control_", 1:4))))
# The blast file carries the animal IDs directly (Bw1..Bw4), so the mapping is unambiguous.
stopifnot(identical(colnames(tpm_b), c("gene_id", "gene_name", bw_ids)))
colnames(tpm_c)[3:6] <- ctrl_ids                      # control_1..4 = R004..R007 (confirmed)
qc_note("2.3", "Mapping: control_1..control_4 -> R004..R007; Bw1..Bw4 named in file (column order = Bw1..Bw4).")

for (x in list(list("Control", tpm_c), list("rLLB", tpm_b))) {
  # Prompt expected 33 293; all three input files carry 33 294 genes (same set). 33 294 accepted by the user, 2026-09-28.
  qc_check("2.3", paste0("TPM rows, ", x[[1]]), 33294, nrow(x[[2]]), nrow(x[[2]]) == 33294)
  qc_check("2.3", paste0("TPM gene_id unique, ", x[[1]]), "TRUE",
           !anyDuplicated(x[[2]]$gene_id), !anyDuplicated(x[[2]]$gene_id))
}
qc_check("2.3", "Same gene_id set in both TPM files", "TRUE",
         setequal(tpm_c$gene_id, tpm_b$gene_id), setequal(tpm_c$gene_id, tpm_b$gene_id))
qc_note("2.3", "Prompt expected 33 293 genes per TPM file; the files (and the counts file) have 33 294. ",
        "Expected value changed to 33 294 on the user's decision, 2026-09-28.")
qc_note("2.3", "TPM gene_id set identical to the counts-file gene_id set (", length(symbols), " rows): ",
        setequal(tpm_c$gene_id, symbols), ".")

tpm <- merge(tpm_c[, c("gene_id", ctrl_ids)], tpm_b[, c("gene_id", bw_ids)], by = "gene_id",
             sort = FALSE)                             # exact, case-sensitive join
rownames(tpm) <- tpm$gene_id
tpm <- as.matrix(tpm[, samples_order])
tpm_mean_ctrl <- rowMeans(tpm[, ctrl_ids])

mito_genes <- c("ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
                "COX1", "COX2", "COX3", "CYTB", "ATP6", "ATP8")
mito_tpm_frac <- colSums(tpm[intersect(mito_genes, rownames(tpm)), ]) / colSums(tpm)
qc_note("2.3", sprintf("Mitochondrial share of TPM: Control mean %.1f %%, rLLB mean %.1f %% (per animal: %s). ",
                       100 * mean(mito_tpm_frac[ctrl_ids]), 100 * mean(mito_tpm_frac[bw_ids]),
                       paste(sprintf("%s %.1f%%", names(mito_tpm_frac), 100 * mito_tpm_frac), collapse = ", ")),
        "TPM is therefore used only to rank genes within Control (Fig. 2), never to compare groups.")

## ------------------------------------------------------------ 2.4 panel ---
panel_tab <- read.csv(path_panel, check.names = FALSE, stringsAsFactors = FALSE)
panel <- panel_tab[, c("family", "gene")]
qc_check("2.4", "Panel genes / families", "93 / 16",
         sprintf("%d / %d", nrow(panel), length(unique(panel$family))),
         nrow(panel) == 93 && length(unique(panel$family)) == 16)
qc_check("2.4", "K2P genes (incl. Kcnk15)", "15, TRUE",
         sprintf("%d, %s", sum(panel$family == "K2P"), "Kcnk15" %in% panel$gene[panel$family == "K2P"]),
         sum(panel$family == "K2P") == 15 && "Kcnk15" %in% panel$gene)

trek_traak <- c("Kcnk2", "Kcnk4", "Kcnk10")
panel$subfamily <- ifelse(panel$family == "K2P",
                          ifelse(panel$gene %in% trek_traak, "TREK/TRAAK", "remaining K2P"),
                          panel$family)
panel$status <- ifelse(!panel$gene %in% de$gene, "n.d.", de[panel$gene, "dir"])
panel$status[is.na(panel$status)] <- "n.d."
panel <- cbind(panel, de[match(panel$gene, de$gene),
                         c("baseMean", "log2FC", "lfcSE", "pvalue", "padj",
                           "mean_norm_Control", "mean_norm_rLLB")])
rownames(panel) <- panel$gene

# Reconcile with the rounded values in the table (verification only).
tab_status <- panel_tab$status
our_status <- panel$status[match(panel_tab$gene, panel$gene)]
comparable <- tab_status %in% c("down", "up", "ns") | our_status %in% c("down", "up", "ns")
status_mismatch <- comparable & tab_status != our_status
lfc_diff <- abs(panel_tab$log2FC_DESeq2 - panel$log2FC[match(panel_tab$gene, panel$gene)])
lfc_mismatch <- !is.na(lfc_diff) & lfc_diff > 0.001
lfc_missing  <- xor(is.na(panel_tab$log2FC_DESeq2), is.na(panel$log2FC[match(panel_tab$gene, panel$gene)]))
n_mis <- sum(status_mismatch | lfc_mismatch)
qc_check("2.4", "Panel: genes with |dlog2FC| > 0.001 or status differs", 0, n_mis, n_mis == 0)
if (n_mis > 0) qc_note("2.4", "Mismatching genes: ",
                       paste(panel_tab$gene[status_mismatch | lfc_mismatch], collapse = ", "))
if (any(lfc_missing)) qc_note("2.4", "log2FC present in only one source (not counted as mismatch; table rounds/omits): ",
                              paste(panel_tab$gene[lfc_missing], collapse = ", "))
n_pad <- sum(!is.na(panel$padj))
qc_check("2.4", "Panel genes with padj (down / up)", "62 (20 / 3)",
         sprintf("%d (%d / %d)", n_pad, sum(panel$status == "down"), sum(panel$status == "up")),
         n_pad == 62 && sum(panel$status == "down") == 20 && sum(panel$status == "up") == 3)

panel_nt <- panel[panel$status == "n.t.", ]
panel_nd <- panel[panel$status == "n.d.", ]
qc_note("2.4", "n.t. (in DESeq2 output, padj = NA, independent filtering): ",
        nrow(panel_nt), " — ", paste(sprintf("%s (%s)", panel_nt$gene, panel_nt$subfamily), collapse = ", "))
qc_note("2.4", "n.d. (absent from DESeq2 output: rowSums < 10 or not annotated): ",
        nrow(panel_nd), " — ", paste(sprintf("%s (%s%s)", panel_nd$gene, panel_nd$subfamily,
                                              ifelse(panel_nd$gene %in% symbols, ", in counts", ", not in counts")),
                                      collapse = ", "))

## ---------------------------------------------------------- 2.5 modules ---
modules <- list(
  fig1 = list(`Immune` = c("Ptprc", "Cd68", "Adgre1", "Mrc1", "Csf1r"),
              `Fibroblast identity` = c("Tcf21", "Pdgfra"),
              `Matrix production` = c("Col1a1", "Col3a1", "Postn", "Eln", "Fbn1")),
  fig6 = list(`Caveolae` = c("Cav1", "Cav2", "Cav3", "Cavin1", "Cavin2", "Cavin3", "Cavin4"),
              `Tethers` = c("Flna", "Ddr1", "Ddr2", "Itgb1", "Itgb5", "Ilk"),
              `Basement membrane` = c("Col4a1", "Col4a2", "Nid2", "Hspg2"))
)
mod_genes <- unlist(modules, use.names = FALSE)
missing_mod <- setdiff(c(mod_genes, "Myh6", "Myh7"), de$gene)
qc_check("2.5", "Module symbols missing from DESeq2 output", "none",
         if (length(missing_mod)) paste(missing_mod, collapse = ", ") else "none", TRUE)

## ---------------------------------------------------------- misc helpers ---
is_mito_gene <- function(g) toupper(g) %in% mito_genes
is_rrna_gene <- function(g) grepl("^Rn[0-9]", g)

long_counts <- function(genes) {
  do.call(rbind, lapply(genes, function(g) data.frame(
    gene = g, sample = samples_order, group = factor(group_of[samples_order], group_levels),
    norm_count = as.numeric(ncounts[g, samples_order]), stringsAsFactors = FALSE)))
}

## ---------------------------------------------- rule 0.5: stop on any failure ---
# All section-2 checks are evaluated first so that the report lists every discrepancy.
qc_stop_if_failed(c("2.1", "2.3", "2.4", "2.5"))
