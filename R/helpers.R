# Shared helpers: data loading, significance labels, p-value formatting, export.
#
# Statistics come exclusively from the existing DESeq2 output
# (results/deseq2_R_results_no_exclusion.csv; DESeq2 1.42.1, R 4.3.3, Bioc 3.18).
# The DESeq2 object built here is used ONLY for normalized counts and VST;
# DESeq() / results() are never called.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(scales)
})

# Unicode labels (−, ×, ↑, ↓) need a UTF-8 locale for plotmath parsing.
invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "C.UTF-8")))

# Arial metrics for devices that consult the PostScript font database (ComplexHeatmap
# measures text on pdf(NULL)); Helvetica metrics are width-compatible with Arial.
if (!"Arial" %in% names(grDevices::pdfFonts()))
  grDevices::pdfFonts(Arial = grDevices::pdfFonts()$Helvetica)

ROOT <- if (exists("ROOT")) ROOT else normalizePath(".")
source(file.path(ROOT, "R", "theme_ijms.R"))

path_counts  <- file.path(ROOT, "data", "salmon.merged.gene_counts_blast2_with_ensembl_ids.csv")
path_deseq   <- file.path(ROOT, "results", "deseq2_R_results_no_exclusion.csv")
path_panel   <- file.path(ROOT, "data", "panel_table1.csv")
dir_fig      <- file.path(ROOT, "figures")
dir_src      <- file.path(ROOT, "source_data")
dir.create(dir_fig, showWarnings = FALSE)
dir.create(dir_src, showWarnings = FALSE)

SAMPLES <- c("R004", "R005", "R006", "R007", "Bw1", "Bw2", "Bw3", "Bw4")
GROUP   <- factor(ifelse(grepl("^R0", SAMPLES), "Control", "rLLB"), levels = c("Control", "rLLB"))
names(GROUP) <- SAMPLES

MITO_GENES <- c("ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
                "COX1", "COX2", "COX3", "ATP6", "ATP8", "CYTB")
is_mito <- function(g) g %in% MITO_GENES | grepl("^mt-", g)
is_rrna <- function(g) grepl("^Rn(5|5-8|18|28|45)s$", g)

## ------------------------------------------------------------ data loading ---
# Everything is computed once per R session and cached in .bw.
if (!exists(".bw")) .bw <- new.env()

get_data <- function() {
  if (!is.null(.bw$d)) return(.bw$d)

  raw <- read.csv(path_counts, check.names = FALSE, fileEncoding = "UTF-8-BOM",
                  stringsAsFactors = FALSE)
  mat <- round(as.matrix(raw[, SAMPLES]))
  storage.mode(mat) <- "integer"
  rownames(mat) <- raw$gene_id
  stopifnot(!anyDuplicated(raw$gene_id))

  keep <- rowSums(mat) >= 10
  cts  <- mat[keep, ]

  coldata <- data.frame(row.names = SAMPLES, group = GROUP)
  dds <- DESeqDataSetFromMatrix(cts, coldata, design = ~ group)
  dds <- estimateSizeFactors(dds)
  ncounts <- counts(dds, normalized = TRUE)

  res <- read.csv(path_deseq, stringsAsFactors = FALSE)
  rownames(res) <- res$gene_id
  res$direction <- ifelse(!is.na(res$padj) & res$padj < 0.05,
                          ifelse(res$log2FoldChange < 0, "Down", "Up"), "NS")

  # baseMean check: mean of normalized counts must equal DESeq2 baseMean.
  common <- intersect(rownames(res), rownames(ncounts))
  bm <- rowMeans(ncounts[common, ])
  bm_rel <- max(abs(bm - res[common, "baseMean"]) / pmax(res[common, "baseMean"], 1e-12))

  vsd <- vst(dds, blind = TRUE)

  panel <- read.csv(path_panel, stringsAsFactors = FALSE)

  .bw$d <- list(raw = mat, counts = cts, dds = dds, ncounts = ncounts, res = res,
                vsd = vsd, panel = panel, bm_rel = bm_rel,
                n_common = length(common), n_res = nrow(res))
  .bw$d
}

## ---------------------------------------------------- significance labels ---
# padj -> stars / status. `present` = gene is in the DESeq2 output.
p_stars <- function(padj, present = TRUE) {
  out <- ifelse(is.na(padj), "n.t.",
         ifelse(padj < 1e-4, "****",
         ifelse(padj < 1e-3, "***",
         ifelse(padj < 1e-2, "**",
         ifelse(padj < 0.05, "*", "ns")))))
  out[!present] <- "n.d."
  out
}

# Gene status relative to the DESeq2 output: "padj", "n.t." or "n.d.".
gene_status <- function(genes, res) {
  ifelse(!genes %in% rownames(res), "n.d.",
         ifelse(is.na(res[match(genes, rownames(res)), "padj"]), "n.t.", "padj"))
}

# Markdown (ggtext/gridtext) sub- and superscripts at an explicit 6 pt, so no text in a
# figure is smaller than 6 pt (plotmath would scale scripts to 0.7 × 7 pt = 4.9 pt).
SCRIPT_PT <- 6
md_sub <- function(x) sprintf("<sub><span style='font-size:%gpt'>%s</span></sub>", SCRIPT_PT, x)
md_sup <- function(x) sprintf("<sup><span style='font-size:%gpt'>%s</span></sup>", SCRIPT_PT, x)
md_it  <- function(x) sprintf("<i>%s</i>", x)
MD_PADJ   <- paste0(md_it("p"), md_sub("adj"))
MD_LOG2FC <- paste0("log", md_sub("2"), " fold change (rLLB vs. control)")

# p-value as a markdown label: "p = 0.036", "p = 1.8 × 10^−4" (Unicode minus/times).
fmt_p_md <- function(p, lab = md_it("p"), digits = 2, sci_below = 1e-3) {
  vapply(p, function(x) {
    if (is.na(x)) return(paste(lab, "= n.t."))
    if (x >= sci_below)
      return(sprintf("%s = %s", lab, formatC(signif(x, digits), format = "fg", digits = digits)))
    e <- floor(log10(x)); m <- signif(x / 10^e, digits)
    if (m >= 10) { m <- m / 10; e <- e + 1 }
    sprintf("%s = %s \u00d7 10%s", lab, formatC(m, format = "f", digits = digits - 1),
            md_sup(sub("-", "\u2212", as.character(e))))
  }, character(1))
}

# Text geom for markdown labels (no box).
geom_md <- function(...) ggtext::geom_richtext(..., fill = NA, label.colour = NA,
                                               label.padding = grid::unit(0, "pt"),
                                               family = "Arial")

# Plain-text p (for reports), "×10^" with Unicode minus.
fmt_p_txt <- function(x, digits = 2) {
  if (is.na(x)) return("NA")
  if (x >= 1e-3) return(formatC(signif(x, digits), format = "fg", digits = digits))
  e <- floor(log10(x)); m <- signif(x / 10^e, digits)
  if (m >= 10) { m <- m / 10; e <- e + 1 }
  sprintf("%s × 10^%s", formatC(m, format = "f", digits = digits - 1), sub("-", "−", e))
}

# Axis labels with a Unicode minus.
lab_minus <- function(...) {
  f <- scales::label_number(...)
  function(x) gsub("-", "−", f(x))
}

## ------------------------------------------------ bar + point summaries ---
# Long table of normalized counts for a set of genes, Control -> rLLB.
long_counts <- function(genes, d = get_data()) {
  genes <- genes[genes %in% rownames(d$ncounts)]
  do.call(rbind, lapply(genes, function(g) data.frame(
    gene = g, sample = SAMPLES, group = GROUP[SAMPLES],
    value = unname(d$ncounts[g, SAMPLES]))))
}

summ_mean_sem <- function(df, by = "gene") {
  keys <- unique(df[, c(by, "group")])
  out <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    sel <- Reduce(`&`, lapply(c(by, "group"), function(k) df[[k]] == keys[[k]][i]))
    v <- df$value[sel]
    cbind(keys[i, , drop = FALSE], n = length(v), mean = mean(v), sd = sd(v),
          sem = sd(v) / sqrt(length(v)))
  }))
  rownames(out) <- NULL
  out
}

# Bar (mean) + SEM whiskers + jittered animal points. `df` = long_counts(),
# `summ` = summ_mean_sem(df). Bars start at zero.
geom_bar_points <- function(df, summ, seed = 1) {
  list(
    geom_col(data = summ, aes(x = group, y = mean, fill = group),
             width = BAR_W, colour = "black", linewidth = 0.18, inherit.aes = FALSE),
    geom_errorbar(data = summ, aes(x = group, ymin = mean - sem, ymax = mean + sem),
                  width = CAP_W, linewidth = ERR_LW, inherit.aes = FALSE),
    geom_point(data = df, aes(x = group, y = value, fill = paste0(group, "_pt")), shape = 21,
               size = PT_SIZE, stroke = PT_STROKE, colour = "black",
               position = position_jitter(width = 0.15 * BAR_W, height = 0, seed = seed),
               inherit.aes = FALSE),
    scale_fill_manual(values = c(grp_fill, setNames(pt_fill, paste0(names(pt_fill), "_pt"))),
                      guide = "none")
  )
}

## ----------------------------------------------------------------- export ---
# Save a figure as vector PDF (cairo_pdf, Arial embedded) and 600-dpi LZW RGB TIFF.
# `plot` is a ggplot/patchwork object, or a function that draws on the current device.
save_fig <- function(plot, name, width_cm, height_cm) {
  w <- width_cm / 2.54; h <- height_cm / 2.54
  draw <- function() if (is.function(plot)) plot() else print(plot)

  pdf_file <- file.path(dir_fig, paste0(name, ".pdf"))
  grDevices::cairo_pdf(pdf_file, width = w, height = h, family = "Arial")
  draw(); invisible(dev.off())

  tif_file <- file.path(dir_fig, paste0(name, ".tiff"))
  ragg::agg_tiff(tif_file, width = w, height = h, units = "in", res = 600,
                 compression = "lzw", background = "white")
  draw(); invisible(dev.off())

  # Lightweight preview for review (not a deliverable).
  png_file <- file.path(dir_fig, "preview", paste0(name, ".png"))
  dir.create(dirname(png_file), showWarnings = FALSE)
  ragg::agg_png(png_file, width = w, height = h, units = "in", res = 200, background = "white")
  draw(); invisible(dev.off())

  .bw$fig_sizes[[name]] <- c(width_cm = width_cm, height_cm = height_cm)
  invisible(pdf_file)
}

write_src <- function(df, name) {
  write.csv(df, file.path(dir_src, paste0(name, ".csv")), row.names = FALSE)
}

# Minimum font size (pt) actually used in a figure — recorded by each script.
note_min_font <- function(name, pt) .bw$min_font[[name]] <- pt
