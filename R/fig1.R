# Fig. 1 — tissue composition (17 cm; A = one third, B–C = two thirds).
#   A: Myh7/Myh6 ratio per animal (normalized counts), bar = mean ± SEM, log10 axis,
#      Welch t-test; mini-panel: DESeq2 log2FC ± lfcSE for Myh6 and Myh7.
#   B–C: log2FC ± lfcSE dot plot for cell-type / matrix marker genes, stars = padj.

if (!exists("ROOT")) ROOT <- normalizePath(".")
source(file.path(ROOT, "R", "helpers.R"))
suppressPackageStartupMessages(library(patchwork))

d <- get_data()
res <- d$res

## ---------------------------------------------------------------- panel A ---
ratio <- data.frame(sample = SAMPLES, group = GROUP[SAMPLES],
                    Myh6 = d$ncounts["Myh6", SAMPLES], Myh7 = d$ncounts["Myh7", SAMPLES])
ratio$value <- ratio$Myh7 / ratio$Myh6
ratio_summ <- summ_mean_sem(transform(ratio, gene = "Myh7/Myh6"))
welch <- t.test(value ~ group, data = ratio, var.equal = FALSE)

write_src(transform(ratio, welch_t = unname(welch$statistic), welch_df = unname(welch$parameter),
                    welch_p = welch$p.value), "Fig1A")
write_src(ratio_summ, "Fig1A_summary")

y_top <- max(ratio$value, ratio_summ$mean + ratio_summ$sem)
br_y  <- y_top * 1.6
pA <- ggplot() +
  geom_bar_points(transform(ratio, gene = "r"), ratio_summ, seed = 11) +
  annotate("segment", x = 1, xend = 2, y = br_y, yend = br_y, linewidth = 0.18) +
  geom_md(data = data.frame(x = 1.5, y = br_y * 1.2, l = fmt_p_md(welch$p.value)),
          aes(x, y, label = l), vjust = 0, size = FONT_MM(7)) +
  scale_y_log10(limits = c(0.01, br_y * 2.6), expand = expansion(mult = c(0, 0)),
                breaks = c(0.01, 0.1, 1), labels = c("0.01", "0.1", "1")) +
  annotation_logticks(sides = "l", linewidth = 0.18, outside = TRUE,
                      short = unit(0.4, "mm"), mid = unit(0.7, "mm"), long = unit(0, "mm")) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "<i>Myh7</i>/<i>Myh6</i> ratio", tag = "A") +
  theme_ijms

# Mini-panel: DESeq2 log2FC ± lfcSE (rLLB vs. Sham).
myh <- res[c("Myh6", "Myh7"), c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "padj", "direction")]
myh$stars <- p_stars(myh$padj)
myh$gene_id <- factor(myh$gene_id, levels = c("Myh7", "Myh6"))
write_src(myh, "Fig1A_log2FC")

lim_m <- c(floor(min(myh$log2FoldChange - myh$lfcSE, 0)), ceiling(max(myh$log2FoldChange + myh$lfcSE, 0)))
pA2 <- ggplot(myh, aes(x = log2FoldChange, y = gene_id)) +
  geom_vline(xintercept = 0, linewidth = 0.18, colour = "grey40") +
  geom_errorbarh(aes(xmin = log2FoldChange - lfcSE, xmax = log2FoldChange + lfcSE),
                 height = 0.25, linewidth = ERR_LW) +
  geom_point(aes(fill = direction), shape = 21, size = PT_SIZE, stroke = PT_STROKE) +
  geom_text(aes(x = Inf, label = paste0("   ", stars)), hjust = 0, size = FONT_MM(7), family = "Arial") +
  scale_fill_manual(values = dir_col, guide = "none") +
  scale_x_continuous(labels = lab_minus(accuracy = 1), breaks = scales::breaks_width(2)) +
  scale_y_discrete(labels = function(x) parse(text = sprintf("italic(\"%s\")", x))) +
  coord_cartesian(xlim = lim_m, clip = "off") +
  labs(x = paste0("log", md_sub("2"), "FC"), y = NULL) +
  theme_ijms + theme(plot.margin = margin(5.5, 16, 5.5, 2))

## ------------------------------------------------------------- panels B–C ---
blocks <- list(
  "Immune"              = c("Ptprc", "Cd68", "Adgre1", "Mrc1", "Csf1r"),
  "Fibroblast identity" = c("Tcf21", "Pdgfra"),
  "Matrix production"   = c("Col1a1", "Col3a1", "Postn", "Eln", "Fbn1"))
mk <- do.call(rbind, lapply(names(blocks), function(b)
  data.frame(block = b, gene = blocks[[b]], stringsAsFactors = FALSE)))
mk$present <- mk$gene %in% rownames(res)
mk <- cbind(mk, res[match(mk$gene, rownames(res)), c("baseMean", "log2FoldChange", "lfcSE", "padj", "direction")])
mk$stars <- p_stars(mk$padj, mk$present)
mk$block <- factor(mk$block, levels = names(blocks))
write_src(mk, "Fig1BC")

x_lim <- c(floor(min(mk$log2FoldChange - mk$lfcSE, 0, na.rm = TRUE) * 2) / 2,
           ceiling(max(mk$log2FoldChange + mk$lfcSE, 0, na.rm = TRUE) * 2) / 2)

dot_plot <- function(dd, tag, show_x_title) {
  dd$gene <- factor(dd$gene, levels = rev(dd$gene))
  ggplot(dd, aes(x = log2FoldChange, y = gene)) +
    geom_vline(xintercept = 0, linewidth = 0.18, colour = "grey40") +
    geom_errorbarh(aes(xmin = log2FoldChange - lfcSE, xmax = log2FoldChange + lfcSE),
                   height = 0.3, linewidth = ERR_LW) +
    geom_point(aes(fill = direction), shape = 21, size = PT_SIZE, stroke = PT_STROKE) +
    facet_grid(block ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_manual(values = dir_col, guide = "none") +
    scale_x_continuous(limits = x_lim, breaks = scales::breaks_width(1),
                       labels = lab_minus(accuracy = 1)) +
    geom_text(aes(x = Inf, label = paste0("   ", stars)), hjust = 0, size = FONT_MM(7),
              family = "Arial") +
    scale_y_discrete(labels = function(x) parse(text = sprintf("italic(\"%s\")", x))) +
    coord_cartesian(clip = "off") +
    labs(x = if (show_x_title) MD_LOG2FC else NULL,
         y = NULL, tag = tag) +
    theme_ijms +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, face = "bold", size = 7, hjust = 1),
          plot.margin = margin(5.5, 16, 5.5, 5.5),
          panel.spacing.y = unit(2, "mm"))
}
pB <- dot_plot(subset(mk, block != "Matrix production"), "B", FALSE)
pC <- dot_plot(subset(mk, block == "Matrix production"), "C", TRUE)

## ---------------------------------------------------------------- layout ---
left  <- (pA | (plot_spacer() / pA2 + plot_layout(heights = c(1, 1.3)))) + plot_layout(widths = c(1.35, 1))
right <- (pB / pC) + plot_layout(heights = c(7, 5))
fig1 <- wrap_plots(left, right, widths = c(1, 2))

save_fig(fig1, "Fig1", width_cm = 17, height_cm = 6.5)
.bw$fig1 <- list(welch = welch, ratio_summ = ratio_summ)
