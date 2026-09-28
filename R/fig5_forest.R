# Fig. 5 — forest plot of the 62 panel genes with padj: log2FC ± lfcSE (DESeq2, no shrinkage).

fam_order <- c("PIEZO", "TMEM63", "TMC", "TREK/TRAAK", "remaining K2P", "TRPC", "TRPV", "TRPM", "TRPML",
               "TRPP", "TMEM16", "ENaC/ASIC", "LRRC8 (VRAC)", "ClC", "Other")
fam5 <- function(sub) {
  x <- ifelse(sub %in% c("TREK/TRAAK", "remaining K2P"), sub,
              ifelse(sub == "LRRC8", "LRRC8 (VRAC)", sub))
  ifelse(x %in% fam_order, x, "Other")                  # KCNMA1, TMEM150C, TRPA
}
fam_strip <- c(`TREK/TRAAK` = "K2P:\nTREK/TRAAK", `remaining K2P` = "K2P:\nremaining")

f5 <- panel[!is.na(panel$padj), ]
f5$fam <- factor(fam5(f5$subfamily), fam_order)
f5$dir <- f5$status                                    # down / up / ns for all 62
f5 <- f5[order(f5$fam, f5$log2FC), ]
f5$gene_f <- factor(f5$gene, rev(f5$gene))
XL5 <- c(-3.5, 3)
f5$clip_hi <- f5$log2FC > XL5[2]; f5$clip_lo <- f5$log2FC < XL5[1]
f5$x <- pmin(pmax(f5$log2FC, XL5[1]), XL5[2])
f5$lo <- pmax(f5$log2FC - f5$lfcSE, XL5[1]); f5$hi <- pmin(f5$log2FC + f5$lfcSE, XL5[2])
f5$whisker_clipped <- (f5$log2FC - f5$lfcSE < XL5[1]) | (f5$log2FC + f5$lfcSE > XL5[2])
f5$low_expr <- f5$gene %in% c("Mcoln2", "Kcnk4")
f5$padj_lab <- fmt_p_math(f5$padj, 2)
f5$stars <- p_stars(f5$padj)
f5$row_lab <- ifelse(f5$low_expr, paste0(f5$gene, "\u2020"), f5$gene)
qc_check("Fig5", "Forest genes (panel with padj)", 62, nrow(f5), nrow(f5) == 62)
qc_note("Fig5", "Whiskers truncated at the axis limits (", XL5[1], ", ", XL5[2], "): ",
        paste(f5$gene[f5$whisker_clipped], collapse = ", "), ". Points beyond the limits (arrow + value): ",
        paste(sprintf("%s %s", f5$gene[f5$clip_hi | f5$clip_lo], fmt_num(f5$log2FC[f5$clip_hi | f5$clip_lo])),
              collapse = ", "), ".")

# Self-check key values.
key <- data.frame(gene = c("Tmc6", "Piezo1", "Piezo2", "Kcnk2", "Kcnk3", "Trpm4"),
                  lfc = c(-1.62, -1.47, -2.82, -1.19, 0.21, -0.50),
                  padj = c("1.5 × 10^-15", "1.4 × 10^-4", NA, NA, "0.56", "0.028"))
for (i in seq_len(nrow(key))) {
  g <- key$gene[i]
  ok_l <- round(de[g, "log2FC"], 2) == key$lfc[i]
  ok_p <- is.na(key$padj[i]) || fmt_p_txt(de[g, "padj"]) == key$padj[i]
  qc_check("Fig5", paste0("Key value ", g), sprintf("%s (%s)", fmt_num(key$lfc[i]), ifelse(is.na(key$padj[i]), "-", key$padj[i])),
           sprintf("%s (%s)", fmt_num(de[g, "log2FC"]), fmt_p_txt(de[g, "padj"])), ok_l && ok_p)
}
lrrc8e <- panel["Lrrc8e", "status"]
qc_note("Fig5", "Lrrc8e status: ", lrrc8e,
        if (lrrc8e %in% c("n.t.", "n.d.")) " -> marked 'not expressed'." else
          sprintf(" (padj %s) -> drawn as a regular row; no 'not expressed' mark needed.", fmt_p_txt(panel["Lrrc8e", "padj"])))

strip_lab <- function(x) ifelse(x %in% names(fam_strip), fam_strip[x], x)
base5 <- theme_ijms() +
  theme(axis.text.y = element_text(face = "italic"), strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 7),
        panel.spacing.y = unit(1.2, "mm"), strip.clip = "off")

pF <- ggplot(f5, aes(x, gene_f)) +
  geom_vline(xintercept = 0, linewidth = lw_axis, colour = "grey40") +
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = dir), height = 0.3, linewidth = lw_whisker) +
  geom_point(data = f5[!f5$low_expr, ], aes(colour = dir), size = 1.5) +
  geom_point(data = f5[f5$low_expr, ], aes(colour = dir), shape = 21, fill = "white", size = 1.5, stroke = 0.6) +
  geom_segment(data = f5[f5$clip_hi, ], aes(x = XL5[2] - 0.45, xend = XL5[2] + 0.05, yend = gene_f, colour = dir),
               arrow = arrow(length = unit(1, "mm"), type = "closed"), linewidth = lw_whisker) +
  geom_text(data = f5[f5$clip_hi, ], aes(x = XL5[2] - 0.5, label = paste0("+", formatC(log2FC, format = "f", digits = 2))),
            hjust = 1, size = 6 / PT, family = FONT) +
  facet_grid(fam ~ ., scales = "free_y", space = "free_y", switch = "y", labeller = as_labeller(strip_lab)) +
  scale_colour_manual(values = dir_cols, guide = "none") +
  scale_y_discrete(labels = function(g) f5$row_lab[match(g, f5$gene)]) +
  scale_x_continuous(limits = c(XL5[1], XL5[2] + 0.1), breaks = seq(-3, 3, 1), labels = lab_minus,
                     expand = c(0.01, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = expression(log[2] ~ "fold change (rLLB vs. Control) \u00B1 lfcSE"), y = NULL) +
  base5

pP <- ggplot(f5, aes(0, gene_f)) +
  geom_text(aes(label = padj_lab), parse = TRUE, hjust = 0, size = 6.5 / PT, family = FONT) +
  geom_text(aes(x = 1.35, label = stars), hjust = 0, size = 6.5 / PT, family = FONT) +
  facet_grid(fam ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0, 2), expand = c(0, 0), position = "top") +
  labs(x = "padj", y = NULL) +
  theme_void(base_family = FONT) +
  theme(strip.text = element_blank(), panel.spacing.y = unit(1.2, "mm"),
        axis.title.x.top = element_text(size = 7, face = "bold", hjust = 0.05, margin = margin(b = 1)))

# n.t. / n.d. listing by family (text only, no points).
miss <- panel[panel$status %in% c("n.t.", "n.d."), ]
miss$fam <- factor(fam5(miss$subfamily), fam_order)
miss_line <- function(st) {
  x <- miss[miss$status == st, ]; x <- x[order(x$fam, x$gene), ]
  parts <- vapply(split(x$gene, droplevels(x$fam)), function(g) paste0("*", g, "*", collapse = ", "), character(1))
  paste(sprintf("**%s:** %s", names(parts), parts), collapse = "; ")
}
txt5 <- paste0("**Not tested** (n.t., padj = NA after independent filtering) — ", miss_line("n.t."), ".<br>",
               "**Not detected** (n.d., absent from the DESeq2 output) — ", miss_line("n.d."), ".<br>",
               "† Low expression: near-absent in Control (mean normalised counts *Mcoln2* ",
               sprintf("%.1f", de["Mcoln2", "mean_norm_Control"]), ", *Kcnk4* ",
               sprintf("%.1f", de["Kcnk4", "mean_norm_Control"]), "); hollow markers, value beyond the axis shown by the arrow.")
pT <- ggplot() +
  geom_textbox(aes(x = 0, y = 1, label = txt5), width = unit(1, "npc"), hjust = 0, vjust = 1, halign = 0,
               box.colour = NA, fill = NA, size = 6.5 / PT, family = FONT, lineheight = 1.25,
               box.padding = unit(c(0, 0, 0, 0), "mm")) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void(base_family = FONT)

fig5 <- (pF + pP + plot_layout(widths = c(4.2, 1))) / pT + plot_layout(heights = c(1, 0.11))
save_fig(fig5, "Fig5", 17, 22)

write_src(data.frame(family = f5$fam, gene = f5$gene, baseMean = f5$baseMean, log2FC = f5$log2FC, lfcSE = f5$lfcSE,
                     pvalue = f5$pvalue, padj = f5$padj, padj_label = fmt_p_txt(f5$padj), stars = f5$stars,
                     status = f5$status, plotted_x = f5$x, whisker_truncated = f5$whisker_clipped,
                     beyond_axis = f5$clip_hi | f5$clip_lo, low_expression = f5$low_expr),
          "Fig5_forest")
write_src(data.frame(family = miss$fam, gene = miss$gene, status = miss$status,
                     log2FC = miss$log2FC, padj = miss$padj), "Fig5_not_tested_not_detected")
