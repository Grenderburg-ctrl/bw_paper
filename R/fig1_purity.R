# Fig. 1 — tissue composition of the preparation.
# A: Myh7/Myh6 ratio per animal (normalised counts) + DESeq2 log2FC of Myh6 and Myh7.
# B–C: log2FC ± lfcSE for immune, fibroblast-identity and matrix-production genes.

## ---------------------------------------------------------------- 1A ---
myh <- data.frame(sample = samples_order,
                  group = factor(group_of[samples_order], group_levels),
                  Myh6 = ncounts["Myh6", samples_order],
                  Myh7 = ncounts["Myh7", samples_order])
myh$ratio <- myh$Myh7 / myh$Myh6
wt <- t.test(ratio ~ group, data = myh, var.equal = FALSE)
wt_log <- t.test(log10(ratio) ~ group, data = myh, var.equal = FALSE)
qc_check("Fig1", "Welch t-test, Myh7/Myh6 ratio: p", "≈ 0.036", sprintf("%.4f", wt$p.value),
         abs(wt$p.value - 0.036) < 0.005)
qc_note("Fig1", sprintf("Myh7/Myh6 ratio (normalised counts): Control %.4f ± %.4f, rLLB %.4f ± %.4f (mean ± SEM); ",
                        mean(myh$ratio[myh$group == "Control"]), sd(myh$ratio[myh$group == "Control"]) / 2,
                        mean(myh$ratio[myh$group == "rLLB"]), sd(myh$ratio[myh$group == "rLLB"]) / 2),
        sprintf("Welch t = %.3f, df = %.2f, p = %.4f. For reference only (not plotted): Welch on log10(ratio) p = %.4f.",
                wt$statistic, wt$parameter, wt$p.value, wt_log$p.value))

sumA <- aggregate(ratio ~ group, myh, function(x) c(mean = mean(x), sem = sd(x) / sqrt(length(x))))
sumA <- data.frame(group = sumA$group, mean = sumA$ratio[, "mean"], sem = sumA$ratio[, "sem"])
yr <- range(c(myh$ratio, sumA$mean - sumA$sem, sumA$mean + sumA$sem))

pA1 <- ggplot(myh, aes(group, ratio)) +
  geom_errorbar(data = sumA, aes(x = group, ymin = mean - sem, ymax = mean + sem), inherit.aes = FALSE,
                width = 0.25, linewidth = lw_whisker) +
  geom_errorbar(data = sumA, aes(x = group, ymin = mean, ymax = mean), inherit.aes = FALSE,
                width = 0.45, linewidth = lw_whisker * 1.5) +
  geom_point(aes(fill = group), shape = 21, size = 1.6, stroke = 0.35, colour = "black",
             position = position_jitter(width = 0.12, height = 0, seed = 1)) +
  scale_fill_manual(values = c(Control = "#000000", rLLB = "#FFFFFF"), guide = "none") +
  scale_y_log10(labels = lab_log10, expand = expansion(mult = c(0.05, 0.25))) +
  annotate("segment", x = 1, xend = 2, y = yr[2] * 1.35, yend = yr[2] * 1.35, linewidth = lw_axis) +
  annotate("text", x = 1.5, y = yr[2] * 1.55, parse = TRUE, family = FONT, size = 7 / PT, vjust = 0,
           label = paste0("italic(p)~", fmt_p_math(wt$p.value, 2, prefix = "= "))) +
  labs(x = NULL, y = "*Myh7*/*Myh6* (normalised counts)") +
  theme_ijms() + theme(axis.title.y = element_markdown(size = 8))

myh_de <- de[c("Myh6", "Myh7"), ]
myh_de$gene <- factor(myh_de$gene, c("Myh6", "Myh7"))
myh_de$stars <- p_stars(myh_de$padj)
pA2 <- ggplot(myh_de, aes(gene, log2FC, colour = dir)) +
  geom_hline(yintercept = 0, linewidth = lw_axis, colour = "grey40") +
  geom_errorbar(aes(ymin = log2FC - lfcSE, ymax = log2FC + lfcSE), width = 0.2, linewidth = lw_whisker) +
  geom_point(size = 1.6) +
  geom_text(aes(y = ifelse(log2FC > 0, log2FC + lfcSE, log2FC - lfcSE), label = stars,
                vjust = ifelse(log2FC > 0, -0.4, 1.3)), size = 7 / PT, family = FONT, colour = "black") +
  scale_colour_manual(values = dir_cols, guide = "none") +
  scale_y_continuous(labels = lab_minus, limits = c(-1.6, 2.4), breaks = -1:2) +
  labs(x = NULL, y = expression(log[2] ~ "FC (rLLB vs. Control)")) +
  theme_ijms() + theme(axis.text.x = element_text(face = "italic", angle = 45, hjust = 1, vjust = 1))

## -------------------------------------------------------------- 1B–C ---
mod1 <- do.call(rbind, lapply(names(modules$fig1), function(b) {
  g <- modules$fig1[[b]]
  data.frame(block = b, gene = g, de[g, c("baseMean", "log2FC", "lfcSE", "pvalue", "padj", "dir",
                                          "mean_norm_Control", "mean_norm_rLLB")], stringsAsFactors = FALSE)
}))
mod1$stars <- p_stars(mod1$padj)
mod1$block <- factor(mod1$block, names(modules$fig1))
xlim1 <- max(abs(c(mod1$log2FC - mod1$lfcSE, mod1$log2FC + mod1$lfcSE))) * 1.05
xstar <- xlim1 * 1.12

dot_block <- function(blocks) {
  d <- mod1[mod1$block %in% blocks, ]
  d$gene <- factor(d$gene, rev(unlist(modules$fig1[blocks])))
  ggplot(d, aes(log2FC, gene, colour = dir)) +
    geom_vline(xintercept = 0, linewidth = lw_axis, colour = "grey40") +
    geom_errorbarh(aes(xmin = log2FC - lfcSE, xmax = log2FC + lfcSE), height = 0.25, linewidth = lw_whisker) +
    geom_point(size = 1.6) +
    geom_text(aes(x = xstar, label = stars), colour = "black", size = 7 / PT, family = FONT, hjust = 0) +
    facet_grid(block ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_colour_manual(values = dir_cols, guide = "none") +
    scale_x_continuous(labels = lab_minus, limits = c(-xlim1, xstar + 0.5), breaks = seq(-4, 4, 1)) +
    coord_cartesian(clip = "off") +
    labs(x = expression(log[2] ~ "fold change (rLLB vs. Control)"), y = NULL) +
    theme_ijms() +
    theme(axis.text.y = element_text(face = "italic"), strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 7),
          panel.spacing.y = unit(1.5, "mm"))
}
pB <- dot_block(c("Immune", "Fibroblast identity")) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.line.x = element_blank())
pC <- dot_block("Matrix production")

imm_ns <- all(mod1$dir[mod1$block %in% c("Immune", "Fibroblast identity")] == "ns")
mat_sig <- all(mod1$dir[mod1$block == "Matrix production"] %in% c("down", "up"))
qc_check("Fig1", "Immune + fibroblast-identity genes all ns", "TRUE", imm_ns, imm_ns)
qc_check("Fig1", "Matrix-production genes all padj < 0.05", "TRUE", mat_sig, mat_sig)
qc_note("Fig1", "Module genes: ", paste(sprintf("%s %s (log2FC %s, padj %s)", mod1$gene, mod1$dir,
                                                 fmt_num(mod1$log2FC), fmt_p_txt(mod1$padj)), collapse = "; "))

fig1 <- (pA1 + pA2 + plot_layout(widths = c(1.7, 1))) - (pB / pC + plot_layout(heights = c(7, 5))) +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(tag_levels = list(c("A", "", "B", "C"))) &
  theme(plot.tag = element_text(size = 10, face = "bold", family = FONT))
save_fig(fig1, "Fig1", 17, 7)

write_src(myh[, c("sample", "group", "Myh6", "Myh7", "ratio")], "Fig1A_Myh7_Myh6_ratio")
write_src(data.frame(test = c("Welch t-test on ratio (plotted)", "Welch t-test on log10(ratio) (not plotted)"),
                     t = c(wt$statistic, wt_log$statistic), df = c(wt$parameter, wt_log$parameter),
                     p = c(wt$p.value, wt_log$p.value)), "Fig1A_welch_test")
write_src(data.frame(gene = rownames(myh_de), myh_de[, c("baseMean", "log2FC", "lfcSE", "pvalue", "padj")]),
          "Fig1A_Myh6_Myh7_deseq2")
write_src(mod1, "Fig1BC_module_log2FC")
