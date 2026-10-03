# Fig. 1 — myosin heavy chain isoform switch. A: Myh7/Myh6 ratio; B: DESeq2 log2FC of Myh6 and Myh7

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
  labs(x = NULL, y = expression(log[2] ~ "fold change (rLLB vs. Control)")) +
  theme_ijms() + theme(axis.text.x = element_text(face = "italic"))

qc_check("Fig1", "Myh6 log2FC (rLLB vs. Control)", "-0.96 ± 0.01", sprintf("%.3f", myh_de["Myh6", "log2FC"]),
         abs(myh_de["Myh6", "log2FC"] - (-0.96)) < 0.01)
qc_check("Fig1", "Myh7 log2FC (rLLB vs. Control)", "+1.70 ± 0.01", sprintf("%.3f", myh_de["Myh7", "log2FC"]),
         abs(myh_de["Myh7", "log2FC"] - 1.70) < 0.01)
qc_check("Fig1", "Myh6 and Myh7 padj < 1e-4", "TRUE", all(myh_de$padj < 1e-4), all(myh_de$padj < 1e-4))

fig1 <- pA1 + pA2 + plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold", family = FONT))
save_fig(fig1, "Fig1", 8.5, 6)

write_src(myh[, c("sample", "group", "Myh6", "Myh7", "ratio")], "Fig1A_Myh7_Myh6_ratio")
write_src(data.frame(test = c("Welch t-test on ratio (plotted)", "Welch t-test on log10(ratio) (not plotted)"),
                     t = c(wt$statistic, wt_log$statistic), df = c(wt$parameter, wt_log$parameter),
                     p = c(wt$p.value, wt_log$p.value)), "Fig1A_welch_test")
write_src(data.frame(gene = rownames(myh_de), myh_de[, c("baseMean", "log2FC", "lfcSE", "pvalue", "padj")]),
          "Fig1B_Myh6_Myh7_deseq2")
