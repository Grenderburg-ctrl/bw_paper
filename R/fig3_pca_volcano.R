# Fig. 3 — global response. A: PCA on VST (top 500). B: volcano of all genes with padj.

## -------------------------------------------------------------- 3A ---
pca <- plotPCA(vsd, intgroup = "condition", ntop = 500, returnData = TRUE)
pv <- 100 * attr(pca, "percentVar")
pca$sample <- rownames(pca)
pca$group <- factor(group_of[pca$sample], group_levels)
qc_check("Fig3", "PCA PC1 % variance (ntop = 500)", "65.7 (± 0.5)", sprintf("%.1f", pv[1]),
         abs(pv[1] - 65.7) <= 0.5)
qc_stop_if_failed("Fig3")

# PC1 vs technical covariates (per sample).
cnt_all <- cnt[, samples_ref]
tech <- data.frame(sample = samples_ref,
                   size_factor = as.numeric(sf[samples_ref]),
                   mito_frac = colSums(cnt_all[is_mito_gene(rownames(cnt_all)), ]) / colSums(cnt_all),
                   rRNA_frac = colSums(cnt_all[is_rrna_gene(rownames(cnt_all)), ]) / colSums(cnt_all),
                   Rn_prefix_frac = colSums(cnt_all[grepl("^Rn", rownames(cnt_all)) & nchar(rownames(cnt_all)) < 8, ]) /
                     colSums(cnt_all),
                   lib_size = colSums(cnt_all))
tech$PC1 <- pca[tech$sample, "PC1"]; tech$PC2 <- pca[tech$sample, "PC2"]
tech$group <- group_of[tech$sample]
cor_tab <- do.call(rbind, lapply(c("size_factor", "mito_frac", "rRNA_frac", "Rn_prefix_frac", "lib_size"), function(v) {
  all <- cor.test(tech$PC1, tech[[v]], method = "pearson")
  within <- sapply(group_levels, function(g) cor(tech$PC1[tech$group == g], tech[[v]][tech$group == g]))
  data.frame(covariate = v, pearson_r_all8 = unname(all$estimate), p_all8 = all$p.value,
             spearman_rho_all8 = cor(tech$PC1, tech[[v]], method = "spearman"),
             r_within_Control = within[["Control"]], r_within_rLLB = within[["rLLB"]])
}))
tech_out <- tech[match(samples_order, tech$sample), ]

pca$lab <- pca$sample
pA <- ggplot(pca, aes(PC1, PC2)) +
  geom_point(aes(shape = group, colour = group, fill = group), size = 2.2, stroke = 0.3) +
  geom_text_repel(aes(label = lab), size = 6.5 / PT, family = FONT, min.segment.length = 0.3,
                  segment.size = 0.2, box.padding = 0.25, point.padding = 0.15, seed = 3) +
  scale_shape_manual(values = c(Control = 21, rLLB = 24), name = NULL) +
  scale_colour_manual(values = c(Control = pal$pca_ctrl, rLLB = pal$pca_bw), name = NULL) +
  scale_fill_manual(values = c(Control = pal$pca_ctrl, rLLB = pal$pca_bw), name = NULL) +
  scale_x_continuous(labels = lab_minus, expand = expansion(mult = 0.12)) +
  scale_y_continuous(labels = lab_minus, expand = expansion(mult = 0.15)) +
  coord_fixed() +
  labs(x = sprintf("PC1 (%.1f%% variance)", pv[1]), y = sprintf("PC2 (%.1f%% variance)", pv[2])) +
  theme_ijms() + theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0),
                       legend.box.spacing = unit(0.5, "mm"))

## -------------------------------------------------------------- 3B ---
XL <- 8
YL <- 30                                               # y cap; genes above drawn as triangles at the top edge
v <- de[!is.na(de$padj), c("gene", "baseMean", "log2FC", "lfcSE", "padj", "dir")]
v$mlog10 <- -log10(v$padj)
v$clipped <- abs(v$log2FC) > XL
v$x <- pmax(pmin(v$log2FC, XL), -XL)
v$clipped_y <- v$mlog10 > YL
v$clipped <- v$clipped | v$clipped_y
v$y <- pmin(v$mlog10, YL)
v$class <- ifelse(is_mito_gene(v$gene), "mito", ifelse(is_rrna_gene(v$gene), "rRNA", "nuclear"))
v$panel <- v$gene %in% panel$gene
v$low <- v$baseMean < 10
lab_genes <- c("Tmc6", "Piezo1", "Piezo2", "Kcnk2", "Kcnk3", "Pkd1", "Tmem63a", "Ano6", "Trpm4", "Mcoln2")
hollow <- c("Mcoln2", "Kcnk4")
qc_check("Fig3", "Volcano genes (padj not NA)", 17362, nrow(v), nrow(v) == 17362)
qc_note("Fig3", "Mean normalised counts (Control / rLLB): ",
        paste(sprintf("%s %.2f / %.2f", hollow, de[hollow, "mean_norm_Control"], de[hollow, "mean_norm_rLLB"]),
              collapse = "; "), " (hollow markers).")
qc_note("Fig3", "Beyond the axes (|log2FC| > 8 and/or −log10 padj > ", YL, "; drawn as edge triangles): ", sum(v$clipped), " genes (x: ", sum(abs(v$log2FC) > XL), ", y: ", sum(v$clipped_y), "); panel genes beyond the axes: ", sum(v$clipped & v$gene %in% panel$gene),
        " genes; mito genes plotted: ", sum(v$class == "mito"), "; rRNA genes plotted: ", sum(v$class == "rRNA"),
        "; genes with baseMean < 10 (semi-transparent): ", sum(v$low), ".")

bg <- v[!v$panel & v$class == "nuclear", ]
bg$alpha <- ifelse(bg$dir == "ns", 0.5, 0.8) * ifelse(bg$low, 0.4, 1)
bg_in <- bg[!bg$clipped, ]; bg_out <- bg[bg$clipped, ]
mt <- v[v$class != "nuclear", ]
pn <- v[v$panel, ]
pn$hollow <- pn$gene %in% hollow
ymax <- YL
# Fixed label positions (down genes stacked on the empty left side, ordered by y; up/ns on the right):
# deterministic, no overlaps, thin leader lines to the points.
labs3 <- pn[match(lab_genes, pn$gene), c("gene", "x", "y", "dir")]
left <- labs3$dir == "down"
labs3$lx <- ifelse(left, -4.6, 5.2)
labs3$hj <- ifelse(left, 1, 0)
ord <- order(-labs3$y[left])
labs3$ly[left][ord] <- seq(26, by = -2.6, length.out = sum(left))
labs3$ly[!left] <- c(Kcnk3 = 14, Mcoln2 = 10)[labs3$gene[!left]]

pB <- ggplot() +
  rasterise_if(geom_point(data = bg_in, aes(x, y, colour = dir, alpha = alpha),
                          size = 0.5, stroke = 0, shape = 16)) +
  geom_point(data = bg_out, aes(x, y, colour = dir, alpha = alpha),
             size = 0.9, stroke = 0, shape = 17) +
  geom_point(data = mt, aes(x, y), shape = 4, size = 1, stroke = 0.35, colour = "grey20") +
  geom_hline(yintercept = -log10(0.05), linetype = "22", linewidth = lw_axis, colour = "grey30") +
  geom_point(data = pn[!pn$hollow, ], aes(x, y, fill = dir), shape = 21, size = 1.5,
             stroke = 0.35, colour = "black") +
  geom_point(data = pn[pn$hollow, ], aes(x, y, colour = dir), shape = 21, fill = "white",
             size = 1.5, stroke = 0.6) +
  geom_segment(data = labs3, aes(x = lx, y = ly, xend = x, yend = y), linewidth = 0.2, colour = "grey30") +
  geom_label(data = labs3, aes(lx, ly, label = gene, hjust = hj), fontface = "italic", family = FONT,
             size = 7 / PT, label.size = 0, label.padding = unit(0.3, "mm"), fill = "white") +
  annotate("text", x = -XL, y = ymax * 1.02, hjust = 0, vjust = 1, label = sprintf("%d ↓", n_down),
           colour = pal$down, size = 7 / PT, family = FONT, fontface = "bold") +
  annotate("text", x = XL, y = ymax * 1.02, hjust = 1, vjust = 1, label = sprintf("%d ↑", n_up),
           colour = pal$up, size = 7 / PT, family = FONT, fontface = "bold") +
  scale_colour_manual(values = dir_cols, guide = "none") +
  scale_fill_manual(values = dir_cols, guide = "none") +
  scale_alpha_identity() +
  scale_x_continuous(limits = c(-XL, XL), breaks = seq(-8, 8, 4), labels = lab_minus) +
  scale_y_continuous(limits = c(0, YL), expand = expansion(mult = c(0.01, 0.04))) +
  labs(x = expression(log[2] ~ "fold change (rLLB vs. Control)"),
       y = expression("−" * log[10] * "(adjusted " * italic(p) * ")")) +
  theme_ijms()

fig3 <- pA + pB + plot_layout(widths = c(1, 1.45)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold", family = FONT))
save_fig(fig3, "Fig3", 17, 8)

write_src(cbind(pca[match(samples_order, pca$sample), c("sample", "group", "PC1", "PC2")],
                percentVar_PC1 = pv[1], percentVar_PC2 = pv[2]), "Fig3A_PCA")
write_src(tech_out, "Fig3A_technical_covariates")
write_src(cor_tab, "Fig3A_PC1_covariate_correlations")
v_out <- v[, c("gene", "baseMean", "log2FC", "lfcSE", "padj", "mlog10", "x", "y", "dir", "class", "panel", "low", "clipped")]
names(v_out)[names(v_out) %in% c("x", "y")] <- c("plotted_x", "plotted_y")
v_out$hollow <- v_out$gene %in% hollow
v_out$labelled <- v_out$gene %in% lab_genes
write_src(v_out, "Fig3B_volcano")
