# Fig. 3 — global response (17 cm; A and B side by side).
#   A: PCA of VST (blind = TRUE), top 500 genes, equal axis scaling.
#   B: volcano of all genes with padj; channel-panel genes highlighted.
# Extra (not drawn): correlation of PC1 with library size, mito and rRNA fractions.

if (!exists("ROOT")) ROOT <- normalizePath(".")
source(file.path(ROOT, "R", "helpers.R"))
suppressPackageStartupMessages({
  library(patchwork)
  library(ggrepel)
  library(ggrastr)
})

d <- get_data()
res <- d$res
panel_ch <- subset(d$panel, module == "channel")

## ---------------------------------------------------------------- panel A ---
pca <- plotPCA(d$vsd, intgroup = "group", ntop = 500, returnData = TRUE)
pv  <- attr(pca, "percentVar")
pca$group <- factor(pca$group, levels = c("Control", "rLLB"))
pca$sample <- rownames(pca)

# Library composition per sample (raw rounded counts, all annotated genes).
lib  <- colSums(d$raw[, SAMPLES])
mito <- colSums(d$raw[is_mito(rownames(d$raw)), SAMPLES, drop = FALSE]) / lib
rrna <- colSums(d$raw[is_rrna(rownames(d$raw)), SAMPLES, drop = FALSE]) / lib
pca$lib_size <- lib[pca$sample]; pca$mito_frac <- mito[pca$sample]; pca$rrna_frac <- rrna[pca$sample]
write_src(pca[, c("sample", "group", "PC1", "PC2", "lib_size", "mito_frac", "rrna_frac")], "Fig3A")

pc1_cor <- do.call(rbind, lapply(c("lib_size", "mito_frac", "rrna_frac"), function(v) {
  ct <- cor.test(pca$PC1, pca[[v]], method = "pearson")
  cs <- suppressWarnings(cor.test(pca$PC1, pca[[v]], method = "spearman"))
  data.frame(variable = v, pearson_r = unname(ct$estimate), pearson_p = ct$p.value,
             spearman_rho = unname(cs$estimate), spearman_p = cs$p.value)
}))
write_src(pc1_cor, "Fig3A_PC1_covariates")

# Square plotting region with identical scale on both axes.
span <- max(diff(range(pca$PC1)), diff(range(pca$PC2))) * 1.2
ctr  <- c(mean(range(pca$PC1)), mean(range(pca$PC2)))
limx <- ctr[1] + c(-1, 1) * span / 2
limy <- ctr[2] + c(-1, 1) * span / 2

pA <- ggplot(pca, aes(PC1, PC2)) +
  geom_point(aes(shape = group, colour = group, fill = group), size = 2, stroke = PT_STROKE) +
  geom_text_repel(aes(label = sample), size = FONT_MM(6), family = "Arial",
                  segment.size = 0.18, segment.colour = "grey50", min.segment.length = 0.3,
                  box.padding = 0.25, point.padding = 0.15, seed = 1) +
  scale_shape_manual(values = c(Control =21, rLLB = 24), name = NULL) +
  scale_colour_manual(values = c(Control ="black", rLLB = "black"), guide = "none") +
  scale_fill_manual(values = pca_col, name = NULL) +
  scale_x_continuous(limits = limx, labels = lab_minus()) +
  scale_y_continuous(limits = limy, labels = lab_minus()) +
  coord_fixed(ratio = 1) +
  labs(x = sprintf("PC1 (%.1f%% variance)", 100 * pv[1]),
       y = sprintf("PC2 (%.1f%% variance)", 100 * pv[2]), tag = "A") +
  theme_ijms +
  theme(legend.position = "inside", legend.position.inside = c(0.99, 0.99),
        legend.justification = c(1, 1), legend.background = element_blank(),
        legend.text = element_text(size = 7))

## ---------------------------------------------------------------- panel B ---
XL <- 8
v <- subset(res, !is.na(padj))
v$y <- -log10(v$padj)
v$clamped <- abs(v$log2FoldChange) > XL
v$x <- pmax(pmin(v$log2FoldChange, XL), -XL)
v$is_panel <- v$gene_id %in% panel_ch$gene
v$is_mtr   <- is_mito(v$gene_id) | is_rrna(v$gene_id)
v$low_expr <- v$baseMean < 5
v$layer <- ifelse(v$is_panel, "panel", ifelse(v$is_mtr, "mito_rRNA", "background"))

label_genes <- c("Piezo1", "Piezo2", "Kcnk2", "Kcnk3", "Tmc6", "Tmem63a", "Pkd1", "Ano6",
                 "Trpm4", "Lrrc8a", "Lrrc8c", "Mcoln2", "Kcnk4")
v$label <- ifelse(v$gene_id %in% label_genes, v$gene_id, "")
write_src(v[, c("gene_id", "baseMean", "log2FoldChange", "padj", "direction", "x", "y",
                "clamped", "layer", "low_expr", "label")], "Fig3B")

n_down <- sum(v$direction == "Down"); n_up <- sum(v$direction == "Up")
bg  <- subset(v, layer == "background")
mtr <- subset(v, layer == "mito_rRNA")
pn  <- subset(v, layer == "panel")
pn$fill_k <- ifelse(pn$low_expr, "hollow", pn$direction)
pn$col_k  <- ifelse(pn$low_expr, pn$direction, "black")
y_max <- max(v$y) * 1.04

pB <- ggplot(mapping = aes(x, y)) +
  geom_hline(yintercept = -log10(0.05), linetype = "22", linewidth = 0.18, colour = "grey30") +
  rasterise(geom_point(data = subset(bg, !clamped), aes(colour = direction),
                       size = 0.5, alpha = 0.5, stroke = 0), dpi = 600) +
  rasterise(geom_point(data = subset(bg, clamped), aes(colour = direction),
                       shape = 17, size = 0.8, alpha = 0.5, stroke = 0), dpi = 600) +
  geom_point(data = mtr, shape = 23, size = 1.3, colour = "grey45", fill = NA, stroke = PT_STROKE) +
  geom_point(data = pn, aes(fill = fill_k, colour = col_k), shape = 21, size = 1.5, stroke = PT_STROKE) +
  # Labels are pushed up into the empty part of the plot; leader lines keep the link.
  geom_text_repel(data = subset(pn, label != ""), aes(label = label),
                  size = FONT_MM(7), family = "Arial", fontface = "italic",
                  segment.size = 0.18, segment.colour = "grey30", min.segment.length = 0,
                  box.padding = 0.3, point.padding = 0.2, max.overlaps = Inf,
                  nudge_y = 22, direction = "both", force = 2, force_pull = 0.5,
                  ylim = c(12, NA), seed = 7) +
  annotate("text", x = -XL, y = y_max, label = paste(n_down, "↓"), hjust = 0, vjust = 1,
           colour = dir_col["Down"], size = FONT_MM(8), family = "Arial") +
  annotate("text", x = XL, y = y_max, label = paste(n_up, "↑"), hjust = 1, vjust = 1,
           colour = dir_col["Up"], size = FONT_MM(8), family = "Arial") +
  scale_colour_manual(values = c(dir_col, black = "black"), guide = "none") +
  scale_fill_manual(values = c(dir_col, hollow = "white"), guide = "none") +
  scale_x_continuous(limits = c(-XL, XL), breaks = seq(-8, 8, 2), labels = lab_minus(),
                     expand = expansion(mult = 0.02)) +
  scale_y_continuous(limits = c(0, y_max), expand = expansion(mult = c(0, 0.01))) +
  labs(x = MD_LOG2FC, y = paste0("\u2212log", md_sub("10"), "(", MD_PADJ, ")"), tag = "B") +
  theme_ijms

fig3 <- (pA | pB) + plot_layout(widths = c(1, 1.35))
save_fig(fig3, "Fig3", width_cm = 17, height_cm = 7.5)

.bw$fig3 <- list(pv = pv, pca = pca, pc1_cor = pc1_cor, n_padj = nrow(v),
                 n_down = n_down, n_up = n_up, n_clamped = sum(v$clamped))
