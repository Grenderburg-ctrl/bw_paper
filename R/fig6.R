# Fig. 6 — force-transmission scaffold (ComplexHeatmap, 17 × ~10 cm).
# Rows: scaffold-module genes in three blocks (Caveolae / Tethers / Basement membrane),
# fixed order, no clustering. Values: row z-score of the same VST as the PCA (Fig. 3A).
# Right annotation: DESeq2 log2FC bar (Down/Up/NS) and padj stars.

if (!exists("ROOT")) ROOT <- normalizePath(".")
source(file.path(ROOT, "R", "helpers.R"))
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
})
ht_opt(message = FALSE)

d <- get_data()
res <- d$res
sc <- subset(d$panel, module == "scaffold")
blocks <- c("Caveolae", "Tethers", "Basement membrane")
sc$family <- factor(sc$family, levels = blocks)
genes <- sc$gene                       # order as listed in the specification
stopifnot(all(genes %in% rownames(d$vsd)))

vm <- SummarizedExperiment::assay(d$vsd)[genes, SAMPLES]
z  <- t(scale(t(vm)))
zc <- pmin(pmax(z, -2), 2)

rr <- res[genes, c("baseMean", "log2FoldChange", "lfcSE", "padj", "direction")]
rr$stars <- p_stars(rr$padj)
write_src(data.frame(gene = genes, block = sc$family, z = round(z, 4), vst = round(vm, 4), rr,
                     check.names = FALSE, row.names = NULL), "Fig6")

col_fun <- colorRamp2(c(-2, 0, 2), heat_cols)
gp7 <- gpar(fontsize = 7, fontfamily = "Arial")

top <- HeatmapAnnotation(
  Group = GROUP[SAMPLES],
  col = list(Group = grp_fill),
  gp = gpar(col = "black", lwd = 0.5),
  simple_anno_size = unit(2.5, "mm"),
  show_annotation_name = FALSE,
  annotation_legend_param = list(Group = list(title = "Group", title_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Arial"),
                                              labels_gp = gp7, border = "black",
                                              grid_width = unit(3, "mm"), grid_height = unit(3, "mm"))))

lfc_lim <- ceiling(max(abs(rr$log2FoldChange)))
right <- rowAnnotation(
  log2FC = anno_barplot(rr$log2FoldChange, baseline = 0, ylim = c(-lfc_lim, lfc_lim),
                        gp = gpar(fill = dir_col[rr$direction], col = NA),
                        bar_width = 0.7, width = unit(18, "mm"),
                        axis_param = list(gp = gp7, at = c(-2, 0, 2),
                                          labels = c("−2", "0", "2"),
                                          labels_rot = 0, side = "bottom")),
  padj = anno_text(rr$stars, gp = gp7, just = "left", location = unit(1, "mm"),
                   width = unit(8, "mm")),
  annotation_name_gp = gpar(fontsize = 7, fontfamily = "Arial"),
  annotation_label = list(log2FC = gt_render(paste0("log", md_sub("2"), "FC"), gp = gp7),
                          padj = gt_render(MD_PADJ, gp = gp7)),
  annotation_name_rot = 0,
  annotation_name_side = "top",
  gap = unit(1.5, "mm"))

ht <- Heatmap(zc, name = "Row z-score", col = col_fun,
              cluster_rows = FALSE, cluster_columns = FALSE,
              row_split = sc$family, row_gap = unit(1.5, "mm"),
              row_title_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Arial"),
              row_title_rot = 0, row_title_side = "left", row_names_side = "left",
              row_names_gp = gpar(fontsize = 7, fontface = "italic", fontfamily = "Arial"),
              column_names_gp = gp7, column_names_rot = 0, column_names_centered = TRUE,
              column_split = GROUP[SAMPLES], column_gap = unit(1, "mm"),
              column_title_gp = gpar(fontsize = 8, fontfamily = "Arial"),
              top_annotation = top, right_annotation = right,
              rect_gp = gpar(col = "white", lwd = 0.5),
              width = unit(8 * 9, "mm"),
              heatmap_legend_param = list(at = -2:2, labels = gsub("-", "−", -2:2),
                                          title_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Arial"),
                                          labels_gp = gp7, legend_height = unit(18, "mm"),
                                          grid_width = unit(3, "mm")))

# Single-panel figure: the panel tag is not drawn.
draw_fig6 <- function() {
  draw(ht, merge_legends = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right",
       padding = unit(c(2, 2, 2, 8), "mm"))
}

save_fig(draw_fig6, "Fig6", width_cm = 17, height_cm = 10)
.bw$fig6 <- list(genes = genes, rr = rr)
