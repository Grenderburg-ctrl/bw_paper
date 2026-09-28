# Fig. 6 — force-transmission module: row z-scores of VST values, 8 animals in fixed order.

g6 <- unlist(modules$fig6, use.names = FALSE)
blk6 <- factor(rep(names(modules$fig6), lengths(modules$fig6)), names(modules$fig6))
blk6_lab <- c(Caveolae = "Caveolae", Tethers = "Tethers", `Basement membrane` = "Basement\nmembrane")
vst_mat <- assay(vsd)[g6, samples_order]
z <- t(scale(t(vst_mat)))
z_clip <- pmin(pmax(z, -2), 2)
d6 <- de[g6, ]
d6$stars <- p_stars(d6$padj)

col_fun <- colorRamp2(c(-2, 0, 2), pal$heat)
top_ann <- HeatmapAnnotation(
  Group = factor(group_of[samples_order], group_levels),
  col = list(Group = c(Control = pal$ctrl_fill, rLLB = pal$bw_fill)),
  gp = gpar(col = "black", lwd = 0.5),
  simple_anno_size = unit(2.5, "mm"),
  annotation_name_gp = gpar(fontsize = 7, fontfamily = FONT),
  annotation_legend_param = list(Group = list(title = NULL, labels_gp = gpar(fontsize = 7, fontfamily = FONT),
                                              grid_height = unit(2.5, "mm"), grid_width = unit(2.5, "mm"),
                                              border = "black")))
lfc_col <- ifelse(d6$dir == "down", pal$down, ifelse(d6$dir == "up", pal$up, pal$ns))
lfc_lim <- ceiling(max(abs(d6$log2FC)) * 2) / 2
right_ann <- rowAnnotation(
  `log2FC` = anno_barplot(d6$log2FC, baseline = 0, gp = gpar(fill = lfc_col, col = NA),
                          ylim = c(-lfc_lim, lfc_lim), width = unit(12, "mm"),
                          axis_param = list(gp = gpar(fontsize = 6, fontfamily = FONT), side = "bottom",
                                            at = c(-lfc_lim, 0, lfc_lim),
                                            labels = lab_minus(c(-lfc_lim, 0, lfc_lim)))),
  padj = anno_text(d6$stars, gp = gpar(fontsize = 7, fontfamily = FONT), just = "left",
                   location = unit(1, "mm")),
  annotation_name_gp = gpar(fontsize = 7, fontfamily = FONT),
  annotation_label = list(log2FC = expression(log[2] * "FC"), padj = "padj"),
  annotation_name_rot = 0, gap = unit(1, "mm"))

ht <- Heatmap(z_clip, name = "z", col = col_fun,
              cluster_rows = FALSE, cluster_columns = FALSE,
              row_split = blk6, row_title = blk6_lab[levels(blk6)], row_title_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = FONT),
              row_title_rot = 0, row_title_side = "left", row_gap = unit(1.2, "mm"),
              row_names_gp = gpar(fontsize = 7, fontface = "italic", fontfamily = FONT),
              row_names_side = "left",
              column_names_gp = gpar(fontsize = 7, fontfamily = FONT), column_names_rot = 90,
              column_split = factor(group_of[samples_order], group_levels), column_title = NULL,
              column_gap = unit(0.8, "mm"),
              top_annotation = top_ann, right_annotation = right_ann,
              rect_gp = gpar(col = "white", lwd = 0.4),
              width = unit(8 * 4, "mm"), height = unit(length(g6) * 3.4, "mm"),
              heatmap_legend_param = list(title = "Row z-score (VST)", title_position = "topcenter", at = c(-2, -1, 0, 1, 2),
                                          labels = lab_minus(c(-2, -1, 0, 1, 2)),
                                          title_gp = gpar(fontsize = 7, fontfamily = FONT),
                                          labels_gp = gpar(fontsize = 7, fontfamily = FONT),
                                          legend_direction = "horizontal", legend_width = unit(18, "mm"),
                                          grid_height = unit(2.5, "mm")))

# Single-panel figure: no panel letter needed.
draw6 <- function() draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
                         merge_legend = TRUE, padding = unit(c(1, 1, 1, 1), "mm"))
save_fig(draw6, "Fig6", 8.5, 10.5)

qc_note("Fig6", "Rows: ", paste(sprintf("%s [%s]", g6, blk6), collapse = ", "),
        ". z-scores clipped at ±2: ", sum(abs(z) > 2), " cells.")
qc_note("Fig6", "Module statistics: ", paste(sprintf("%s %s (log2FC %s, padj %s)", g6, d6$dir, fmt_num(d6$log2FC),
                                                      fmt_p_txt(d6$padj)), collapse = "; "))

src6 <- data.frame(block = blk6, gene = g6, vst_mat, check.names = FALSE)
names(src6)[3:10] <- paste0("VST_", samples_order)
src6 <- cbind(src6, setNames(as.data.frame(z), paste0("z_", samples_order)),
              d6[, c("baseMean", "log2FC", "lfcSE", "pvalue", "padj", "dir")], stars = d6$stars)
write_src(src6, "Fig6_heatmap")
