# Fig. 2 — channelome of Control atria (TPM, mean of R004–R007).
# TPM ranks different genes within Control only; group comparisons are never made on TPM.
# A: summed TPM per family. B: ranked panel genes (n.d. genes not drawn). C: reserved (not built).

panel$tpm_ctrl <- unname(tpm_mean_ctrl[panel$gene])       # NA if the symbol is absent from TPM
panel$in_tpm <- !is.na(panel$tpm_ctrl)

# TPM verification against the (rounded) table: 1 % relative tolerance.
chk <- merge(panel[, c("gene", "tpm_ctrl")], panel_tab[, c("gene", "mean_TPM_sham")], by = "gene")
chk$rel_diff <- ifelse(chk$mean_TPM_sham == 0 & chk$tpm_ctrl == 0, 0,
                       abs(chk$tpm_ctrl - chk$mean_TPM_sham) / pmax(chk$mean_TPM_sham, 1e-12))
chk_ok <- chk[!is.na(chk$mean_TPM_sham) | !is.na(chk$tpm_ctrl), ]
bad_tpm <- chk_ok[is.na(chk_ok$rel_diff) | chk_ok$rel_diff > 0.01, ]
qc_check("Fig2", "Mean Control TPM vs panel table (≤ 1 %)", "0 mismatches", nrow(bad_tpm), nrow(bad_tpm) == 0)
if (nrow(bad_tpm)) qc_note("Fig2", "TPM mismatches: ",
                           paste(sprintf("%s table %s vs %s", bad_tpm$gene, bad_tpm$mean_TPM_sham,
                                         signif(bad_tpm$tpm_ctrl, 4)), collapse = "; "))
qc_note("Fig2", "Max relative TPM difference vs table: ", sprintf("%.2e", max(chk_ok$rel_diff, na.rm = TRUE)),
        ". Panel genes absent from the TPM files: ",
        paste(panel$gene[!panel$in_tpm], collapse = ", "), ".")

## -------------------------------------------------------------- 2A ---
fam_disp <- function(f) ifelse(f == "LRRC8", "LRRC8 (VRAC)", f)
famA <- aggregate(tpm_ctrl ~ family, data = transform(panel[panel$in_tpm, ], family = fam_disp(family)), FUN = sum)
famA$n_genes <- as.integer(table(fam_disp(panel$family[panel$in_tpm]))[famA$family])
famA <- famA[order(-famA$tpm_ctrl), ]
top_fams <- head(famA$family, 6)
fam_cols <- c(setNames(ggsci::pal_npg()(6), top_fams),
              setNames(rep(pal$other_family, nrow(famA) - 6), famA$family[-(1:6)]))
famA$colour_group <- ifelse(famA$family %in% top_fams, famA$family, "Other families")
famA$family <- factor(famA$family, rev(famA$family))
x_floor <- 1e-4
pA <- ggplot(famA) +
  geom_rect(aes(xmin = x_floor, xmax = tpm_ctrl, ymin = as.numeric(family) - 0.38,
                ymax = as.numeric(family) + 0.38, fill = family), colour = NA) +
  geom_text(aes(x = tpm_ctrl * 1.25, y = as.numeric(family), label = paste0("n = ", n_genes)),
            hjust = 0, size = 6 / PT, family = FONT) +
  scale_fill_manual(values = fam_cols, guide = "none") +
  scale_x_log10(limits = c(x_floor, 500), breaks = 10^(-4:2), labels = lab_log10, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(levels(famA$family)), labels = levels(famA$family),
                     expand = expansion(add = 0.6)) +
  labs(x = expression("TPM (" * log[10] * " scale)"), y = NULL) +
  theme_ijms()

## -------------------------------------------------------------- 2C placeholder ---
pCph <- ggplot() +
  annotate("text", x = 0, y = 0, label = "Reserved: comparison with sinoatrial node\n(Turner et al., 2021) — source data pending",
           size = 6 / PT, family = FONT, colour = "grey55") +
  theme_void() + theme(plot.tag = element_text(size = 10, face = "bold", family = FONT))
qc_note("Fig2", "Fig. 2C (SAN comparison, Turner 2021) NOT built: source data are not in the repository; ",
        "an empty placeholder panel keeps the layout.")

## -------------------------------------------------------------- 2B ---
b <- panel[panel$status != "n.d." & panel$in_tpm, ]
floor_b <- 3e-4                                        # plotting position for TPM = 0
b$zero <- b$tpm_ctrl == 0
b$y <- ifelse(b$zero, floor_b, b$tpm_ctrl)
b <- b[order(-b$tpm_ctrl, b$gene), ]
b$gene_f <- factor(b$gene, b$gene)
b$fam_col <- ifelse(fam_disp(b$family) %in% top_fams, fam_disp(b$family), "Other families")
leg_cols <- c(setNames(ggsci::pal_npg()(6), top_fams), `Other families` = pal$other_family)
other_lab <- "Other families (grey bars in A)"

pB <- ggplot(b, aes(gene_f, y)) +
  geom_segment(aes(xend = gene_f, y = floor_b * 0.6, yend = y), linewidth = lw_axis, colour = "grey55") +
  geom_point(aes(colour = fam_col, shape = zero), size = 1.3, stroke = 0.4) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
  scale_colour_manual(values = leg_cols, breaks = c(top_fams, "Other families"),
                      labels = c(top_fams, other_lab), name = NULL) +
  scale_y_log10(limits = c(floor_b * 0.6, 60), breaks = c(floor_b, 10^(-3:1)),
                labels = c("0", lab_log10(10^(-3:1))), expand = c(0, 0)) +
  labs(x = NULL, y = expression("Mean TPM, Control (" * log[10] * " scale)")) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.6))) +
  theme_ijms() +
  theme(axis.text.x = element_text(face = "italic", angle = 90, hjust = 1, vjust = 0.5, size = 6),
        legend.position = "top", legend.justification = "left",
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"))

fig2 <- (pA + pCph + plot_layout(widths = c(1, 1))) / pB + plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = list(c("A", "C", "B"))) &
  theme(plot.tag = element_text(size = 10, face = "bold", family = FONT))
save_fig(fig2, "Fig2", 17, 15)

write_src(data.frame(family = as.character(famA$family), sum_mean_TPM_Control = famA$tpm_ctrl,
                     n_genes_in_TPM = famA$n_genes, colour_group = famA$colour_group),
          "Fig2A_family_TPM")
b_out <- panel[, c("family", "gene", "status", "tpm_ctrl")]
b_out <- cbind(b_out, tpm[match(panel$gene, rownames(tpm)), ctrl_ids])
b_out$plotted_in_2B <- panel$status != "n.d." & panel$in_tpm
b_out$note <- ifelse(panel$status == "n.d.", "n.d. (absent from DESeq2 output) - not drawn",
                     ifelse(!panel$in_tpm, "absent from TPM files",
                            ifelse(panel$tpm_ctrl == 0, "TPM = 0 in Control, drawn at the '0' position", "")))
b_out <- b_out[order(-b_out$tpm_ctrl, na.last = TRUE), ]
names(b_out)[names(b_out) == "tpm_ctrl"] <- "mean_TPM_Control"
write_src(b_out, "Fig2B_panel_gene_TPM")
