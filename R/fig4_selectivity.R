# Fig. 4 — selectivity of the channelome response.
# A: direction of DEGs, genome vs panel (Fisher). B: expression-matched permutation control.
# C: family set tests (limma::fry) + descriptive rLLB/Control ratio. D: Kcnk2 vs Kcnk3.

## -------------------------------------------------------------- 4A ---
in_panel <- de$gene %in% panel$gene
tab4 <- matrix(c(sum(de$dir == "down" & !in_panel), sum(de$dir == "up" & !in_panel),
                 sum(de$dir == "down" & in_panel),  sum(de$dir == "up" & in_panel)),
               nrow = 2, byrow = TRUE, dimnames = list(c("Genome", "Panel"), c("down", "up")))
ft <- fisher.test(tab4[c("Panel", "Genome"), ])                 # OR = odds(down | panel) / odds(down | genome)
qc_check("Fig4", "Genome (non-panel) down / up", "2543 / 2824",
         sprintf("%d / %d", tab4["Genome", "down"], tab4["Genome", "up"]),
         tab4["Genome", "down"] == 2543 && tab4["Genome", "up"] == 2824)
qc_check("Fig4", "Panel down / up", "20 / 3", sprintf("%d / %d", tab4["Panel", "down"], tab4["Panel", "up"]),
         tab4["Panel", "down"] == 20 && tab4["Panel", "up"] == 3)
qc_check("Fig4", "Fisher OR (two-sided)", "7.40", sprintf("%.2f", ft$estimate), round(ft$estimate, 2) == 7.40)
qc_check("Fig4", "Fisher p", "1.8 × 10^-4", fmt_p_txt(ft$p.value), fmt_p_txt(ft$p.value) == "1.8 × 10^-4")

dA <- data.frame(set = factor(rep(c("Genome", "Panel"), each = 2), c("Genome", "Panel")),
                 dir = factor(rep(c("down", "up"), 2), c("up", "down")),
                 n = c(tab4["Genome", ], tab4["Panel", ]))
dA$frac <- dA$n / ave(dA$n, dA$set, FUN = sum)
dA$ypos <- ifelse(dA$dir == "down", dA$frac / 2, 1 - dA$frac / 2)
pA <- ggplot(dA, aes(set, frac, fill = dir)) +
  geom_col(width = 0.65, colour = "black", linewidth = lw_axis) +
  geom_text(aes(y = ypos, label = n), colour = "white", size = 7 / PT, family = FONT) +
  annotate("text", x = 1.5, y = 1.06, vjust = 0, size = 7 / PT, family = FONT, parse = TRUE,
           label = sprintf('"OR = %.2f,"~italic(p)~%s', ft$estimate, fmt_p_math(ft$p.value, 2, prefix = "= "))) +
  scale_fill_manual(values = c(down = pal$down, up = pal$up), breaks = c("down", "up"),
                    labels = c("Down", "Up"), name = NULL) +
  scale_y_continuous(labels = function(x) paste0(100 * x), breaks = seq(0, 1, 0.25),
                     expand = c(0, 0), limits = c(0, 1.16)) +
  labs(x = NULL, y = "DEGs (padj < 0.05), %") +
  theme_ijms() + theme(legend.position = "right")

## -------------------------------------------------------------- 4B ---
tested <- de[!is.na(de$padj), ]
tested$bin <- cut(log10(tested$baseMean), breaks = quantile(log10(tested$baseMean), probs = seq(0, 1, 0.05)),
                  include.lowest = TRUE, labels = FALSE)
pan_t <- tested[tested$gene %in% panel$gene, ]
bg_t  <- tested[!tested$gene %in% panel$gene, ]
need <- table(factor(pan_t$bin, levels = 1:20))
bg_by_bin <- split(seq_len(nrow(bg_t)), factor(bg_t$bin, levels = 1:20))
obs_stat <- sum(pan_t$dir == "down") - sum(pan_t$dir == "up")
qc_check("Fig4", "Observed panel (down − up)", 17, obs_stat, obs_stat == 17)

set.seed(20260928)
NPERM <- 2000
perm <- t(vapply(seq_len(NPERM), function(i) {
  idx <- unlist(lapply(names(need)[need > 0], function(b) {
    pool <- bg_by_bin[[b]]; pool[sample.int(length(pool), need[[b]])]
  }))
  d <- bg_t$dir[idx]
  c(down = sum(d == "down"), up = sum(d == "up"))
}, numeric(2)))
perm <- data.frame(perm, diff = perm[, "down"] - perm[, "up"])
n_ge <- sum(perm$diff >= obs_stat)
p_perm <- n_ge / NPERM
p_lab <- if (n_ge == 0) sprintf('italic(p) < "5 %s 10"^"%s4"~"(0 of %d permutations)"', TIMES, MINUS, NPERM) else
  sprintf('italic(p)~"= %s (%d of %d permutations)"', formatC(p_perm, format = "g", digits = 2), n_ge, NPERM)
ms <- sapply(perm, function(x) c(mean = mean(x), sd = sd(x)))
qc_note("Fig4", sprintf("Permutations (seed 20260928, %d, 20 quantile bins of log10 baseMean over %d genes): ",
                        NPERM, nrow(tested)),
        sprintf("background down %.1f ± %.1f, up %.1f ± %.1f, down − up %.1f ± %.1f (manuscript: 7.7 ± 2.5, 9.2 ± 2.7, −1.5 ± 4.1). ",
                ms["mean", "down"], ms["sd", "down"], ms["mean", "up"], ms["sd", "up"], ms["mean", "diff"], ms["sd", "diff"]),
        sprintf("Permutations with (down − up) ≥ %d: %d.", obs_stat, n_ge))
for (k in c("down", "up", "diff")) {
  ref_mean <- c(down = 7.7, up = 9.2, diff = -1.5)[[k]]
  qc_check("Fig4", paste0("Permutation mean ", k, " vs manuscript (|Δ| ≤ 1)"), ref_mean,
           sprintf("%.2f", ms["mean", k]), abs(ms["mean", k] - ref_mean) <= 1)
}

pB <- ggplot(perm, aes(diff)) +
  geom_histogram(binwidth = 1, fill = "grey70", colour = "white", linewidth = 0.1) +
  geom_vline(xintercept = obs_stat, colour = pal$up, linewidth = lw_whisker * 1.5) +
  annotate("text", x = obs_stat - 0.6, y = Inf, vjust = 1.3, hjust = 1, colour = pal$up, size = 7 / PT,
           family = FONT, label = sprintf("Panel: %d", obs_stat)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = 6.5 / PT, family = FONT,
           parse = TRUE, label = p_lab) +
  scale_x_continuous(labels = lab_minus, breaks = seq(-20, 20, 5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Down − up among 62 expression-matched genes", y = "Permutations") +
  theme_ijms()

## -------------------------------------------------------------- 4C ---
y <- DGEList(counts = counts(dds)[, samples_ref])
y <- calcNormFactors(y)
design <- model.matrix(~ condition, data = as.data.frame(colData(dds)))
stopifnot("conditionblast" %in% colnames(design))
vo <- voom(y, design)

fam_label <- function(f) ifelse(f == "LRRC8", "VRAC/LRRC8", f)
panel$set <- fam_label(panel$subfamily)
sets_all <- split(panel$gene, panel$set)
sets_in  <- lapply(sets_all, function(g) intersect(g, rownames(vo)))
sets_in  <- sets_in[lengths(sets_in) > 0]
multi    <- sets_in[lengths(sets_in) >= 2]
single   <- sets_in[lengths(sets_in) == 1]
absent   <- setdiff(names(sets_all), names(sets_in))

fr <- fry(vo, index = multi, design = design, contrast = which(colnames(design) == "conditionblast"),
          sort = "none")
fr$FDR_BH <- p.adjust(fr$PValue, method = "BH")

famC <- do.call(rbind, lapply(names(sets_in), function(s) {
  g <- sets_in[[s]]
  sc <- de[g, "mean_norm_Control"]; sb <- de[g, "mean_norm_rLLB"]
  dom <- g[which.max(sc + sb)]
  data.frame(set = s, n_genes = length(g), genes = paste(g, collapse = " "),
             sum_mean_norm_Control = sum(sc), sum_mean_norm_rLLB = sum(sb),
             ratio = sum(sb) / sum(sc),
             dominant_gene = dom, dominant_share = max(sc + sb) / sum(sc + sb),
             fry_Direction = if (s %in% rownames(fr)) fr[s, "Direction"] else NA,
             fry_PValue = if (s %in% rownames(fr)) fr[s, "PValue"] else NA,
             FDR = if (s %in% rownames(fr)) fr[s, "FDR_BH"] else NA,
             tested = s %in% rownames(fr), stringsAsFactors = FALSE)
}))
famC$colour <- ifelse(!is.na(famC$FDR) & famC$FDR < 0.05,
                      ifelse(famC$fry_Direction == "Down", "down", "up"), "ns")
famC <- famC[order(famC$ratio), ]

ref_ratio <- c(PIEZO = 0.34, TMC = 0.40, `TREK/TRAAK` = 0.44, TRPP = 0.51, TMEM63 = 0.64, `VRAC/LRRC8` = 0.68,
               `ENaC/ASIC` = 0.70, TRPC = 0.83, TRPV = 0.86, TRPM = 0.95, ClC = 1.03, TRPML = 1.04,
               TMEM16 = 1.06, `remaining K2P` = 1.14, TMEM150C = 1.21)
famC$ratio_manuscript <- unname(ref_ratio[famC$set])
famC$ratio_diff <- famC$ratio - famC$ratio_manuscript
bad_ratio <- famC[!is.na(famC$ratio_diff) & abs(famC$ratio_diff) > 0.01, ]
not_in_ms <- famC$set[is.na(famC$ratio_manuscript)]
if (length(not_in_ms)) qc_note("Fig4", "Sets without a manuscript ratio (shown, not compared): ",
                               paste(sprintf("%s %.3f", not_in_ms, famC$ratio[match(not_in_ms, famC$set)]), collapse = ", "), ".")
qc_check("Fig4", "Family ratios vs manuscript (|Δ| ≤ 0.01; 15 sets)", "0 differences", nrow(bad_ratio), nrow(bad_ratio) == 0)
if (nrow(bad_ratio)) qc_note("Fig4", "Ratio differences > 0.01: ",
                             paste(sprintf("%s %.3f vs %s", bad_ratio$set, bad_ratio$ratio, bad_ratio$ratio_manuscript),
                                   collapse = "; "))
qc_check("Fig4", "Single-gene sets (shown, not tested)", "KCNMA1, TMEM150C",
         paste(sort(names(single)), collapse = ", "), identical(sort(names(single)), c("KCNMA1", "TMEM150C")))
qc_check("Fig4", "Families absent from the 18 114 genes", "TRPA", paste(absent, collapse = ", "),
         identical(absent, "TRPA"))

famC$set_f <- factor(famC$set, famC$set)
famC$right <- ifelse(famC$tested, fmt_p_math(famC$FDR, 2), '"single gene"')
xr <- range(famC$ratio)
x_txt <- 2^(log2(xr[2]) + 0.55)
pC <- ggplot(famC, aes(ratio, set_f)) +
  geom_vline(xintercept = 1, linewidth = lw_axis, colour = "grey40", linetype = "22") +
  geom_point(aes(colour = colour), size = 2) +
  geom_text(aes(x = x_txt, label = right), parse = TRUE, hjust = 0, size = 6.5 / PT, family = FONT) +
  annotate("text", x = x_txt, y = nrow(famC) + 0.9, label = "FDR (fry)", hjust = 0, size = 6.5 / PT,
           family = FONT, fontface = "bold") +
  scale_colour_manual(values = dir_cols, guide = "none") +
  scale_x_continuous(trans = "log2", breaks = c(0.25, 0.5, 1, 2), labels = c("0.25", "0.5", "1", "2"),
                     limits = c(0.25, 2^(log2(x_txt) + 1.2))) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.2))) +
  coord_cartesian(clip = "off") +
  labs(x = expression("rLLB / Control (sum of normalised counts, " * log[2] * " scale)"), y = NULL) +
  theme_ijms()

## -------------------------------------------------------------- 4D ---
d4 <- long_counts(c("Kcnk2", "Kcnk3"))
d4$gene_lab <- factor(sprintf("*%s* (%s)", d4$gene, ifelse(d4$gene == "Kcnk2", "TREK-1", "TASK-1")))
s4 <- do.call(rbind, lapply(split(d4, list(d4$gene, d4$group)), function(x)
  data.frame(gene = x$gene[1], gene_lab = x$gene_lab[1], group = x$group[1],
             mean = mean(x$norm_count), sem = sd(x$norm_count) / sqrt(nrow(x)))))
top4 <- do.call(rbind, lapply(split(d4, d4$gene), function(x) {
  s <- s4[s4$gene == x$gene[1], ]
  data.frame(gene = x$gene[1], gene_lab = x$gene_lab[1], top = max(x$norm_count, s$mean + s$sem))
}))
top4$padj <- de[top4$gene, "padj"]
top4$lab <- paste0("padj~", fmt_p_math(top4$padj, 2, prefix = "= "))
# Jitter ≤ 0.15 of the bar width (0.6) -> ±0.09.
jit <- position_jitter(width = 0.09, height = 0, seed = 11)
pD <- ggplot(d4, aes(group, norm_count)) +
  geom_col(data = s4, aes(group, mean, fill = group), width = 0.6, colour = "black", linewidth = lw_axis) +
  geom_errorbar(data = s4, aes(group, ymin = mean - sem, ymax = mean + sem), inherit.aes = FALSE,
                width = 0.22, linewidth = lw_whisker) +
  geom_point(data = d4[d4$group == "Control", ], shape = 21, fill = "#000000", colour = "black",
             size = 1.3, stroke = 0.35, position = jit) +
  geom_point(data = d4[d4$group == "rLLB", ], shape = 21, fill = "#FFFFFF", colour = "black",
             size = 1.3, stroke = 0.35, position = jit) +
  geom_segment(data = top4, aes(x = 1, xend = 2, y = top * 1.08, yend = top * 1.08), inherit.aes = FALSE,
               linewidth = lw_axis) +
  geom_text(data = top4, aes(x = 1.5, y = top * 1.11, label = lab), inherit.aes = FALSE, parse = TRUE,
            vjust = 0, size = 6.5 / PT, family = FONT) +
  geom_blank(data = top4, aes(x = 1.5, y = top * 1.3), inherit.aes = FALSE) +
  facet_wrap(~ gene_lab, scales = "free_y") +
  scale_fill_manual(values = c(Control = pal$ctrl_fill, rLLB = pal$bw_fill), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)), limits = c(0, NA)) +
  labs(x = NULL, y = "Normalised counts") +
  theme_ijms() + theme(strip.text = element_markdown(face = "plain", size = 7))

fig4 <- (pA + pB + plot_layout(widths = c(1, 1.3))) / (pC + pD + plot_layout(widths = c(1.35, 1))) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold", family = FONT))
save_fig(fig4, "Fig4", 17, 15)

write_src(data.frame(set = rownames(tab4), down = tab4[, "down"], up = tab4[, "up"],
                     fisher_OR = unname(ft$estimate), fisher_p = ft$p.value,
                     OR_definition = "odds(down | panel) / odds(down | genome), two-sided fisher.test"),
          "Fig4A_direction_counts")
write_src(data.frame(permutation = seq_len(NPERM), perm, observed_diff = obs_stat), "Fig4B_permutations")
bm <- rbind(data.frame(set = "panel (62)", gene = pan_t$gene, baseMean = pan_t$baseMean, bin = pan_t$bin),
            data.frame(set = "background pool (non-panel, padj not NA)", gene = bg_t$gene,
                       baseMean = bg_t$baseMean, bin = bg_t$bin))
write_src(bm, "Fig4B_baseMean_panel_vs_background")
write_src(data.frame(bin = 1:20, panel_genes = as.integer(need),
                     background_pool = lengths(bg_by_bin),
                     panel_median_baseMean = sapply(1:20, function(b) median(pan_t$baseMean[pan_t$bin == b])),
                     background_median_baseMean = sapply(1:20, function(b) median(bg_t$baseMean[bg_t$bin == b]))),
          "Fig4B_bin_matching")
write_src(famC[, c("set", "n_genes", "genes", "sum_mean_norm_Control", "sum_mean_norm_rLLB", "ratio",
                   "ratio_manuscript", "ratio_diff", "fry_Direction", "fry_PValue", "FDR", "tested",
                   "dominant_gene", "dominant_share", "colour")], "Fig4C_family_fry")
write_src(merge(d4[, c("gene", "sample", "group", "norm_count")],
                data.frame(gene = top4$gene, padj = top4$padj), by = "gene"), "Fig4D_Kcnk2_Kcnk3")
