# fig7_western.R — Fig. 7: western blot densitometry (PIEZO1, TRPM4) and mRNA vs protein.
# Needs theme_ijms.R and the QC helpers from 00_setup.R (sourced by R/run_all.R).
# Normalisation: each protein's 8 raw values / mean of its own 4 Control values (nothing else).
# Welch t-test (two-sided, unpaired) is the second t-test exception (first: Myh7/Myh6, Fig. 1A).
# Mann–Whitney is computed as a robustness check only; it is not drawn and not used to pick the test.

path_wb_piezo1 <- "data/westerm_blot_control_piezo1.csv"   # sic: "westerm"
path_wb_trpm4  <- "data/westerm_blot_control_trpm4.csv"
near <- function(x, target, tol) isTRUE(abs(x - target) <= tol)

rd <- function(f) read.csv(f, fileEncoding = "UTF-8-BOM", check.names = FALSE, stringsAsFactors = FALSE)
wb_raw <- list(PIEZO1 = rd(path_wb_piezo1), TRPM4 = rd(path_wb_trpm4))
stopifnot(identical(colnames(wb_raw$PIEZO1), c("group", "control", "piezo1-rLLB")),
          identical(colnames(wb_raw$TRPM4),  c("group", "control", "trpm4")))

## ------------------------------------------------ densitometry (long) ---
wb_long <- do.call(rbind, lapply(names(wb_raw), function(p) {
  d <- wb_raw[[p]]
  data.frame(protein = p, group = rep(group_levels, each = nrow(d)), rat_id = rep(d$group, 2),
             raw = c(d[[2]], d[[3]]), stringsAsFactors = FALSE)
}))
wb_long$group <- factor(wb_long$group, group_levels)
wb_long$normalised <- NA_real_
for (p in names(wb_raw)) {
  i <- wb_long$protein == p
  wb_long$normalised[i] <- wb_long$raw[i] / mean(wb_long$raw[i & wb_long$group == "Control"])
}
wb_long$log2_normalised <- log2(wb_long$normalised)

## ------------------------------------------------------------ tests ---
wb_test <- do.call(rbind, lapply(names(wb_raw), function(p) {
  d <- wb_long[wb_long$protein == p, ]
  x <- d$normalised[d$group == "Control"]; y <- d$normalised[d$group == "rLLB"]
  lx <- log2(x); ly <- log2(y)
  tt <- t.test(y, x, var.equal = FALSE)
  mw <- suppressWarnings(wilcox.test(y, x, exact = TRUE))
  tl <- t.test(ly, lx, var.equal = FALSE)
  sem <- function(v) sd(v) / sqrt(length(v))
  data.frame(protein = p, n_control = length(x), n_rLLB = length(y),
             mean_control = mean(x), mean_rLLB = mean(y), SEM_control = sem(x), SEM_rLLB = sem(y),
             ratio = mean(y) / mean(x), t = unname(tt$statistic), df = unname(tt$parameter),
             welch_p = tt$p.value, mwu_p = mw$p.value,
             log2_diff = mean(ly) - mean(lx), log2_SE = sqrt(var(ly) / length(ly) + var(lx) / length(lx)),
             welch_p_log2 = tl$p.value, stringsAsFactors = FALSE)
}))
rownames(wb_test) <- wb_test$protein
raw_ctrl_mean <- sapply(wb_raw, function(d) mean(d$control))

## --------------------------------------------------------------- mRNA ---
ref7 <- read.csv(path_ref_res, check.names = FALSE, stringsAsFactors = FALSE)
mrna <- ref7[ref7$gene_id %in% c("Piezo1", "Trpm4"), c("gene_id", "log2FoldChange", "lfcSE", "padj")]
mrna <- mrna[match(c("Piezo1", "Trpm4"), mrna$gene_id), ]
rownames(mrna) <- mrna$gene_id
stopifnot(nrow(mrna) == 2)

## ---------------------------------------------------------- checks ---
Tt <- wb_test
qc_check("Fig7", "n per group, both proteins (4 + 4)", "4 + 4, 4 + 4",
         paste(sprintf("%d + %d", Tt$n_control, Tt$n_rLLB), collapse = ", "),
         all(Tt$n_control == 4 & Tt$n_rLLB == 4))
ctrl_differ <- !isTRUE(all.equal(wb_raw$PIEZO1$control, wb_raw$TRPM4$control))
qc_check("Fig7", "`control` columns differ between the two files", "TRUE", ctrl_differ, ctrl_differ)
qc_check("Fig7", "Raw Control mean, PIEZO1 / TRPM4 (tol 0.001)", "4303.735 / 9617.051",
         sprintf("%.3f / %.3f", raw_ctrl_mean["PIEZO1"], raw_ctrl_mean["TRPM4"]),
         near(raw_ctrl_mean["PIEZO1"], 4303.735, 0.001) && near(raw_ctrl_mean["TRPM4"], 9617.051, 0.001))
qc_check("Fig7", "SEM of normalised Control, PIEZO1 / TRPM4 (tol 0.001)", "0.070 / 0.265",
         sprintf("%.3f / %.3f", Tt["PIEZO1", "SEM_control"], Tt["TRPM4", "SEM_control"]),
         near(Tt["PIEZO1", "SEM_control"], 0.070, 0.001) && near(Tt["TRPM4", "SEM_control"], 0.265, 0.001))
qc_check("Fig7", "PIEZO1 rLLB mean ± SEM (tol 0.001)", "0.574 ± 0.095",
         sprintf("%.3f ± %.3f", Tt["PIEZO1", "mean_rLLB"], Tt["PIEZO1", "SEM_rLLB"]),
         near(Tt["PIEZO1", "mean_rLLB"], 0.574, 0.001) && near(Tt["PIEZO1", "SEM_rLLB"], 0.095, 0.001))
qc_check("Fig7", "PIEZO1 Welch t / df / p (tol 0.001 on p)", "−3.61 / 5.55 / 0.0128",
         sprintf("%s / %.2f / %.4f", fmt_num(Tt["PIEZO1", "t"], 2), Tt["PIEZO1", "df"], Tt["PIEZO1", "welch_p"]),
         near(Tt["PIEZO1", "welch_p"], 0.0128, 0.001) && near(Tt["PIEZO1", "t"], -3.61, 0.01) &&
           near(Tt["PIEZO1", "df"], 5.55, 0.01))
qc_check("Fig7", "PIEZO1 Mann–Whitney exact p (tol 0.001)", "0.0286", sprintf("%.4f", Tt["PIEZO1", "mwu_p"]),
         near(Tt["PIEZO1", "mwu_p"], 0.0286, 0.001))
qc_check("Fig7", "TRPM4 rLLB mean ± SEM (tol 0.001)", "1.054 ± 0.246",
         sprintf("%.3f ± %.3f", Tt["TRPM4", "mean_rLLB"], Tt["TRPM4", "SEM_rLLB"]),
         near(Tt["TRPM4", "mean_rLLB"], 1.054, 0.001) && near(Tt["TRPM4", "SEM_rLLB"], 0.246, 0.001))
qc_check("Fig7", "TRPM4 Welch t / df / p (tol 0.002 on p)", "0.15 / 5.97 / 0.887",
         sprintf("%.2f / %.2f / %.3f", Tt["TRPM4", "t"], Tt["TRPM4", "df"], Tt["TRPM4", "welch_p"]),
         near(Tt["TRPM4", "welch_p"], 0.887, 0.002) && near(Tt["TRPM4", "t"], 0.15, 0.01) &&
           near(Tt["TRPM4", "df"], 5.97, 0.01))
qc_check("Fig7", "TRPM4 Mann–Whitney exact p (tol 0.001)", "0.886", sprintf("%.3f", Tt["TRPM4", "mwu_p"]),
         near(Tt["TRPM4", "mwu_p"], 0.886, 0.001))
qc_check("Fig7", "log2 group difference PIEZO1 ± SE (tol 0.002)", "−0.851 ± 0.265",
         sprintf("%s ± %.3f", fmt_num(Tt["PIEZO1", "log2_diff"], 3), Tt["PIEZO1", "log2_SE"]),
         near(Tt["PIEZO1", "log2_diff"], -0.851, 0.002) && near(Tt["PIEZO1", "log2_SE"], 0.265, 0.002))
qc_check("Fig7", "log2 group difference TRPM4 ± SE (tol 0.002)", "+0.093 ± 0.479",
         sprintf("%+.3f ± %.3f", Tt["TRPM4", "log2_diff"], Tt["TRPM4", "log2_SE"]),
         near(Tt["TRPM4", "log2_diff"], 0.093, 0.002) && near(Tt["TRPM4", "log2_SE"], 0.479, 0.002))
for (x in list(list("Piezo1", -1.473, 0.334, 1.4e-4), list("Trpm4", -0.497, 0.186, 0.028))) {
  m <- mrna[x[[1]], ]
  qc_check("Fig7", sprintf("mRNA %s: log2FC / lfcSE / padj", x[[1]]),
           sprintf("%s / %.3f / %s", fmt_num(x[[2]], 3), x[[3]], fmt_p_txt(x[[4]])),
           sprintf("%s / %.3f / %s", fmt_num(m$log2FoldChange, 3), m$lfcSE, fmt_p_txt(m$padj)),
           near(m$log2FoldChange, x[[2]], 0.001) && near(m$lfcSE, x[[3]], 0.001) && signif(m$padj, 2) == x[[4]])
}
qc_stop_if_failed("Fig7")   # rule 0.7: no figure when a section-2 check fails

## ---------------------------------------------------- extra QC notes ---
cors <- do.call(rbind, lapply(names(wb_raw), function(p) {
  d <- wb_long[wb_long$protein == p, ]
  ct <- cor.test(d$normalised[d$group == "Control"], d$normalised[d$group == "rLLB"])
  data.frame(protein = p, r = unname(ct$estimate), p = ct$p.value)
}))
cv <- function(p, g) { v <- wb_long$normalised[wb_long$protein == p & wb_long$group == g]; 100 * sd(v) / mean(v) }
qc_note("Fig7", "**Data and assumptions.** Values are raw densitometry (ImageJ integrated band density, a.u.), ",
        "normalised only to the mean of the four Control values of the same file; normalisation to a loading control ",
        "and total protein is absent. Columns = groups (Control, rLLB; 4 animals each); `rat_N` is a within-group ",
        "ordinal, not a pair, so no pairing, no paired test, no connecting lines, no animal labels. ",
        "Files `data/westerm_blot_control_piezo1.csv` and `data/westerm_blot_control_trpm4.csv` (corrected PIEZO1 ",
        "version; Control columns differ between files: ", ctrl_differ, ").")
qc_note("Fig7", "**Statistics.** Welch t-test (two-sided, unpaired) on normalised densitometry, exact p on the figure; ",
        "this is the second exception to the 'no t-tests for stars' rule (first: Welch on *Myh7*/*Myh6*, Fig. 1A). ",
        "No stars. Mann–Whitney (exact) is recorded below as a robustness check only: not drawn, not used to choose the test. ",
        "Power was not computed.")
qc_note("Fig7", paste(c(
  "| protein | mean Control | mean rLLB | SEM Control | SEM rLLB | rLLB/Control | Welch t | df | Welch p | MWU p | log2 diff ± SE | Welch p (log2) |",
  "|---|---|---|---|---|---|---|---|---|---|---|---|",
  sprintf("| %s | %.3f | %.3f | %.3f | %.3f | %.3f | %s | %.2f | %.4f | %.4f | %s ± %.3f | %.3f |",
          Tt$protein, Tt$mean_control, Tt$mean_rLLB, Tt$SEM_control, Tt$SEM_rLLB, Tt$ratio,
          fmt_num(Tt$t, 2), Tt$df, Tt$welch_p, Tt$mwu_p, fmt_num(Tt$log2_diff, 3), Tt$log2_SE, Tt$welch_p_log2)),
  collapse = "\n"))
qc_note("Fig7", sprintf("Correlation within `rat_N` between Control and rLLB columns (normalised; reference only, not evidence of pairing; no paired test): PIEZO1 r = %.2f, p = %.2f; TRPM4 r = %.2f, p = %.2f.",
                        cors$r[1], cors$p[1], cors$r[2], cors$p[2]))
qc_note("Fig7", sprintf("CV of normalised values: PIEZO1 Control %.0f %%, rLLB %.0f %%; TRPM4 Control %.0f %%, rLLB %.0f %%.",
                        cv("PIEZO1", "Control"), cv("PIEZO1", "rLLB"), cv("TRPM4", "Control"), cv("TRPM4", "rLLB")))
qc_note("Fig7", sprintf("Welch on log2 normalised values (check): PIEZO1 p = %.3f, TRPM4 p = %.3f. With 4 vs 4 the exact Mann–Whitney p cannot go below 0.0286; PIEZO1 sits exactly at this minimum (all rLLB values below all Control values: %s). Reference only.",
                        Tt["PIEZO1", "welch_p_log2"], Tt["TRPM4", "welch_p_log2"],
                        max(wb_long$normalised[wb_long$protein == "PIEZO1" & wb_long$group == "rLLB"]) <
                          min(wb_long$normalised[wb_long$protein == "PIEZO1" & wb_long$group == "Control"])))
qc_note("Fig7", sprintf("Ratio of means rLLB/Control: PIEZO1 %.3f, TRPM4 %.3f.", Tt["PIEZO1", "ratio"], Tt["TRPM4", "ratio"]))
qc_note("Fig7", "Values are raw densitometry normalised only to the mean Control; normalisation to a loading control and total protein is absent.")

## ----------------------------------------------------------- panels ---
ymax <- 2.1
sem_tab <- function(p) {
  d <- wb_long[wb_long$protein == p, ]
  do.call(rbind, lapply(group_levels, function(g) {
    v <- d$normalised[d$group == g]
    data.frame(group = factor(g, group_levels), mean = mean(v), sem = sd(v) / sqrt(length(v)))
  }))
}
y_bracket <- 1.86     # common to 7A and 7B (maximum of both proteins is 1.76)
bar_panel <- function(p, title) {
  d <- wb_long[wb_long$protein == p, ]
  s <- sem_tab(p)
  jit <- position_jitter(width = 0.09, height = 0, seed = 7)   # <= 0.15 of bar width (0.6)
  lab <- paste0("p~", fmt_p_math(Tt[p, "welch_p"], 2, prefix = "= "))
  ggplot(d, aes(group, normalised)) +
    geom_col(data = s, aes(group, mean, fill = group), width = 0.6, colour = "black", linewidth = lw_axis) +
    geom_errorbar(data = s, aes(group, ymin = mean - sem, ymax = mean + sem), inherit.aes = FALSE,
                  width = 0.22, linewidth = lw_whisker) +
    geom_point(data = d[d$group == "Control", ], shape = 21, fill = "#000000", colour = "black",
               size = 1.3, stroke = 0.35, position = jit) +
    geom_point(data = d[d$group == "rLLB", ], shape = 21, fill = "#FFFFFF", colour = "black",
               size = 1.3, stroke = 0.35, position = jit) +
    annotate("segment", x = 1, xend = 2, y = y_bracket, yend = y_bracket, linewidth = lw_axis) +
    annotate("text", x = 1.5, y = y_bracket + 0.04, label = lab, parse = TRUE, vjust = 0,
             size = 6.5 / PT, family = FONT) +
    scale_fill_manual(values = c(Control = pal$ctrl_fill, rLLB = pal$bw_fill), guide = "none") +
    scale_y_continuous(limits = c(0, ymax), breaks = seq(0, 2, 0.5), labels = function(x) format(x, nsmall = 1),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = "Band density (a.u.,\nControl mean = 1)", title = title, subtitle = "n = 4 + 4") +
    coord_cartesian(clip = "off") +
    theme_ijms() +
    theme(plot.subtitle = element_text(size = 6.5, hjust = 0), plot.title = element_text(size = 8, face = "bold"))
}
pA <- bar_panel("PIEZO1", "PIEZO1")
pB <- bar_panel("TRPM4", "TRPM4")

## --- 7C: mRNA vs protein ---
dirf <- function(est, p) ifelse(p < 0.05, ifelse(est < 0, "down", "up"), "ns")
c7 <- data.frame(
  gene = rep(c("Piezo1", "Trpm4"), 2), protein = rep(c("PIEZO1", "TRPM4"), 2),
  assay = rep(c("mRNA (RNA-seq, DESeq2)", "Protein (Western blot)"), each = 2),
  estimate = c(mrna$log2FoldChange, Tt$log2_diff),
  error = c(mrna$lfcSE, Tt$log2_SE),
  error_type = rep(c("lfcSE", "SE (Welch, log2 densities)"), each = 2),
  p_or_padj = c(mrna$padj, Tt$welch_p), stringsAsFactors = FALSE)
c7$direction <- dirf(c7$estimate, c7$p_or_padj)
c7$kind <- factor(rep(c("mRNA (RNA-seq)", "Protein (Western blot)"), each = 2),
                  c("mRNA (RNA-seq)", "Protein (Western blot)"))
c7$row <- ifelse(c7$protein == "PIEZO1", 2, 1)
c7$y <- c7$row + ifelse(c7$kind == "mRNA (RNA-seq)", 0.17, -0.17)
c7$lab <- fmt_p_math(c7$p_or_padj, 2)
x_txt <- 1.75
pC <- ggplot(c7, aes(estimate, y)) +
  geom_vline(xintercept = 0, linewidth = lw_axis, colour = "#808080", linetype = "22") +
  geom_segment(aes(x = estimate - error, xend = estimate + error, yend = y, colour = direction),
               linewidth = lw_whisker) +
  geom_point(aes(colour = direction, shape = kind), size = 2.2) +
  geom_text(aes(x = x_txt, label = lab), parse = TRUE, hjust = 0, size = 6.5 / PT, family = FONT) +
  annotate("text", x = x_txt, y = 2.6, label = "mRNA padj\nprotein p", hjust = 0, vjust = 0.5, lineheight = 0.95, size = 6.5 / PT,
           family = FONT, fontface = "bold") +
  scale_colour_manual(values = dir_cols, guide = "none") +
  scale_shape_manual(values = c(`mRNA (RNA-seq)` = 16, `Protein (Western blot)` = 15), name = NULL,
                     guide = guide_legend(override.aes = list(colour = "black", size = 2.2))) +
  scale_x_continuous(limits = c(-2.5, 1.5), breaks = seq(-2, 1, 1), labels = lab_minus,
                     expand = expansion(mult = 0), oob = scales::oob_keep) +   # keep the p column right of the axis
  scale_y_continuous(breaks = c(2, 1), labels = c("PIEZO1", "TRPM4"), limits = c(0.55, 2.75),
                     expand = expansion(mult = 0)) +
  coord_cartesian(clip = "off") +
  labs(x = expression("log"[2] * " fold change (rLLB vs. Control)"), y = NULL) +
  theme_ijms() +
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-2, 0, 0, 0),
        plot.margin = margin(5.5, 50, 5.5, 5.5))

fig7 <- (pA + pB + pC + plot_layout(widths = c(1, 1, 1.5))) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold", family = FONT))
save_fig(fig7, "Fig7", 17, 6)

## ------------------------------------------------------- source data ---
sd_cols <- c("group", "rat_id", "raw", "normalised", "log2_normalised")
write_src(wb_long[wb_long$protein == "PIEZO1", sd_cols], "Fig7A_PIEZO1_densitometry")
write_src(wb_long[wb_long$protein == "TRPM4", sd_cols], "Fig7B_TRPM4_densitometry")
write_src(Tt[, c("protein", "mean_control", "mean_rLLB", "SEM_control", "SEM_rLLB", "ratio", "t", "df",
                 "welch_p", "mwu_p", "log2_diff", "log2_SE", "welch_p_log2")], "Fig7_western_tests")
write_src(c7[, c("gene", "protein", "assay", "estimate", "error", "error_type", "p_or_padj", "direction")],
          "Fig7C_mRNA_vs_protein")
