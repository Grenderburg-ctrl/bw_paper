# Fig. 4 — selectivity (A, B, D drawn; C left as an empty labelled slot).
#   A: 100 % stacked bars, Genome (panel removed) vs. Panel, Down/Up; Fisher exact test.
#   B: expression-matched permutation null for (down − up) in the panel.
#   D: Kcnk2 (TREK-1) and Kcnk3 (TASK-1) normalized counts; padj from DESeq2 output.

if (!exists("ROOT")) ROOT <- normalizePath(".")
source(file.path(ROOT, "R", "helpers.R"))
suppressPackageStartupMessages(library(patchwork))

d <- get_data()
res <- d$res
panel_ch <- subset(d$panel, module == "channel")
in_panel <- res$gene_id %in% panel_ch$gene

## ---------------------------------------------------------------- panel A ---
cnt <- function(sel) c(Down = sum(res$direction[sel] == "Down"), Up = sum(res$direction[sel] == "Up"))
g_cnt <- cnt(!in_panel); p_cnt <- cnt(in_panel)
tab <- rbind(Panel = p_cnt, Genome = g_cnt)
fisher <- fisher.test(tab, alternative = "two.sided")

fa <- data.frame(set = rep(c("Genome", "Panel"), each = 2), dir = rep(c("Down", "Up"), 2),
                 n = c(g_cnt, p_cnt))
fa$frac <- fa$n / ave(fa$n, fa$set, FUN = sum)
fa$set <- factor(fa$set, levels = c("Genome", "Panel"))
fa$dir <- factor(fa$dir, levels = c("Up", "Down"))   # Down at the bottom of the stack
write_src(transform(fa, fisher_OR = unname(fisher$estimate), fisher_p = fisher$p.value), "Fig4A")

or_lab <- sprintf("OR = %.2f, %s", unname(fisher$estimate), fmt_p_md(fisher$p.value))
pA <- ggplot(fa, aes(set, frac, fill = dir)) +
  geom_col(width = BAR_W, colour = "black", linewidth = 0.18) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), colour = "white",
            size = FONT_MM(7), family = "Arial") +
  annotate("segment", x = 1, xend = 2, y = 1.05, yend = 1.05, linewidth = 0.18) +
  geom_md(data = data.frame(x = 1.5, y = 1.08, l = or_lab), aes(x, y, label = l),
          inherit.aes = FALSE, vjust = 0, size = FONT_MM(7)) +
  scale_fill_manual(values = dir_col[c("Down", "Up")], breaks = c("Down", "Up"), name = NULL) +
  scale_y_continuous(labels = scales::label_percent(), breaks = seq(0, 1, 0.25),
                     limits = c(0, 1.2), expand = expansion(mult = 0)) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = paste0("DEGs (", MD_PADJ, " &lt; 0.05)"), tag = "A") +
  theme_ijms + theme(legend.position = "right")

## ---------------------------------------------------------------- panel B ---
# For every panel gene with padj, draw one gene at random among the 200 non-panel
# genes with padj that are nearest in log10(baseMean). 2000 permutations, seed 1.
N_PERM <- 2000; K_NN <- 200
tested <- subset(res, !is.na(padj))
pg <- subset(tested, gene_id %in% panel_ch$gene)
bgp <- subset(tested, !gene_id %in% panel_ch$gene)
lbm_bg <- log10(bgp$baseMean)
nn <- lapply(log10(pg$baseMean), function(x) order(abs(lbm_bg - x))[seq_len(K_NN)])

set.seed(1)
null <- t(vapply(seq_len(N_PERM), function(i) {
  pick <- vapply(nn, function(ix) ix[sample.int(K_NN, 1)], integer(1))
  dd <- bgp$direction[pick]
  c(down = sum(dd == "Down"), up = sum(dd == "Up"))
}, numeric(2)))
null <- data.frame(perm = seq_len(N_PERM), down = null[, "down"], up = null[, "up"])
null$diff <- null$down - null$up
write.csv(null, file.path(dir_src, "Fig4B_null.csv"), row.names = FALSE)

obs <- sum(pg$direction == "Down") - sum(pg$direction == "Up")
n_ge <- sum(null$diff >= obs)
write_src(data.frame(observed_down = sum(pg$direction == "Down"), observed_up = sum(pg$direction == "Up"),
                     observed_diff = obs, n_perm = N_PERM, n_perm_ge_obs = n_ge,
                     mean_down = mean(null$down), sd_down = sd(null$down),
                     mean_up = mean(null$up), sd_up = sd(null$up),
                     mean_diff = mean(null$diff), sd_diff = sd(null$diff)), "Fig4B")

y_lab <- max(table(null$diff)) * 1.08
pB <- ggplot(null, aes(diff)) +
  geom_histogram(binwidth = 1, fill = "grey70", colour = "white", linewidth = 0.1, boundary = 0.5) +
  geom_vline(xintercept = obs, colour = dir_col["Down"], linewidth = 0.4) +
  annotate("text", x = obs - 0.6, y = y_lab, vjust = 0, hjust = 1, colour = dir_col["Down"],
           label = sprintf("observed = %d;\n%d/%d permutations ≥ observed", obs, n_ge, N_PERM),
           size = FONT_MM(7), family = "Arial", lineheight = 0.95) +
  scale_x_continuous(labels = lab_minus(), breaks = scales::breaks_width(5)) +
  scale_y_continuous(limits = c(0, y_lab * 1.22), expand = expansion(mult = c(0, 0))) +
  labs(x = "Down − up (expression-matched null)", y = "Permutations", tag = "B") +
  theme_ijms

## ---------------------------------------------------------------- panel C ---
pC <- ggplot() + theme_void() + labs(tag = "C") +
  theme(plot.tag = element_text(size = 10, face = "bold", family = "Arial"))

## ---------------------------------------------------------------- panel D ---
genes_d <- c("Kcnk2", "Kcnk3")
prot    <- c(Kcnk2 = "TREK-1", Kcnk3 = "TASK-1")
dfd <- long_counts(genes_d, d)
smd <- summ_mean_sem(dfd)
padj_d <- res[genes_d, "padj"]
ann <- data.frame(gene = genes_d, padj = padj_d,
                  top = tapply(pmax(dfd$value, 0), dfd$gene, max)[genes_d])
ann$top <- pmax(ann$top, (smd$mean + smd$sem)[match(genes_d, smd$gene)])
ann$lab <- fmt_p_md(ann$padj, lab = MD_PADJ)
write_src(merge(smd, ann[, c("gene", "padj")]), "Fig4D")
write_src(dfd, "Fig4D_points")

strip_lab <- setNames(sprintf('italic("%s") ~ "(%s)"', genes_d, prot[genes_d]), genes_d)
for (df_name in c("dfd", "smd", "ann")) {
  x <- get(df_name); x$gene <- factor(strip_lab[x$gene], levels = strip_lab); assign(df_name, x)
}

pD <- ggplot() +
  geom_bar_points(dfd, smd, seed = 3) +
  geom_segment(data = ann, aes(x = 1, xend = 2, y = top * 1.1, yend = top * 1.1), linewidth = 0.18) +
  geom_md(data = ann, aes(x = 1.5, y = top * 1.13, label = lab), vjust = 0, size = FONT_MM(7)) +
  geom_blank(data = ann, aes(x = 1.5, y = top * 1.35)) +
  facet_wrap(~ gene, scales = "free_y", labeller = label_parsed) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, NA),
                     labels = scales::label_number(big.mark = "")) +
  labs(x = NULL, y = "Normalized counts", tag = "D") +
  theme_ijms + theme(strip.text = element_text(face = "plain", size = 7))

## ---------------------------------------------------------------- layout ---
fig4 <- wrap_plots(pA, pB, pC, pD, ncol = 2, widths = c(1, 1.3))
save_fig(fig4, "Fig4", width_cm = 17, height_cm = 12)

.bw$fig4 <- list(tab = tab, fisher = fisher, null = null, obs = obs, n_ge = n_ge)
