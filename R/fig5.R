# Fig. 5 — channel genes by subfamily (17 cm wide, height by content, <= 22 cm).
# Bar = mean ± SEM of normalized counts with animal points, one free-y facet per gene,
# stars = DESeq2 padj. Only genes with padj get bars; n.t. (padj = NA) and n.d.
# (absent from the DESeq2 output) genes are listed in a text line under each panel.

if (!exists("ROOT")) ROOT <- normalizePath(".")
source(file.path(ROOT, "R", "helpers.R"))
suppressPackageStartupMessages(library(patchwork))

d <- get_data()
res <- d$res
panel_ch <- subset(d$panel, module == "channel")

NCOL   <- 10      # facets per full 17-cm row
ROW_CM <- 2.15    # height of one facet row (cm)
TXT_CM <- 0.45    # n.t./n.d. text line

panels <- list(
  A = c("PIEZO"),
  B = c("TMEM63", "TMC", "TMEM150C"),
  C = c("K2P", "KCNMA1"),
  D = c("TRPC", "TRPV", "TRPM", "TRPA", "TRPML", "TRPP"),
  E = c("TMEM16/ANO", "ENaC/ASIC"),
  F = c("LRRC8", "ClC"))

panel_genes <- function(fams) {
  g <- unlist(lapply(fams, function(f) panel_ch$gene[panel_ch$family == f]))
  st <- gene_status(g, res)
  list(all = g, bars = g[st == "padj"], nt = g[st == "n.t."], nd = g[st == "n.d."])
}

status_expr <- function(nt, nd) {
  it <- function(x) paste(sprintf("italic(\"%s\")", x), collapse = ' * ", " * ')
  parts <- c(if (length(nt)) paste0('"n.t.: " * ', it(nt)),
             if (length(nd)) paste0('"n.d.: " * ', it(nd)))
  if (!length(parts)) return(NULL)
  parse(text = paste(parts, collapse = ' * "; " * '))[[1]]
}

src_all <- list()
build_panel <- function(tag, fams, ncol, y_title = TRUE) {
  pg <- panel_genes(fams)
  df <- long_counts(pg$bars, d)
  sm <- summ_mean_sem(df)
  ann <- data.frame(gene = pg$bars, padj = res[pg$bars, "padj"])
  ann$stars <- p_stars(ann$padj)
  ann$top <- vapply(pg$bars, function(g) max(df$value[df$gene == g], sm$mean[sm$gene == g] + sm$sem[sm$gene == g]), 0)
  lv <- pg$bars
  df$gene <- factor(df$gene, lv); sm$gene <- factor(sm$gene, lv); ann$gene <- factor(ann$gene, lv)

  fam_of <- setNames(panel_ch$family, panel_ch$gene)
  mm <- merge(sm[, c("gene", "group", "n", "mean", "sd", "sem")],
              ann[, c("gene", "padj", "stars")], by = "gene", sort = FALSE)
  mm <- mm[order(match(mm$gene, lv), mm$group), ]
  miss <- c(pg$nt, pg$nd)
  src_all[[tag]] <<- rbind(
    data.frame(panel = tag, gene = as.character(mm$gene), family = unname(fam_of[as.character(mm$gene)]),
               status = "padj", group = as.character(mm$group), n = mm$n, mean = mm$mean, sd = mm$sd,
               sem = mm$sem, padj = mm$padj, stars = mm$stars),
    if (length(miss)) data.frame(panel = tag, gene = miss, family = unname(fam_of[miss]),
      status = rep(c("n.t.", "n.d."), c(length(pg$nt), length(pg$nd))), group = NA, n = NA, mean = NA,
      sd = NA, sem = NA, padj = NA, stars = rep(c("n.t.", "n.d."), c(length(pg$nt), length(pg$nd)))))

  p <- ggplot() +
    geom_bar_points(df, sm, seed = 5) +
    geom_segment(data = ann, aes(x = 1, xend = 2, y = top * 1.08, yend = top * 1.08), linewidth = 0.18) +
    geom_text(data = ann, aes(x = 1.5, y = top * ifelse(stars == "ns", 1.12, 1.06), label = stars,
                              vjust = ifelse(stars == "ns", 0, 0.1)),
              size = FONT_MM(7), family = "Arial") +
    geom_blank(data = ann, aes(x = 1.5, y = top * 1.3)) +
    facet_wrap(~ gene, scales = "free_y", ncol = ncol) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0)),
                       breaks = scales::breaks_extended(n = 3),
                       labels = scales::label_number(big.mark = "", scale_cut = scales::cut_short_scale())) +
    labs(x = NULL, y = if (y_title) "Normalized counts" else NULL, tag = tag, caption = status_expr(pg$nt, pg$nd)) +
    theme_ijms +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.y = ggtext::element_markdown(size = 7),
          strip.clip = "off",
          panel.spacing.x = unit(1.2, "mm"), panel.spacing.y = unit(1.5, "mm"),
          plot.caption = element_text(hjust = 0, size = 7, margin = margin(t = 1)),
          plot.caption.position = "plot",
          plot.margin = margin(2, 3, 2, 2))
  list(plot = p, nrow = ceiling(length(pg$bars) / ncol), has_txt = !is.null(status_expr(pg$nt, pg$nd)),
       n = length(pg$bars))
}

pa <- build_panel("A", panels$A, 2)
pb <- build_panel("B", panels$B, 7, y_title = FALSE)
pc <- build_panel("C", panels$C, NCOL)
pd <- build_panel("D", panels$D, NCOL)
pe <- build_panel("E", panels$E, NCOL)
pf <- build_panel("F", panels$F, NCOL)
write_src(do.call(rbind, src_all), "Fig5")

# Shared group key (drawn once, right of panel B).
key <- ggplot(data.frame(group = factor(c("Control", "rLLB"), c("Control", "rLLB")), y = c(2, 1))) +
  geom_tile(aes(x = 0, y = y, fill = group), width = 0.6, height = 0.55, colour = "black", linewidth = 0.18) +
  geom_text(aes(x = 0.5, y = y, label = group), hjust = 0, size = FONT_MM(7), family = "Arial") +
  scale_fill_manual(values = grp_fill, guide = "none") +
  scale_x_continuous(limits = c(-0.5, 3.5)) + scale_y_continuous(limits = c(-2, 3.5)) +
  theme_void()
key_legend <- key

row_h <- function(p) p$nrow * ROW_CM + if (p$has_txt) TXT_CM else 0
w_col <- 17 / NCOL
row1 <- (pa$plot | pb$plot | key_legend) + plot_layout(widths = c(pa$n, pb$n, NCOL - pa$n - pb$n))
rows <- list(row1, pc$plot, pd$plot, pe$plot, pf$plot)
heights <- c(max(row_h(pa), row_h(pb)), row_h(pc), row_h(pd), row_h(pe), row_h(pf))
fig5 <- wrap_plots(rows, ncol = 1, heights = heights)
height_cm <- sum(heights) + 0.3
if (height_cm > 22) warning("Fig. 5 is taller than 22 cm: ", round(height_cm, 1))

save_fig(fig5, "Fig5", width_cm = 17, height_cm = height_cm)
.bw$fig5 <- list(height_cm = height_cm, split_D = FALSE,
                 n_bars = sum(pa$n, pb$n, pc$n, pd$n, pe$n, pf$n))
