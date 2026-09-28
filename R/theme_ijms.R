# theme_ijms.R — IJMS/MDPI figure style: theme, palette, text helpers, export.

FONT <- "Arial"
PT <- ggplot2::.pt                     # mm per pt for geom_text sizes
lw_axis    <- 0.18                     # ≈ 0.5 pt
lw_whisker <- 0.26                     # ≈ 0.75 pt

pal <- list(
  ctrl_fill = "#FFFFFF", ctrl_line = "#000000",
  bw_fill   = "#595959",
  pca_ctrl  = "#A6A6A6", pca_bw = "#000000",
  down = "#2166AC", up = "#B2182B", ns = "#BFBFBF",
  heat = c("#2166AC", "#F7F7F7", "#B2182B"),
  other_family = "#8C8C8C"
)
dir_cols <- c(down = pal$down, up = pal$up, ns = pal$ns)

theme_ijms <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = FONT) %+replace%
    theme(
      line = element_line(linewidth = lw_axis, colour = "black"),
      axis.line = element_line(linewidth = lw_axis, colour = "black"),
      axis.ticks = element_line(linewidth = lw_axis, colour = "black"),
      axis.ticks.length = unit(1.2, "mm"),
      axis.text = element_text(size = 7, colour = "black"),
      axis.title = element_text(size = 8, colour = "black"),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 7),
      legend.key.size = unit(3, "mm"),
      legend.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 7, face = "bold"),
      plot.tag = element_text(size = 10, face = "bold", family = FONT),
      plot.title = element_text(size = 8, hjust = 0),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      complete = TRUE
    )
}

## ---------------------------------------------------------- text helpers ---
MINUS <- "−"
TIMES <- "×"

p_stars <- function(p) {
  ifelse(is.na(p), "n.t.",
         ifelse(p < 1e-4, "****", ifelse(p < 1e-3, "***",
                ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))))
}

# Signed number with a true minus sign.
fmt_num <- function(x, digits = 2) {
  s <- formatC(x, format = "f", digits = digits)
  sub("^-", MINUS, s)
}

# p value as a plotmath string: "0.036", "1.5 × 10^-15" (all glyphs from the text font;
# the minus is U+2212 inside a string, so no Symbol font is pulled into the PDF).
fmt_p_math <- function(p, digits = 2, prefix = "") {
  vapply(p, function(v) {
    if (is.na(v)) return(sprintf('"%sn.t."', prefix))
    if (v >= 0.001) return(sprintf('"%s%s"', prefix, formatC(signif(v, digits), format = "fg", digits = digits,
                                                              flag = "#")))
    e <- floor(log10(v)); mant <- v / 10^e
    mant_s <- formatC(round(mant, digits - 1), format = "f", digits = digits - 1)
    if (mant_s == "10.0" || mant_s == "10") { e <- e + 1; mant_s <- formatC(1, format = "f", digits = digits - 1) }
    sprintf('"%s%s %s 10"^"%s%d"', prefix, mant_s, TIMES, MINUS, abs(e))
  }, character(1))
}

# Same, as a plain string for CSV / markdown ("1.5 × 10^-15").
fmt_p_txt <- function(p, digits = 2) {
  vapply(p, function(v) {
    if (is.na(v)) return("NA")
    if (v >= 0.001) return(formatC(signif(v, digits), format = "fg", digits = digits, flag = "#"))
    e <- floor(log10(v)); mant <- round(v / 10^e, digits - 1)
    if (mant >= 10) { mant <- 1; e <- e + 1 }
    sprintf("%s × 10^%d", formatC(mant, format = "f", digits = digits - 1), e)
  }, character(1))
}

# Axis labels with a true minus sign.
lab_minus <- function(x) sub("^-", MINUS, format(x, trim = TRUE, drop0trailing = TRUE))
lab_log10 <- function(x) {                   # 0.001, 0.01, 0.1, 1, 10, 100 without sci notation
  vapply(x, function(v) if (is.na(v)) NA_character_ else
    format(v, scientific = FALSE, drop0trailing = TRUE, trim = TRUE), character(1))
}

## ---------------------------------------------------------------- export ---
# `plot`: a ggplot/patchwork object, or a function that draws on the current device.
save_fig <- function(plot, name, width_cm, height_cm) {
  pdf_file  <- file.path(dir_fig, paste0(name, ".pdf"))
  tiff_file <- file.path(dir_fig, paste0(name, ".tiff"))
  draw <- function() if (is.function(plot)) plot() else print(plot)
  grDevices::cairo_pdf(pdf_file, width = width_cm / 2.54, height = height_cm / 2.54, family = FONT)
  draw(); invisible(dev.off())
  grDevices::tiff(tiff_file, width = width_cm, height = height_cm, units = "cm", res = 300,
                  compression = "lzw", type = "cairo", bg = "white", family = FONT)
  draw(); invisible(dev.off())
  message("saved ", pdf_file, " and ", tiff_file)
  invisible(c(pdf_file, tiff_file))
}

rasterise_if <- function(layer, dpi = 600) {
  if (requireNamespace("ggrastr", quietly = TRUE)) ggrastr::rasterise(layer, dpi = dpi, dev = "ragg") else layer
}

write_src <- function(df, name) {
  f <- file.path(dir_src, paste0(name, ".csv"))
  write.csv(df, f, row.names = FALSE)
  invisible(f)
}
