# Theme and palettes for all article figures (IJMS style). Do not edit values:
# they are fixed by the figure specification.

library(ggplot2)

theme_ijms <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.title = element_text(size = 8),
        axis.line = element_line(linewidth = 0.18),
        axis.ticks = element_line(linewidth = 0.18),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 7),
        legend.key.size = unit(3, "mm"),
        plot.tag = element_text(size = 10, face = "bold"))

grp_fill  <- c(Sham = "#FFFFFF", rLLB = "#595959")   # bars: white with black outline / dark grey
pt_fill   <- c(Sham = "#000000", rLLB = "#FFFFFF")   # animal points: black / white with black outline
pca_col   <- c(Sham = "#A6A6A6", rLLB = "#000000")   # PCA: light-grey circle / black triangle
dir_col   <- c(Down = "#2166AC", Up = "#B2182B", NS = "#BFBFBF")
heat_cols <- c("#2166AC", "#F7F7F7", "#B2182B")      # z-score from -2 to 2, values clipped

# Shared geometry constants (ggplot units).
# ggplot point stroke -> pt: lwd = stroke * .stroke / 2, 1 lwd = 1/96 in = 0.75 pt.
PT_STROKE <- 0.5 / 0.75 * 2 / .stroke   # 0.5 pt outline
PT_SIZE   <- 1.5                        # ~1.5 mm animal points
ERR_LW    <- 0.26                       # error-bar linewidth
BAR_W     <- 0.7                        # bar width
CAP_W     <- 0.35 * BAR_W               # cap = 35 % of bar width
JIT_W     <- 0.15 * BAR_W               # jitter <= 0.15 bar width
FONT_MM   <- function(pt) pt / .pt      # font size in pt -> geom_text size (mm)

# General rules (section "Общие правила"): tick labels and legend text are 7 pt.
# theme_classic() would otherwise draw them at rel(0.8) = 5.6 pt.
# Markdown titles (ggtext) let sub/superscripts be set to an explicit 6 pt.
theme_ijms <- theme_ijms +
  theme(axis.text = element_text(size = 7, colour = "black"),
        legend.text = element_text(size = 7),
        axis.title.x = ggtext::element_markdown(size = 8),
        axis.title.y = ggtext::element_markdown(size = 8))
