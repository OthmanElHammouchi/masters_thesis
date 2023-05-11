library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(latex2exp)

plot.dir <- "plots"

theme_set(theme_bw())

# in mm
width <- 410
height <- 630

mm.per.pt <- 0.3528
width.mm <- mm.per.pt * width
height.mm <- mm.per.pt * height

theme_update(
  plot.title = element_text(hjust = 0.5, size = 10),
  axis.title = element_text(size = 8),
  axis.text = element_text(size = 6),
  legend.text = element_text(size = 8),
  legend.key.size = unit(0.6, "cm"),
  strip.text.x = element_text(size = 6)
)

single.res <- readRDS("results/glm_single.RDS")

factors <- unique(single.res$factor)
i <- 3
j <- 3

plot.list <- list()
for (f in factors) {
  contaminated <- single.res[(excl.rowidx != outlier.rowidx | excl.colidx != outlier.colidx) & factor == f]
  clean <- single.res[excl.rowidx == outlier.rowidx & excl.colidx == outlier.colidx & factor == f]

  n <- with(
    contaminated[outlier.rowidx == i & outlier.colidx == j],
    length(unique(interaction(excl.rowidx, excl.colidx))))
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

  p <- ggplot() +
    geom_density(
      aes(reserve, colour = interaction(excl.rowidx, excl.colidx)),
      contaminated[outlier.rowidx == i & outlier.colidx == j],
      key_glyph = draw_key_path
    ) +
    geom_density(aes(reserve),
      clean[outlier.rowidx == i & outlier.colidx == j],
      colour = "black",
      key_glyph = draw_key_path) +
    scale_color_manual(values = getPalette(n)) +
    ggtitle(TeX(sprintf(r"($c_\lambda$ = %3.2f)", f))) +
    guides(colour = guide_legend(title = "Excluded point")) +
    xlab("Reserve") +
    ylab("Density")

  plot.list <- c(plot.list, list(p))
}

p <- do.call(wrap_plots, plot.list) + plot_layout(guides = "collect")

ggsave(
  file.path(plot.dir, "glm_sim_by_factor.eps"),
  p,
  units = "mm",
  height = width.mm,
  width = height.mm
)
