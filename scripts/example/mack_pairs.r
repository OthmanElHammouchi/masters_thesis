if (opt$recompute) {
  res.param <- mackPairsBoot(UKMotor, nboot, nsim, sim.type = "parametric")
  res.semiparam <- mackPairsBoot(UKMotor, nboot, nsim, sim.type = "semiparametric")
  saveRDS(res.param, file.path(res.dir, "mack_pairs_param.RDS"))
  saveRDS(res.semiparam, file.path(res.dir, "mack_pairs_semiparam.RDS"))
} else {
  res.param <- readRDS(file.path(res.dir, "mack_pairs_param.RDS"))
  res.semiparam <- readRDS(file.path(res.dir, "mack_pairs_semiparam.RDS"))

}

write.csv(
  round(res.param$est, 2),
  file = file.path(res.dir, "mack_pairs_param.csv"),
  quote = FALSE,
  row.names = FALSE
)

write.csv(
  round(res.semiparam$est, 2),
  file = file.path(res.dir, "mack_pairs_semiparam.csv"),
  quote = FALSE,
  row.names = FALSE
)

theme_update(
  legend.position = "top",
  axis.title = element_text(size = 8),
  axis.text = element_text(size = 6),
  legend.text = element_text(size = 8),
  legend.key.size = unit(0.6, "cm")
)

plotdf <- melt(res.param$dist, value.name = "dist", variable.name = "type")

p <- ggplot(plotdf) +
  geom_density(
    aes(x = dist, colour = type),
    key_glyph = draw_key_path
  ) +
  scale_color_manual(
    labels = c("Fitted", "Predictive"),
    values = c("red", "blue"),
    name = NULL
  ) +
  xlab("Reserve") +
  ylab("Density")

ggsave(
  file.path(plot.dir, "mack_pairs_param.eps"),
  p,
  units = "mm",
  height = height.mm / 3,
  width = width.mm / 2.2
)

plotdf <- melt(res.semiparam$dist, value.name = "dist", variable.name = "type")

p <- ggplot(plotdf) +
  geom_density(
    aes(x = dist, colour = type),
    key_glyph = draw_key_path
  ) +
  scale_color_manual(
    labels = c("Fitted", "Predictive"),
    values = c("red", "blue"),
    name = NULL
  ) +
  xlab("Reserve") +
  ylab("Density")

ggsave(
  file.path(plot.dir, "mack_pairs_semiparam.eps"),
  p,
  units = "mm",
  height = height.mm / 3,
  width = width.mm / 2.2
)

# In the later columns it's very possible to draw the same pair
# for all resamples and then you have no variability.
