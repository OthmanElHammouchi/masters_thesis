if (opt$recompute) {
  res.norm <- glmParamBoot(triangle, nboot, nsim, "normal")
  res.gamma <- glmParamBoot(triangle, nboot, nsim, "gamma")
  res.pois <- glmParamBoot(triangle, nboot, nsim, "pois")
  saveRDS(res.norm, file.path(res.dir, "pois_param_norm.RDS"))
  saveRDS(res.norm, file.path(res.dir, "pois_param_gamma.RDS"))
  saveRDS(res.norm, file.path(res.dir, "pois_param_pois.RDS"))
} else {
  res.norm <- readRDS(file.path(res.dir, "pois_param_norm.RDS"))
  res.gamma <- readRDS(file.path(res.dir, "pois_param_gamma.RDS"))
  res.pois <- readRDS(file.path(res.dir, "pois_param_pois.RDS"))
}

res.list <- list(res.norm, res.gamma, res.pois)
names(res.list) <- c("norm", "gamma", "pois")

for (i in seq_along(res.list)) {
  res <- res.list[[i]]
  write.csv(
    round(res$table, 2),
    file = file.path(res.dir, paste0("pois_param_", names(res.list)[i], "_table.csv")),
    quote = FALSE,
    row.names = FALSE
  )
  write.csv(
    round(res$point, 2),
    file = file.path(res.dir, paste0("pois_param_", names(res.list)[i], "_point.csv")),
    quote = FALSE,
    row.names = FALSE
  )
}

theme_update(
  legend.position = c(0.8, 0.8),
  legend.background = element_rect(fill = NA),
  legend.title = element_text(size = 10),
  axis.title = element_text(size = 8),
  axis.text = element_text(size = 6),
  legend.text = element_text(size = 8),
  legend.key.size = unit(0.3, "cm"),
  strip.text.x = element_text(size = 8)
)

plotdf <- data.frame(
  Normal = res.norm$pred.dist,
  Gamma = res.gamma$pred.dist,
  Poisson = res.pois$pred.dist
)
plotdf <- melt(plotdf, variable.name = "dist", value.name = "pred.dist")

p <- ggplot(plotdf) +
  geom_density(
    aes(x = pred.dist, colour = dist),
    key_glyph = draw_key_path
  ) +
  scale_color_brewer(palette = "Dark2", name = "Distribution") +
  xlab("Reserve") +
  ylab("Density")

ggsave(
  file.path(plot.dir, "pois_param.eps"),
  p,
  units = "mm",
  height = height.mm / 3,
  width = width.mm
)
