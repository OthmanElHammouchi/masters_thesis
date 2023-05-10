if (opt$recompute) {
  res.pears <- glmSemiParamBoot(triangle, nboot, nsim, resids.type = "pearson")
  res.dev <- glmSemiParamBoot(triangle, nboot, nsim, resids.type = "deviance")
  # res.quant <- semiParamBoot(triangle, nboot, nsim, resids.type = "quantile") # not quasi
  saveRDS(res.pears, file.path(res.dir, "pois_semiparam_pears.RDS"))
  saveRDS(res.dev, file.path(res.dir, "pois_semiparam_dev.RDS"))
  # saveRDS(res.quant, file.path(res.dir, "pois_semiparam_quant.RDS"))

} else {
  res.pears <- readRDS(file.path(res.dir, "pois_semiparam_pears.RDS"))
  res.dev <- readRDS(file.path(res.dir, "pois_semiparam_dev.RDS"))
  # res.quant <- readRDS(file.path(res.dir, "pois_semiparam_quant.RDS")) # not quasi
}

res.list <- list(res.pears, res.dev)
names(res.list) <- c("pears", "dev")

for (i in seq_along(res.list)) {
  res <- res.list[[i]]
  write.csv(
    round(res$table, 2),
    file = file.path(res.dir, paste0("pois_semiparam_", names(res.list)[i], "_table.csv")),
    quote = FALSE,
    row.names = FALSE
  )
  write.csv(
    round(res$point, 2),
    file = file.path(res.dir, paste0("pois_semiparam_", names(res.list)[i], "_point.csv")),
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
  Pearson = res.pears$pred.dist,
  Deviance = res.dev$pred.dist
  )
plotdf <- melt(plotdf, variable.name = "resids.type", value.name = "pred.dist")

p <- ggplot(plotdf) +
  geom_density(
    aes(x = pred.dist, colour = resids.type),
    key_glyph = draw_key_path
  ) +
  scale_color_brewer(palette = "Dark2", name = "Residuals type") +
  xlab("Reserve") +
  ylab("Density")

ggsave(
  file.path(plot.dir, "pois_semiparam.eps"),
  p,
  units = "mm",
  height = height.mm / 2,
  width = width.mm
)
