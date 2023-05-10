if (opt$recompute) {
  res.cond.norm <- mackParamBoot(UKMotor, nboot, nsim, cond = TRUE, dist = "normal")
  res.uncond.norm <- mackParamBoot(UKMotor, nboot, nsim, cond = FALSE, dist = "normal") 

  res.cond.gamma <- mackParamBoot(UKMotor, nboot, nsim, cond = TRUE, dist = "gamma")
  res.uncond.gamma <- mackParamBoot(UKMotor, nboot, nsim, cond = FALSE, dist = "gamma") 

  saveRDS(res.cond.norm, file.path(res.dir, "mack_param_cond_norm.RDS"))
  saveRDS(res.uncond.norm, file.path(res.dir, "mack_param_uncond_norm.RDS"))
  saveRDS(res.cond.gamma, file.path(res.dir, "mack_param_cond_gamma.RDS"))
  saveRDS(res.uncond.gamma, file.path(res.dir, "mack_param_uncond_gamma.RDS"))
} else {
  res.cond.norm <- readRDS(file.path(res.dir, "mack_param_cond_norm.RDS"))
  res.uncond.norm <- readRDS(file.path(res.dir, "mack_param_uncond_norm.RDS"))
  res.cond.gamma <- readRDS(file.path(res.dir, "mack_param_cond_gamma.RDS"))
  res.uncond.gamma <- readRDS(file.path(res.dir, "mack_param_uncond_gamma.RDS"))
}

res.list <- list(res.cond.norm, res.uncond.norm, res.cond.gamma, res.uncond.gamma)
names(res.list) <- c("cond_norm", "uncond_norm", "cond_gamma", "uncond_gamma")

for (i in seq_along(res.list)) {
  res <- res.list[[i]]
  write.csv(
  round(res$est, 2),
  file = file.path(res.dir, paste0("mack_param_", names(res.list)[i], ".csv")),
  quote = FALSE,
  row.names = FALSE
)
}

theme_update(
  legend.position = "top",
  axis.title = element_text(size = 8),
  axis.text = element_text(size = 6),
  legend.text = element_text(size = 8),
  legend.key.size = unit(0.6, "cm"),
  strip.text.x = element_text(size = 8)
)

plotdf <- rbind(
  cbind(
    melt(res.cond.norm$dist, value.name = "reserve", variable.name = "type"),
    data.frame(dist = rep("Normal", nrow(res.cond.norm$dist)))
  ),
  cbind(
    melt(res.cond.gamma$dist, value.name = "reserve", variable.name = "type"),
    data.frame(dist = rep("Gamma", nrow(res.cond.norm$dist)))
  )
)

p <- ggplot(plotdf) +
  geom_density(
    aes(x = reserve, colour = type),
    key_glyph = draw_key_path
  ) +
  facet_wrap(vars(dist), scales = "free") +
  scale_colour_manual(
    labels = c(fit.dist = "Fitted", pred.dist = "Predictive"), 

    values = c(fit.dist = "red", pred.dist = "blue"),
    name = NULL
    ) +
  xlab("x") +
  ylab("Density")

ggsave(
  file.path(plot.dir, "mack_param_cond.eps"),
  p,
  units = "mm",
  height = width.mm / 2.5, # landscape
  width = height.mm
)

plotdf <- rbind(
  cbind(
    melt(res.uncond.norm$dist, value.name = "reserve", variable.name = "type"),
    data.frame(dist = rep("Normal", nrow(res.uncond.norm$dist)))
  ),
  cbind(
    melt(res.uncond.gamma$dist, value.name = "reserve", variable.name = "type"),
    data.frame(dist = rep("Gamma", nrow(res.uncond.norm$dist)))
  )
)

p <- ggplot(plotdf) +
  geom_density(
    aes(x = reserve, colour = type),
    key_glyph = draw_key_path
  ) +
  facet_wrap(vars(dist), scales = "free") +
  scale_colour_manual(
    labels = c(fit.dist = "Fitted", pred.dist = "Predictive"), 

    values = c(fit.dist  = "red", pred.dist = "blue"),
    name = NULL
    ) +
  xlab("x") +
  ylab("Density")

ggsave(
  file.path(plot.dir, "mack_param_uncond.eps"),
  p,
  units = "mm",
  height = width.mm / 2.5, # landscape
  width = height.mm
)
