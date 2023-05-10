if (opt$recompute) {
  res.cond.standard <- mackSemiParamBoot(UKMotor, nboot, nsim, cond = TRUE, resids.type = "standardised")
  res.uncond.standard <- mackSemiParamBoot(UKMotor, nboot, nsim, cond = FALSE, resids.type = "standardised")
  res.cond.student <- mackSemiParamBoot(UKMotor, nboot, nsim, cond = TRUE, resids.type = "studentised")
  res.uncond.student <- mackSemiParamBoot(UKMotor, nboot, nsim, cond = FALSE, resids.type = "studentised")
  res.cond.lognorm <- mackSemiParamBoot(UKMotor, nboot, nsim, cond = TRUE, resids.type = "log-normal")
  res.uncond.lognorm <- mackSemiParamBoot(UKMotor, nboot, nsim, cond = FALSE, resids.type = "log-normal")
  saveRDS(res.cond.standard, file.path(res.dir, "mack_semiparam_cond_standard.RDS"))
  saveRDS(res.uncond.standard, file.path(res.dir, "mack_semiparam_uncond_standard.RDS"))
  saveRDS(res.cond.student, file.path(res.dir, "mack_semiparam_cond_student.RDS"))
  saveRDS(res.uncond.student, file.path(res.dir, "mack_semiparam_uncond_student.RDS"))
  saveRDS(res.cond.lognorm, file.path(res.dir, "mack_semiparam_cond_log_normal.RDS"))
  saveRDS(res.uncond.lognorm, file.path(res.dir, "mack_semiparam_uncond_log_normal.RDS"))
} else {
  res.cond.standard <- readRDS(file.path(res.dir, "mack_semiparam_cond_standard.RDS"))
  res.uncond.standard <- readRDS(file.path(res.dir, "mack_semiparam_uncond_standard.RDS"))
  res.cond.student <- readRDS(file.path(res.dir, "mack_semiparam_cond_student.RDS"))
  res.uncond.student <- readRDS(file.path(res.dir, "mack_semiparam_uncond_student.RDS"))
  res.cond.lognorm <- readRDS(file.path(res.dir, "mack_semiparam_cond_log_normal.RDS"))
  res.uncond.lognorm <- readRDS(file.path(res.dir, "mack_semiparam_uncond_log_normal.RDS"))
}

res.list <- list(res.cond.standard,
  res.uncond.standard,
  res.cond.student,
  res.uncond.student,
  res.cond.lognorm,
  res.uncond.lognorm
)

names(res.list) <- c(
  "cond_standard",
  "uncond_standard",
  "cond_student",
  "uncond_student",
  "cond_log_normal",
  "uncond_log_normal")

for (i in seq_len(length(res.list))) {
  write.csv(
    round(res.list[[i]]$est, 2),
    file = file.path(res.dir, paste0("mack_semiparam_", names(res.list)[i], ".csv")),
    quote = FALSE,
    row.names = FALSE
  )
}

theme_update(
  legend.position = "top",
  axis.title = element_text(size = 8),
  axis.text = element_text(size = 6),
  legend.text = element_text(size = 8),
  legend.key.size = unit(0.6, "cm")
)

for (i in seq_along(res.list)) {
  res <- res.list[[i]]
  plotdf <- melt(res$dist, value.name = "reserve", variable.name = "type")

  p <- ggplot(plotdf) +
    geom_density(
      aes(x = reserve, colour = type),
      key_glyph = draw_key_path
    ) +
    scale_colour_manual(
      labels = c(fit.dist = "Fitted", pred.dist = "Predictive"),

      values = c(fit.dist = "red", pred.dist = "blue"),
      name = NULL
    ) +
    xlab("Reserve") +
    ylab("Density")

  ggsave(
    file.path(plot.dir, paste0("mack_semiparam_", names(res.list)[i], ".eps")),
    p,
    units = "mm",
    height = height.mm / 3,
    width = width.mm / 2.2
  )
}

# What about the penultimate column?
