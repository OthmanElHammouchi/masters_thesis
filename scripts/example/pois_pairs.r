if (opt$recompute) {
  res.randomised <- glmPairsBoot(triangle, nboot, nsim, keep = "randomised")
  res.rowcol <- glmPairsBoot(triangle, nboot, nsim, keep = "rowcol")
  res.corners <- glmPairsBoot(triangle, nboot, nsim, keep = "corners")
  saveRDS(res.randomised, file.path(res.dir, "pois_pairs_randomised.RDS"))
  saveRDS(res.rowcol, file.path(res.dir, "pois_pairs_rowcol.RDS"))
  saveRDS(res.corners, file.path(res.dir, "pois_pairs_corners.RDS"))
} else {
  res.randomised <- readRDS(file.path(res.dir, "pois_pairs_randomised.RDS"))
  res.rowcol <- readRDS(file.path(res.dir, "pois_pairs_rowcol.RDS"))
  res.corners <- readRDS(file.path(res.dir, "pois_pairs_corners.RDS"))
}

table <- cbind(res.randomised$table, res.rowcol$table$prederror, res.corners$table$prederror)
colnames(table) <- c(names(res.randomised)[-5], "randomised", "rowcol", "corners")
write.csv(
  round(table, 2),
  file = file.path(res.dir, "pois_pairs_table.csv"),
  quote = FALSE,
  row.names = FALSE
)

write.csv(
  round(res.randomised$point, 2),
  file = file.path(res.dir, "pois_pairs_point.csv"),
  quote = FALSE,
  row.names = FALSE
)

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