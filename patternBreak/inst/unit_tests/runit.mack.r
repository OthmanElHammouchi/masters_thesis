test.mackBoot <- function() {
  config <- expand.grid(
    list(
      c("parametric", "residuals", "pairs"),
      c("normal", "gamma"),
      c(TRUE, FALSE),
      c("standardised", "modified", "studentised", "lognormal")
    ),
    stringsAsFactors = FALSE
  )

  for (i in nrow(config)) {
    boot.type <- config[i, 1]
    proc.dist <- config[i, 2]
    conditional <- config[i, 3]
    resids.type <- config[i, 4]
    reserve <- mackBoot(
      test.triangle,
      1e3,
      boot.type,
      proc.dist,
      conditional,
      resids.type
    )
  }

  checkTrue(!any(is.na(reserve)))
}


test.mackSim <- function() {
  factors <- seq(0.5, 1.5, by = 0.25)
  boot_types <- c("parametric", "residuals", "pairs")
  proc_dists <- c("normal", "gamma")
  conds <- c(TRUE, FALSE)
  resids_types <- c("standardised", "modified", "studentised", "lognormal")

  res <- mackSim(test.triangle,
    "single",
    1e3,
    factors,
    boot_types,
    proc_dists,
    conds,
    resids_types,
    show_progress = FALSE
  )

  checkTrue(!any(is.na(res$reserve)))
}
