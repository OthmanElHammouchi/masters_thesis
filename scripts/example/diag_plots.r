# bench <- MackChainLadder(UKMotor)

# bench.std.resids <- as.triangle(residuals(bench),
#   origin = "origin.period",
#   dev = "dev.period",
#   value = "standard.residuals")

# bench.resids <- as.triangle(residuals(bench),
#   origin = "origin.period",
#   dev = "dev.period",
#   value = "residuals")

# bench.fitted <- as.triangle(residuals(bench),
#   origin = "origin.period",
#   dev = "dev.period",
#   value = "fitted.value")

diagPlotsMack <- function(triangle) {
  ndev <- ncol(triangle)
  npts <- (ndev**2 + ndev) / 2
  nresids <- (ndev**2 - ndev) / 2
  bench <- suppressWarnings(MackChainLadder(triangle))

  std.resids <- matrix(rep(NA, (ndev - 1)**2), ncol = ndev - 1)
  for (j in seq_len(ndev - 1)) {
    nrows <- ndev - j
    std.resids[1:nrows, j] <- (triangle[1:nrows, j + 1] - bench$f[j] *
      triangle[1:nrows, j]) / (bench$sigma[j] * sqrt(triangle[1:nrows, j]) *
      sqrt(1 - triangle[1:nrows, j] / sum(triangle[1:nrows, j])))
  }

  resids <- matrix(rep(NA, (ndev - 1)**2), ncol = ndev - 1)
  for (j in seq_len(ndev - 1)) {
    nrows <- ndev - j
    resids[1:nrows, j] <- triangle[1:nrows, j + 1] - bench$f[j] *
      triangle[1:nrows, j]
  }

  fitted <- matrix(rep(NA, ndev**2), ncol = ndev)
  fitted[, 1] <- triangle[, 1]
  for (j in seq_len(ndev - 1)) {
    nrows <- ndev - j
    fitted[1:nrows, j + 1] <- bench$f[j] * triangle[1:nrows, j]
  }

  plot.data <- data.frame(
    origin = integer(nresids - 3),
    dev = integer(nresids - 3),
    fitted = numeric(nresids - 3),
    resid = numeric(nresids - 3),
    std.resid = numeric(nresids - 3)
  )

  k <- 1
  for (j in seq(2, ndev - 2)) {
    for (i in seq_len(ndev + 1 - j)) {
      plot.data$origin[k] <- i
      plot.data$dev[k] <- j
      plot.data$fitted[k] <- fitted[i, j]
      plot.data$resid[k] <- resids[i, j - 1]
      plot.data$std.resid[k] <- std.resids[i, j - 1]
      k <- k + 1
    }
  }

  res <- list(
    plot.data = plot.data,
    resids = resids,
    std.resids = std.resids,
    fitted = fitted
  )

  return(res)
}

singleOutlier <- function(outlier.rowidx, outlier.colidx, factor, ...,
                          initcol,
                          devfacs,
                          sigmas,
                          dist = "normal") {

  ndev <- length(initcol)

  if (outlier.colidx == 1) {
    stop("Outlier column index must be greater than 1.")
  }

  if (dist == "normal") {
    triangle <- matrix(ncol = ndev, nrow = ndev)
    triangle[, 1] <- initcol

    for (colidx in 2:ndev) {
      for (rowidx in setdiff(1:(ndev + 1 - colidx), outlier.rowidx)) {
        prevc <- triangle[rowidx, colidx - 1]
        triangle[rowidx, colidx] <-
          rnorm(1,
            devfacs[colidx - 1] * prevc,
            sigmas[colidx - 1] * sqrt(prevc))
      }
    }

    if (outlier.colidx > 2) {
      for (colidx in 2:(outlier.colidx - 1)) {
        prevc <- triangle[outlier.rowidx, colidx - 1]
        triangle[outlier.rowidx, colidx] <- rnorm(1,
          devfacs[colidx - 1] * prevc,
          sigmas[colidx - 1] * sqrt(prevc))
      }
    }

    prevc <- triangle[outlier.rowidx, outlier.colidx - 1]
    triangle[outlier.rowidx, outlier.colidx] <-
      rnorm(
        1,
        factor * devfacs[outlier.colidx - 1] * prevc,
        sigmas[outlier.colidx - 1] * sqrt(prevc)
      )

    if (outlier.colidx < ndev) {
      begin <- outlier.colidx + 1
      end <- ndev + 1 - outlier.rowidx
      for (colidx in seq(begin, end, length = max(0, end - begin + 1))) {
        prevc <- triangle[outlier.rowidx, colidx - 1]
        triangle[outlier.rowidx, colidx] <- rnorm(1,
          devfacs[colidx - 1] * prevc,
          sigmas[colidx - 1] * sqrt(prevc))
      }
    }
    return(triangle)
  } else {
    triangle <- matrix(ncol = ndev, nrow = ndev)
    triangle[, 1] <- initcol

    for (colidx in 2:ndev) {
      for (rowidx in setdiff(1:(ndev + 1 - colidx), outlier.rowidx)) {
        prevc <- triangle[rowidx, colidx - 1]
        alpha <- devfacs[colidx - 1]**2 * prevc / sigmas[colidx - 1]**2
        beta <- devfacs[colidx - 1] / (sigmas[colidx - 1]**2)

        triangle[rowidx, colidx] <-
          rgamma(1, shape = alpha, rate = beta)
      }
    }

    if (outlier.colidx > 2) {
      for (colidx in 2:(outlier.colidx - 1)) {
        prevc <- triangle[outlier.rowidx, colidx - 1]
        alpha <- devfacs[colidx - 1]**2 * prevc / sigmas[colidx - 1]**2
        beta <- devfacs[colidx - 1] / (sigmas[colidx - 1]**2)

        triangle[outlier.rowidx, colidx] <-
          rgamma(1, shape = alpha, rate = beta)
      }
    }

    prevc <- triangle[outlier.rowidx, outlier.colidx - 1]
    alpha <- factor * devfacs[outlier.colidx - 1]**2 * prevc /
      sigmas[outlier.colidx - 1]**2
    beta <- devfacs[outlier.colidx - 1] / (sigmas[outlier.colidx - 1]**2)

    triangle[outlier.rowidx, outlier.colidx] <-
      rgamma(1, shape = alpha, rate = beta)

    if (outlier.colidx < ndev) {
      for (colidx in (outlier.colidx + 1):(ndev + 1 - outlier.rowidx)) {
        prevc <- triangle[outlier.rowidx, colidx - 1]
        alpha <- devfacs[colidx - 1]**2 * prevc / sigmas[colidx - 1]**2
        beta <- devfacs[colidx - 1] / (sigmas[colidx - 1]**2)

        triangle[outlier.rowidx, colidx] <-
          rgamma(1, shape = alpha, rate = beta)
      }
    }

    return(triangle)
  }
}

bench <- suppressWarnings(MackChainLadder(UKMotor))
plots.pert <- list()
plots.mod.pert <- list()

# Examples which best illustrate the point:
points <- list(c(2, 5), c(4, 4))

for (k in seq_along(points)) {

  i <- points[[k]][1]
  j <- points[[k]][2]

  mod.pert <- singleOutlier(i, j, 1.5,
    initcol = UKMotor[, 1],
    devfacs = bench$f[-length(bench$f)],
    sigmas = bench$sigma
  )

  pert <- UKMotor
  pert[i, j] <- 1.5 * pert[i, j]

  plot.data.pert <- diagPlotsMack(pert)$plot.data
  plot.data.mod.pert <- diagPlotsMack(pert)$plot.data

  point.pert <- plot.data.pert[
    plot.data.pert$origin == i & plot.data.pert$dev == j,
    c("fitted", "std.resid")
  ]

  point.mod.pert <- plot.data.mod.pert[
    plot.data.mod.pert$origin == i &
      plot.data.mod.pert$dev == j,
    c("fitted", "std.resid")
  ]

  plot.pert <- ggplot(data = plot.data.pert, aes(fitted, std.resid)) +
    geom_point(size = 1) +
    geom_point(size = 1,
      data = point.pert,
      mapping = aes(fitted, std.resid),
      colour = "red") +
    geom_text_repel(
      size = 2.5,
      mapping = aes(label = paste0("(", origin, ",", dev, ")"))
    ) +
    ylim(c(-3, 3)) +
    ylab("Standardised residuals") +
    xlab("Fitted values")

  plot.mod.pert <- ggplot(data = plot.data.mod.pert, aes(fitted, std.resid)) +
    geom_point(size = 1) +
    geom_point(size = 1,
      data = point.mod.pert,
      aes(fitted, std.resid),
      colour = "red") +
    geom_text_repel(
      size = 2.5,
      mapping = aes(label = paste0("(", origin, ",", dev, ")"))
    ) +
    ylim(c(-3, 3)) +
    ylab("Standardised residuals") +
    xlab("Fitted values")

  plots.pert <- c(plots.pert, list(plot.pert))
  plots.mod.pert <- c(plots.mod.pert, list(plot.mod.pert))
}

plots.pert <- wrap_plots(plots.pert, ncol = 2)
plots.mod.pert <- wrap_plots(plots.mod.pert, ncol = 2)

theme_update(
  plot.title = element_text(size = 9),
  axis.title = element_text(size = 8),
  axis.text = element_text(size = 6))

# in cm
width <- 14.69785
height <- 20.78696

ggsave(
  file.path(plot.dir, "perturbed_resids.eps"),
  plots.pert,
  units = "cm",
  height = height / 3,
  width = width
)

ggsave(
  file.path(plot.dir, "model_perturbed_resids.eps"),
  plots.mod.pert,
  units = "cm",
  height = height / 3,
  width = width
)

plot.data <- diagPlotsMack(UKMotor)$plot.data

plot.original <- ggplot(data = plot.data, aes(fitted, std.resid)) +
  geom_point(size = 1) +
  geom_text_repel(
    size = 2.5,
    mapping = aes(label = paste0("(", origin, ",", dev, ")"))
  ) +
  ylim(c(-3, 3)) +
  ylab("Standardised residuals") +
  xlab("Fitted values")

theme_update(
  plot.title = element_text(size = 11),
  axis.title = element_text(size = 9),
  axis.text = element_text(size = 7))

ggsave(
  file.path(plot.dir, "original_resids.eps"),
  plot.original,
  units = "cm",
  height = height / 3,
  width = width / 1.5
)
