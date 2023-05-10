poisFit <- function(triangle) {
  ndev <- ncol(triangle)

  long <- as.data.frame(triangle)
  long <- transform(long, origin = factor(origin, levels = dimnames(triangle)[[1]]))
  long <- long[order(long$origin, long$dev), ]

  long.lower <- long[is.na(long$value), ]
  long.upper <- long[!is.na(long$value), ]

  model <- glm(value ~ factor(origin) + factor(dev), quasipoisson(), long.upper, maxit = 1e3)
  a <- model$coefficients[2:ndev]
  b <- model$coefficients[(ndev + 1):(2 * ndev - 1)]
  intercept <- model$coefficients[1]
  disp <- summary(model)$dispersion

  long.lower$value <- predict(model, long.lower, type = "response")
  triangle.proj <- incr2cum(as.triangle(rbind(long.upper, long.lower)))

  latest <- triangle.proj[
    row(triangle.proj) + col(triangle.proj) == ndev + 1][1:(ndev - 1)
  ]
  reserve <- triangle.proj[2:ndev, ndev] - rev(latest)

  resids.pears <- resid(model, type = "pearson")
  resids.dev <- sign(model$y - model$fitted.values) *
    sqrt(zapsmall(2 * (model$y * log(model$y / model$fitted.values) + model$fitted.values - model$y)))
  resids.quant <- ppois(model$y, model$fitted.values)

  param.var.quasi <- summary(model)$cov.scaled # (X^T X)^-1 * \phi
  param.var.pois <- summary(model)$cov.unscaled # (X^T X)^-1
  model.terms <- delete.response(terms(model))
  X <- model.matrix(model.terms, long.lower, xlev = model$xlevels)

  eta.vars.pois <- list()
  eta.vars.quasi <- list()
  idx <- 1
  for (i in 2:ndev) {
    ncols <- i - 1
    eta.vars.pois[[i - 1]] <- X[idx:(idx + ncols - 1), , drop = FALSE] %*% param.var.pois %*% t(X[idx:(idx + ncols - 1), , drop = FALSE]) # nolint
    eta.vars.quasi[[i - 1]] <- X[idx:(idx + ncols - 1), , drop = FALSE] %*% param.var.quasi %*% t(X[idx:(idx + ncols - 1), , drop = FALSE]) # nolint
    idx <- idx + ncols
  }

  triangle.lower <- as.triangle(long.lower)
  pred.error.pois <- rep(0, ndev - 1)
  pred.error.quasi <- rep(0, ndev - 1)
  for (i in 1:(ndev - 1)) {
    mu <- triangle.lower[i, !is.na(triangle.lower[i, ]), drop = FALSE]
    pred.error.pois[i] <- sqrt(sum(triangle.lower[i, ], na.rm = TRUE) + mu %*% eta.vars.pois[[i]] %*% t(mu))
    pred.error.quasi[i] <- sqrt(disp * sum(triangle.lower[i, ], na.rm = TRUE) + mu %*% eta.vars.quasi[[i]] %*% t(mu))
  }

  res <- list(
    table = data.frame(
      idx = 2:ndev,
      a = a,
      b = b,
      reserve = reserve,
      prederrorpois = pred.error.pois,
      prederrorquasi = pred.error.quasi
    ),
    point = data.frame(
      intercept = intercept,
      disp = disp
    ),
    resids = data.frame(
      resids.pears = resids.pears,
      resids.dev = resids.dev,
      resids.quant = resids.quant
    ),
    fitted = model$fitted.values,
    model = model
  )

  return(res)
}

long <- as.data.frame(triangle)
long <- long[!is.na(long$value), ]
res <- poisFit(triangle)

write.csv(
  round(res$table, 2),
  file = file.path(res.dir, "pois_bench_table.csv"),
  quote = FALSE,
  row.names = FALSE
)

write.csv(
  round(res$point, 2),
  file = file.path(res.dir, "pois_bench_point.csv"),
  quote = FALSE,
  row.names = FALSE
)

write.csv(
  long,
  file = file.path(res.dir, "UKMotor_long.csv"),
  quote = FALSE,
  row.names = FALSE
)

plotdf <- cbind(
  melt(res$resids, variable.name = "type", value.name = "resid"),
  fitted = rep(res$fitted, 3)
)
nbins <- ceiling(2 * nrow(plotdf)**(1/3))
p <- ggplot(plotdf) +
  geom_histogram(aes(x = resid), bins = 8) +
  facet_wrap(vars(type), scales = "free")

ggsave(
  file.path(plot.dir, "pois_bench_resids.eps"),
  p,
  units = "mm",
  height = height.mm / 3,
  width = width.mm
)

# diffs <- rep(0, 100)
# for (iexp in 1:100) {
#   triangle <- triangle
#   ndev <- ncol(triangle)
#   for (j in 1:ndev) {
#     col.total <- sum(triangle[, j], na.rm = TRUE)
#     row.total <- sum(triangle[ndev + 1 - j, ], na.rm = TRUE)
#     max.total <- min(col.total, row.total, min(triangle[, 1]) / 3)
#     nmod <- min(3, ndev + 1 - j)
#     subtract <- round(runif(1, 0, max.total))
#     idxs <- sample(1:(ndev + 1 - j), 3, replace = TRUE)
#     triangle[idxs, j] <- triangle[idxs, j] - subtract / nmod
#   }
#   suppressWarnings({
#   bench <- summary(MackChainLadder(incr2cum(triangle)))$ByOrigin$IBNR[-1]
#   })
#   res <- poisFit(triangle)
#   if (any(is.na(res$resids$resids.dev))) { stop("Problem!") }
#   reserve <- res$table$reserve
#   diffs[iexp] <- norm(reserve - bench, type = "2")
# }
