semiParamBoot <- function(X, y, nboot, x.new, nsim) {
  n <- nrow(X)
  p <- ncol(X)
  model <- lm(y ~ X - 1)
  betas <- model$coefficients
  sigma <- summary(model)$sigma

  leverages <- hatvalues(model)
  resids <- model$residuals / sqrt(1 - leverages)
  resids <- resids - mean(resids)

  pgb <- txtProgressBar(max = nboot)
  progress <- function(n) setTxtProgressBar(pgb, n)
  opts <- list(progress = progress)

  combine <- function(list1, list2) {
    res <- list()
    res$betas.boot <- rbind(list1$betas.boot, list2$betas.boot)
    res$dist <- rbind(list1$dist, list2$dist)
    return(res)
  }

  res <- foreach(
    iboot = 1:nboot,
    .combine = combine,
    .options.snow = opts) %dopar% {
    resids.boot <- sample(resids, n, replace = TRUE)
    y.boot <- X %*% betas + resids.boot

    model.boot <- lm(y.boot ~ X - 1)
    betas.boot <- model.boot$coefficients

    resids.sim <- sample(resids, nsim, replace = TRUE)
    pred.dist <- x.new %*% betas.boot + resids.sim
    fit.dist <- x.new %*% betas + resids.sim

    list(
      betas.boot = betas.boot,
      dist = data.frame(
        fit.dist = fit.dist,
        pred.dist = pred.dist)
    )
  }
  close(pgb)

  return(res)
}

paramBoot <- function(X, y, nboot, x.new, nsim) {
  n <- nrow(X)
  p <- ncol(X)
  model <- lm(y ~ X - 1)
  betas <- model$coefficients
  sigma <- summary(model)$sigma

  pgb <- txtProgressBar(max = nboot)
  progress <- function(n) setTxtProgressBar(pgb, n)
  opts <- list(progress = progress)

  combine <- function(list1, list2) {
    res <- list()
    res$betas.boot <- rbind(list1$betas.boot, list2$betas.boot)
    res$dist <- rbind(list1$dist, list2$dist)
    return(res)
  }

  res <- foreach(
    iboot = 1:nboot,
    .combine = combine,
    .options.snow = opts,
    .packages = "MASS") %dopar% {
    y.boot <- mvrnorm(1, X %*% betas, diag(sigma**2, nrow = n))

    model.boot <- lm(y.boot ~ X - 1)
    betas.boot <- model.boot$coefficients

    pred.dist <- rnorm(n, x.new %*% betas.boot, sigma)
    fit.dist <- rnorm(n, x.new %*% betas, sigma)

    list(
      betas.boot = betas.boot,
      dist = data.frame(
        fit.dist = fit.dist,
        pred.dist = pred.dist)
    )
  }
  close(pgb)

  return(res)
}

p <- 10
n <- 1e3
betas <- c(4.1, 4.4, -2.1, 3.3, 1.4, 0.2, 2.4, -3.7, 1.6, 2.1)
sigma <- 5.7

X <- matrix(
  rep(seq(0, 10, length.out = n), p),
  ncol = p
)
X <- X + matrix(rnorm(n * p), ncol = p) # avoid collinearity
y <- X %*% betas + rnorm(n, 0, sigma)
x.new <- runif(p, 10, 14)

if (opt$recompute) {
  res.semiparam <- semiParamBoot(X, y, nboot, x.new, nsim)
  res.param <- paramBoot(X, y, nboot, x.new, nsim)
  saveRDS(res.semiparam, file.path(res.dir, "boot_reg_semiparam_res.RDS"))
  saveRDS(res.param, file.path(res.dir, "boot_reg_param_res.RDS"))
} else {
  res.semiparam <- readRDS(file.path(res.dir, "boot_reg_semiparam_res.RDS"))
  res.param <- readRDS(file.path(res.dir, "boot_reg_param_res.RDS"))
}

theme_update(
  legend.position = "top",
  axis.title = element_text(size = 10),
  axis.text = element_text(size = 8),
  legend.text = element_text(size = 10),
  legend.key.size = unit(0.6, "cm")
)

plotdf <- melt(res.param$dist, value.name = "dist", variable.name = "type")

p <- ggplot(plotdf) +
  geom_density(
    mapping = aes(x = dist, y = after_stat(scaled), colour = type),
    key_glyph = draw_key_path,
    bw = 2
  ) +
  scale_color_manual(
    labels = c("Fitted distribution", "Predictive distribution"),
    values = c("red", "blue"),
    name = NULL
  ) +
  xlab(TeX("$y_+$")) +
  ylab("Density")

ggsave(
  file.path(plot.dir, "boot_reg_param.eps"),
  p,
  units = "mm",
  height = height.mm / 2,
  width = width.mm
)

plotdf <- melt(res.semiparam$dist, value.name = "dist", variable.name = "type")

p <- ggplot(plotdf) +
  geom_density(
    mapping = aes(x = dist, after_stat(scaled), colour = type),
    key_glyph = draw_key_path,
    bw = 2
  ) +
  scale_color_manual(
    labels = c("Fitted distribution", "Predictive distribution"),
    values = c("red", "blue"),
    name = NULL) +
  xlab(TeX("$y_+$")) +
  ylab("Density")

ggsave(
  file.path(plot.dir, "boot_reg_semiparam.eps"),
  p,
  units = "mm",
  height = height.mm / 2,
  width = width.mm
)
