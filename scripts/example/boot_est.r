nonParamBoot <- function(x, nboot) {
  n <- length(x)

  pgb <- txtProgressBar(max = nboot)
  progress <- function(n) setTxtProgressBar(pgb, n)
  opts <- list(progress = progress)

  res <- foreach(
    iboot = 1:nboot,
    .combine = rbind,
    .packages = "fitdistrplus",
    .options.snow = opts) %dopar% {
    x.boot <- sample(x, n, replace = TRUE)
    model <- fitdist(x.boot, "gamma")
    alpha.boot <- model$estimate["shape"]
    beta.boot <- model$estimate["rate"]

    data.frame(alpha = alpha.boot, beta = beta.boot)
  }
  close(pgb)

  return(res)
}

paramBoot <- function(x, nboot) {
  n <- length(x)
  model <- fitdist(x, "gamma")
  alpha <- model$estimate["shape"]
  beta <- model$estimate["rate"]

  pgb <- txtProgressBar(max = nboot)
  progress <- function(n) setTxtProgressBar(pgb, n)
  opts <- list(progress = progress)
  res <- foreach(
    iboot = 1:nboot,
    .combine = rbind,
    .packages = "fitdistrplus",
    .options.snow = opts) %dopar% {
    x.boot <- rgamma(n, alpha, beta)
    model <- fitdist(x.boot, "gamma")
    alpha.boot <- model$estimate["shape"]
    beta.boot <- model$estimate["rate"]

    data.frame(alpha = alpha.boot, beta = beta.boot)
  }
  close(pgb)

  return(res)
}

n <- 1e4
alpha <- 2
beta <- 0.5
x <- rgamma(n, alpha, beta)

# Asymptotic distribution
fisher.inf <- matrix(rep(0, 4), nrow = 2)
fisher.inf[2, 2] <- alpha / beta**2
fisher.inf[1, 1] <- trigamma(alpha)
fisher.inf[1, 2] <- fisher.inf[2, 1] <- -1 / beta
fisher.inf <- n * fisher.inf

var.cov <- solve(fisher.inf)
alpha.sd <- sqrt(var.cov[1, 1])
beta.sd <- sqrt(var.cov[2, 2])

# Bootstrap estimates
if (opt$recompute) {
  res.nonparam <- nonParamBoot(x, nboot)
  res.param <- paramBoot(x, nboot)
  saveRDS(res.nonparam, file.path(res.dir, "boot_est_nonparam.RDS"))
  saveRDS(res.param, file.path(res.dir, "boot_est_param.RDS"))
} else {
  res.nonparam <- readRDS(file.path(res.dir, "boot_est_nonparam.RDS"))
  res.param <- readRDS(file.path(res.dir, "boot_est_param.RDS"))
}

# Comparison plots
alpha.lbound <- alpha - 3 * alpha.sd
alpha.ubound <- alpha + 4.5 * alpha.sd
beta.lbound <- beta - 3 * beta.sd
beta.ubound <- beta + 4.5 * beta.sd

nboot <- nrow(res.param)

alpha.x <- seq(alpha.lbound, alpha.ubound, length.out = nboot)
beta.x <- seq(beta.lbound, beta.ubound, length.out = nboot)

theme_update(
  legend.position = "top",
  axis.title = element_text(size = 8),
  axis.text = element_text(size = 6),
  legend.text = element_text(size = 8),
  legend.key.size = unit(0.6, "cm"),
  strip.text.x = element_text(size = 8)
)

colours <- c(exact = "blue", bootstrap = "red")

# Parametric
alpha.density.param <- density(
  res.param$alpha,
  n = nboot,
  from = alpha.lbound,
  to = alpha.ubound
)
beta.density.param <- density(
  res.param$beta,
  n = nboot,
  from = beta.lbound,
  to = beta.ubound
)

alpha.y.param <- alpha.density.param$y
beta.y.param <- beta.density.param$y

plotdf <- cbind(
  data.frame(x = c(alpha.x, beta.x)),
  melt(data.frame(
    alpha = alpha.y.param,
    beta = beta.y.param
  ),
  variable.name = "param",
  value.name = "boot",
  variable.factor = TRUE
  ),
  melt(
    data.frame(
      alpha = dnorm(alpha.x, alpha, alpha.sd),
      beta = dnorm(beta.x, beta, beta.sd)
    ),
    value.name = "exact",
    variable.factor = TRUE
  )[, "exact", drop = FALSE]
)

p <- ggplot(plotdf) +
  geom_line(aes(x, boot, colour = "bootstrap")) +
  geom_line(aes(x, exact, colour = "exact")) +
  facet_wrap(vars(param), scales = "free", labeller = label_parsed) +
  scale_colour_manual(
    labels = c(exact = "Exact", bootstrap = "Bootstrap"),
    values = colours,
    name = NULL
  ) +
  ylab("Density") +
  xlab("x")

ggsave(
  file.path(plot.dir, "boot_est_param.eps"),
  p,
  units = "mm",
  height = width.mm / 2.5, # landscape
  width = height.mm
)

# Nonparameteric
alpha.density.nonparam <- density(
  res.nonparam$alpha,
  n = nboot,
  from = alpha.lbound,
  to = alpha.ubound
)
beta.density.nonparam <- density(
  res.nonparam$beta,
  n = nboot,
  from = beta.lbound,
  to = beta.ubound
)

alpha.y.nonparam <- alpha.density.nonparam$y
beta.y.nonparam <- beta.density.nonparam$y

plotdf <- cbind(
  data.frame(x = c(alpha.x, beta.x)),
  melt(data.frame(
    alpha = alpha.y.nonparam,
    beta = beta.y.nonparam
  ),
  variable.name = "param",
  value.name = "boot",
  variable.factor = TRUE
  ),
  melt(
    data.frame(
      alpha = dnorm(alpha.x, alpha, alpha.sd),
      beta = dnorm(beta.x, beta, beta.sd)
    ),
    value.name = "exact",
    variable.factor = TRUE
  )[, "exact", drop = FALSE]
)

p <- ggplot(plotdf) +
  geom_line(aes(x, boot, colour = "bootstrap")) +
  geom_line(aes(x, exact, colour = "exact")) +
  facet_wrap(vars(param), scales = "free", labeller = label_parsed) +
  scale_colour_manual(
    values = colours,
    labels = c(exact = "Exact", bootstrap = "Bootstrap"),
    name = NULL
  ) +
  ylab("Density") +
  xlab("x")

ggsave(
  file.path(plot.dir, "boot_est_nonparam.eps"),
  p,
  units = "mm",
  height = width.mm / 2.5, # landscape
  width = height.mm
)
