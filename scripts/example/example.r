#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(ggplot2); library(latex2exp);
  library(patchwork); library(Hmisc);
  library(fitdistrplus); library(reshape2); library(MASS);
  library(ChainLadder); library(Rmpfr); library(RColorBrewer);
  library(ggrepel); library(gridExtra); library(claimsBoot)
})

option_list <- list(
  make_option(
    c("-r", "--recompute"),
    type = "logical",
    action = "store_true",
    default = FALSE,
    help = "recompute results"),
  make_option(
    c("--nboot"),
    type = "integer",
    default = 1e2,
    help = "number of bootstrap iterations"
  ),
  make_option(
    c("--nsim"),
    type = "integer",
    default = 10,
    help = "number of simulations for the predictive distribution"
  )
)

opt.parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt.parser)

nboot <- opt$nboot
nsim <- opt$nsim

plot.dir <- file.path("plots", "example")
res.dir <- file.path("results", "example")

if (opt$recompute) {
  set.seed(42)
  UKMotor.incr <- cum2incr(UKMotor)
  ## Mack benchmark
  res <- mackFit(UKMotor)

  ## GLM benchmark
  res <- glmFit(UKMotor.incr)

  ## Mack
  mack.names <- c("$j$",
    "$\\widehat{f}^B_j$",
    "$\\widehat{\\sigma}^B_j$",
    "$\\widehat{R}^B_j$",
    "\\resizebox{4em}{!}{$\\widehat{\\mathrm{MSEP}}(\\widehat{R}^B_j)$}"
  )
  # semiparametric
  options <- expand.grid(cond = c(TRUE, FALSE), resids.type = c("standardised", "studentised", "log-normal"))
  mack.semiparam <- list()
  for (i in seq_len(nrow(options))) {
    option <- options[i, ]
    name <- paste(ifelse(option$cond, "cond", "uncond"), option$resids.type, sep = "_")

    res <- mackSemiParamBoot(UKMotor, nboot, nsim, cond = option$cond, resids.type = option$resids.type)
    colnames(res$table) <- mack.names

    mack.semiparam[[name]] <- res
    saveRDS(res,
      file.path(res.dir, paste0(paste("mack_semiparam", name, sep = "_"), ".RDS"))
    )
    latex(round(res$table, 2),
      title = "",
      file = file.path(res.dir, paste0(paste("mack_semiparam", name, sep = "_"), ".tex")),
      booktabs = TRUE,
      rowname = NULL,
      table.env = FALSE
    )
  }

  # parametric
  options <- expand.grid(cond = c(TRUE, FALSE), dist = c("normal", "gamma"))
  mack.param <- list()
  for (i in seq_len(nrow(options))) {
    option <- options[i, ]
    name <- paste(ifelse(option$cond, "cond", "uncond"), option$dist, sep = "_")

    res <- mackParamBoot(UKMotor, nboot, nsim, cond = option$cond, dist = option$dist)
    colnames(res$table) <- mack.names

    mack.param[[name]] <- res
    saveRDS(res,
      file.path(res.dir, paste0(paste("mack_param", name, sep = "_"), ".RDS"))
    )
    latex(round(res$table, 2),
      title = "",
      file = file.path(res.dir, paste0(paste("mack_param", name, sep = "_"), ".tex")),
      booktabs = TRUE,
      rowname = NULL,
      table.env = FALSE
    )
  }

  # pairs
  options <- c("parametric", "semiparametric")
  mack.pairs <- list()
  for (option in options) {
    res <- mackPairsBoot(UKMotor, nboot, nsim, sim.type = option)
    colnames(res$table) <- mack.names

    mack.pairs[[option]] <- res
    saveRDS(res,
      file.path(res.dir, paste0(paste("mack_pairs", option, sep = "_"), ".RDS"))
    )
    latex(round(res$table, 2),
      title = "",
      file = file.path(res.dir, paste0(paste("mack_pairs", option, sep = "_"), ".tex")),
      booktabs = TRUE,
      rowname = NULL,
      table.env = FALSE
    )
  }

  ## GLM
  glm.names <- c("$i/j$",
    "$\\widehat{a}_i$",
    "$\\widehat{b}_j$",
    "$\\widehat{R}^B_j$",
    "\\resizebox{4em}{!}{$\\widehat{\\mathrm{MSEP}}(\\widehat{R}^B_j)$}"
  )
  # semiparametric
  options <- c("pearson", "deviance")
  glm.semiparam <- list()
  for (option in options) {
    res <- glmSemiParamBoot(UKMotor.incr, nboot, nsim, resids.type = option)
    colnames(res$table) <- glm.names

    glm.semiparam[[option]] <- res
    saveRDS(res,
      file.path(res.dir, paste0(paste("glm_semiparam", option, sep = "_"), ".RDS"))
    )
    latex(round(res$table, 2),
      title = "",
      file = file.path(res.dir, paste0(paste("glm_semiparam", option, sep = "_"), ".tex")),
      booktabs = TRUE,
      rowname = NULL,
      table.env = FALSE
    )
  }
  res.pears <- glmSemiParamBoot(UKMotor.incr, nboot, nsim, resids.type = "pearson")
  res.dev <- glmSemiParamBoot(UKMotor.incr, nboot, nsim, resids.type = "deviance")
  saveRDS(res.pears, file.path(res.dir, "pois_semiparam_pears.RDS"))
  saveRDS(res.dev, file.path(res.dir, "pois_semiparam_dev.RDS"))

  # parametric
  options <- c("normal", "gamma", "poisson")
  glm.param <- list()
  for (option in options) {
    res <- glmParamBoot(UKMotor.incr, nboot, nsim, dist = option)
    colnames(res$table) <- glm.names

    glm.param[[option]] <- res
    saveRDS(res,
      file.path(res.dir, paste0(paste("glm_param", option, sep = "_"), ".RDS"))
    )
    latex(round(res$table, 2),
      title = "",
      file = file.path(res.dir, paste0(paste("glm_param", option, sep = "_"), ".tex")),
      booktabs = TRUE,
      rowname = NULL,
      table.env = FALSE
    )
  }

  # pairs
  options <- c("rowcol", "corners", "randomised")
  glm.pairs <- list()
  for (option in options) {
    res <- glmPairsBoot(UKMotor.incr, nboot, nsim, keep = option)
    colnames(res$table) <- glm.names

    glm.pairs[[option]] <- res
    saveRDS(res,
      file.path(res.dir, paste0(paste("glm_pairs", option, sep = "_"), ".RDS"))
    )
  }

  pairs.tab <- as.data.frame(lapply(glm.pairs, function(res) { res$table[, 5] }))
  colnames(pairs.tab) <- c("First row and column", "Corner points", "Randomised")
  latex(round(pairs.tab, 2),
    title = "",
    file = file.path(res.dir, "glm_pairs.tex"),
    booktabs = TRUE,
    rowlabel = "Origin",
    rowname = 2:7,
    table.env = FALSE,
    cgroup = "\\resizebox{4em}{!}{$\\bm{\\widehat{\\mathrm{MSEP}}(\\widehat{R}^B_j)$}}",
    col.just = rep("c", ncol(pairs.tab))
  )
  file.conn <- file("results/example/glm_pairs.tex")
  str <- readLines(file.conn)
  str <- sub("cline", "cmidrule", str)
  writeLines(str, file.conn)
  close(file.conn)


} else { # don't recompute
  ## Mack
  # semiparametric
  options <- expand.grid(cond = c(TRUE, FALSE), resids.type = c("standardised", "studentised", "log-normal"))
  mack.semiparam <- list()
  for (i in seq_len(nrow(options))) {
    option <- options[i, ]
    name <- paste(ifelse(option$cond, "cond", "uncond"), option$resids.type, sep = "_")
    mack.semiparam[[name]] <- readRDS(file.path(res.dir,
      paste0(paste("mack", "semiparam", name, sep = "_"), ".RDS")
    ))
  }

  # parametric
  options <- expand.grid(cond = c(TRUE, FALSE), dist = c("normal", "gamma"))
  mack.param <- list()
  for (i in seq_len(nrow(options))) {
    option <- options[i, ]
    name <- paste(ifelse(option$cond, "cond", "uncond"), option$dist, sep = "_")
    mack.param[[name]] <- readRDS(file.path(res.dir,
      paste0(paste("mack", "param", name, sep = "_"), ".RDS")
    ))
  }

  # pairs
  options <- c("parametric", "semiparametric")
  mack.pairs <- list()
  for (option in options) {
    mack.pairs[[option]] <- readRDS(file.path(res.dir, paste0(paste("mack_pairs", option, sep = "_"), ".RDS")))
  }

  ## GLM
  # semiparametric
  options <- c("pearson", "deviance")
  glm.semiparam <- list()
  for (option in options) {
    glm.semiparam[[option]] <- readRDS(file.path(res.dir, paste0(paste("glm_semiparam", option, sep = "_"), ".RDS")))
  }

  # parametric
  options <- c("normal", "gamma", "poisson")
  glm.param <- list()
  for (option in options) {
    glm.param[[option]] <- readRDS(file.path(res.dir, paste0(paste("glm_param", option, sep = "_"), ".RDS")))
  }

  # pairs
  options <- c("rowcol", "corners", "randomised")
  glm.pairs <- list()
  for (option in options) {
    glm.pairs[[option]] <- readRDS(file.path(res.dir, paste0(paste("glm_pairs", option, sep = "_"), ".RDS")))
  }
}

### Plots
width <- 418
height <- 591
mm.per.pt <- 0.3528
width.mm <- mm.per.pt * width
height.mm <- mm.per.pt * height

theme_set(theme_bw())

theme_update(
  legend.position = "top",
  axis.title = element_text(size = 8),
  axis.text = element_text(size = 6),
  legend.text = element_text(size = 8),
  legend.key.size = unit(0.6, "cm"),
  legend.background = element_rect(fill = NA),
  legend.title = element_text(size = 10)
)

## Mack
# semiparametric
for (i in seq_along(mack.semiparam)) {
  res <- mack.semiparam[[i]]

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

  name <- names(mack.semiparam)[i]

  if (grepl("standardised", name)) {
    ggsave(
      file.path(plot.dir, paste0(paste("mack_semiparam", name, sep = "_"), ".eps")),
      p,
      units = "mm",
      height = height.mm / 3,
      width = width.mm / 2.2
    )
  } else {
    ggsave(
      file.path(plot.dir, paste0(paste("mack_semiparam", name, sep = "_"), ".eps")),
      p,
      units = "mm",
      height = width.mm / 2.2, # landscape
      width = height.mm / 2.2
    )
  }
}

# parametric
for (i in seq_along(mack.param)) {
  res <- mack.param[[i]]

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
    xlab("x") +
    ylab("Density")

  name <- names(mack.param)[i]

  ggsave(
    file.path(plot.dir, paste0(paste("mack_param", name, sep = "_"), ".eps")),
    p,
    units = "mm",
    height = width.mm / 2.2, # landscape
    width = height.mm / 2.2
  )
}

# pairs
for (i in seq_along(mack.pairs)) {
  plotdf <- melt(res$dist, value.name = "dist", variable.name = "type")
  p <- ggplot(plotdf) +
    geom_density(
      aes(x = dist, colour = type),
      key_glyph = draw_key_path
    ) +
    scale_color_manual(
      labels = c("Fitted", "Predictive"),
      values = c("red", "blue"),
      name = NULL
    ) +
    xlab("Reserve") +
    ylab("Density")

  name <- names(mack.pairs)[i]

  ggsave(
    file.path(plot.dir, paste0(paste("mack_pairs", name, sep = "_"), ".eps")),
    p,
    units = "mm",
    height = height.mm / 3,
    width = width.mm / 2.2
  )
}

## GLM
# semiparametric
plotdf <- data.frame(reserve = glm.semiparam[["pearson"]]$pred.dist)
p <- ggplot(plotdf) +
  geom_density(
    aes(x = reserve),
    colour = "blue",
    key_glyph = draw_key_path
  ) +
  xlab("Reserve") +
  ylab("Density")

ggsave(
  file.path(plot.dir, "glm_semiparam.eps"),
  p,
  units = "mm",
  height = height.mm / 2.5,
  width = width.mm
)

# parametric
plotdf <- data.frame(
  Normal = glm.param[["normal"]]$pred.dist,
  Gamma = glm.param[["gamma"]]$pred.dist,
  Poisson = glm.param[["poisson"]]$pred.dist
)
plotdf <- melt(plotdf, variable.name = "dist", value.name = "reserve")

p <- ggplot(plotdf) +
  geom_density(
    aes(x = reserve, colour = dist),
    key_glyph = draw_key_path
  ) +
  scale_color_manual(
    values = c("red", "blue", "green"),
    name = NULL
  ) +
  xlab("Reserve") +
  ylab("Density")

ggsave(
  file.path(plot.dir, "glm_param.eps"),
  p,
  units = "mm",
  height = height.mm / 2.5,
  width = width.mm
)