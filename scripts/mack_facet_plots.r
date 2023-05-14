library(ggplot2)
library(ggforce)
library(claimsBoot)
library(data.table)
suppressPackageStartupMessages(library(ChainLadder))

triangle <- UKMotor
ndev <- ncol(triangle)

# in mm
width <- 410
height <- 630

mm.per.pt <- 0.3528
width.mm <- mm.per.pt * width
height.mm <- mm.per.pt * height

plot.dir <- "plots"

theme_set(theme_bw())

theme_update(
  legend.position = "top",
  axis.title = element_text(size = 8),
  axis.text = element_text(size = 6),
  legend.text = element_text(size = 8),
  legend.key.size = unit(0.6, "cm"),
  strip.text.x = element_text(size = 6)
)

### Single outlier:
###############################################################################
single.res <- readRDS("results/mack_single.RDS")

## Semiparametric
# Conditional with standardised residuals
bt <- "residuals"
o <- "standardised"
c <- TRUE
mf <- 0.5
sf <- 0.5

subset.idxs <- data.frame(outlier.rowidx = c(1, 2, 4, 5), outlier.colidx = c(2, 3, 2, 3))
subset <- do.call(rbind, lapply(seq_len(nrow(subset.idxs)), function(i) {
  single.res[outlier.rowidx == subset.idxs[i, 1] & outlier.colidx == subset.idxs[i, 2]]
}))
contaminated <- subset[
  boot.type == bt &
    opt == o &
    (excl.rowidx != outlier.rowidx | excl.colidx != outlier.colidx) &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]
cleaned <- subset[
  boot.type == bt &
    opt == o &
    excl.rowidx == outlier.rowidx &
    excl.colidx == outlier.colidx &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]

p <- ggplot() +
  geom_density(aes(reserve, group = interaction(excl.colidx, excl.rowidx)), contaminated) +
  geom_density(mapping = aes(reserve), cleaned, colour = "red") +
  facet_wrap(
    vars(factor(outlier.rowidx), factor(outlier.colidx)),
    scales = "free",
    labeller = label_wrap_gen(multi_line = FALSE),
    ncol = 2,
    nrow = 2
  ) +
  xlab("Reserve") +
  ylab("Density")

path <- file.path(plot.dir, paste0(paste(
  "mack_single_densities_cond", bt, o, sep = "_"),
".eps")
)

ggsave(
  path,
  p,
  units = "mm",
  height = width.mm, # landscape
  width = height.mm
)

## Semiparametric
# Unconditional with standardised residuals
bt <- "residuals"
o <- "standardised"
c <- FALSE
mf <- 2
sf <- 0.5

subset.idxs <- data.frame(outlier.rowidx = c(1, 2, 4, 5), outlier.colidx = c(2, 3, 2, 3))
subset <- do.call(rbind, lapply(seq_len(nrow(subset.idxs)), function(i) {
  single.res[outlier.rowidx == subset.idxs[i, 1] & outlier.colidx == subset.idxs[i, 2]]
}))
contaminated <- subset[
  boot.type == bt &
    opt == o &
    (excl.rowidx != outlier.rowidx | excl.colidx != outlier.colidx) &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]
cleaned <- subset[
  boot.type == bt &
    opt == o &
    excl.rowidx == outlier.rowidx &
    excl.colidx == outlier.colidx &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]

p <- ggplot() +
  geom_density(aes(reserve, group = interaction(excl.colidx, excl.rowidx)), contaminated) +
  geom_density(mapping = aes(reserve), cleaned, colour = "red") +
  facet_wrap(
    vars(factor(outlier.rowidx), factor(outlier.colidx)),
    scales = "free",
    labeller = label_wrap_gen(multi_line = FALSE),
    ncol = 2,
    nrow = 2
  ) +
  xlab("Reserve") +
  ylab("Density")

path <- file.path(plot.dir, paste0(paste(
  "mack_single_densities_uncond", bt, o, sep = "_"),
".eps")
)

ggsave(
  path,
  p,
  units = "mm",
  height = width.mm, # landscape
  width = height.mm
)

# Conditional with log-normal residuals
bt <- "residuals"
o <- "log-normal"
c <- TRUE
mf <- 0.5
sf <- 2

subset.idxs <- data.frame(outlier.rowidx = c(1, 2, 4, 5), outlier.colidx = c(2, 3, 2, 3))
subset <- do.call(rbind, lapply(seq_len(nrow(subset.idxs)), function(i) {
  single.res[outlier.rowidx == subset.idxs[i, 1] & outlier.colidx == subset.idxs[i, 2]]
}))
contaminated <- subset[
  boot.type == bt &
    opt == o &
    (excl.rowidx != outlier.rowidx | excl.colidx != outlier.colidx) &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]
cleaned <- subset[
  boot.type == bt &
    opt == o &
    excl.rowidx == outlier.rowidx &
    excl.colidx == outlier.colidx &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]

p <- ggplot() +
  geom_density(aes(reserve, group = interaction(excl.colidx, excl.rowidx)), contaminated) +
  geom_density(mapping = aes(reserve), cleaned, colour = "red") +
  facet_wrap(
    vars(factor(outlier.rowidx), factor(outlier.colidx)),
    scales = "free",
    labeller = label_wrap_gen(multi_line = FALSE),
    ncol = 2,
    nrow = 2
  ) +
  xlab("Reserve") +
  ylab("Density")

path <- file.path(plot.dir, paste0(paste(
  "mack_single_densities_cond", bt, o, sep = "_"),
".eps")
)

ggsave(
  path,
  p,
  units = "mm",
  height = width.mm, # landscape
  width = height.mm
)

# Unconditional with log-normal residuals
bt <- "residuals"
o <- "log-normal"
c <- FALSE
mf <- 2
sf <- 2

subset.idxs <- data.frame(outlier.rowidx = c(1, 2, 4, 5), outlier.colidx = c(2, 3, 2, 3))
subset <- do.call(rbind, lapply(seq_len(nrow(subset.idxs)), function(i) {
  single.res[outlier.rowidx == subset.idxs[i, 1] & outlier.colidx == subset.idxs[i, 2]]
}))
contaminated <- subset[
  boot.type == bt &
    opt == o &
    (excl.rowidx != outlier.rowidx | excl.colidx != outlier.colidx) &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]
cleaned <- subset[
  boot.type == bt &
    opt == o &
    excl.rowidx == outlier.rowidx &
    excl.colidx == outlier.colidx &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]

p <- ggplot() +
  geom_density(aes(reserve, group = interaction(excl.colidx, excl.rowidx)), contaminated) +
  geom_density(mapping = aes(reserve), cleaned, colour = "red") +
  facet_wrap(
    vars(factor(outlier.rowidx), factor(outlier.colidx)),
    scales = "free",
    labeller = label_wrap_gen(multi_line = FALSE),
    ncol = 2,
    nrow = 2
  ) +
  xlab("Reserve") +
  ylab("Density")

path <- file.path(plot.dir, paste0(paste(
  "mack_single_densities_uncond", bt, o, sep = "_"),
".eps")
)

ggsave(
  path,
  p,
  units = "mm",
  height = width.mm, # landscape
  width = height.mm
)

## Parametric
# Conditional with normal distribution
bt <- "parametric"
o <- "normal"
c <- TRUE
mf <- 2
sf <- 1

subset.idxs <- data.frame(outlier.rowidx = 4, outlier.colidx = 2)
subset <- do.call(rbind, lapply(seq_len(nrow(subset.idxs)), function(i) {
  single.res[outlier.rowidx == subset.idxs[i, 1] & outlier.colidx == subset.idxs[i, 2]]
}))
contaminated <- subset[
  boot.type == bt &
    opt == o &
    (excl.rowidx != outlier.rowidx | excl.colidx != outlier.colidx) &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]
cleaned <- subset[
  boot.type == bt &
    opt == o &
    excl.rowidx == outlier.rowidx &
    excl.colidx == outlier.colidx &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]

p <- ggplot() +
  geom_density(aes(reserve, group = interaction(excl.colidx, excl.rowidx)), contaminated) +
  geom_density(mapping = aes(reserve), cleaned, colour = "red") +
  xlab("Reserve") +
  ylab("Density")

path <- file.path(plot.dir, paste0(paste(
  "mack_single_densities_cond", bt, o, sep = "_"),
".eps")
)

ggsave(
  path,
  p,
  units = "mm",
  height = width.mm / 2.5, # landscape
  width = height.mm / 2
)

# Unconditional with normal distribution
bt <- "parametric"
o <- "normal"
c <- FALSE
mf <- 0.5
sf <- 1

subset.idxs <- subset.idxs <- data.frame(outlier.rowidx = 4, outlier.colidx = 2)
subset <- do.call(rbind, lapply(seq_len(nrow(subset.idxs)), function(i) {
  single.res[outlier.rowidx == subset.idxs[i, 1] & outlier.colidx == subset.idxs[i, 2]]
}))
contaminated <- subset[
  boot.type == bt &
    opt == o &
    (excl.rowidx != outlier.rowidx | excl.colidx != outlier.colidx) &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]
cleaned <- subset[
  boot.type == bt &
    opt == o &
    excl.rowidx == outlier.rowidx &
    excl.colidx == outlier.colidx &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]

p <- ggplot() +
  geom_density(aes(reserve, group = interaction(excl.colidx, excl.rowidx)), contaminated) +
  geom_density(mapping = aes(reserve), cleaned, colour = "red") +
  xlab("Reserve") +
  ylab("Density")

path <- file.path(plot.dir, paste0(paste(
  "mack_single_densities_uncond", bt, o, sep = "_"),
".eps")
)

ggsave(
  path,
  p,
  units = "mm",
  height = width.mm / 2.5, # landscape
  width = height.mm / 2
)

# Conditional with gamma distribution
bt <- "parametric"
o <- "gamma"
c <- TRUE
mf <- 2
sf <- 1

subset.idxs <- data.frame(outlier.rowidx = 4, outlier.colidx = 2)
subset <- do.call(rbind, lapply(seq_len(nrow(subset.idxs)), function(i) {
  single.res[outlier.rowidx == subset.idxs[i, 1] & outlier.colidx == subset.idxs[i, 2]]
}))
contaminated <- subset[
  boot.type == bt &
    opt == o &
    (excl.rowidx != outlier.rowidx | excl.colidx != outlier.colidx) &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]
cleaned <- subset[
  boot.type == bt &
    opt == o &
    excl.rowidx == outlier.rowidx &
    excl.colidx == outlier.colidx &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]

p <- ggplot() +
  geom_density(aes(reserve, group = interaction(excl.colidx, excl.rowidx)), contaminated) +
  geom_density(mapping = aes(reserve), cleaned, colour = "red") +
  xlab("Reserve") +
  ylab("Density")

path <- file.path(plot.dir, paste0(paste(
  "mack_single_densities_cond", bt, o, sep = "_"),
".eps")
)

ggsave(
  path,
  p,
  units = "mm",
  height = width.mm / 2.5, # landscape
  width = height.mm / 2
)

# Unconditional with gamma distribution
bt <- "parametric"
o <- "gamma"
c <- FALSE
mf <- 2
sf <- 1

subset.idxs <- data.frame(outlier.rowidx = 4, outlier.colidx = 2)
subset <- do.call(rbind, lapply(seq_len(nrow(subset.idxs)), function(i) {
  single.res[outlier.rowidx == subset.idxs[i, 1] & outlier.colidx == subset.idxs[i, 2]]
}))
contaminated <- subset[
  boot.type == bt &
    opt == o &
    (excl.rowidx != outlier.rowidx | excl.colidx != outlier.colidx) &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]
cleaned <- subset[
  boot.type == bt &
    opt == o &
    excl.rowidx == outlier.rowidx &
    excl.colidx == outlier.colidx &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]

p <- ggplot() +
  geom_density(aes(reserve, group = interaction(excl.colidx, excl.rowidx)), contaminated) +
  geom_density(mapping = aes(reserve), cleaned, colour = "red") +
  xlab("Reserve") +
  ylab("Density")

path <- file.path(plot.dir, paste0(paste(
  "mack_single_densities_uncond", bt, o, sep = "_"),
".eps")
)

ggsave(
  path,
  p,
  units = "mm",
  height = width.mm / 2.5, # landscape
  width = height.mm / 2
)

## Pairs
bt <- "pairs"
mf <- 2
sf <- 1

subset.idxs <- data.frame(outlier.rowidx = c(1, 2, 4, 5), outlier.colidx = c(2, 3, 2, 3))
subset <- do.call(rbind, lapply(seq_len(nrow(subset.idxs)), function(i) {
  single.res[outlier.rowidx == subset.idxs[i, 1] & outlier.colidx == subset.idxs[i, 2]]
}))

contaminated <- subset[
  boot.type == bt &
    (excl.rowidx != outlier.rowidx | excl.colidx != outlier.colidx) &
    mean.factor == mf &
    sd.factor == sf
]
cleaned <- subset[
  boot.type == bt &
    excl.rowidx == outlier.rowidx &
    excl.colidx == outlier.colidx &
    mean.factor == mf &
    sd.factor == sf
]

p <- ggplot() +
  geom_density(aes(reserve, group = interaction(excl.colidx, excl.rowidx)), contaminated) +
  geom_density(mapping = aes(reserve), cleaned, colour = "red") +
  facet_wrap(
    vars(factor(outlier.rowidx), factor(outlier.colidx)),
    scales = "free",
    labeller = label_wrap_gen(multi_line = FALSE),
    ncol = 2,
    nrow = 2
  ) +
  xlab("Reserve") +
  ylab("Density")

path <- file.path(plot.dir, paste0(paste(
  "mack_single_densities", bt, sep = "_"),
".eps")
)

ggsave(
  path,
  p,
  units = "mm",
  height = width.mm, # landscape
  width = height.mm
)

### Calendar outlier:
###############################################################################
calendar.res <- readRDS("results/mack_calendar.RDS")

## Semiparametric
# Conditional with standardised residuals
bt <- "residuals"
o <- "standardised"
c <- TRUE
mf <- 0.5
sf <- 0.5

subset.idxs <- c(1, 6)
subset <- do.call(rbind, lapply(seq_len(length(subset.idxs)), function(i) {
  calendar.res[outlier.diagidx == subset.idxs[i]]
}))

contaminated <- subset[
  boot.type == bt &
    opt == o &
    excl.diagidx != outlier.diagidx &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]

cleaned <- subset[
  boot.type == bt &
    opt == o &
    excl.diagidx == outlier.diagidx &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]

p <- ggplot() +
  geom_density(aes(reserve, group = factor(excl.diagidx)), contaminated) +
  geom_density(mapping = aes(reserve), cleaned, colour = "red") +
  facet_wrap(
    vars(factor(outlier.diagidx)),
    scales = "free",
    labeller = label_wrap_gen(multi_line = FALSE),
    ncol = 1,
    nrow = 2
  ) +
  xlab("Reserve") +
  ylab("Density")

path <- file.path(plot.dir, paste0(paste(
  "mack_calendar_densities_cond", bt, o, sep = "_"),
".eps")
)

ggsave(
  path,
  p,
  units = "mm",
  height = width.mm,
  width = height.mm / 2
)

# Unconditional with standardised residuals
bt <- "residuals"
o <- "standardised"
c <- TRUE
mf <- 2
sf <- 0.5

subset.idxs <- c(1, 6)
subset <- do.call(rbind, lapply(seq_len(length(subset.idxs)), function(i) {
  calendar.res[outlier.diagidx == subset.idxs[i]]
}))

contaminated <- subset[
  boot.type == bt &
    opt == o &
    excl.diagidx != outlier.diagidx &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]

cleaned <- subset[
  boot.type == bt &
    opt == o &
    excl.diagidx == outlier.diagidx &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]

p <- ggplot() +
  geom_density(aes(reserve, group = factor(excl.diagidx)), contaminated) +
  geom_density(mapping = aes(reserve), cleaned, colour = "red") +
  facet_wrap(
    vars(factor(outlier.diagidx)),
    scales = "free",
    labeller = label_wrap_gen(multi_line = FALSE),
    ncol = 1,
    nrow = 2
  ) +
  xlab("Reserve") +
  ylab("Density")

path <- file.path(plot.dir, paste0(paste(
  "mack_calendar_densities_uncond", bt, o, sep = "_"),
".eps")
)

ggsave(
  path,
  p,
  units = "mm",
  height = width.mm,
  width = height.mm / 2
)

# Conditional with log-normal residuals
bt <- "residuals"
o <- "log-normal"
c <- TRUE
mf <- 0.5
sf <- 2

subset.idxs <- c(1, 6)
subset <- do.call(rbind, lapply(seq_len(length(subset.idxs)), function(i) {
  calendar.res[outlier.diagidx == subset.idxs[i]]
}))

contaminated <- subset[
  boot.type == bt &
    opt == o &
    excl.diagidx != outlier.diagidx &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]

cleaned <- subset[
  boot.type == bt &
    opt == o &
    excl.diagidx == outlier.diagidx &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]

p <- ggplot() +
  geom_density(aes(reserve, group = factor(excl.diagidx)), contaminated) +
  geom_density(mapping = aes(reserve), cleaned, colour = "red") +
  facet_wrap(
    vars(factor(outlier.diagidx)),
    scales = "free",
    labeller = label_wrap_gen(multi_line = FALSE),
    ncol = 1,
    nrow = 2
  ) +
  xlab("Reserve") +
  ylab("Density")

path <- file.path(plot.dir, paste0(paste(
  "mack_calendar_densities_cond", bt, o, sep = "_"),
".eps")
)

ggsave(
  path,
  p,
  units = "mm",
  height = width.mm,
  width = height.mm / 2
)

# Unconditional with log-normal residuals
bt <- "residuals"
o <- "log-normal"
c <- FALSE
mf <- 2
sf <- 2

subset.idxs <- c(1, 6)
subset <- do.call(rbind, lapply(seq_len(length(subset.idxs)), function(i) {
  calendar.res[outlier.diagidx == subset.idxs[i]]
}))

contaminated <- subset[
  boot.type == bt &
    opt == o &
    excl.diagidx != outlier.diagidx &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]

cleaned <- subset[
  boot.type == bt &
    opt == o &
    excl.diagidx == outlier.diagidx &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]

p <- ggplot() +
  geom_density(aes(reserve, group = factor(excl.diagidx)), contaminated) +
  geom_density(mapping = aes(reserve), cleaned, colour = "red") +
  facet_wrap(
    vars(factor(outlier.diagidx)),
    scales = "free",
    labeller = label_wrap_gen(multi_line = FALSE),
    ncol = 1,
    nrow = 2
  ) +
  xlab("Reserve") +
  ylab("Density")

path <- file.path(plot.dir, paste0(paste(
  "mack_calendar_densities_uncond", bt, o, sep = "_"),
".eps")
)

ggsave(
  path,
  p,
  units = "mm",
  height = width.mm,
  width = height.mm / 2
)

## Parametric
# Condition with normal distribution
bt <- "parametric"
o <- "normal"
c <- TRUE
mf <- 0.5
sf <- 1

subset.idxs <- c(1, 6)
subset <- do.call(rbind, lapply(seq_len(length(subset.idxs)), function(i) {
  calendar.res[outlier.diagidx == subset.idxs[i]]
}))

contaminated <- subset[
  boot.type == bt &
    opt == o &
    excl.diagidx != outlier.diagidx &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]

cleaned <- subset[
  boot.type == bt &
    opt == o &
    excl.diagidx == outlier.diagidx &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]

p <- ggplot() +
  geom_density(aes(reserve, group = factor(excl.diagidx)), contaminated) +
  geom_density(mapping = aes(reserve), cleaned, colour = "red") +
  facet_wrap(
    vars(factor(outlier.diagidx)),
    scales = "free",
    labeller = label_wrap_gen(multi_line = FALSE),
    ncol = 1,
    nrow = 2
  ) +
  xlab("Reserve") +
  ylab("Density")

path <- file.path(plot.dir, paste0(paste(
  "mack_calendar_densities_cond", bt, o, sep = "_"),
".eps")
)

ggsave(
  path,
  p,
  units = "mm",
  height = width.mm,
  width = height.mm / 2
)

# Unconditional with gamma distribution
bt <- "parametric"
o <- "gamma"
c <- FALSE
mf <- 2
sf <- 1

subset.idxs <- c(1, 6)
subset <- do.call(rbind, lapply(seq_len(length(subset.idxs)), function(i) {
  calendar.res[outlier.diagidx == subset.idxs[i]]
}))

contaminated <- subset[
  boot.type == bt &
    opt == o &
    excl.diagidx != outlier.diagidx &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]

cleaned <- subset[
  boot.type == bt &
    opt == o &
    excl.diagidx == outlier.diagidx &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]

p <- ggplot() +
  geom_density(aes(reserve, group = factor(excl.diagidx)), contaminated) +
  geom_density(mapping = aes(reserve), cleaned, colour = "red") +
  facet_wrap(
    vars(factor(outlier.diagidx)),
    scales = "free",
    labeller = label_wrap_gen(multi_line = FALSE),
    ncol = 1,
    nrow = 2
  ) +
  xlab("Reserve") +
  ylab("Density")

path <- file.path(plot.dir, paste0(paste(
  "mack_calendar_densities_uncond", bt, o, sep = "_"),
".eps")
)

ggsave(
  path,
  p,
  units = "mm",
  height = width.mm,
  width = height.mm / 2
)

## Pairs
bt <- "pairs"
mf <- 2
sf <- 1

subset.idxs <- c(1, 6)
subset <- do.call(rbind, lapply(seq_len(length(subset.idxs)), function(i) {
  calendar.res[outlier.diagidx == subset.idxs[i]]
}))

contaminated <- subset[
  boot.type == bt &
    excl.diagidx != outlier.diagidx &
    mean.factor == mf &
    sd.factor == sf
]

uncontaminated <- subset[
  boot.type == bt &
    excl.diagidx == outlier.diagidx &
    mean.factor == mf &
    sd.factor == sf
]

p <- ggplot() +
  geom_density(aes(reserve, group = factor(excl.diagidx)), contaminated) +
  geom_density(mapping = aes(reserve), uncontaminated, colour = "red") +
  facet_wrap(
    vars(outlier.diagidx),
    scales = "free",
    labeller = label_wrap_gen(multi_line = FALSE),
    ncol = 2,
    nrow = 2
  ) +
  xlab("Reserve") +
  ylab("Density")

path <- file.path(plot.dir, paste0(paste(
  "mack_calendar_densities", bt, sep = "_"),
".eps")
)

ggsave(
  path,
  p,
  units = "mm",
  height = width.mm, # landscape
  width = height.mm
)