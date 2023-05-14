library(ggplot2)
library(ggforce)
library(claimsBoot)
library(data.table)
suppressPackageStartupMessages(library(ChainLadder))

triangle <- UKMotor
ndev <- ncol(triangle)

plot.dir <- "plots"

theme_set(theme_bw())

# in mm
width <- 410
height <- 630

mm.per.pt <- 0.3528
width.mm <- mm.per.pt * width
height.mm <- mm.per.pt * height

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
single.res <- as.data.table(readRDS("results/glm_single.RDS"))

## Parametric
# normal
bt <- "parametric"
o <- "normal"
f <- 2

subset.idxs <- data.frame(outlier.rowidx = c(1, 2, 4, 5), outlier.colidx = c(2, 3, 2, 3))
subset <- do.call(rbind, lapply(seq_len(nrow(subset.idxs)), function(i) {
  single.res[outlier.rowidx == subset.idxs[i, 1] & outlier.colidx == subset.idxs[i, 2]]
}))

contaminated <- subset[
  boot.type == bt &
    opt == o &
    (excl.rowidx != outlier.rowidx | excl.colidx != outlier.colidx) &
    factor == f
]

cleaned <- subset[
  boot.type == bt &
    opt == o &
    excl.rowidx == outlier.rowidx &
    excl.colidx == outlier.colidx &
    factor == f
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
  "glm_single_densities", bt, o, sep = "_"),
".eps")
)

ggsave(
  path,
  p,
  units = "mm",
  height = width.mm, # landscape
  width = height.mm
)

# gamma
bt <- "parametric"
o <- "gamma"
f <- 2

subset.idxs <- data.frame(outlier.rowidx = c(1, 2, 4, 5), outlier.colidx = c(2, 3, 2, 3))
subset <- do.call(rbind, lapply(seq_len(nrow(subset.idxs)), function(i) {
  single.res[outlier.rowidx == subset.idxs[i, 1] & outlier.colidx == subset.idxs[i, 2]]
}))

contaminated <- subset[
  boot.type == bt &
    opt == o &
    (excl.rowidx != outlier.rowidx | excl.colidx != outlier.colidx) &
    factor == f
]

cleaned <- subset[
  boot.type == bt &
    opt == o &
    excl.rowidx == outlier.rowidx &
    excl.colidx == outlier.colidx &
    factor == f
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
  "glm_single_densities", bt, o, sep = "_"),
".eps")
)

ggsave(
  path,
  p,
  units = "mm",
  height = width.mm, # landscape
  width = height.mm
)

# poisson
bt <- "parametric"
o <- "poisson"
f <- 2

subset.idxs <- data.frame(outlier.rowidx = c(1, 2, 4, 5), outlier.colidx = c(2, 3, 2, 3))
subset <- do.call(rbind, lapply(seq_len(nrow(subset.idxs)), function(i) {
  single.res[outlier.rowidx == subset.idxs[i, 1] & outlier.colidx == subset.idxs[i, 2]]
}))

contaminated <- subset[
  boot.type == bt &
    opt == o &
    (excl.rowidx != outlier.rowidx | excl.colidx != outlier.colidx) &
    factor == f
]

cleaned <- subset[
  boot.type == bt &
    opt == o &
    excl.rowidx == outlier.rowidx &
    excl.colidx == outlier.colidx &
    factor == f
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
  "glm_single_densities", bt, o, sep = "_"),
".eps")
)

ggsave(
  path,
  p,
  units = "mm",
  height = width.mm, # landscape
  width = height.mm
)

## Residuals
bt <- "residuals"
f <- 2

subset.idxs <- data.frame(outlier.rowidx = c(1, 2, 4, 5), outlier.colidx = c(2, 3, 2, 3))
subset <- do.call(rbind, lapply(seq_len(nrow(subset.idxs)), function(i) {
  single.res[outlier.rowidx == subset.idxs[i, 1] & outlier.colidx == subset.idxs[i, 2]]
}))

contaminated <- subset[
  boot.type == bt &
    (excl.rowidx != outlier.rowidx | excl.colidx != outlier.colidx) &
    factor == f
]
cleaned <- subset[
  boot.type == bt &
    excl.rowidx == outlier.rowidx &
    excl.colidx == outlier.colidx &
    factor == f
]

p <- ggplot() +
  geom_density(aes(reserve, group = interaction(excl.colidx, excl.rowidx)), contaminated) +
  geom_density(mapping = aes(reserve), cleaned, colour = "red") +
  facet_wrap_paginate(
    vars(factor(outlier.rowidx), factor(outlier.colidx)),
    scales = "free",
    labeller = label_wrap_gen(multi_line = FALSE),
    ncol = 2,
    nrow = 2
  ) +
  xlab("Reserve") +
  ylab("Density")

path <- file.path(plot.dir, paste0(paste(
  "glm_single_densities", bt, sep = "_"),
".eps")
)

ggsave(
  path,
  p,
  units = "mm",
  height = width.mm, # landscape
  width = height.mm
)

## Calendar outlier:
###############################################################################
calendar.res <- as.data.table(readRDS("results/glm_calendar.RDS"))
f <- 2

bt <- "parametric"
opts <- unique(calendar.res[boot.type == bt]$opt)
for (o in opts[!is.na(opts)]) {
  contaminated <- calendar.res[
    boot.type == bt &
      opt == o &
      excl.diagidx != outlier.diagidx &
      factor == f
  ]

  uncontaminated <- calendar.res[
    boot.type == bt &
      opt == o &
      excl.diagidx == outlier.diagidx &
      factor == f
  ]

  contaminated <- contaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]
  uncontaminated <- uncontaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]

  p <- ggplot() +
    geom_density(aes(reserve, group = factor(excl.diagidx)), contaminated) +
    geom_density(mapping = aes(reserve), uncontaminated, colour = "red") +
    facet_wrap(
      vars(factor(outlier.diagidx)),
      scales = "free",
      labeller = label_wrap_gen(multi_line = FALSE),
      ncol = 4,
      nrow = 2
    ) +
    xlab("Reserve") +
    ylab("Density")

  path <- file.path(plot.dir, paste0(paste(
    "glm_calendar_densities", bt, o, sep = "_"),
  ".eps")
  )

  ggsave(
    path,
    p,
    units = "mm",
    height = width.mm, # landscape
    width = height.mm
  )
}


bt <- "residuals"
contaminated <- calendar.res[
  boot.type == bt &
    excl.diagidx != outlier.diagidx &
    factor == f
]

uncontaminated <- calendar.res[
  boot.type == bt &
    excl.diagidx == outlier.diagidx &
    factor == f
]

contaminated <- contaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]
uncontaminated <- uncontaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]

p <- ggplot() +
  geom_density(aes(reserve, group = factor(excl.diagidx)), contaminated) +
  geom_density(mapping = aes(reserve), uncontaminated, colour = "red") +
  facet_wrap(
    vars(factor(outlier.diagidx)),
    scales = "free",
    labeller = label_wrap_gen(multi_line = FALSE),
    ncol = 4,
    nrow = 2
  ) +
  xlab("Reserve") +
  ylab("Density")

path <- file.path(plot.dir, paste0(paste(
  "glm_calendar_densities", bt, sep = "_"),
".eps")
)

ggsave(
  path,
  p,
  units = "mm",
  height = width.mm, # landscape
  width = height.mm
)