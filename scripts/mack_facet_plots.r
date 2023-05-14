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
outlier.pts <- as.matrix(unique(single.res[opt == "standardised", .(outlier.rowidx, outlier.colidx)]))
outlier.pts <- rbind(outlier.pts, matrix(c(c(1, 7), c(1, 6), c(2, 6)), byrow = TRUE, nrow = 3))

## Semiparametric
# Conditional with standardised residuals
bt <- "residuals"
o <- "standardised"
c <- TRUE
mf <- 0.5
sf <- 0.5

subset.idxs <- outlier.pts[1:4, ]
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
    vars(i = factor(outlier.rowidx), j = factor(outlier.colidx)),
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

subset.idxs <- outlier.pts[5:8, ]
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
    vars(i = factor(outlier.rowidx), j = factor(outlier.colidx)),
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

subset.idxs <- outlier.pts[9:12, ]
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
    vars(i = factor(outlier.rowidx), j = factor(outlier.colidx)),
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

subset.idxs <- outlier.pts[18:21, ]
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
    vars(i = factor(outlier.rowidx), j = factor(outlier.colidx)),
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

subset.idxs <- outlier.pts[13:16, ]
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
    vars(i = factor(outlier.rowidx), j = factor(outlier.colidx)),
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

# Unconditional with gamma distribution
bt <- "parametric"
o <- "gamma"
c <- FALSE
mf <- 0.5
sf <- 1

subset.idxs <- outlier.pts[13:16, ]
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
    vars(i = factor(outlier.rowidx), j = factor(outlier.colidx)),
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

## Pairs
bt <- "pairs"
contaminated <- single.res[
  boot.type == bt &
    (excl.rowidx != outlier.rowidx | excl.colidx != outlier.colidx)
]

uncontaminated <- single.res[
  boot.type == bt &
    excl.rowidx == outlier.rowidx &
    excl.colidx == outlier.colidx
]

contaminated <- contaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]
uncontaminated <- uncontaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]

for (i_page in 1:2) {

  p <- ggplot() +
    geom_density(aes(reserve, group = interaction(excl.colidx, excl.rowidx)), contaminated[cond == TRUE]) +
    geom_density(mapping = aes(reserve), uncontaminated[cond == TRUE], colour = "red") +
    facet_wrap_paginate(
      vars(i = factor(outlier.rowidx), j = factor(outlier.colidx)),
      scales = "free",
      labeller = label_wrap_gen(multi_line = FALSE),
      ncol = 4,
      nrow = 3,
      page = i_page
    ) +
    xlab("Reserve") +
    ylab("Density")

  path <- file.path(plot.dir, paste0(paste(
    "mack_single_densities_cond", bt, i_page, sep = "_"),
  ".eps")
  )

  ggsave(
    path,
    p,
    units = "mm",
    height = width.mm, # landscape
    width = height.mm
  )

  p <- ggplot() +
    geom_density(aes(reserve, group = interaction(excl.colidx, excl.rowidx)), contaminated[cond == FALSE]) +
    geom_density(mapping = aes(reserve), uncontaminated[cond == FALSE], colour = "red") +
    facet_wrap_paginate(
      vars(i = factor(outlier.rowidx), j = factor(outlier.colidx)),
      scales = "free",
      labeller = label_wrap_gen(multi_line = FALSE),
      ncol = 4,
      nrow = 3,
      page = i_page
    ) +
    xlab("Reserve") +
    ylab("Density")

  path <- file.path(plot.dir, paste0(paste(
    "mack_single_densities_uncond", bt, i_page, sep = "_"),
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

subset.idxs <- 1:2
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

subset.idxs <- 3:4
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

subset.idxs <- 5:6
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

subset.idxs <- 3:4
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

subset.idxs <- 1:2
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

subset.idxs <- 3:4
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
contaminated <- calendar.res[
  boot.type == bt &
    excl.diagidx != outlier.diagidx &
    cond == c &
    mean.factor == mf &
    sd.factor == sf
]

uncontaminated <- calendar.res[
  boot.type == bt &
    excl.diagidx == outlier.diagidx &
    mean.factor == mf &
    sd.factor == sf
]

contaminated <- contaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]
uncontaminated <- uncontaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]

p <- ggplot() +
  geom_density(aes(reserve, group = factor(excl.diagidx)), contaminated) +
  geom_density(mapping = aes(reserve), uncontaminated, colour = "red") +
  facet_wrap(
    vars(outlier.diagidx),
    scales = "free",
    labeller = label_wrap_gen(multi_line = FALSE),
    ncol = 3,
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

## Origin outlier:
###############################################################################
origin.res <- readRDS("results/mack_origin.RDS")
mean.factors <- unique(origin.res$mean.factor)
sd.factors <- unique(origin.res$sd.factor)
boot.types <- unique(origin.res$boot.type)

mf <- 0.5
sf <- 0.5

# Semiparametric and parametric
for (bt in c("residuals", "parametric")) {
  opts <- unique(origin.res[boot.type == bt]$opt)
  for (o in opts) {
    contaminated <- origin.res[
      boot.type == bt &
        opt == o &
        excl.rowidx != outlier.rowidx &
        mean.factor == mf &
        sd.factor == sf
    ]

    uncontaminated <- origin.res[
      boot.type == bt &
        opt == o &
        excl.rowidx == outlier.rowidx &
        mean.factor == mf &
        sd.factor == sf
    ]

    contaminated <- contaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]
    uncontaminated <- uncontaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]

    p <- ggplot() +
      geom_density(aes(reserve, group = factor(excl.rowidx)), contaminated[cond == TRUE]) +
      geom_density(mapping = aes(reserve), uncontaminated[cond == TRUE], colour = "red") +
      facet_wrap(
        vars(factor(outlier.rowidx)),
        scales = "free",
        labeller = label_wrap_gen(multi_line = FALSE),
        ncol = 4,
        nrow = 2
      ) +
      xlab("Reserve") +
      ylab("Density")

    path <- file.path(plot.dir, paste0(paste(
      "mack_origin_densities_cond", bt, o, sep = "_"),
    ".eps")
    )

    ggsave(
      path,
      p,
      units = "mm",
      height = width.mm, # landscape
      width = height.mm
    )

    p <- ggplot() +
      geom_density(aes(reserve, group = factor(excl.rowidx)), contaminated[cond == FALSE]) +
      geom_density(mapping = aes(reserve), uncontaminated[cond == FALSE], colour = "red") +
      facet_wrap(
        vars(factor(outlier.rowidx)),
        scales = "free",
        labeller = label_wrap_gen(multi_line = FALSE),
        ncol = 4,
        nrow = 2
      ) +
      xlab("Reserve") +
      ylab("Density")

    path <- file.path(plot.dir, paste0(paste(
      "mack_origin_densities_uncond", bt, o, sep = "_"),
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
}

# Pairs
bt <- "pairs"
contaminated <- origin.res[
  boot.type == bt &
    excl.rowidx != outlier.rowidx &
    mean.factor == mf &
    sd.factor == sf
]

uncontaminated <- origin.res[
  boot.type == bt &
    excl.rowidx == outlier.rowidx &
    mean.factor == mf &
    sd.factor == sf
]

contaminated <- contaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]
uncontaminated <- uncontaminated[reserve < quantile(reserve, 0.975, na.rm = TRUE)]

p <- ggplot() +
  geom_density(aes(reserve, group = factor(excl.rowidx)), contaminated) +
  geom_density(mapping = aes(reserve), uncontaminated, colour = "red") +
  facet_wrap(
    vars(factor(outlier.rowidx)),
    scales = "free",
    labeller = label_wrap_gen(multi_line = FALSE),
    ncol = 4,
    nrow = 2
  ) +
  xlab("Reserve") +
  ylab("Density")

path <- file.path(plot.dir, paste0(paste(
  "mack_origin_densities", bt, sep = "_"),
".eps")
)

ggsave(
  path,
  p,
  units = "mm",
  height = width.mm, # landscape
  width = height.mm
)
