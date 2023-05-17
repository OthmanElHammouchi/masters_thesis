library(claimsBoot)
suppressPackageStartupMessages(library(ChainLadder))
library(ggplot2)

calendar.res <- mackSim(UKMotor, "calendar", 1e3, c(2), c(1), "pairs")

single.res <- readRDS("results/mack_single.RDS")
calendar.res <- readRDS("results/mack_calendar.RDS")
origin.res <- readRDS("results/mack_origin.RDS")

contaminated <- single.res[
  (outlier.rowidx != excl.rowidx | outlier.colidx != excl.colidx) &
    boot.type == "pairs"
]
clean <- single.res[
  (outlier.rowidx == excl.rowidx & outlier.colidx == excl.colidx) &
    boot.type == "pairs"
]

ggplot() +
  geom_density(aes(reserve, group = interaction(excl.colidx, excl.rowidx)), contaminated) +
  geom_density(aes(reserve), clean, colour = "red") +
  facet_wrap(
    vars(factor(outlier.rowidx), factor(outlier.colidx)),
    scales = "free",
    labeller = label_wrap_gen(multi_line = FALSE)
  )

#######################################################



library(formattable)
min_formatter <- formatter("span", style = x ~ style(color = ifelse(x == min(x), "red", "black")))
formattable(mean.table, lapply(seq_len(nrow(mean.table)), function(row) {
  area(row) ~ min_formatter
}))

########################################################


quasipoisson <- function (link = "log")
## Amended by David Firth, 2003.01.16, at points labelled ###
## to cope with negative y values
##
## Computes Pearson X^2 rather than Poisson deviance
##
## Starting values are all equal to the global mean
{
     linktemp <- substitute(link)
     if (!is.character(linktemp)) {
         linktemp <- deparse(linktemp)
         if (linktemp == "link")
             linktemp <- eval(link)
     }
     if (any(linktemp == c("log", "identity", "sqrt")))
         stats <- make.link(linktemp)
     else stop(paste(linktemp, "link not available for poisson",
         "family; available links are", "\"identity\", \"log\" and 
\"sqrt\""))
     variance <- function(mu) mu
     validmu <- function(mu) all(mu > 0)
     dev.resids <- function(y, mu, wt) wt*(y-mu)^2/mu   ###
     aic <- function(y, n, mu, wt, dev) NA
     initialize <- expression({
         n <- rep(1, nobs)
         mustart <- rep(mean(y), length(y))             ###
     })
     structure(list(family = "quasipoisson", link = linktemp,
         linkfun = stats$linkfun, linkinv = stats$linkinv, variance = 
variance,
         dev.resids = dev.resids, aic = aic, mu.eta = stats$mu.eta,
         initialize = initialize, validmu = validmu, valideta = 
stats$valideta),
         class = "family")
}

suppressPackageStartupMessages(library(ChainLadder))
triangle <- cum2incr(UKMotor)
triangle[1, 1] <- -triangle[1, 1]

long <- as.data.frame(triangle)
long <- transform(long, origin = factor(origin, levels = dimnames(triangle)[[1]]))
long <- long[order(long$origin, long$dev), ]

long.lower <- long[is.na(long$value), ]
long.upper <- long[!is.na(long$value), ]

model <- glm(value ~ factor(origin) + factor(dev), quasipoisson(), long.upper, maxit = 1e3)
long.lower$value <- predict(model, long.lower, type = "response")
triangle.proj <- incr2cum(as.triangle(rbind(long.upper, long.lower)))

#############################################################

suppressPackageStartupMessages(library(ChainLadder))
library(knitr)

long <- as.data.frame(triangle)
long[, c("origin", "dev")] <- lapply(long[, c("origin", "dev")], as.factor)
long.upper <- long[!is.na(long$value), ]
long.lower <- long[is.na(long$value), ]
model <- glm(value ~ origin + dev, quasipoisson(), long.upper)
long.lower$value <- predict(model, long.lower, type = "response")

lower <- as.triangle(as.data.frame(long.lower))
pois.res <- rowSums(lower, na.rm = TRUE)
suppressWarnings({
  mack.res <- summary(MackChainLadder(UKMotor))$ByOrigin$IBNR[-1]
})
res.df <- data.frame(Origin = rownames(lower), Mack = mack.res, Pois = pois.res)
kable(res.df)

#############################################################

# What about the penultimate column?

# In the later columns it's very possible to draw the same pair
# for all resamples and then you have no variability.
