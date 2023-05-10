suppressPackageStartupMessages({
  library(optparse); library(ggplot2); library(latex2exp);
  library(patchwork); library(parallel); library(doSNOW);
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
    default = 10,
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

# in mm
width <- 418
height <- 591

mm.per.pt <- 0.3528
width.mm <- mm.per.pt * width
height.mm <- mm.per.pt * height

theme_set(theme_bw())

plot.dir <- file.path("plots", "example")
res.dir <- file.path("results", "example")
example.dir <- file.path("scripts", "example")

triangle <- cum2incr(UKMotor)

if (!exists("par.init")) {
  par.init <- TRUE
  cores <- detectCores()
  cl <- makeSOCKcluster(round(cores / 2))
  registerDoSNOW(cl)
}


# ## Amended by David Firth, 2003.01.16, at points labelled ###
# ## to cope with negative y values
# ##
# ## Computes Pearson X^2 rather than Poisson deviance
# ##
# ## Starting values are all equal to the global mean
# quasipoisson <- function(link = "log") {
#   linktemp <- substitute(link)
#   if (!is.character(linktemp)) {
#     linktemp <- deparse(linktemp)
#     if (linktemp == "link")
#       linktemp <- eval(link)
#   }
#   if (any(linktemp == c("log", "identity", "sqrt")))
#     stats <- make.link(linktemp)
#   else stop(paste(linktemp, "link not available for poisson",
#     "family; available links are", "\"identity\", \"log\" and
# \"sqrt\""))
#   variance <- function(mu) mu
#   validmu <- function(mu) all(mu > 0)
#   dev.resids <- function(y, mu, wt) wt * (y - mu)^2 / mu ###
#   aic <- function(y, n, mu, wt, dev) NA
#   initialize <- expression({
#     n <- rep(1, nobs)
#     mustart <- rep(mean(y), length(y)) ###
#   })
#   structure(list(family = "quasipoisson", link = linktemp,
#     linkfun = stats$linkfun, linkinv = stats$linkinv, variance =
#       variance,
#     dev.resids = dev.resids, aic = aic, mu.eta = stats$mu.eta,
#     initialize = initialize, validmu = validmu, valideta =
#       stats$valideta),
#   class = "family")
# }