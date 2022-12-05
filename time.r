library(ChainLadder)
library(ggplot2)
source("pattern_break.r")

triangle <- unname(UKMotor)
ndev <- ncol(triangle)
nboot <- 1e3

exclude.resids <- list(c(1, 1), c(1, 2))


res <- reserveBoot(gamma.triangle, 1e3, distribution = "gamma", bootstrap.type = "conditional", resids.type = "scaled")

times <- rep(0, 10)
times.fortran <- rep(0, 10)

for (i in seq_len(length(times))) {

start <- Sys.time()
reserve <- reserveBoot(triangle, 1e3, bootstrap.type = "unconditional")
end <- Sys.time()

times[i] <- end - start

start <- Sys.time()
reserveFortran <- reserveBootFortran(triangle, 1e3, distribution = "gamma")
end <- Sys.time()

times.fortran[i] <- end - start
}

result <- mean(times)
result.fortran <- mean(times.fortran)