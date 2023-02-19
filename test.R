library(patternBreak)
suppressPackageStartupMessages(library(ChainLadder))

triangle <- UKMotor
# datasets <- data(package="ChainLadder")$results[, "Item"]

# for (dataset in datasets) {

#     triangle <- eval(parse(text = dataset))

#     if (!is.list(triangle)) {

#     triangle[is.na(triangle)] <- 0
#     ndev <- nrow(triangle)
#     nboot <- 1e3

#     reserve <- glmBoot(triangle, nboot)

#     }
# }

triangle <- cum2incr(triangle)
triangle[is.na(triangle)] <- 0
ndev <- nrow(triangle)
nboot <- 10

reserve <- glmBoot(triangle, nboot)
# reserve <- mackBoot(triangle, nboot, 1, 1, 1)
