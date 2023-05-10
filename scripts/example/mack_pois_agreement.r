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