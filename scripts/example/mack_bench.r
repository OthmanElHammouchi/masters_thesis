mackFit <- function(triangle) {
  ndev <- ncol(triangle)

  dev.facs <- rep(0, ndev - 1)
  sigmas <- rep(0, ndev - 1)

  for (j in 1:(ndev - 1)) {
    nrows <- ndev - j

    model <- lm(
      y ~ x - 1,
      weights = 1 / triangle[1:nrows, j],
      data = data.frame(
        x = triangle[1:nrows, j],
        y = triangle[1:nrows, j + 1]
      )
    )

    dev.facs[j] <- unname(model$coefficients)
    if (j < ndev - 1) {
      sigmas[j] <- summary(model)$sigma
    } else {
      sigmas[j] <- sqrt(min(
        sigmas[(j - 1)]**2,
        sigmas[j - 2]**2,
        sigmas[j - 1]**4 / sigmas[j - 2]**2
      ))
    }
  }

  triangle.proj <- triangle

  for (i in seq(2, ndev)) {
    latest.colidx <- ndev + 1 - i
    triangle.proj[i, (latest.colidx + 1):ndev] <- triangle[i, latest.colidx] *
      cumprod(dev.facs[latest.colidx:(ndev - 1)])
  }

  latest <- triangle[row(triangle) + col(triangle) == ndev + 1][1:(ndev - 1)]
  reserve <- triangle.proj[2:ndev, ndev] - rev(latest)

  sums <- rep(0, ndev)
  for (j in 1:ndev) {
    nrows <- ndev - j
    sums[j] <- sum(triangle[1:nrows, j])
  }

  rmse <- rep(0, ndev - 1)
  for (i in 2:ndev) {
    j <- ndev + 1 - i
    rmse[i - 1] <- sqrt(triangle.proj[i, ndev]**2 *
      sum((sigmas[j:(ndev - 1)]**2 / dev.facs[j:(ndev - 1)]**2) *
        (1 / triangle.proj[i, j:(ndev - 1)] + 1 / sums[j:(ndev - 1)])))
  }

  res <- list(
    dev.facs = dev.facs,
    sigmas = sigmas,
    reserve = reserve,
    rmse = rmse
    )
  return(res)
}

ndev <- ncol(UKMotor)

res <- mackFit(UKMotor)
res <- data.frame(
  idx = 2:ndev, 2,
  devfacs = res$dev.facs,
  sigmas = res$sigmas,
  reserve = res$reserve,
  prederror = res$rmse
)

write.csv(
  round(res, 2),
  file = file.path(res.dir, "mack_bench.csv"),
  quote = FALSE,
  row.names = FALSE
)