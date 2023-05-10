#' Perform pairs bootstrap for the Mack chain ladder
#'
#' FUNCTION_DESCRIPTION
#'
#' @param triangle Cumulative claims triangle.
#' @param nboot Number of bootstrap iterations.
#' @param nsim Number of iterations for the predictive distribution.
#' @param sim.type Type of simulation for the predictive distribution.
#'
#' @return RETURN_DESCRIPTION
#' @export
mackPairsBoot <- function(triangle, nboot, nsim, sim.type) {
  ndev <- ncol(triangle)
  npts <- (ndev**2 + ndev) / 2
  nresids <- (ndev**2 - ndev) / 2

  dev.facs <- rep(0, ndev - 1)
  sigmas <- rep(0, ndev - 1)
  resids <- matrix(rep(0, (ndev - 1)**2), nrow = ndev - 1)

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

    resids.temp <- rstandard(model) # leverage-adjusted and standardised
    resids.temp <- resids.temp - mean(resids.temp)
    resids[1:nrows, j] <- resids.temp
  }

  combine <- function(list1, list2) {
    res <- list()
    res$dev.facs <- rbind(list1$dev.facs, list2$dev.facs)
    res$sigmas <- rbind(list1$sigmas, list2$sigmas)
    res$reserve <- rbind(list1$reserve, list2$reserve)
    res$pred.error <- rbind(list1$pred.error, list2$pred.error)
    res$pred.dist <- c(list1$pred.dist, list2$pred.dist)
    res$fit.dist <- c(list1$fit.dist, list2$fit.dist)
    return(res)
  }

  pgb <- txtProgressBar(max = nboot)
  progress <- function(n) setTxtProgressBar(pgb, n)
  opts <- list(progress = progress)

  res.list <- foreach(
    iboot = 1:nboot,
    .combine = combine,
    .options.snow = opts) %dopar% {
    dev.facs.boot <- rep(0, ndev - 1)
    sigmas.boot <- rep(0, ndev - 1)

    for (j in seq_len(ndev - 1)) {
      nrows <- ndev - j

      cols <- triangle[1:nrows, c(j, j + 1), drop = FALSE]
      cols <- cols[sample(1:nrows, nrows, replace = TRUE), , drop = FALSE]

      indiv <- cols[, 2] / cols[, 1]
      dev.facs.boot[j] <- sum(cols[, 2]) / sum(cols[, 1])
      if (j < ndev - 1) {
        sigmas.boot[j] <- sum(cols[, 1] * (indiv -
          dev.facs.boot[j])**2) / (nrows - 1)
      } else {
        sigmas.boot[j] <-
          sqrt(min(sigmas.boot[j - 1]**4 / sigmas.boot[j - 2]**2,
            sigmas.boot[j - 1]**2, sigmas.boot[j - 2]**2, na.rm = TRUE))
      }
    }

    triangle.proj <- triangle

    for (i in seq(2, ndev)) {
      latest.colidx <- ndev + 1 - i
      triangle.proj[i, (latest.colidx + 1):ndev] <- triangle[i, latest.colidx] *
        cumprod(dev.facs.boot[latest.colidx:(ndev - 1)])
    }

    latest <- triangle[row(triangle) + col(triangle) == ndev + 1][1:(ndev - 1)]
    reserve.proj <- triangle.proj[2:ndev, ndev] - rev(latest)

    pred.error <- matrix(rep(0, nsim * (ndev - 1)), nrow = nsim)
    pred.dist <- rep(0, nsim)
    fit.dist <- rep(0, nsim)
    if (sim.type == "parametric") {
      for (isim in 1:nsim) {
        triangle.pred <- triangle
        triangle.fit <- triangle
        for (diagidx in 1:(ndev - 1)) {
          for (rowidx in (diagidx + 1):ndev) {
            colidx <- ndev + diagidx + 1 - rowidx
            # off-diagonal elements satisfy i + j = I + 1 + (diagonal number)
            triangle.pred[rowidx, colidx] <-
              rnorm(
                1,
                triangle.pred[rowidx, colidx - 1] * dev.facs.boot[colidx - 1],
                sqrt(triangle.pred[rowidx, colidx - 1]) *
                  sigmas.boot[colidx - 1]
              )
            triangle.fit[rowidx, colidx] <-
              rnorm(
                1,
                triangle.fit[rowidx, colidx - 1] * dev.facs[colidx - 1],
                sqrt(triangle.fit[rowidx, colidx - 1]) *
                  sigmas[colidx - 1]
              )
          }
        }
        reserve.fit <- triangle.fit[2:ndev, ndev] - rev(latest)
        reserve.pred <- triangle.pred[2:ndev, ndev] - rev(latest)

        pred.error[isim, ] <- reserve.proj - reserve.fit
        pred.dist[isim] <- sum(reserve.pred)
        fit.dist[isim] <- sum(reserve.fit)
      }
    } else {
      for (isim in 1:nsim) {
        resids.sim <- matrix(rep(0, ndev**2), nrow = ndev)
        resids.sim[row(resids.sim) + col(resids.sim) > ndev + 1] <- sample(
          resids[row(resids) + col(resids) <= ndev & col(resids) != ndev - 1],
          replace = TRUE,
          size = nresids - ndev + 1
        )
        resids.sim <- resids.sim[, -1]

        triangle.fit <- triangle
        triangle.pred <- triangle
        for (diagidx in 1:(ndev - 1)) {
          for (rowidx in (diagidx + 1):ndev) {
            colidx <- ndev + diagidx + 1 - rowidx
            # off-diagonal elements satisfy i + j = I + 1 + (diagonal number)
            mean <- triangle.pred[rowidx, colidx - 1] * dev.facs.boot[colidx - 1] 

            sd <- sigmas.boot[colidx - 1] * sqrt(triangle.pred[rowidx, colidx - 1]) 

            triangle.pred[rowidx, colidx] <- mean + sd * resids.sim[rowidx, colidx - 1] 


            mean <- triangle.fit[rowidx, colidx - 1] * dev.facs[colidx - 1]
            sd <- sigmas[colidx - 1] * sqrt(triangle.fit[rowidx, colidx - 1]) 

            triangle.fit[rowidx, colidx] <- mean + sd * resids.sim[rowidx, colidx - 1] 

          }
        }
        reserve.pred <- triangle.pred[2:ndev, ndev] - rev(latest)
        reserve.fit <- triangle.fit[2:ndev, ndev] - rev(latest)

        pred.error[isim, ] <- reserve.proj - reserve.fit
        pred.dist[isim] <- sum(reserve.pred)
        fit.dist[isim] <- sum(reserve.fit)
      }
    }

    list(
      dev.facs = dev.facs.boot,
      sigmas = sigmas.boot,
      reserve = reserve.proj,
      pred.error = pred.error,
      pred.dist = pred.dist,
      fit.dist = fit.dist
    )
  }
  close(pgb)

  res <- list(
    est = data.frame(
      idx = 2:ndev,
      devfacs = colMeans(res.list$dev.facs),
      sigmas = colMeans(res.list$sigmas),
      reserve = colMeans(res.list$reserve),
      prederror = sqrt(colMeans(res.list$pred.error**2))
    ),
    dist = data.frame(
      fit.dist = res.list$fit.dist,
      pred.dist = res.list$pred.dist
    )
  )
  return(res)
}

#' Perform semiparametric bootstrap for the Mack chain ladder
#'
#' FUNCTION_DESCRIPTION
#'
#' @param triangle DESCRIPTION.
#' @param nboot DESCRIPTION.
#' @param nsim DESCRIPTION.
#' @param cond DESCRIPTION.
#' @param resids.type DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @export
mackSemiParamBoot <- function(triangle, nboot, nsim, cond, resids.type) {
  ndev <- ncol(triangle)
  nresids <- (ndev**2 - ndev) / 2

  dev.facs <- rep(0, ndev - 1)
  sigmas <- rep(0, ndev - 1)
  resids <- matrix(rep(0, (ndev - 1)**2), nrow = ndev - 1)

  if (resids.type == "log-normal") {
    shifts <- matrix(rep(0, (ndev - 1)**2), nrow = ndev - 1)
    sds <- matrix(rep(0, (ndev - 1)**2), nrow = ndev - 1)
    means <- matrix(rep(0, (ndev - 1)**2), nrow = ndev - 1)
  }

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

    if (resids.type == "standardised") {
      resids[1:nrows, j] <- rstandard(model) # leverage-adjusted and standardised

    } else if (resids.type == "studentised") {
      resids[1:nrows, j] <- rstudent(model)
    } else {

      shifts[1:nrows, j] <- dev.facs[j] * sqrt(triangle[1:nrows, j]) / sigmas[j]
      sds[1:nrows, j] <- sqrt(log(1 + 1 / shifts[1:nrows, j]**2))
      means[1:nrows, j] <- log(shifts[1:nrows, j]) - sds[1:nrows, j]**2 / 2
      eps <- (triangle[1:nrows, j + 1] - dev.facs[j] * triangle[1:nrows, j]) / (sigmas[j] * sqrt(triangle[1:nrows, j]))
      resids[1:nrows, j] <- (log(eps + shifts[1:nrows, j]) - means[1:nrows, j]) / sds[1:nrows, j]
    }
  }

  resids[, c(ndev - 2, ndev - 1)] <- 0
  resids <- resids - mean(resids[row(resids) + col(resids) <= ndev & !(col(resids) %in% c(ndev - 2, ndev - 1))])
  resids[row(resids) + col(resids) > ndev | col(resids) %in% c(ndev - 2, ndev - 1)] <- 0

  combine <- function(list1, list2) {
    res <- list()
    res$dev.facs <- rbind(list1$dev.facs, list2$dev.facs)
    res$sigmas <- rbind(list1$sigmas, list2$sigmas)
    res$reserve <- rbind(list1$reserve, list2$reserve)
    res$pred.error <- rbind(list1$pred.error, list2$pred.error)
    res$pred.dist <- c(list1$pred.dist, list2$pred.dist)
    res$fit.dist <- c(list1$fit.dist, list2$fit.dist)
    return(res)
  }

  pgb <- txtProgressBar(max = nboot)
  progress <- function(n) setTxtProgressBar(pgb, n)
  opts <- list(progress = progress)

  res.list <- foreach(
    iboot = 1:nboot,
    .combine = combine,
    .options.snow = opts) %dopar% {
    resids.boot <- resids
    vals <- sample(
      resids[row(resids) + col(resids) <= ndev & !(col(resids) %in% c(ndev - 2, ndev - 1))],
      replace = TRUE,
      size = nresids
    )
    resids.boot[row(resids.boot) + col(resids.boot) <= ndev] <- vals

    if (resids.type == "log-normal") {

      idxs <- matrix(rep(0, 2 * nresids), ncol = 2)
      for (i in 1:nresids) {
        idxs[i, ] <- which(resids == vals[i], arr.ind = TRUE)
      }
      resids.boot[row(resids.boot) + col(resids.boot) <= ndev] <- exp(means[idxs] + sds[idxs] *
        resids.boot[row(resids.boot) + col(resids.boot) <= ndev]) - shifts[idxs]
    }

    dev.facs.boot <- rep(0, ndev - 1)
    sigmas.boot <- rep(0, ndev - 1)

    triangle.boot <- matrix(rep(0, ndev**2), nrow = ndev)
    triangle.boot[, 1] <- triangle[, 1]

    if (cond) {
      for (j in 1:(ndev - 1)) {
        nrows <- ndev - j
        triangle.boot[1:nrows, j + 1] <- triangle[1:nrows, j] * dev.facs[j] +
          resids.boot[1:nrows, j] * sigmas[j] * sqrt(triangle[1:nrows, j])

        model <- lm(
          y ~ x - 1,
          weights = 1 / triangle.boot[1:nrows, j],
          data = data.frame(
            x = triangle.boot[1:nrows, j],
            y = triangle.boot[1:nrows, j + 1]
          )
        )

        dev.facs.boot[j] <- unname(model$coefficients)
        if (j < ndev - 1) {
          sigmas.boot[j] <- summary(model)$sigma
        } else {
          sigmas.boot[j] <- sqrt(min(
            sigmas.boot[(j - 1)]**2,
            sigmas.boot[j - 2]**2,
            sigmas.boot[j - 1]**4 / sigmas.boot[j - 2]**2
          ))
        }
      }

    } else {
      for (j in 1:(ndev - 1)) {
        nrows <- ndev - j
        triangle.boot[1:nrows, j + 1] <- triangle.boot[1:nrows, j] * dev.facs[j] +
          resids.boot[1:nrows, j] * sigmas[j] * sqrt(triangle.boot[1:nrows, j])

        model <- lm(
          y ~ x - 1,
          weights = 1 / triangle.boot[1:nrows, j],
          data = data.frame(
            x = triangle.boot[1:nrows, j],
            y = triangle.boot[1:nrows, j + 1]
          )
        )

        dev.facs.boot[j] <- unname(model$coefficients)
        if (j < ndev - 1) {
          sigmas.boot[j] <- summary(model)$sigma
        } else {
          sigmas.boot[j] <- sqrt(min(
            sigmas.boot[(j - 1)]**2,
            sigmas.boot[j - 2]**2,
            sigmas.boot[j - 1]**4 / sigmas.boot[j - 2]**2
          ))
        }
      }
    }

    triangle.proj <- triangle

    for (i in seq(2, ndev)) {
      latest.colidx <- ndev + 1 - i
      triangle.proj[i, (latest.colidx + 1):ndev] <- triangle[i, latest.colidx] *
        cumprod(dev.facs.boot[latest.colidx:(ndev - 1)])
    }

    latest <- triangle[row(triangle) + col(triangle) == ndev + 1][1:(ndev - 1)]
    reserve.proj <- triangle.proj[2:ndev, ndev] - rev(latest)

    pred.error <- matrix(rep(0, nsim * (ndev - 1)), nrow = nsim)
    pred.dist <- rep(0, nsim)
    fit.dist <- rep(0, nsim)

    for (isim in 1:nsim) {
      resids.sim <- matrix(rep(0, ndev**2), nrow = ndev)
      vals <- sample(
        resids[row(resids) + col(resids) <= ndev & !(col(resids) %in% c(ndev - 2, ndev - 1))],
        replace = TRUE,
        size = nresids
      )
      resids.sim[row(resids.sim) + col(resids.sim) > ndev + 1] <- vals

      if (resids.type == "log-normal") {
        idxs <- matrix(rep(0, 2 * nresids), ncol = 2)
        for (i in 1:nresids) {
          idxs[i, ] <- which(resids == vals[i], arr.ind = TRUE)
        }
        resids.sim[row(resids.sim) + col(resids.sim) > ndev + 1] <- exp(means[idxs] + sds[idxs] *
          resids.sim[row(resids.sim) + col(resids.sim) > ndev + 1]) - shifts[idxs]
      }
      resids.sim <- resids.sim[, -1]

      triangle.fit <- triangle
      triangle.pred <- triangle
      for (diagidx in 1:(ndev - 1)) {
        for (rowidx in (diagidx + 1):ndev) {
          colidx <- ndev + diagidx + 1 - rowidx
          # off-diagonal elements satisfy i + j = I + 1 + (diagonal number)

          mean <- triangle.pred[rowidx, colidx - 1] * dev.facs.boot[colidx - 1]
          sd <- sigmas.boot[colidx - 1] * sqrt(triangle.pred[rowidx, colidx - 1])
          triangle.pred[rowidx, colidx] <- mean + sd * resids.sim[rowidx, colidx - 1]

          mean <- triangle.fit[rowidx, colidx - 1] * dev.facs[colidx - 1]
          sd <- sigmas[colidx - 1] * sqrt(triangle.fit[rowidx, colidx - 1])
          triangle.fit[rowidx, colidx] <- mean + sd * resids.sim[rowidx, colidx - 1]

        }
      }
      reserve.pred <- triangle.pred[2:ndev, ndev] - rev(latest)
      reserve.fit <- triangle.fit[2:ndev, ndev] - rev(latest)

      pred.error[isim, ] <- reserve.proj - reserve.fit
      pred.dist[isim] <- sum(reserve.pred)
      fit.dist[isim] <- sum(reserve.fit)
    }

    list(
      dev.facs = dev.facs.boot,
      sigmas = sigmas.boot,
      reserve = reserve.proj,
      pred.error = pred.error,
      pred.dist = pred.dist,
      fit.dist = fit.dist
    )

  }
  close(pgb)

  res <- list(
    est = data.frame(
      idx = 2:ndev,
      devfacs = colMeans(res.list$dev.facs),
      sigmas = colMeans(res.list$sigmas),
      reserve = colMeans(res.list$reserve),
      prederror = sqrt(colMeans(res.list$pred.error**2))
    ),
    dist = data.frame(
      pred.dist = res.list$pred.dist,
      fit.dist = res.list$fit.dist
    )
  )

  return(res)
}

#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param triangle DESCRIPTION.
#' @param nboot DESCRIPTION.
#' @param nsim DESCRIPTION.
#' @param cond DESCRIPTION.
#' @param dist DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @export
mackParamBoot <- function(triangle, nboot, nsim, cond, dist) {
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

  combine <- function(list1, list2) {
    res <- list()
    res$dev.facs <- rbind(list1$dev.facs, list2$dev.facs)
    res$sigmas <- rbind(list1$sigmas, list2$sigmas)
    res$reserve <- rbind(list1$reserve, list2$reserve)
    res$pred.error <- rbind(list1$pred.error, list2$pred.error)
    res$pred.dist <- c(list1$pred.dist, list2$pred.dist)
    res$fit.dist <- c(list1$fit.dist, list2$fit.dist)
    return(res)
  }

  pgb <- txtProgressBar(max = nboot)
  progress <- function(n) setTxtProgressBar(pgb, n)
  opts <- list(progress = progress)

  res.list <- foreach(
    iboot = 1:nboot,
    .combine = combine,
    .packages = "MASS",
    .options.snow = opts) %dopar% {
    dev.facs.boot <- rep(0, ndev - 1)
    sigmas.boot <- rep(0, ndev - 1)

    triangle.boot <- matrix(rep(0, ndev**2), nrow = ndev)
    triangle.boot[, 1] <- triangle[, 1]

    if (cond) {
      for (j in 1:(ndev - 1)) {
        nrows <- ndev - j

        triangle.boot[1:nrows, j + 1] <- mvrnorm(
          1,
          mu = triangle[1:nrows, j] * dev.facs[j],
          Sigma = diag(triangle[1:nrows, j] * sigmas[j]**2, nrow = nrows)
        )

        model <- lm(
          y ~ x - 1,
          weights = 1 / triangle.boot[1:nrows, j],
          data = data.frame(
            x = triangle.boot[1:nrows, j],
            y = triangle.boot[1:nrows, j + 1]
          )
        )

        dev.facs.boot[j] <- unname(model$coefficients)
        if (j < ndev - 1) {
          sigmas.boot[j] <- summary(model)$sigma
        } else {
          sigmas.boot[j] <- sqrt(min(
            sigmas.boot[(j - 1)]**2,
            sigmas.boot[j - 2]**2,
            sigmas.boot[j - 1]**4 / sigmas.boot[j - 2]**2
          ))
        }
      }
    } else {
      for (j in 1:(ndev - 1)) {
        nrows <- ndev - j

        triangle.boot[1:nrows, j + 1] <- mvrnorm(
          1,
          mu = triangle.boot[1:nrows, j] * dev.facs[j],
          Sigma = diag(triangle.boot[1:nrows, j] * sigmas[j]**2, nrow = nrows)
        )

        model <- lm(
          y ~ x - 1,
          weights = 1 / triangle.boot[1:nrows, j],
          data = data.frame(
            x = triangle.boot[1:nrows, j],
            y = triangle.boot[1:nrows, j + 1]
          )
        )

        dev.facs.boot[j] <- unname(model$coefficients)
        if (j < ndev - 1) {
          sigmas.boot[j] <- summary(model)$sigma
        } else {
          sigmas.boot[j] <- sqrt(min(
            sigmas.boot[(j - 1)]**2,
            sigmas.boot[j - 2]**2,
            sigmas.boot[j - 1]**4 / sigmas.boot[j - 2]**2
          ))
        }
      }
    }

    triangle.proj <- triangle

    for (i in seq(2, ndev)) {
      latest.colidx <- ndev + 1 - i
      triangle.proj[i, (latest.colidx + 1):ndev] <- triangle[i, latest.colidx] *
        cumprod(dev.facs.boot[latest.colidx:(ndev - 1)])
    }

    latest <- triangle[row(triangle) + col(triangle) == ndev + 1][1:(ndev - 1)]
    reserve.proj <- triangle.proj[2:ndev, ndev] - rev(latest)

    pred.error <- matrix(rep(0, nsim * (ndev - 1)), nrow = nsim)
    pred.dist <- rep(0, nsim)
    fit.dist <- rep(0, nsim)
    for (isim in 1:nsim) {
      triangle.pred <- triangle
      triangle.fit <- triangle
      for (diagidx in 1:(ndev - 1)) {
        for (rowidx in (diagidx + 1):ndev) {
          colidx <- ndev + diagidx + 1 - rowidx
          # off-diagonal elements satisfy i + j = I + 1 + (diagonal number)
          if (dist == "normal") {
            triangle.pred[rowidx, colidx] <-
              rnorm(
                1,
                triangle.pred[rowidx, colidx - 1] * dev.facs.boot[colidx - 1], 

                sqrt(triangle.pred[rowidx, colidx - 1]) *
                  sigmas.boot[colidx - 1]
              )
            triangle.fit[rowidx, colidx] <-
              rnorm(
                1,
                triangle.fit[rowidx, colidx - 1] * dev.facs[colidx - 1],
                sqrt(triangle.fit[rowidx, colidx - 1]) *
                  sigmas[colidx - 1]
              )
          } else {
            alpha <- dev.facs.boot[colidx - 1]**2 *
              triangle.pred[rowidx, colidx - 1] / sigmas.boot[colidx - 1]**2
            beta <- dev.facs.boot[colidx - 1] / (sigmas.boot[colidx - 1]**2)
            triangle.pred[rowidx, colidx] <- rgamma(1, alpha, beta)

            alpha <- dev.facs[colidx - 1]**2 *
              triangle.fit[rowidx, colidx - 1] / sigmas[colidx - 1]**2
            beta <- dev.facs[colidx - 1] / (sigmas[colidx - 1]**2)
            triangle.fit[rowidx, colidx] <- rgamma(1, alpha, beta)
          }
        }
        reserve.pred <- triangle.pred[2:ndev, ndev] - rev(latest)
        reserve.fit <- triangle.fit[2:ndev, ndev] - rev(latest)

        pred.error[isim, ] <- reserve.proj - reserve.pred
        pred.dist[isim] <- sum(reserve.pred)
        fit.dist[isim] <- sum(reserve.fit)
      }
    }

    list(
      dev.facs = dev.facs.boot,
      sigmas = sigmas.boot,
      reserve = reserve.proj,
      pred.error = pred.error,
      pred.dist = pred.dist,
      fit.dist = fit.dist
    )

  }
  close(pgb)

  res <- list(
    est = data.frame(
      idx = 2:ndev,
      devfacs = colMeans(res.list$dev.facs),
      sigmas = colMeans(res.list$sigmas),
      reserve = colMeans(res.list$reserve),
      prederror = sqrt(colMeans(res.list$pred.error**2))
    ),
    dist = data.frame(
      pred.dist = res.list$pred.dist,
      fit.dist = res.list$fit.dist
    )
  )
  return(res)
}