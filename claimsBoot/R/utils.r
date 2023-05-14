mackConfig <- function(ndev,
                       mean.factors,
                       sd.factors,
                       sim.type,
                       boot.type
) {
  key <- claimsBoot:::.global$key
  if (boot.type != "pairs") {
    conds <- c(TRUE, FALSE)
    if (boot.type == "residuals") {
      opt <- as.double(key[c("standardised", "studentised", "lognormal")])
    } else if (boot.type == "parametric") {
      opt <- as.double(key[c("gamma", "normal")])
    }
  }

  if (sim.type == "single") {
    npts <- (ndev**2 - ndev) / 2
    outlier.points <- matrix(rep(0, 2 * npts), ncol = 2)
    k <- 1
    for (i in seq_len(ndev - 1)) {
      for (j in seq(2, ndev + 1 - i)) {
        outlier.points[k, ] <- c(i, j)
        k <- k + 1
      }
    }
    excl.points <- outlier.points

    if (boot.type == "parametric") {
      outlier.points <- outlier.points[!(outlier.points[, 1] == 1 & outlier.points[, 2] == 7), ]
      excl.points <- outlier.points
      idxs <- do.call(expand.grid,
        list(
          seq_len(nrow(outlier.points)),
          seq_along(mean.factors),
          seq_along(sd.factors),
          seq_len(nrow(excl.points)),
          seq_along(opt),
          seq_along(conds)
      ))
      config <- cbind(
        outlier.points[idxs[, 1], ],
        mean.factors[idxs[, 2]],
        sd.factors[idxs[, 3]],
        excl.points[idxs[, 4], ],
        opt[idxs[, 5]],
        as.integer(conds[idxs[, 6]])
      )

    } else if (boot.type == "residuals") {
      opt <- as.double(key[c("lognormal")])
      idxs <- do.call(expand.grid,
        list(
          seq_len(nrow(outlier.points)),
          seq_along(mean.factors),
          seq_along(sd.factors),
          seq_len(nrow(excl.points)),
          seq_along(opt),
          seq_along(conds)
      ))
      df1 <- cbind(
        outlier.points[idxs[, 1], ],
        mean.factors[idxs[, 2]],
        sd.factors[idxs[, 3]],
        excl.points[idxs[, 4], ],
        opt[idxs[, 5]],
        as.integer(conds[idxs[, 6]])
      )

      outlier.points <- outlier.points[
        !(outlier.points[, 1] == 1 & outlier.points[, 2] == 7) &
        !(outlier.points[, 1] == 1 & outlier.points[, 2] == 6) &
        !(outlier.points[, 1] == 2 & outlier.points[, 2] == 6),
      ]
      excl.points <- outlier.points
      opt <- as.double(key[c("standardised", "studentised")])
      idxs <- do.call(expand.grid,
        list(
          seq_len(nrow(outlier.points)),
          seq_along(mean.factors),
          seq_along(sd.factors),
          seq_len(nrow(outlier.points)),
          seq_along(opt),
          seq_along(conds)
      ))
      df2 <- cbind(
        outlier.points[idxs[, 1], ],
        mean.factors[idxs[, 2]],
        sd.factors[idxs[, 3]],
        excl.points[idxs[, 4], ],
        opt[idxs[, 5]],
        as.integer(conds[idxs[, 6]])
      )
      config <- rbind(df1, df2)

    } else {
      outlier.points <- outlier.points[!(outlier.points[, 1] == 1 & outlier.points[, 2] == 7), ]
      excl.points <- outlier.points
      idxs <- do.call(expand.grid,
        list(
          seq_len(nrow(outlier.points)),
          seq_along(mean.factors),
          seq_along(sd.factors),
          seq_len(nrow(excl.points))
      ))
      config <- cbind(
        outlier.points[idxs[, 1], ],
        mean.factors[idxs[, 2]],
        sd.factors[idxs[, 3]],
        excl.points[idxs[, 4], ]
      )
    }

  } else {
    excl.diags <- outlier.diags <- seq_len(ndev - 1)
    if (boot.type != "pairs") {
      idxs <- do.call(expand.grid,
        list(
          seq_along(outlier.diags),
          seq_along(mean.factors),
          seq_along(sd.factors),
          seq_along(excl.diags),
          seq_along(opt),
          seq_along(conds)
      ))
      config <- cbind(
        outlier.diags[idxs[, 1]],
        mean.factors[idxs[, 2]],
        sd.factors[idxs[, 3]],
        excl.diags[idxs[, 4]],
        opt[idxs[, 5]],
        as.integer(conds[idxs[, 6]])
      )

    } else {
      idxs <- do.call(expand.grid,
        list(
          seq_along(outlier.diags),
          seq_along(mean.factors),
          seq_along(sd.factors),
          seq_along(excl.diags)
      ))
      config <- cbind(
        outlier.diags[idxs[, 1]],
        mean.factors[idxs[, 2]],
        sd.factors[idxs[, 3]],
        excl.diags[idxs[, 4]]
      )
    }
  }
  return(config)
}

mackPost <- function(res.list, res.names, sim.type) {
  if (sim.type == "single") {
    out.names <- c(
      "boot.type",
      "outlier.rowidx",
      "outlier.colidx",
      "mean.factor",
      "sd.factor",
      "excl.rowidx",
      "excl.colidx",
      "opt",
      "cond",
      "reserve"
    )
  } else if (sim.type == "calendar") {
    out.names <- c(
      "boot.type",
      "outlier.diagidx",
      "mean.factor",
      "sd.factor",
      "excl.diagidx",
      "opt",
      "cond",
      "reserve"
    )
  } else {
    out.names <- c(
      "boot.type",
      "outlier.rowidx",
      "mean.factor",
      "sd.factor",
      "excl.rowidx",
      "opt",
      "cond",
      "reserve"
    )
  }

  if ("pairs" %in% res.names) {
    idx <- which(res.names == "pairs")
    df.pairs <- res.list[[idx]]
    nrows <- nrow(df.pairs)
    res.list[[idx]] <- cbind(
      df.pairs[, -ncol(df.pairs)],
      matrix(rep(NA, nrows * 2), ncol = 2),
      df.pairs[, ncol(df.pairs)]
    )
  }
  res.list <- lapply(seq_len(length(res.list)),
    function(idx) {
      df <- as.data.table(res.list[[idx]])
      df <- cbind(rep(res.names[idx], nrow(df)), df)
      colnames(df) <- out.names

      df[, cond := ifelse(cond == 1, TRUE, FALSE)]

      if (res.names[idx] == "residuals") {
        df[, opt := ifelse(opt == 1, "standardised", opt)]
        df[, opt := ifelse(opt == 2, "studentised", opt)]
        df[, opt := ifelse(opt == 3, "log-normal", opt)]
      } else if (res.names[idx] == "parametric") {
        df[, opt := ifelse(opt == 1, "normal", opt)]
        df[, opt := ifelse(opt == 2, "gamma", opt)]
      }
      return(df)
    }
  )
  out <- do.call(rbind, res.list)
  return(out)
}

glmConfig <- function(ndev,
                      factors,
                      sim.type,
                      boot.type) {
  key <- claimsBoot:::.global$key
  if (boot.type == "parametric") {
    opt <- as.double(key[c("normal", "gamma", "poisson")])
  }

  if (sim.type == "single") {
    if (boot.type == "parametric") {

      npts <- (ndev**2 + ndev) / 2 - 2 # exclude corner points
      outlier.points <- matrix(rep(0, 2 * npts), ncol = 2)
      k <- 1
      for (i in seq_len(ndev)) {
        for (j in seq_len(ndev + 1 - i)) {
          if (!(list(c(i, j)) %in% list(c(1, ndev), c(ndev, 1)))) {
            outlier.points[k, ] <- c(i, j)
            k <- k + 1
          }
        }
      }

      excl.points <- outlier.points

      indices <- do.call(expand.grid,
        list(
          seq_len(npts),
          seq_along(factors),
          seq_len(npts),
          seq_along(opt)
      ))

      config <- do.call(cbind,
        list(
          outlier.points[indices[, 1], ],
          factors[indices[, 2]],
          excl.points[indices[, 3], ],
          opt[indices[, 4]]
      ))

    } else {

      npts <- (ndev**2 + ndev) / 2
      outlier.points <- matrix(rep(0, 2 * npts), ncol = 2)
      k <- 1
      for (i in seq_len(ndev)) {
        for (j in seq_len(ndev + 1 - i)) {
          outlier.points[k, ] <- c(i, j)
          k <- k + 1
        }
      }

      excl.points <- outlier.points

      indices <- do.call(expand.grid,
        list(
          seq_len(npts),
          seq_along(factors),
          seq_len(npts)
      ))

      config <- do.call(cbind,
        list(
          outlier.points[indices[, 1], ],
          factors[indices[, 2]],
          excl.points[indices[, 3], ]
      ))
    }
    return(config)

  } else {

    if (sim.type == "calendar") {
      excl.idxs <- outlier.idxs <- seq_len(ndev)
    } else {
      excl.idxs <- outlier.idxs <- 2:(ndev - 1)
    }

    if (boot.type == "parametric") {
      indices <- do.call(expand.grid,
        list(
          seq_along(outlier.idxs),
          seq_along(factors),
          seq_along(excl.idxs),
          seq_along(opt)
      ))

      config <- do.call(cbind,
        list(
          outlier.idxs[indices[, 1]],
          factors[indices[, 2]],
          excl.idxs[indices[, 3]],
          opt[indices[, 4]]
      ))
    } else {
      indices <- do.call(expand.grid,
        list(
          seq_along(outlier.idxs),
          seq_along(factors),
          seq_along(excl.idxs)
      ))

      config <- do.call(cbind,
        list(
          outlier.idxs[indices[, 1]],
          factors[indices[, 2]],
          excl.idxs[indices[, 3]]
      ))
    }

    return(config)
  }
}

glmPost <- function(res.list, res.names, sim.type) {
  if (sim.type == "single") {
    out.names <- c(
      "boot.type",
      "outlier.rowidx",
      "outlier.colidx",
      "factor",
      "excl.rowidx",
      "excl.colidx",
      "opt",
      "reserve"
    )
  } else if (sim.type == "calendar") {
    out.names <- c(
      "boot.type",
      "outlier.diagidx",
      "factor",
      "excl.diagidx",
      "opt",
      "reserve"
    )
  } else {
    out.names <- c(
      "boot.type",
      "outlier.rowidx",
      "factor",
      "excl.rowidx",
      "opt",
      "reserve"
    )
  }

  if ("residuals" %in% res.names) {
    idx <- which(res.names == "residuals")
    df.residuals <- res.list[[idx]]
    nrows <- nrow(df.residuals)
    res.list[[idx]] <- cbind(
      df.residuals[, -ncol(df.residuals)],
      matrix(rep(NA, nrows), ncol = 1),
      df.residuals[, ncol(df.residuals)]
    )
  }

  res.list <- lapply(seq_len(length(res.list)),
    function(idx) {
      df <- as.data.table(res.list[[idx]])
      df <- cbind(rep(res.names[idx], nrow(df)), df)
      colnames(df) <- out.names

      if (res.names[idx] == "parametric") {
        df[, opt := ifelse(opt == 1, "normal", opt)]
        df[, opt := ifelse(opt == 2, "gamma", opt)]
        df[, opt := ifelse(opt == 3, "poisson", opt)]
      }
      return(df)
    }
  )
  out <- do.call(rbind, res.list)
  return(out)
}

#' Convert cumulative triangle to incremental one.
#'
#' @param triangle.cum Cumulative claims triangle.
#'
#' @return Incremental claims triangle.
#' @export
cum2incr <- function(triangle.cum) {

  ndev <- ncol(triangle.cum)

  for (j in seq(2, ndev)) {
    for (i in seq(1, ndev + 1 - j)) {
      triangle.cum[i, j] <- triangle.cum[i, j] - sum(triangle.cum[i, 1:(j - 1)])
    }
  }

  return(triangle.cum)
}

#' Convert incremental triangle to cumulative one.
#'
#' @param triangle.inrc Incremental claims triangle.
#'
#' @return Cumulative claims triangle.
#' @export
incr2cum <- function(triangle.incr) {

  ndev <- ncol(triangle.incr)

  for (j in seq(ndev, 2)) {
    for (i in seq(1, ndev + 1 - j)) {
      triangle.incr[i, j] <- triangle.incr[i, j] + sum(triangle.incr[i, 1:(j - 1)])
    }
  }

  return(triangle.incr)
}
