mackConfig <- function(ndev,
                       mean.factors,
                       sd.factors,
                       sim.type = c("single", "calendar", "origin"),
                       boot.type = c("residuals, parametric, pairs")
) {
  key <- patternBreak:::.global$key
  if (boot.type != "pairs") {
    conds <- c(TRUE, FALSE)
    if (boot.type == "residuals") {
      opt <- as.double(key[c("standardised", "studentised")])
    } else if (boot.type == "parametric") {
      opt <- as.double(key[c("normal", "gamma")])
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

    if (boot.type != "pairs") {
      idxs <- do.call(expand.grid,
        list(
          seq_len(npts),
          seq_along(mean.factors),
          seq_along(sd.factors),
          seq_len(npts),
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

      temp <- cbind(
        outlier.points[idxs[, 1], ],
        mean.factors[idxs[, 2]],
        sd.factors[idxs[, 3]],
        excl.points[idxs[, 4], ],
        opt[idxs[, 5]],
        as.integer(conds[idxs[, 6]])
      )

    } else {
      idxs <- do.call(expand.grid,
        list(
          seq_len(npts),
          seq_along(mean.factors),
          seq_along(sd.factors),
          seq_len(npts)
      ))
      config <- cbind(
        outlier.points[idxs[, 1], ],
        mean.factors[idxs[, 2]],
        sd.factors[idxs[, 3]],
        excl.points[idxs[, 4], ]
      )
      conf.int <- config
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
  } else {
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
  }

  if ("pairs" %in% res.names) {
    idx <- which(res.names == "pairs")
    df.pairs <- res.list[[idx]]
    n_row <- nrow(df.pairs)
    res.list[[idx]] <- cbind(
      df.pairs[, -ncol(df.pairs)],
      matrix(rep(NA, n_row * 2), ncol = 2),
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
                      type = "single") {

  if (type == "single") {

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
        seq_len(nrow(outlier.points)),
        seq_along(factors),
        seq_len(nrow(excl.points))
    ))

    config <- do.call(cbind,
      list(
        outlier.points[indices[, 1], ],
        factors[indices[, 2]],
        excl.points[indices[, 3], ]
    ))

    config <- as.data.table(config)

    names(config) <- c("outlier.rowidx", "outlier.colidx", "factor", "excl.rowidx", "excl.colidx")

    return(config)

  } else if (type == "calendar" || type == "origin") {

    excl.diags <- outlier.diags <- seq_len(ndev - 1)

    indices <- do.call(expand.grid,
      list(
        seq_along(outlier.diags),
        seq_along(excl.diags),
        seq_along(factors)
    ))

    config <- do.call(cbind,
      list(
        outlier.diags[indices[, 1]],
        factors[indices[, 2]],
        excl.diags[indices[, 3]]
    ))

    config <- as.data.table(config)

    names(config) <- c("outlier.diagidx", "factor", "excl.diagidx")

    return(config)
  }
}

#' Convert cumulative triangle to incremental one.
#'
#' FUNCTION_DESCRIPTION
#'
#' @param triangle DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
#' @export
cum2incr <- function(triangle) {

  ndev <- ncol(triangle)

  for (j in seq(2, ndev)) {
    for (i in seq(1, ndev + 1 - j)) {
      triangle[i, j] <- triangle[i, j] - sum(triangle[i, 1:(j - 1)])
    }
  }

  return(triangle)
}


#' Convert incremental triangle to cumulative one.
#'
#' FUNCTION_DESCRIPTION
#'
#' @param triangle DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
#' @export
incr2cum <- function(triangle) {

  ndev <- ncol(triangle)

  for (j in seq(ndev, 2)) {
    for (i in seq(1, ndev + 1 - j)) {
      triangle[i, j] <- triangle[i, j] + sum(triangle[i, 1:(j - 1)])
    }
  }

  return(triangle)
}
