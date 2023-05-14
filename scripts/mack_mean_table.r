library(Hmisc)
single.res <- readRDS("results/mack_single.RDS")

### Studentised residuals:
###############################################################################
res.student <- single.res[
  boot.type == "residuals" &
    opt == "studentised" &
    cond == TRUE &
    mean.factor == 2 &
    sd.factor == 1
]

outlier.pts <- unique(res.student[, c("outlier.rowidx", "outlier.colidx")])
outlier.pts <- outlier.pts[order(outlier.pts$outlier.colidx), ]
outlier.pts <- paste0(
  rep("(", nrow(outlier.pts)),
  outlier.pts$outlier.rowidx,
  ", ",
  outlier.pts$outlier.colidx,
  rep(")", nrow(outlier.pts))
)
by.outliers <- split(res.student, list(res.student$outlier.rowidx, res.student$outlier.colidx), drop = TRUE)
means <- lapply(by.outliers, function(df) {
  sapply(split(df, list(df$excl.rowidx, df$excl.colidx), drop = TRUE), function(df) { mean(df$reserve) })
})

mean.table <- as.data.frame(do.call(rbind, means))
markup <- matrix(rep("", prod(dim(mean.table))), nrow = nrow(mean.table))
for (i in seq_len(nrow(markup))) {
  row <- mean.table[i, ]
  markup[i, ][row == min(row)] <- "cellcolor{red}"
}
mean.table <- mean.table / 1e3 # in thousands
formatted <- as.data.frame(lapply(mean.table,
  function(col) {
    sapply(col,
      function(val) {
        ifelse(val > 1000, formatC(val, format = "e", digits = 0), formatC(val, format = "f", digits = 2))
      })
  })
)
rownames(formatted) <- colnames(formatted) <- outlier.pts

latex(formatted,
  file = "results/mean_table_studentised.tex",
  title = "",
  cgroup = c("Excluded"),
  n.cgroup = c(18),
  rowlabel = "Outlier",
  rowlabel.just = "r",
  booktabs = TRUE,
  caption = "Predictive reserve mean for different contaminated and excluded points (thousands)",
  caption.loc = "bottom",
  dec = 2,
  cellTexCmds = markup,
  label = "tab:resids-student-res-mean"
)

file.conn <- file("results/mean_table_studentised.tex")
str <- readLines(file.conn)
str <- sub("cline", "cmidrule", str)
writeLines(str, file.conn)
close(file.conn)

### Pairs bootstrap:
###############################################################################
res.pairs <- single.res[
  boot.type == "pairs" &
    mean.factor == 2 &
    sd.factor == 1]

outlier.pts <- unique(res.pairs[, c("outlier.rowidx", "outlier.colidx")])
outlier.pts <- outlier.pts[order(outlier.pts$outlier.colidx), ]
outlier.pts <- paste0(
  rep("(", nrow(outlier.pts)),
  outlier.pts$outlier.rowidx,
  ", ",
  outlier.pts$outlier.colidx,
  rep(")", nrow(outlier.pts))
)
by.outliers <- split(res.pairs, list(res.pairs$outlier.rowidx, res.pairs$outlier.colidx), drop = TRUE)
means <- lapply(by.outliers, function(df) {
  sapply(split(df, list(df$excl.rowidx, df$excl.colidx), drop = TRUE), function(df) { mean(df$reserve) })
})

mean.table <- as.data.frame(do.call(rbind, means))
markup <- matrix(rep("", prod(dim(mean.table))), nrow = nrow(mean.table))
for (i in seq_len(nrow(markup))) {
  row <- mean.table[i, ]
  markup[i, ][row == min(row)] <- "cellcolor{red}"
}
mean.table <- mean.table / 1e3 # in thousands
formatted <- as.data.frame(lapply(mean.table,
  function(col) {
    sapply(col,
      function(val) {
        ifelse(val > 1000, formatC(val, format = "e", digits = 0), formatC(val, format = "f", digits = 2))
      })
  })
)
rownames(formatted) <- colnames(formatted) <- outlier.pts

latex(formatted,
  file = "results/mean_table_pairs.tex",
  title = "",
  cgroup = c("Excluded"),
  n.cgroup = c(20),
  rowlabel = "Outlier",
  rowlabel.just = "r",
  booktabs = TRUE,
  caption = "Predictive reserve mean for different contaminated and excluded points (thousands)",
  caption.loc = "bottom",
  dec = 2,
  cellTexCmds = markup,
  label = "tab:resids-pairs-res-mean"
)

file.conn <- file("results/mean_table_pairs.tex")
str <- readLines(file.conn)
str <- sub("cline", "cmidrule", str)
writeLines(str, file.conn)
close(file.conn)
