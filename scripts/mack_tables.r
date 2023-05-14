single.res <- readRDS("results/mack_single.RDS")

res.student <- single.res[
    boot.type == "residuals" &
    opt == "log-normal" &
    cond == TRUE
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

library(Hmisc)

formatted1 <- formatted[, 1:9]
formatted2 <- formatted[, 10:18]

latex(formatted1,
  file = "results/mean_table_1.tex",
  title = "",
  cgroup = c("Excluded"),
  n.cgroup = c(9),
  rowlabel = "Outlier",
  rowlabel.just = "r",
  booktabs = TRUE,
  caption = "Predictive reserve mean for different contaminated and excluded points (thousands)",
  caption.loc = "bottom",
  dec = 2,
  cellTexCmds = markup[, 1:9]
)

latex(formatted2,
  file = "results/mean_table_2.tex",
  title = "",
  cgroup = c("Excluded"),
  n.cgroup = c(9),
  rowlabel = "Outlier",
  rowlabel.just = "r",
  booktabs = TRUE,
  caption = "Predictive reserve mean for different contaminated and excluded points (thousands)",
  caption.loc = "bottom",
  dec = 2,
  cellTexCmds = markup[, 10:18]
)

library(formattable)
min_formatter <- formatter("span", style = x ~ style(color = ifelse(x == min(x), "red", "black")))
formattable(mean.table, lapply(seq_len(nrow(mean.table)), function(row) {
  area(row) ~ min_formatter
}))