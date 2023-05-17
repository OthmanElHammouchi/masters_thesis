library(claimsBoot)
suppressPackageStartupMessages(library(ChainLadder))
library(Hmisc)

triangle <- UKMotor
boot.types <- c("residuals")

# Single outlier
# mean.factors <- c(100)
# sd.factors <- c(1)
# single.res <- mackSim(triangle, "single", 10, mean.factors, sd.factors, boot.types)

# Example of simulated triangle from Fortran routine
# The way we obtain it for the moment is very primitive, we just interrupt the routine
triangle <- matrix(c(
  c(3511, 4001, 4355, 4295, 4150, 5102, 6283),
  c(663172.2816, 7455.209576, 8291.400468, 7850.695796, 7842.247430, 9389.750391, 0),
  c(847468.650625, 9760.09601149, 10770.1962569, 9992.00281948, 9869.05642957, 0, 0),
  c(972420.52285004, 11230.454541154, 12376.929359844, 10968.636365075, 0, 0, 0),
  c(1069036.050766686, 12303.744241880328, 13605.408490054315, 0, 0, 0, 0),
  c(1123357.236118498, 12949.43337276258, 0, 0, 0, 0, 0),
  c(1154289.5436912179, 0, 0, 0, 0, 0, 0)), nrow = 7)
ndev <- ncol(triangle)

dev.facs <- rep(0, ndev - 1)
sigmas <- rep(0, ndev - 1)
shifts <- matrix(rep(0, (ndev - 1)**2), nrow = ndev - 1)
sds <- matrix(rep(0, (ndev - 1)**2), nrow = ndev - 1)
means <- matrix(rep(0, (ndev - 1)**2), nrow = ndev - 1)

resids.ln <- matrix(rep(0, (ndev - 1)**2), nrow = ndev - 1)
resids.standard <- matrix(rep(0, (ndev - 1)**2), nrow = ndev - 1)
resids.student <- matrix(rep(0, (ndev - 1)**2), nrow = ndev - 1)

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

  shifts[1:nrows, j] <- dev.facs[j] * sqrt(triangle[1:nrows, j]) / sigmas[j]
  sds[1:nrows, j] <- sqrt(log(1 + 1 / shifts[1:nrows, j]**2))
  means[1:nrows, j] <- log(shifts[1:nrows, j]) - sds[1:nrows, j]**2 / 2
  eps <- (triangle[1:nrows, j + 1] - dev.facs[j] * triangle[1:nrows, j]) / (sigmas[j] * sqrt(triangle[1:nrows, j]))

  resids.ln[1:nrows, j] <- (log(eps + shifts[1:nrows, j]) - means[1:nrows, j]) / sds[1:nrows, j]
  resids.standard[1:nrows, j] <- rstandard(model)
  resids.student[1:nrows, j] <- rstudent(model)

}

rownames(triangle) <- colnames(triangle) <- seq_len(ndev)
rownames(resids.ln) <- colnames(resids.ln) <- 2:ndev
rownames(resids.standard) <- colnames(resids.standard) <- 2:ndev
rownames(resids.student) <- colnames(resids.student) <- 2:ndev

triangle <- round(triangle, 2)
resids.ln <- round(resids.ln, 2)
resids.standard <- round(resids.standard, 2)
resids.student <- round(resids.student, 2)

triangle[triangle == 0] <- NA
resids.ln[resids.ln == 0] <- NA
resids.standard[resids.standard == 0] <- NA
resids.student[resids.student == 0] <- NA

resids.student[is.nan(resids.student)] <- NA
resids.standard[is.nan(resids.standard)] <- NA

latex(triangle,
  title = "",
  file = "results/sim_triangle_example.tex",
  cgroup = c("Dev"),
  n.cgroup = c(ndev),
  booktabs = TRUE,
  rowlabel.just = "r",
  rowlabel = "Origin",
  table.env = FALSE
)

latex(resids.standard,
  title = "",
  file = "results/resids_standardised_example.tex",
  cgroup = c("Dev"),
  n.cgroup = c(ndev - 1),
  booktabs = TRUE,
  rowlabel.just = "r",
  rowlabel = "Origin",
  table.env = FALSE
)

latex(resids.ln,
  title = "",
  file = "results/resids_log_normal_example.tex",
  cgroup = c("Dev"),
  n.cgroup = c(ndev - 1),
  booktabs = TRUE,
  rowlabel.just = "r",
  rowlabel = "Origin",
  table.env = FALSE
)

latex(resids.student,
  title = "",
  file = "results/resids_studentised_example.tex",
  cgroup = c("Dev"),
  n.cgroup = c(ndev - 1),
  booktabs = TRUE,
  rowlabel.just = "r",
  rowlabel = "Origin",
  table.env = FALSE
)

for (path in c("results/sim_triangle_example.tex",
  "results/resids_standardised_example.tex",
  "results/resids_log_normal_example.tex",
  "results/resids_studentised_example.tex")) {

  file.conn <- file(path)
  str <- readLines(file.conn)
  str <- sub("cline", "cmidrule", str)
  writeLines(str, file.conn)
  close(file.conn)
}
