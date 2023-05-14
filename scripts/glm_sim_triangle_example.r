library(claimsBoot)
suppressPackageStartupMessages(library(ChainLadder))
library(Hmisc)

# boot.types <- c("residuals")
# factors <- c(100)
# single.res <- glmSim(UKMotor, "single", 10, factors, boot.types)

# Example of simulated triangle from Fortran routine
# The way we obtain it for the moment is very primitive, we just interrupt the routine

ndev <- ncol(UKMotor)
rownames(UKMotor) <- colnames(UKMotor) <- seq_len(ndev)
long <- as.data.frame(UKMotor)
long <- as.data.frame(lapply(long, as.numeric))

long.upper <- long[long$origin + long$dev <= ndev + 1, ]
long.lower <- long[long$origin + long$dev > ndev + 1, ]

model <- glm(value ~ factor(origin) + factor(dev), quasipoisson, long.upper)
phi <- summary(model)$dispersion
fitted <- as.triangle(cbind(long.upper[, 1:2], data.frame(value = predict(model, long.upper, type = "response"))))
resids <- as.triangle(cbind(long.upper[, 1:2], data.frame(value = resid(model, type = "pearson"))))

rownames(resids) <- colnames(resids) <- seq_len(ndev)
resids <- round(zapsmall(resids), 2)

latex(resids,
  title = "",
  file = "results/glm_resids_original_example.tex",
  cgroup = c("Dev"),
  n.cgroup = c(ndev),
  booktabs = TRUE,
  rowlabel.just = "r",
  rowlabel = "Origin",
  table.env = FALSE,
  where = "!htb"
)

ndev <- 7
triangle <- as.triangle(matrix(c(
  c(376313, 4001, 4355, 4295, 4150, 5102, 6283),
  c(6726, 7703, 8287, 7750, 7897, 9650, 0),
  c(8992, 9981, 10233, 9773, 10217, 0, 0),
  c(10704, 11161, 11755, 11093, 0, 0, 0),
  c(11763, 12117, 12993, 0, 0, 0, 0),
  c(12350, 12746, 0, 0, 0, 0, 0),
  c(12690, 0, 0, 0, 0, 0, 0)), nrow = ndev))

long <- as.data.frame(triangle)
long <- as.data.frame(lapply(long, as.numeric))

long.upper <- long[long$origin + long$dev <= ndev + 1, ]
long.lower <- long[long$origin + long$dev > ndev + 1, ]

model <- glm(value ~ factor(origin) + factor(dev), quasipoisson, long.upper)
phi <- summary(model)$dispersion
fitted <- as.triangle(cbind(long.upper[, 1:2], data.frame(value = predict(model, long.upper, type = "response"))))
resids <- as.triangle(cbind(long.upper[, 1:2], data.frame(value = resid(model, type = "pearson"))))

rownames(triangle) <- colnames(triangle) <- seq_len(ndev)
rownames(resids) <- colnames(triangle) <- seq_len(ndev)
triangle <- round(triangle / 1e3, 2) # in thousands

resids <- round(zapsmall(resids), 2)

latex(triangle,
  title = "",
  file = "results/glm_triangle_example.tex",
  cgroup = c("Dev"),
  n.cgroup = c(ndev),
  caption = "Simulated triangle where observation $X_{11}$ has been perturbed, with $c_\\lambda = 100$ (thousands)",
  caption.loc = "bottom",
  label = "tab:glm-triangle-example",
  booktabs = TRUE,
  rowlabel.just = "r",
  rowlabel = "Origin",
  where = "!htb"
)

latex(resids,
  title = "",
  file = "results/glm_resids_perturbed_example.tex",
  cgroup = c("Dev"),
  n.cgroup = c(ndev),
  booktabs = TRUE,
  rowlabel.just = "r",
  rowlabel = "Origin",
  table.env = FALSE,
  where = "!htb"
)
