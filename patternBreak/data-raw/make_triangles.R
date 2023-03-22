file.names <- c("comauto_pos.csv", "medmal_pos.csv", "othliab_pos.csv", "ppauto_pos.csv", "prodliab_pos.csv", "wkcomp_pos.csv")

for (name in file.names) {
  if (!file.exists(paste0("data-raw/", name))) {
    url <- paste0("https://www.casact.org/sites/default/files/2021-04/", name)
    download.file(url, paste0("data-raw/", name))
  }
}

# It appears that the wrong company name is given in some of the datasets,
# leading to duplicates. The GRCODE does appear to be unique though.
# Need to contact CAS about this, in the meantime exclude the affected data.

datasets <- list()

for (name in file.names) {
  dataset <- read.csv(paste0("data-raw/", name))
  if (nrow(dataset) / length(unique(dataset$GRNAME)) != 100) {
    occurrences <- table(dataset$GRNAME)
    duplicates <- names(occurrences)[occurrences > 100]
    dataset <- dataset[!(dataset$GRNAME %in% duplicates), ]
  }
  datasets <- c(datasets, list(dataset))
}

# Check if things are ok now.
for (dataset in datasets) {
  stopifnot(nrow(dataset) / length(unique(dataset$GRNAME)) == 100)
}

# Put data in triangle format.
as.triangle <- function(long.triangle) {
  start <- min(long.triangle[, 1])
  ndev <- max(long.triangle[long.triangle[, 1] == long.triangle[1, 1], 2]) - min(long.triangle[long.triangle[, 1] == long.triangle[1, 1], 2]) + 1
  nacc <- max(long.triangle[, 1]) - min(long.triangle[, 1]) + 1

  triangle <- matrix(rep(NA, ndev * nacc), nrow = nacc)
  for (i in seq_len(nacc)) {
    accident.data <- long.triangle[long.triangle[, 1] == start + i - 1, ]
    accident.start <- min(accident.data[, 1])
    for (j in seq_len(ndev)) {
      triangle[i, j] <- accident.data[accident.data[, 2] == accident.start + j - 1, 3]
    }
  }
  return(triangle)
}

triangles <- list()
triangle <- matrix(rep(0, 100), nrow = 10)
for (dataset in datasets) {
  for (code in unique(dataset$GRCODE)) {
    paid.triangle <- cbind(dataset[dataset$GRCODE == code, c("AccidentYear", "DevelopmentYear")], dataset[dataset$GRCODE == code, 7])
    names(paid.triangle) <- c("accident.year", "development.year", "loss")
    paid.triangle <- as.triangle(paid.triangle)
    triangles <- c(triangles, list(paid.triangle))
  }
}


test.triangle <- dataset[dataset$GRCODE == unique(dataset$GRCODE)[1], c("AccidentYear", "DevelopmentYear", "CumPaidLoss_D")]
names(test.triangle) <- c("accident.year", "development.year", "loss")

res <- as.triangle(test.triangle)
