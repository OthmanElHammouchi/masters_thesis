#!/usr/bin/env Rscript
source(file.path("scripts", "example", "common.r"))
set.seed(42)
files <- setdiff(list.files(example.dir, full.names = TRUE), file.path(example.dir, c("main.r", "common.r")))
for (file in files) {
  print(paste0("File: ", file))
  source(file)
}
stopCluster(cl)