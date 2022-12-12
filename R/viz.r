library(ggplot2)
library(data.table)
source("R/pattern_break.r")

library(tidyr)

results <- readRDS("results/data_objects/single_outlier.RDS")

factor <- seq(0.5, 1.5, by = 0.25)
resids.type <- c("parametric", "raw", "scaled")
boot.type <- c("conditional", "unconditional")
dist <- c("normal", "gamma")

plot.config <- genConfig(factor, resids.type, boot.type, dist)
names(plot.config) <- c("factor", "resids.type", "boot.type", "dist")

progress.bar <- txtProgressBar(min = 0, max = nconfig, initial = 0, style = 3)

for (rowidx in seq_len(nrow(plot.config))) {

    setTxtProgressBar(progress.bar, rowidx)

    factor <- plot.config$factor[rowidx]
    resids.type <- plot.config$resids.type[rowidx]
    boot.type <- plot.config$boot.type[rowidx]
    dist <- plot.config$dist[rowidx]

    plot.dt <- as.data.table(results)[
            resids.type == resids.type &
            boot.type == boot.type &
            dist == dist &
            factor == factor &
            outlier.colidx != excl.colidx &
            outlier.rowidx != excl.rowidx,
            .(reserve = unlist(reserve)),
            by = setdiff(names(results), "reserve")]

    point.dt <- as.data.table(results)[
        resids.type == resids.type &
        boot.type == boot.type &
        dist == dist &
        factor == factor &
        outlier.colidx == excl.colidx &
        outlier.rowidx == excl.rowidx,
        .(reserve = unlist(reserve)),
        by = setdiff(names(results), "reserve")]


    # density plot
    ggplot() +
        geom_density(
            data = plot.dt,
            aes(x = reserve)) +
        geom_density(
            data = point.dt,
            aes(x = reserve),
            colour = "red") +
        facet_grid(outlier.colidx ~ outlier.rowidx, scales = "free") +
        labs(
            title = "Reserve distributions for different outlier points",
            subtitle = sprintf("Perturbation factor: %.2f", factor),
            x_lab = "Reserve",
            y_lab = "Density") +
        theme(axis.text.y = element_blank())

    ggsave(sprintf("results/graphs/single_outlier/densities_%s_%s_%s_factor_%.2f.svg", dist, resids.type, boot.type, factor))

}

close(progress.bar)


densityPlot <- function(results, resids.type, boot.type, dist, factor) {

    plot.dt <- results[,
            c("lower", "upper") := list(quantile(reserve, 0.01, na.rm = TRUE), quantile(reserve, 0.999, na.rm = TRUE)), by = setdiff(names(results), "reserve")][
            resids.type == resids.type &
            boot.type == boot.type &
            dist == dist &
            factor == factor &
            outlier.colidx != excl.colidx &
            outlier.rowidx != excl.rowidx]

    point.dt <- results[,
            c("lower", "upper") := list(quantile(reserve, 0.01, na.rm = TRUE), quantile(reserve, 0.999, na.rm = TRUE)), by = setdiff(names(results), "reserve")][
            resids.type == resids.type &
            boot.type == boot.type &
            dist == dist &
            factor == factor &
            outlier.colidx == excl.colidx &
            outlier.rowidx == excl.rowidx]

    ggplot() +
        geom_density(
            data = plot.dt,
            aes(x = reserve, group = interaction(excl.colidx, excl.rowidx))) +
        geom_density(
            data = point.dt,
            aes(x = reserve),
            colour = "red") +
        facet_grid(outlier.colidx ~ outlier.rowidx, scales = "free_x") +
        labs(
            title = "Reserve distributions for different outlier points",
            subtitle = sprintf("Perturbation factor: %.2f", factor),
            x_lab = "Reserve",
            y_lab = "Density") +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1)) +
        xlim(min(plot.dt$lower, point.dt$lower), min(plot.dt$upper, point.dt$upper))

}

densityPlot(results, "scaled", "conditional", "normal", 1.25)
