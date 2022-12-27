library(ggplot2)
library(data.table)
source("R/helpers.r")

single.outlier.results <- readRDS("results/data_objects/single_outlier.RDS")
calendar.year.results <- readRDS("results/data_objects/calendar_outlier.RDS")
origin.year.results <- readRDS("results/data_objects/origin_outlier.RDS")

ndev <- max(single.outlier.results$outlier.colidx)

factor <- seq(0.5, 1.5, by = 0.25)
resids.type <- c("parametric", "raw", "scaled")
boot.type <- c("conditional", "unconditional")
dist <- c("normal", "gamma")

plot.config <- genConfig(factor, resids.type, boot.type, dist)
names(plot.config) <- c("factor", "resids.type", "boot.type", "dist")
plot.config <- plot.config[!(resids.type == "parametric" & dist == "gamma")]

nconfig <- nrow(plot.config)

# meanPlot <- function(contaminated, uncontaminated) {

#     labels <- paste(rep("Outlier Row:"), 1:ndev)
#     names(labels) <- 1:ndev

#     p <- ggplot() +
#         geom_violin(
#             data = contaminated,
#             aes(x = as.factor(outlier.rowidx), y = reserve.mean)
#         ) +
#         geom_point(
#             data = uncontaminated,
#             aes(x = as.factor(outlier.rowidx), y = reserve.mean),
#             colour = "red"
#         ) +
#         facet_wrap(vars(outlier.colidx),
#             scales = "free",
#             labeller = as_labeller(labels)
#         ) +
#         labs(
#             title = "Reserve mean for different outlier points",
#             subtitle = sprintf("Perturbation factor: %.2f", factor),
#             x = "Outlier Column",
#             y = "Reserve Mean"
#         ) +
#         theme(
#             axis.text.y = element_blank(),
#             strip.text.x = element_text(size = 12)
#         )

#     return(p)
# }

progress.bar <- txtProgressBar(min = 0, max = nconfig, initial = 0, style = 3)

for (rowidx in seq_len(nconfig)) {

    setTxtProgressBar(progress.bar, rowidx)

    #density plots
    contaminated <- single.outlier.results[
            resids.type == plot.config$resids.type[rowidx] &
            boot.type == plot.config$boot.type[rowidx] &
            dist == plot.config$dist[rowidx] &
            factor == plot.config$factor[rowidx] &
            outlier.colidx != excl.colidx &
            outlier.rowidx != excl.rowidx]

    uncontaminated <- single.outlier.results[
        resids.type == plot.config$resids.type[rowidx] &
        boot.type == plot.config$boot.type[rowidx] &
        dist == plot.config$dist[rowidx] &
        factor == plot.config$factor[rowidx] &
        outlier.colidx == excl.colidx &
        outlier.rowidx == excl.rowidx]

    row.labels <- paste(rep("Outlier row:"), 1:ndev)
    names(row.labels) <- 1:ndev
    col.labels <- paste(rep("Outlier column:"), 1:ndev)
    names(col.labels) <- 1:ndev


    p <- ggplot() +
        geom_density(
            data = contaminated,
            aes(reserve, after_stat(scaled), group = interaction(excl.colidx, excl.rowidx))) +
        geom_density(
            data = uncontaminated,
            aes(reserve, after_stat(scaled)),
            colour = "red") +
        facet_wrap(vars(outlier.rowidx, outlier.colidx),
            scales = "free_x",
            labeller = labeller(
                outlier.rowidx = as_labeller(row.labels),
                outlier.colidx = as_labeller(col.labels))
          )  +
        labs(
            title = "Reserve distributions for different outlier points",
            subtitle = sprintf("Perturbation factor: %.2f", factor),
            x = "Reserve",
            y = "Density")

    suppressMessages(
        ggsave(
            sprintf(
                "results/graphs/single_outlier/densities_%s_%s_%s_factor_%.2f.svg",
                plot.config$resids.type[rowidx],
                plot.config$boot.type[rowidx],
                plot.config$dist[rowidx],
                plot.config$factor[rowidx]
            ),
            plot = p
        )
    )

    contaminated <- calendar.year.results[
        resids.type == plot.config$resids.type[rowidx] &
        boot.type == plot.config$boot.type[rowidx] &
        dist == plot.config$dist[rowidx] &
        factor == plot.config$factor[rowidx] &
        outlier.diagidx != excl.diagidx]

    uncontaminated <- calendar.year.results[
        resids.type == plot.config$resids.type[rowidx] &
        boot.type == plot.config$boot.type[rowidx] &
        dist == plot.config$dist[rowidx] &
        factor == plot.config$factor[rowidx] &
        outlier.diagidx == excl.diagidx]

    labels <- paste(rep("Outlier diag:"), 1:(ndev - 1))
    names(row.labels) <- 1:(ndev - 1)

    p <- ggplot() +
        geom_density(
            data = contaminated,
            aes(reserve, after_stat(scaled), group = excl.diagidx)) +
        geom_density(
            data = uncontaminated,
            aes(reserve, after_stat(scaled)),
            colour = "red") +
        facet_wrap(vars(outlier.diagidx),
            scales = "free_x",
            labeller = as_labeller(labels)
          )  +
        labs(
            title = "Reserve distributions for different outlier points",
            subtitle = sprintf("Perturbation factor: %.2f", factor),
            x = "Reserve",
            y = "Density")

    suppressMessages(
        ggsave(
            sprintf(
                "results/graphs/calendar_outlier/densities_%s_%s_%s_factor_%.2f.svg",
                plot.config$resids.type[rowidx],
                plot.config$boot.type[rowidx],
                plot.config$dist[rowidx],
                plot.config$factor[rowidx]
            ),
            plot = p
        )
    )

    contaminated <- single.outlier.results[
            resids.type == plot.config$resids.type[rowidx] &
            boot.type == plot.config$boot.type[rowidx] &
            dist == plot.config$dist[rowidx] &
            factor == plot.config$factor[rowidx] &
            outlier.rowidx != excl.rowidx]

    uncontaminated <- single.outlier.results[
        resids.type == plot.config$resids.type[rowidx] &
        boot.type == plot.config$boot.type[rowidx] &
        dist == plot.config$dist[rowidx] &
        factor == plot.config$factor[rowidx] &
        outlier.rowidx == excl.rowidx]

    labels <- paste(rep("Outlier row:"), 1:ndev)
    names(labels) <- 1:ndev

    p <- ggplot() +
        geom_density(
            data = contaminated,
            aes(reserve, after_stat(scaled), group = excl.rowidx)) +
        geom_density(
            data = uncontaminated,
            aes(reserve, after_stat(scaled)),
            colour = "red") +
        facet_wrap(vars(outlier.rowidx),
            scales = "free_x",
            labeller = as_labeller(labels)
          )  +
        labs(
            title = "Reserve distributions for different outlier points",
            subtitle = sprintf("Perturbation factor: %.2f", factor),
            x = "Reserve",
            y = "Density")

    suppressMessages(
        ggsave(
            sprintf(
                "results/graphs/origin_outlier/densities_%s_%s_%s_factor_%.2f.svg",
                plot.config$resids.type[rowidx],
                plot.config$boot.type[rowidx],
                plot.config$dist[rowidx],
                plot.config$factor[rowidx]
            ),
            plot = p
        )
    )

    #mean plots
    # contaminated.mean <- results[
    #         resids.type == plot.config$resids.type[rowidx] &
    #         boot.type == plot.config$boot.type[rowidx] &
    #         dist == plot.config$dist[rowidx] &
    #         factor == plot.config$factor[rowidx] &
    #         outlier.colidx != excl.colidx &
    #         outlier.rowidx != excl.rowidx,
    #         .(reserve.mean = mean(reserve)),
    #         by = .(outlier.rowidx, outlier.colidx, excl.rowidx, excl.colidx)]

    # uncontaminated.mean <- results[
    #     resids.type == plot.config$resids.type[rowidx] &
    #     boot.type == plot.config$boot.type[rowidx] &
    #     dist == plot.config$dist[rowidx] &
    #     factor == plot.config$factor[rowidx] &
    #     outlier.colidx == excl.colidx &
    #     outlier.rowidx == excl.rowidx,
    #     .(reserve.mean = mean(reserve)),
    #     by = .(outlier.rowidx, outlier.colidx, excl.rowidx, excl.colidx)]

    # p <- meanPlot(contaminated.mean, uncontaminated.mean)

    # suppressMessages(
    #     ggsave(
    #         sprintf(
    #             "results/graphs/single_outlier/means_%s_%s_%s_factor_%.2f.svg",
    #             plot.config$resids.type[rowidx],
    #             plot.config$boot.type[rowidx],
    #             plot.config$dist[rowidx],
    #             plot.config$factor[rowidx]
    #         ),
    #     plot = p
    #     )
    # )

}

close(progress.bar)