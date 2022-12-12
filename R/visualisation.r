library(shiny)
library(data.table)
library(ggplot2)

results <- readRDS("../results/data_objects/single_outlier.RDS")

resids.type <- c("parametric", "raw", "scaled")
boot.type <- c("conditional", "unconditional")
dist <- c("normal", "gamma")

ui <- fluidPage(
    sidebarLayout(
        sidebarPanel(
        selectInput("dist", "Distribution", dist),
        selectInput("boot.type", "Bootstrap type", boot.type),
        selectInput("resids.type", "Residuals type", resids.type)
        ),
        mainPanel(
            plotOutput("plot")
        )
    )
)


server <- function(input, output, session) {
    output$plot <- renderPlot({

        plot.dt <- as.data.table(results)[
            resids.type == input$resids.type & 
            boot.type == input$boot.type & 
            dist == input$dist &
            factor == factor &
            outlier.colidx != excl.colidx &
            outlier.rowidx != excl.rowidx,
            .(reserve = unlist(reserve)), 
            by = setdiff(names(results), "reserve")]

        point.dt <- as.data.table(results)[
            resids.type == input$resids.type &
            boot.type == input$boot.type &
            dist == input$dist & 
            factor == factor &
            outlier.colidx == excl.colidx &
            outlier.rowidx == excl.rowidx,
            .(reserve = unlist(reserve)), 
            by = setdiff(names(results), "reserve")]

        ggplot() +
            geom_density(
                data = plot.dt,
                aes(x = reserve, group = excl.colidx), alpha = 0.7,
            ) +
            geom_density(
                data = point.dt,
                aes(x = reserve), colour = "red"
            ) +
            facet_wrap(vars(outlier.rowidx, outlier.colidx), scales = "free") 
            # +
            # ggtitle("Reserve distributions for different outlier points",
            #     subtitle = sprintf("Perturbation factor: %f", factor))


})
}

shinyApp(ui, server)