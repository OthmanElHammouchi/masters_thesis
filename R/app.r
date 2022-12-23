library(shiny)
library(data.table)
library(ggplot2)

results <- readRDS("../results/data_objects/single_outlier.RDS")

resids.type <- c("parametric", "raw", "scaled")
boot.type <- c("conditional", "unconditional")
dist <- c("normal", "gamma")
factor <- seq(0.5, 1.5, by = 0.25)

dist.options <- list(normal = 1L, gamma = 2L)
resids.type.options <- list(raw = 1L, scaled = 2L, parametric = 3L)
boot.type.options <- list(conditional = 1L, unconditional = 2L)

ndev <- max(results$outlier.colidx)

ui <- fluidPage(
    fluidRow(),
    wellPanel(
    fluidRow(
        column(3, selectInput("dist", "Distribution", dist)),
        column(3, selectInput("boot.type", "Bootstrap type", boot.type)),
        column(3, selectInput("resids.type", "Residuals type", resids.type)),
        column(3, selectInput("factor", "Perturbation factor", factor))
    ),
    fluidRow(
        column(6, selectInput("outlier.colidx", "Outlier column", 2:ndev, selected = 2)),
        column(6, selectInput("outlier.rowidx", "Outlier row", 1))
    )),
    hr(),
    plotOutput("plot")
)

server <- function(input, output, session) {

     observe({
          updateSelectInput(session, "outlier.rowidx", choices = 1:(ndev + 1 - as.numeric(input$outlier.colidx)))
    })

    other.dt <- reactive({

        results[
            outlier.rowidx == as.numeric(input$outlier.rowidx) &
            outlier.colidx == as.numeric(input$outlier.colidx) &
            resids.type == resids.type.options[[input$resids.type]] &
            boot.type == boot.type.options[[input$boot.type]] &
            dist == dist.options[[input$dist]] &
            factor == as.numeric(input$factor) &
            excl.colidx != outlier.colidx &
            excl.rowidx != outlier.rowidx, ]

    })

    same.dt <- reactive({

        results[
            outlier.rowidx == as.numeric(input$outlier.rowidx) &
            outlier.colidx == as.numeric(input$outlier.colidx) &
            resids.type == resids.type.options[[input$resids.type]] &
            boot.type == boot.type.options[[input$boot.type]] &
            dist == dist.options[[input$dist]] &
            factor == as.numeric(input$factor) &
            excl.colidx == outlier.colidx &
            excl.rowidx == outlier.rowidx, ]

    })

    output$plot <- renderPlot({

        p <- ggplot() +
            geom_density(
                data = other.dt(),
                aes(x = reserve, group = interaction(excl.colidx, excl.rowidx))) +
            geom_density(
                data = same.dt(),
                aes(x = reserve),
                colour = "red")
        print(p)
    })
}

shinyApp(ui, server)