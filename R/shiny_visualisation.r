library(shiny)
library(data.table)
library(ggplot2)

single.outlier.results <- readRDS("../results/data_objects/single_outlier.RDS")
calendar.year.results <- readRDS("../results/data_objects/calendar_outlier.RDS")
origin.year.results <- readRDS("../results/data_objects/origin_outlier.RDS")

ndev <- max(single.outlier.results$outlier.colidx)

ui <- navbarPage(
    "Simulation type",
    tabPanel(
        "Single outlier",
        fluidRow(),
        wellPanel(
            fluidRow(
                column(
                    6,
                    selectInput("single.outlier.dist", "Distribution", c("normal", "gamma")),
                    selectInput("single.outlier.boot.type", "Bootstrap type", c("conditional", "unconditional")),
                    selectInput("single.outlier.resids.type", "Residuals type", c("parametric", "raw", "scaled"))
                ),
                fluidRow(
                    column(
                        6,
                        numericInput("single.outlier.factor", "Perturbation factor", 1, min = 0.5, max = 1.5, step = 0.25),
                        numericInput("single.outlier.colidx", "Outlier column", 2, min = 2, max = ndev, step = 1),
                        numericInput("single.outlier.rowidx", "Outlier row", 1, min = 1, max = 1)
                    )
                )
            )
        ),
        hr(),
        plotOutput("single.outlier.plot")
    ),
    tabPanel(
        "Calendar year outlier",
        fluidRow(),
        wellPanel(
            fluidRow(
                column(
                    6,
                    selectInput("calendar.year.dist", "Distribution", c("normal", "gamma")),
                    selectInput("calendar.year.boot.type", "Bootstrap type", c("conditional", "unconditional")),
                    selectInput("calendar.year.resids.type", "Residuals type", c("parametric", "raw", "scaled"))
                ),
                fluidRow(
                    column(
                        6,
                        numericInput("calendar.year.factor", "Perturbation factor", 1, min = 0.5, max = 1.5, step = 0.25),
                        numericInput("calendar.year.outlier.diagidx", "Outlier antidiagonal", 1, min = 1, max = ndev - 1, step = 1)
                    )
                )
            )
        ),
        hr(),
        plotOutput("calendar.year.plot")
    ),
    tabPanel(
        "Origin year outlier",
        fluidRow(),
        wellPanel(
            fluidRow(
                column(
                    6,
                    selectInput("origin.year.dist", "Distribution", c("normal", "gamma")),
                    selectInput("origin.year.boot.type", "Bootstrap type", c("conditional", "unconditional")),
                    selectInput("origin.year.resids.type", "Residuals type", c("parametric", "raw", "scaled"))
                ),
                fluidRow(
                    column(
                        6,
                        numericInput("origin.year.factor", "Perturbation factor", 1, min = 0.5, max = 1.5, step = 0.25),
                        numericInput("origin.year.outlier.rowidx", "Outlier row", 1, min = 1, max = ndev, step = 1),
                    )
                )
            )
        ),
        hr(),
        plotOutput("origin.year.plot")
    )
)

server <- function(input, output, session) {
    observe({
        updateNumericInput(session, "single.outlier.rowidx", max = ndev + 1 - as.numeric(input$single.outlier.colidx))
    })

    observe({
        if (input$single.outlier.resids.type == "parametric") {
            choices <- c("normal")
        } else {
            choices <- c("normal", "gamma")
        }

        updateSelectInput(session, "single.outlier.dist", choices = choices)
    })

    observe({
        if (input$calendar.year.resids.type == "parametric") {
            choices <- c("normal")
        } else {
            choices <- c("normal", "gamma")
        }

        updateSelectInput(session, "calendar.year.dist", choices = choices)
    })

    observe({
            if (input$origin.year.resids.type == "parametric") {
                choices <- c("normal")
            } else {
                choices <- c("normal", "gamma")
            }

            updateSelectInput(session, "origin.year.dist", choices = choices)
        })


    single.outlier.contaminated <- reactive({
        single.outlier.results[
            outlier.rowidx == input$single.outlier.rowidx &
                outlier.colidx == input$single.outlier.colidx &
                resids.type == input$single.outlier.resids.type &
                boot.type == input$single.outlier.boot.type &
                dist == input$single.outlier.dist &
                factor == as.numeric(input$single.outlier.factor) &
                excl.colidx != outlier.colidx &
                excl.rowidx != outlier.rowidx
        ]
    })

    single.outlier.uncontaminated <- reactive({
        single.outlier.results[
            outlier.rowidx == input$single.outlier.rowidx &
                outlier.colidx == input$single.outlier.colidx &
                resids.type == input$single.outlier.resids.type &
                boot.type == input$single.outlier.boot.type &
                dist == input$single.outlier.dist &
                factor == as.numeric(input$single.outlier.factor) &
                excl.colidx == outlier.colidx &
                excl.rowidx == outlier.rowidx
        ]
    })


    output$single.outlier.plot <- renderPlot({
        p <- ggplot() +
            geom_density(
                data = single.outlier.contaminated(),
                aes(x = reserve, group = interaction(excl.colidx, excl.rowidx))
            ) +
            geom_density(
                data = single.outlier.uncontaminated(),
                aes(x = reserve),
                colour = "red"
            )
        print(p)
    })

    calendar.year.contaminated <- reactive({
        calendar.year.results[
            outlier.diagidx == input$calendar.year.outlier.diagidx &
                resids.type == input$calendar.year.resids.type &
                boot.type == input$calendar.year.boot.type &
                dist == input$calendar.year.dist &
                factor == as.numeric(input$calendar.year.factor) &
                excl.diagidx != outlier.diagidx
        ]
    })

    calendar.year.uncontaminated <- reactive({
        calendar.year.results[
            outlier.diagidx == input$calendar.year.outlier.diagidx &
                resids.type == input$calendar.year.resids.type &
                boot.type == input$calendar.year.boot.type &
                dist == input$calendar.year.dist &
                factor == as.numeric(input$calendar.year.factor) &
                excl.diagidx == outlier.diagidx
        ]
    })

    output$calendar.year.plot <- renderPlot({

        p <- ggplot() +
            geom_density(
                data = calendar.year.contaminated(),
                aes(x = reserve, group = excl.diagidx)
            ) +
            geom_density(
                data = calendar.year.uncontaminated(),
                aes(x = reserve),
                colour = "red"
            )

        print(p)
    })

    origin.year.contaminated <- reactive({
        origin.year.results[
            outlier.rowidx == input$origin.year.outlier.rowidx &
                resids.type == input$origin.year.resids.type &
                boot.type == input$origin.year.boot.type &
                dist == input$origin.year.dist &
                factor == as.numeric(input$origin.year.factor) &
                excl.rowidx != outlier.rowidx
        ]
    })

    origin.year.uncontaminated <- reactive({
        origin.year.results[
            outlier.rowidx == input$origin.year.outlier.rowidx &
                resids.type == input$origin.year.resids.type &
                boot.type == input$origin.year.boot.type &
                dist == input$origin.year.dist &
                factor == as.numeric(input$origin.year.factor) &
                excl.rowidx == outlier.rowidx
        ]
    })

    output$origin.year.plot <- renderPlot({
        p <- ggplot() +
            geom_density(
                data = origin.year.contaminated(),
                aes(x = reserve, group = excl.rowidx)
            ) +
            geom_density(
                data = origin.year.uncontaminated(),
                aes(x = reserve),
                colour = "red"
            )
        print(p)
    })
}


shinyApp(ui, server)
