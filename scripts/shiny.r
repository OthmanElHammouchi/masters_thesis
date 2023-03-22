library(shiny)
library(data.table)
library(ggplot2)
library(patternBreak)
suppressPackageStartupMessages(library(ChainLadder))

triangle <- UKMotor
ndev <- nrow(triangle)

load("../results/single.rda")
load("../results/calendar.rda")
load("../results/origin.rda")

single.res <- as.data.table(single.res)
calendar.res <- as.data.table(calendar.res)
origin.res <- as.data.table(origin.res)

ui <- navbarPage(
  "Simulation type",
  tabPanel(
    "Single outlier",
    fluidRow(),
    wellPanel(
      fluidRow(
        column(
          6,
          selectInput("single.boot.type", "Bootstrap type", c("parametric", "residuals", "pairs")),
          selectInput("single.dist", "Process distribution", c("normal", "gamma")),
          selectInput("single.cond", "Conditional", c(TRUE, FALSE)),
          selectInput("single.resids.type", "Residuals type", c("standardised", "modified", "studentised", "log-normal")),
        ),
        fluidRow(
          column(
            6,
            numericInput("single.factor", "Perturbation factor", 1, min = 0.5, max = 1.5, step = 0.25),
            numericInput("single.colidx", "Outlier column", 2, min = 2, max = ndev, step = 1),
            numericInput("single.rowidx", "Outlier row", 1, min = 1, max = 1)
          )
        )
      )
    ),
    hr(),
    plotOutput("single.plot")
  ),
  tabPanel(
    "Calendar year outlier",
    fluidRow(),
    wellPanel(
      fluidRow(
        column(
          6,
          selectInput("calendar.boot.type", "Bootstrap type", c("parametric", "residuals", "pairs")),
          selectInput("calendar.dist", "Process distribution", c("normal", "gamma")),
          selectInput("calendar.cond", "Conditional", c(TRUE, FALSE)),
          selectInput("calendar.resids.type", "Residuals type", c("standardised", "modified", "studentised", "log-normal"))
        ),
        fluidRow(
          column(
            6,
            numericInput("calendar.factor", "Perturbation factor", 1, min = 0.5, max = 1.5, step = 0.25),
            numericInput("calendar.diagidx", "Outlier antidiagonal", 1, min = 1, max = ndev - 1, step = 1)
          )
        )
      )
    ),
    hr(),
    plotOutput("calendar.plot")
  ),
  tabPanel(
    "Origin year outlier",
    fluidRow(),
    wellPanel(
      fluidRow(
        column(
          6,
          selectInput("origin.boot.type", "Bootstrap type", c("parametric", "residuals", "pairs")),
          selectInput("origin.dist", "Process distribution", c("normal", "gamma")),
          selectInput("origin.cond", "Conditional", c(TRUE, FALSE)),
          selectInput("origin.resids.type", "Residuals type", c("standardised", "modified", "studentised", "log-normal"))
        ),
        fluidRow(
          column(
            6,
            numericInput("origin.factor", "Perturbation factor", 1, min = 0.5, max = 1.5, step = 0.25),
            numericInput("origin.rowidx", "Outlier row", 1, min = 1, max = ndev, step = 1),
          )
        )
      )
    ),
    hr(),
    plotOutput("origin.plot")
  )
)
server <- function(input, output, session) {

  observeEvent(input$single.colidx, {
    updateNumericInput(inputId = "single.rowidx",
      max = ndev + 1 - as.numeric(input$single.colidx))
  })

  observeEvent(input$single.boot.type, {
    if (input$single.boot.type == "parametric") {
      updateSelectInput(inputId = "single.resids.type",
        choices = NA
      )
      updateSelectInput(inputId = "single.cond",
        choices = c(TRUE, FALSE))

    } else if (input$single.boot.type == "pairs") {
      updateSelectInput(inputId = "single.resids.type",
        choices = NA)
      updateSelectInput(inputId = "single.cond",
        choices = NA)
        
    } else {
        updateSelectInput(inputId = "single.resids.type",
        choices = NA)
      updateSelectInput(inputId = "single.cond",
        choices = c(TRUE, FALSE))
    }
  })

  single.contaminated <- reactive({
    single.res[
      outlier.rowidx == input$single.rowidx &
      outlier.colidx == input$single.colidx &
      factor == as.numeric(input$single.factor) &
      boot.type == input$single.boot.type &
      proc.dist == input$single.dist &
      cond == as.logical(input$single.cond) &
      resids.type == input$single.resids.type &
      excl.colidx != outlier.colidx &
      excl.rowidx != outlier.rowidx
    ]
  })

  single.uncontaminated <- reactive({
    single.res[
      outlier.rowidx == input$single.rowidx &
      outlier.colidx == input$single.colidx &
      factor == as.numeric(input$single.factor) &
      boot.type == input$single.boot.type &
      proc.dist == input$single.dist &
      cond == as.logical(input$single.cond) &
      resids.type == input$single.resids.type &
      excl.colidx == outlier.colidx &
      excl.rowidx == outlier.rowidx
    ]
  })

  output$single.plot <- renderPlot({
    p <- ggplot() +
      geom_density(
        data = single.contaminated(),
        aes(x = reserve, group = interaction(excl.colidx, excl.rowidx))
      ) +
      geom_density(
        data = single.uncontaminated(),
        aes(x = reserve),
        colour = "red"
      )
    print(p)
  })

  observeEvent(input$calendar.colidx, {
    updateNumericInput(inputId = "calendar.rowidx",
      max = ndev + 1 - as.numeric(input$calendar.colidx))
  })

  observeEvent(input$calendar.boot.type, {
    if (input$calendar.boot.type == "parametric") {
      updateSelectInput(inputId = "calendar.resids.type",
        choices = NA
      )
      updateSelectInput(inputId = "calendar.cond",
        choices = c(TRUE, FALSE))

    } else if (input$calendar.boot.type == "pairs") {
      updateSelectInput(inputId = "calendar.resids.type",
        choices = NA)
      updateSelectInput(inputId = "calendar.cond",
        choices = NA)

    } else {
        updateSelectInput(inputId = "calendar.resids.type",
        choices = c("standardised", "modified", "studentised", "log-normal"))
      updateSelectInput(inputId = "calendar.cond",
        choices = c(TRUE, FALSE))
    }
  })

  calendar.contaminated <- reactive({
    calendar.res[
      outlier.diagidx == input$calendar.diagidx &
      factor == as.numeric(input$calendar.factor) &
      excl.diagidx != outlier.diagidx &
      boot.type == input$calendar.boot.type &
      proc.dist == input$calendar.dist &
      cond == as.logical(input$calendar.cond) &
      resids.type == input$calendar.resids.type
    ]
  })

  calendar.uncontaminated <- reactive({
    calendar.res[
      outlier.diagidx == input$calendar.diagidx &
      factor == input$calendar.factor &
      excl.diagidx == outlier.diagidx &
      resids.type == input$calendar.resids.type &
      boot.type == input$calendar.boot.type &
      cond == as.logical(input$calendar.cond) &
      proc.dist == input$calendar.dist
    ]
  })

  output$calendar.plot <- renderPlot({
    p <- ggplot() +
      geom_density(
        data = calendar.contaminated(),
        aes(x = reserve, group = excl.diagidx)
      ) +
      geom_density(
        data = calendar.uncontaminated(),
        aes(x = reserve),
        colour = "red"
      )

    print(p)
  })

  observeEvent(input$origin.colidx, {
    updateNumericInput(inputId = "origin.rowidx",
      max = ndev + 1 - as.numeric(input$origin.colidx))
  })

  observeEvent(input$origin.boot.type, {
    if (input$origin.boot.type == "parametric") {
      updateSelectInput(inputId = "origin.resids.type",
        choices = NA
      )
      updateSelectInput(inputId = "origin.cond",
        choices = c(TRUE, FALSE))

    } else if (input$origin.boot.type == "pairs") {
      updateSelectInput(inputId = "origin.resids.type",
        choices = NA)
      updateSelectInput(inputId = "origin.cond",
        choices = NA)

    } else {
        updateSelectInput(inputId = "origin.resids.type",
        choices = c("standardised", "modified", "studentised", "log-normal"))
      updateSelectInput(inputId = "origin.cond",
        choices = c(TRUE, FALSE))
    }
  })

  origin.contaminated <- reactive({
    origin.res[
      outlier.rowidx == input$origin.rowidx &
      factor == as.numeric(input$origin.factor) &
      excl.rowidx != outlier.rowidx &
      boot.type == input$origin.boot.type &
      proc.dist == input$origin.dist &
      cond == as.logical(input$origin.cond) &
      resids.type == input$origin.resids.type
    ]
  })

  origin.uncontaminated <- reactive({
    origin.res[
      outlier.rowidx == input$origin.rowidx &
      factor == as.numeric(input$origin.factor) &
      excl.rowidx == outlier.rowidx &
      boot.type == input$origin.boot.type &
      proc.dist == input$origin.dist &
      cond == as.logical(input$origin.cond) &
      resids.type == input$origin.resids.type
    ]
  })

  output$origin.plot <- renderPlot({
    p <- ggplot() +
      geom_density(
        data = origin.contaminated(),
        aes(x = reserve, group = excl.rowidx)
      ) +
      geom_density(
        data = origin.uncontaminated(),
        aes(x = reserve),
        colour = "red"
      )
    print(p)
  })
}

shinyApp(ui, server)
