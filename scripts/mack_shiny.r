library(shiny)
library(data.table)
library(ggplot2)
library(claimsBoot)
library(arrow)

options(shiny.autoreload = TRUE)
options(shiny.reactlog = TRUE)

single.res <- read_feather("../results/mack_single.feather")
calendar.res <- read_feather("../results/mack_calendar.feather")
origin.res <- read_feather("../results/mack_origin.feather")

ndev <- max(unique(single.res$outlier.colidx))
mean.factors <- unique(single.res$mean.factor)
sd.factors <- unique(single.res$sd.factor)

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      selectInput("boot.type", "Bootstrap type",
        c("parametric", "semiparametric" = "residuals", "nonparametric" = "pairs")
      ),
      selectInput("opt", "Option", c("normal", "gamma")),
      selectInput("cond", "Conditional", c(TRUE, FALSE))
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Single outlier",
          wellPanel(
            fluidRow(
              selectInput("single.meanfac", "Mean perturbation factor", mean.factors),
              selectInput("single.sdfac", "Stdev perturbation factor", sd.factors),
              selectInput("single.colidx", "Outlier column", 2:ndev),
              selectInput("single.rowidx", "Outlier row", 1)
            ),
            plotOutput("single.plot")
          )
        ),
        tabPanel(
          "Calendar year outlier",
          wellPanel(
            fluidRow(
              selectInput("calendar.meanfac", "Mean perturbation factor", mean.factors),
              selectInput("calendar.sdfac", "Stdev perturbation factor", sd.factors),
              selectInput("calendar.diagidx", "Outlier antidiagonal", seq_len(ndev - 1))
            ),
            plotOutput("calendar.plot")
          )
        ),
        tabPanel(
          "Origin year outlier",
          wellPanel(
            fluidRow(
              selectInput("origin.meanfac", "Mean perturbation factor", mean.factors),
              selectInput("origin.sdfac", "Stdev perturbation factor", sd.factors),
              selectInput("origin.rowidx", "Outlier row", seq_len(ndev - 1)),
            ),
            plotOutput("origin.plot")
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  observeEvent(input$boot.type, {
    if (input$boot.type == "parametric") {
      updateSelectInput(
        inputId = "opt",
        choices = c("normal", "gamma")
      )
      updateSelectInput(
        inputId = "cond",
        choices = c(TRUE, FALSE)
      )
    } else if (input$boot.type == "pairs") {
      updateSelectInput(
        inputId = "opt",
        choices = "NA"
      )
      updateSelectInput(
        inputId = "cond",
        choices = TRUE
      )
    } else {
      updateSelectInput(
        inputId = "opt",
        choices = c("standardised", "studentised", "log-normal")
      )
      updateSelectInput(
        inputId = "cond",
        choices = c(TRUE, FALSE)
      )
    }
  })

  observeEvent(input$single.colidx, {
    updateSelectInput(
      inputId = "single.rowidx",
      choices = 1:(ndev + 1 - as.numeric(input$single.colidx))
    )
  })

  single.contaminated <- reactive({
    if (input$boot.type != "pairs") {
      single.res[
        outlier.rowidx == as.numeric(input$single.rowidx) &
          outlier.colidx == as.numeric(input$single.colidx) &
          mean.factor == as.numeric(input$single.meanfac) &
          sd.factor == as.numeric(input$single.sdfac) &
          boot.type == input$boot.type &
          opt == input$opt &
          cond == as.logical(input$cond) &
          (excl.colidx != outlier.colidx | excl.rowidx != outlier.rowidx)
      ]
    } else {
      single.res[
        outlier.rowidx == as.numeric(input$single.rowidx) &
          outlier.colidx == as.numeric(input$single.colidx) &
          mean.factor == as.numeric(input$single.meanfac) &
          sd.factor == as.numeric(input$single.sdfac) &
          boot.type == input$boot.type &
          (excl.colidx != outlier.colidx | excl.rowidx != outlier.rowidx)
      ]
    }
  })

  single.uncontaminated <- reactive({
    if (input$boot.type != "pairs") {
      single.res[
        outlier.rowidx == as.numeric(input$single.rowidx) &
          outlier.colidx == as.numeric(input$single.colidx) &
          mean.factor == as.numeric(input$single.meanfac) &
          sd.factor == as.numeric(input$single.sdfac) &
          boot.type == input$boot.type &
          opt == input$opt &
          cond == as.logical(input$cond) &
          excl.colidx == outlier.colidx &
          excl.rowidx == outlier.rowidx
      ]
    } else {
      single.res[
        outlier.rowidx == as.numeric(input$single.rowidx) &
          outlier.colidx == as.numeric(input$single.colidx) &
          mean.factor == as.numeric(input$single.meanfac) &
          sd.factor == as.numeric(input$single.sdfac) &
          boot.type == input$boot.type &
          excl.colidx == outlier.colidx &
          excl.rowidx == outlier.rowidx
      ]
    }
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

  calendar.contaminated <- reactive({
    if (input$boot.type != "pairs") {
      calendar.res[
        outlier.diagidx == as.numeric(input$calendar.diagidx) &
          mean.factor == as.numeric(input$calendar.meanfac) &
          sd.factor == as.numeric(input$calendar.sdfac) &
          excl.diagidx != outlier.diagidx &
          boot.type == input$boot.type &
          opt == input$opt &
          cond == as.logical(input$cond)
      ]
    } else {
      calendar.res[
        outlier.diagidx == as.numeric(input$calendar.diagidx) &
          mean.factor == as.numeric(input$calendar.meanfac) &
          sd.factor == as.numeric(input$calendar.sdfac) &
          excl.diagidx != outlier.diagidx &
          boot.type == input$boot.type
      ]
    }
  })

  calendar.uncontaminated <- reactive({
    if (input$boot.type != "pairs") {
      calendar.res[
        outlier.diagidx == as.numeric(input$calendar.diagidx) &
          mean.factor == as.numeric(input$calendar.meanfac) &
          sd.factor == as.numeric(input$calendar.sdfac) &
          excl.diagidx == outlier.diagidx &
          boot.type == input$boot.type &
          cond == as.logical(input$cond) &
          opt == input$opt
      ]
    } else {
      calendar.res[
        outlier.diagidx == as.numeric(input$calendar.diagidx) &
          mean.factor == as.numeric(input$calendar.meanfac) &
          sd.factor == as.numeric(input$calendar.sdfac) &
          excl.diagidx == outlier.diagidx &
          boot.type == input$boot.type
      ]
    }
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

  origin.contaminated <- reactive({
    if (input$boot.type != "pairs") {
      origin.res[
        outlier.rowidx == as.numeric(input$origin.rowidx) &
          mean.factor == as.numeric(input$origin.meanfac) &
          sd.factor == as.numeric(input$origin.sdfac) &
          excl.rowidx != outlier.rowidx &
          boot.type == input$boot.type &
          opt == input$opt &
          cond == as.logical(input$cond)
      ]
    } else {
      origin.res[
        outlier.rowidx == as.numeric(input$origin.rowidx) &
          mean.factor == as.numeric(input$origin.meanfac) &
          sd.factor == as.numeric(input$origin.sdfac) &
          excl.rowidx != outlier.rowidx &
          boot.type == input$boot.type
      ]
    }
  })

  origin.uncontaminated <- reactive({
    if (input$boot.type != "pairs") {
      origin.res[
        outlier.rowidx == as.numeric(input$origin.rowidx) &
          mean.factor == as.numeric(input$origin.meanfac) &
          sd.factor == as.numeric(input$origin.sdfac) &
          excl.rowidx == outlier.rowidx &
          boot.type == input$boot.type &
          opt == input$opt &
          cond == as.logical(input$cond)
      ]
    } else {
      origin.res[
        outlier.rowidx == as.numeric(input$origin.rowidx) &
          mean.factor == as.numeric(input$origin.meanfac) &
          sd.factor == as.numeric(input$origin.sdfac) &
          excl.rowidx == outlier.rowidx &
          boot.type == input$boot.type
      ]
    }
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
