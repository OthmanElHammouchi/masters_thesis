library(shiny)
library(data.table)
library(ggplot2)
library(claimsBoot)
library(arrow)

options(shiny.autoreload = TRUE)
options(shiny.reactlog = TRUE)

single.res <- read_feather("../results/glm_single.feather")
calendar.res <- read_feather("../results/glm_calendar.feather")
origin.res <- read_feather("../results/glm_origin.feather")

ndev <- max(unique(single.res$outlier.colidx))
factors <- unique(single.res$factor)

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      selectInput("boot.type", "Bootstrap type", c("parametric", "semiparametric" = "residuals")),
      selectInput("opt", "Option", c("normal", "gamma", "poisson"))
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Single outlier",
          wellPanel(
            fluidRow(
              selectInput("single.factor", "Perturbation factor", factors),
              selectInput("single.colidx", "Outlier column", 1:ndev),
              selectInput("single.rowidx", "Outlier row", 1)
            ),
            plotOutput("single.plot")
          )
        ),
        tabPanel(
          "Calendar year outlier",
          wellPanel(
            fluidRow(
              selectInput("calendar.factor", "Perturbation factor", factors),
              selectInput("calendar.diagidx", "Outlier antidiagonal", seq_len(ndev - 1))
            ),
            plotOutput("calendar.plot")
          )
        ),
        tabPanel(
          "Origin year outlier",
          wellPanel(
            fluidRow(
              selectInput("origin.factor", "Perturbation factor", factors),
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
        choices = c("normal", "gamma", "poisson")
      )
    } else {
      updateSelectInput(
        inputId = "opt",
        choices = "pearson"
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
    if (input$opt == "pearson") {
      single.res[
        outlier.rowidx == as.numeric(input$single.rowidx) &
          outlier.colidx == as.numeric(input$single.colidx) &
          factor == as.numeric(input$single.factor) &
          boot.type == input$boot.type &
          (excl.colidx != outlier.colidx | excl.rowidx != outlier.rowidx)
      ]
    } else {
      single.res[
        outlier.rowidx == as.numeric(input$single.rowidx) &
          outlier.colidx == as.numeric(input$single.colidx) &
          factor == as.numeric(input$single.factor) &
          boot.type == input$boot.type &
          opt == input$opt &
          (excl.colidx != outlier.colidx | excl.rowidx != outlier.rowidx)
      ]
    }
  })

  single.uncontaminated <- reactive({
    if (input$opt == "pearson") {
      single.res[
        outlier.rowidx == as.numeric(input$single.rowidx) &
          outlier.colidx == as.numeric(input$single.colidx) &
          factor == as.numeric(input$single.factor) &
          boot.type == input$boot.type &
          excl.colidx == outlier.colidx &
          excl.rowidx == outlier.rowidx
      ]
    } else {
      single.res[
        outlier.rowidx == as.numeric(input$single.rowidx) &
          outlier.colidx == as.numeric(input$single.colidx) &
          factor == as.numeric(input$single.factor) &
          boot.type == input$boot.type &
          opt == input$opt &
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
    if (input$opt == "pearson") {
      calendar.res[
        outlier.diagidx == as.numeric(input$calendar.diagidx) &
          factor == as.numeric(input$calendar.factor) &
          excl.diagidx != outlier.diagidx &
          boot.type == input$boot.type
      ]
    } else {
      calendar.res[
        outlier.diagidx == as.numeric(input$calendar.diagidx) &
          factor == as.numeric(input$calendar.factor) &
          excl.diagidx != outlier.diagidx &
          boot.type == input$boot.type &
          opt == input$opt
      ]
    }
  })

  calendar.uncontaminated <- reactive({
    if (input$opt == "pearson") {
      calendar.res[
        outlier.diagidx == as.numeric(input$calendar.diagidx) &
          factor == as.numeric(input$calendar.factor) &
          excl.diagidx == outlier.diagidx &
          boot.type == input$boot.type &
          is.na(opt)
      ]
    } else {
      calendar.res[
        outlier.diagidx == as.numeric(input$calendar.diagidx) &
          factor == as.numeric(input$calendar.factor) &
          excl.diagidx == outlier.diagidx &
          boot.type == input$boot.type &
          opt == input$opt
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
    if (input$opt == "pearson") {
      origin.res[
        outlier.rowidx == as.numeric(input$origin.rowidx) &
          factor == as.numeric(input$origin.factor) &
          excl.rowidx != outlier.rowidx &
          boot.type == input$boot.type &
          is.na(opt)
      ]
    } else {
      origin.res[
        outlier.rowidx == as.numeric(input$origin.rowidx) &
          factor == as.numeric(input$origin.factor) &
          excl.rowidx != outlier.rowidx &
          boot.type == input$boot.type &
          opt == input$opt
      ]
    }
  })

  origin.uncontaminated <- reactive({
    if (input$opt == "pearson") {
      origin.res[
        outlier.rowidx == as.numeric(input$origin.rowidx) &
          factor == as.numeric(input$origin.factor) &
          excl.rowidx == outlier.rowidx &
          boot.type == input$boot.type &
          is.na(opt)
      ]
    } else {
      origin.res[
        outlier.rowidx == as.numeric(input$origin.rowidx) &
          factor == as.numeric(input$origin.factor) &
          excl.rowidx == outlier.rowidx &
          boot.type == input$boot.type &
          opt == input$opt
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
