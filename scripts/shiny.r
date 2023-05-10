library(shiny)
library(data.table)
library(ggplot2)
library(claimsBoot)

load("../results/single.rda")
load("../results/calendar.rda")
load("../results/origin.rda")

single.res <- as.data.table(single.res)
calendar.res <- as.data.table(calendar.res)
origin.res <- as.data.table(origin.res)

ndev <- max(unique(single.res$outlier.colidx))
factors <- unique(single.res$factor)

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      selectInput("boot.type", "Bootstrap type", c("parametric", "residuals", "pairs")),
      selectInput("dist", "Process distribution", c("normal", "gamma")),
      selectInput("cond", "Conditional", c(TRUE, FALSE)),
      selectInput("resids.type", "Residuals type", c("standardised", "modified", "studentised", "log-normal")),
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Single outlier",
          wellPanel(
            fluidRow(
              selectInput("single.factor", "Perturbation factor", factors),
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
        inputId = "resids.type",
        choices = NA
      )
      updateSelectInput(
        inputId = "cond",
        choices = c(TRUE, FALSE)
      )
    } else if (input$boot.type == "pairs") {
      updateSelectInput(
        inputId = "resids.type",
        choices = NA
      )
      updateSelectInput(
        inputId = "cond",
        choices = NA
      )
    } else {
      updateSelectInput(
        inputId = "resids.type",
        choices = c("standardised", "modified", "studentised", "log-normal")
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
    single.res[
      outlier.rowidx == as.numeric(input$single.rowidx) &
        outlier.colidx == as.numeric(input$single.colidx) &
        factor == as.numeric(input$single.factor) &
        boot.type == input$boot.type &
        proc.dist == input$dist &
        cond == as.logical(input$cond) &
        resids.type == input$resids.type &
        excl.colidx != outlier.colidx &
        excl.rowidx != outlier.rowidx
    ]
  })

  single.uncontaminated <- reactive({
    single.res[
      outlier.rowidx == as.numeric(input$single.rowidx) &
        outlier.colidx == as.numeric(input$single.colidx) &
        factor == as.numeric(input$single.factor) &
        boot.type == input$boot.type &
        proc.dist == input$dist &
        cond == as.logical(input$cond) &
        resids.type == input$resids.type &
        excl.colidx == outlier.colidx &
        excl.rowidx == outlier.rowidx
    ]
  })

  single.lims <- quantile(single.res$reserve, c(0.005, 0.995))

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
      ) +
      xlim(single.lims)

    print(p)
  })

  calendar.contaminated <- reactive({
    calendar.res[
      outlier.diagidx == as.numeric(input$calendar.diagidx) &
        factor == as.numeric(input$calendar.factor) &
        excl.diagidx != outlier.diagidx &
        boot.type == input$boot.type &
        proc.dist == input$dist &
        cond == as.logical(input$cond) &
        resids.type == input$resids.type
    ]
  })

  calendar.uncontaminated <- reactive({
    calendar.res[
      outlier.diagidx == as.numeric(input$calendar.diagidx) &
        factor == as.numeric(input$calendar.factor) &
        excl.diagidx == outlier.diagidx &
        resids.type == input$resids.type &
        boot.type == input$boot.type &
        cond == as.logical(input$cond) &
        proc.dist == input$dist
    ]
  })


  calendar.lims <- quantile(calendar.res$reserve, c(0.005, 0.995))

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
      ) +
      xlim(calendar.lims)

    print(p)
  })

  origin.contaminated <- reactive({
    origin.res[
      outlier.rowidx == as.numeric(input$origin.rowidx) &
        factor == as.numeric(input$origin.factor) &
        excl.rowidx != outlier.rowidx &
        boot.type == input$boot.type &
        proc.dist == input$dist &
        cond == as.logical(input$cond) &
        resids.type == input$resids.type
    ]
  })

  origin.uncontaminated <- reactive({
    origin.res[
      outlier.rowidx == as.numeric(input$origin.rowidx) &
        factor == as.numeric(input$origin.factor) &
        excl.rowidx == outlier.rowidx &
        boot.type == input$boot.type &
        proc.dist == input$dist &
        cond == as.logical(input$cond) &
        resids.type == input$resids.type
    ]
  })

  origin.lims <- quantile(origin.res$reserve, c(0.005, 0.995))

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
      ) +
      xlim(origin.lims)
    print(p)
  })
}

shinyApp(ui, server)
