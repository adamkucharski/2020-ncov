library(shiny)
library(rootSolve)

# Define UI for slider demo application
shinyUI(fluidPage(
  
  #  Application title
  titlePanel("Probability of a large nCoV outbreak following introduction of cases"),
  
  p("This analysis uses a model that incorporates randomness and individual-level variation in transmission
    (i.e. potential for 'superspreading') to calculate the probability that a given number of independently
    introduced cases will eventually lead to a large outbreak. (Source: Lloyd-Smith et al, Nature, 2005)"),

  p("Infections with more individual-level variation (such as SARS) lead to more fragile initial transmission chains, 
    and hence are less likely to spark a large outbreak following an introduced case."),
  
  hr(),
  
  # Sidebar with sliders that demonstrate various available
  # options
  sidebarLayout(
    sidebarPanel(
      # Decimal interval with step value

      sliderInput("R0value", "Assumed local reproduction number:", 
                  min = 1, max = 4, value = 2, step= 0.1),
      selectInput("k", 
                  label = "Assumed individual-level variation in transmission",
                  choices = list("SARS-like", "MERS-like",
                                 "Random-mixing"),
                  selected = "SARS-like")

    ),
    
    # Show a table summarizing the values entered
    mainPanel(
      plotOutput('plot')
    )
  ),
  
  hr(),
  
  h4('About this app'),
  
  p("Contributors: Adam Kucharski, Tim Russell, Charlie Diamond, CMMID nCoV working group, Sebastian Funk, Rosalind Eggo."),
  p("Note: this is preliminary analysis and has not yet been peer-reviewed. Parameters are taken from Kucharski & Althaus (Eurosurveillance, 2015)."),
  p("More real-time nCoV analysis by CMMID is ", a(href = 'https://cmmid.github.io/ncov','available here')),
  
  hr()
  
))