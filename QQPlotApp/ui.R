library(shiny)

# Define UI for application that draws a histogram
shinyUI
(
  fluidPage
  (
    titlePanel(h1(strong("Bursty Processes"))),
    
    
    sidebarLayout
    (
      sidebarPanel
      (
        helpText("Demonstrate Bursty Process Model Fit on Datasets"),
        
        selectInput("dataSet",
                    label = "Choose a dataset to model",
                    choices = c("Earthquakes", "Brain Neurons",
                                "Stock Prices"),
                    selected = "Earthquakes"
                    ),
        sliderInput
        (
          "range",
          label = "Number of Observations above threshold:",
          5,750,200
        )
      ),
      
      mainPanel
      (
        plotOutput("QQplot")
      )
    )

  )
)