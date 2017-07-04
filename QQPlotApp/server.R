#Code below this line runs when the app is launched

library(shiny)
library(maps)
library(mapproj)

source("helpers.R")

Vipaka <- read.csv("Data/Vipaka.csv")
Vipaka <- Vipaka[-(1:27), ]
parse_iso_8601(Vipaka$time) -> Vipaka$time
n=dim(Vipaka)[1]

shinyServer
(
  #Code below this line runs for every user
  function(input, output) 
  {
    #Code below this line runs everytime an input changes
    output$QQplot <- renderPlot({
    data = switch(input$dataSet,
                  "Earthquakes" = Vipaka,
                  "Brain Neurons" = Vipaka,
                  "Stock Prices" = Vipaka)
    ##TT Times, JJ magnitudes
    TT=as.vector(data$time)
    JJ=data$mag
    
    # set lowest threshold at this many values
    k_max = n
    idxJ <- order(JJ, decreasing = T)[1:k_max]
    # set highest threshold at this many values
    k_min = 5
    
    #set K here
    k=input$range
    x=get_durations(TT,idxJ,k)
    estimates=ml.est(x)
    
    x=sort(x)
    y=qml(ppoints(k-1),estimates$tail, estimates$scale)
    
    qqplot(x,y,xlab="Sample Quantiles",ylab="Theoretical Quantiles")
    })
  }
)