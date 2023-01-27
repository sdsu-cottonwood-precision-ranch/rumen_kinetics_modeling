library(shiny)
library(deSolve)
library(ggplot2)
require(gridExtra)
library(scales)
library(readr)
library(data.table)
library(purrr)
library(DT)
library(shinythemes)

dmi = read_csv('https://raw.githubusercontent.com/rhensen/modnut/main/dmi_example.csv')

shiny_app_function= function(dmi,num_feeds,interval,cow,ka,kd,kp){
  
  #call the first cow
  dmi_2=dmi[1:num_feeds,cow] #can subset by [row1:20,column]
  
  #get the starting stock for s1. 
  dmi_0=as.numeric(dmi_2[1,1])
  
  #this is the function that runs our daily step
  full_model= function(stocks){
    
    START<-0; FINISH<-(24/interval)+1; STEP<-1
    
    # Create time vector
    simtime <- seq(START, FINISH, by=STEP)
    
    # Create stocks vector, with initial values
    # Create auxiliaries vector, with values
    auxs  <- c(aKa=ka, aKd=kd, aKp=kp)
    
    # Write callback function (model equations)
    model <- function(time, stocks, auxs){
      with(as.list(c(stocks, auxs)),{ 
        
        
        KaR1 <- sS1 * aKa
        KdR1 <- sS1 * aKd
        KpR1 <- sS1 * aKp
        
        KaR2 <- sS2 * aKa
        KdR2 <- sS2 * aKd
        KpR2 <- sS2 * aKp
        
        KaR3 <- sS3 * aKa
        KdR3 <- sS3 * aKd
        KpR3 <- sS3 * aKp
        
        KaR4 <- sS4 * aKa
        KdR4 <- sS4 * aKd
        KpR4 <- sS4 * aKp
        
        d_sS1_dt  <- -(KaR1+KdR1+KpR1)
        
        d_sS2_dt  <- (KdR1)-(KaR2+KdR2+KpR2)
        
        d_sS3_dt  <- (KdR2)-(KaR3+KdR3+KpR3)
        
        d_sS4_dt  <- (KdR3)-(KaR4+KdR4+KpR4)
        
        
        return (list(c(d_sS1_dt, d_sS2_dt, d_sS3_dt, d_sS4_dt),
                     fKaR1=KaR1, fKdR1=KdR1, fKpR1=KpR1, 
                     fKaR2=KaR2, fKdR2=KdR2, fKpR2=KpR2,
                     fKaR3=KaR3, fKdR3=KdR3, fKpR3=KpR3,
                     fKaR4=KaR4, fKdR4=KdR4, fKpR4=KpR4
        ))   
      })
    }
    
    # Call Solver, and store results in a data frame
    o<-data.frame(ode(y=stocks, times=simtime, func = model, 
                      parms=auxs, method="euler"))
    return(o)
  }
  
  
  #grab the initial days data. Set all other stocks to 0 to start the model.
  stocks_initial=c(sS1=dmi_0,sS2=0,sS3=0,sS4=0)
  
  #f_0 grabs the initial day's model run
  f_0=full_model(stocks_initial)
  
  
  results=list()
  
  # grab the last stocks at the end of the day in order to use as the starting point for the next day. We also set day 1 with these results. 
  stocks=c(f_0[interval,]$sS1,f_0[interval,]$sS2,f_0[interval,]$sS3,f_0[interval,]$sS4)
  f_0$day=1
  results[[1]]=f_0
  
  #the following loops over the dmi file and adds in those daily intakes to the model in addition to whatever remains in each stock at the end of the day.
  for(i in 2:nrow(dmi_2)) {
    dmi_daily=as.numeric(dmi_2[i,1])
    stocks=stocks+(c(dmi_daily,0,0,0))#add dmi to end of the previous day's stocks
    names(stocks)=c('sS1','sS2','sS3','sS4')
    output=full_model(stocks)
    stocks=c(output[interval,]$sS1,output[interval,]$sS2,output[interval,]$sS3,output[interval,]$sS4)
    day=i 
    output$day=rep(i,nrow(output))
    results[[i]]=output
  }
  
  #bind the lists for each day together
  o=rbindlist(results)
  
  #add an x axis variable
  o$Time=seq.int(nrow(o))
  
  
  # Plot output
  p3<-ggplot()+
    geom_line(data=o,aes(Time,sS1,color="1. S1"))+
    geom_line(data=o,aes(Time,sS2,color="2. S2"))+
    geom_line(data=o,aes(Time,sS3,color="2. S3"))+
    geom_line(data=o,aes(Time,sS4,color="2. S4"))+
    
    
    scale_y_continuous(labels = comma)+
    ylab("Rumen Content") +
    xlab("time (hours)") +
    labs(color="")+
    theme(legend.position="bottom")
  
  return(p3)
}


# Define UI for application that draws a histogram
ui <- fluidPage(theme = shinytheme("paper"),h1("Precision Rumen Kinetics Model", align = "center"),
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          fileInput(inputId = "dmi_file", label = "Import  DMI file", multiple = F,placeholder = " No file selected", accept = c(".csv")),
          sliderInput("feeds",
                      "Feeding Period (Days)",
                      min = 0,
                      max = NROW(dmi),
                      value = 3),
          sliderInput("interval",
                      "Number of Feedings Per Day",
                      min = 1,
                      max = 24,
                      value = 1),
          selectInput("cow", "Select unique animal to display",
                      choices=  names(dmi[,-c(1)]),
                      multiple=F),
          sliderInput("ka",
                      "aKa",
                      min = 0,
                      max = 1,
                      value =0.02),
          sliderInput("kd",
                      "Kd",
                      min = 0,
                      max = 1,
                      value =0.25),
          sliderInput("kp",
                      "Kp",
                      min = 0,
                      max = 1,
                      value =0.05)),

        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(type = "tabs",
            tabPanel("Instructions",
            h3('Jaimie to add', align='center')),
            tabPanel("Model Results",
           h3('Model Output Plot', align='center'),
           plotOutput("distPlot"),
           h3('DMI Data', align='center'),
           dataTableOutput('data')
            ))
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output,session) {
    
  DMI <- reactive({
    if(is_empty(input$dmi_file)){
      return(dmi)
    }
    else({
      print(input$dmi_file)
      df=read_csv(input$dmi_file$datapath)
      updateSelectInput(session, "cow", "Choose the cow based on their tag number",
                      choices=  names(df[,-c(1)]))
      updateSliderInput(session, "feeds",
                        "Number of Feedings",
                        min = 0,
                        max = NROW(df),
                        value = 3)
      return(df)
    })
  })

    output$distPlot <- renderPlot({
        shiny_app_function(DMI(),input$feeds,input$interval,input$cow,input$ka,input$kd,input$kp)
    })
    output$data <- renderDataTable(DT::datatable(DMI()))
}

# Run the application 
shinyApp(ui = ui, server = server)
