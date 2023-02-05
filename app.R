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

dmi = read_csv('data/dmi_example.csv')

shiny_app_function= function(dmi,num_feeds,interval,cow,ka,kd,kp){
  interval=round(24/interval)
  #call the first cow
  dmi_2=dmi[num_feeds[[1]]:num_feeds[[2]],cow] #can subset by [row1:20,column]
  #get the starting stock for s1. 
  dmi_0=as.numeric(dmi_2[1,1])
  
  #this is the function that runs our daily step
  full_model= function(stocks){
    START<-0; FINISH<-interval+1; STEP<-1
    
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
ui <- fluidPage( 
  #tags$head(includeHTML(("google-analytics.html"))),
                  titlePanel(tags$head(tags$link(rel = "icon", type = "image/png", href = "favicon.png"),
                                        tags$title("Rumen Kinetics")),title=div(img(src="SDSU_logo2.png", height = '100px',
                                           width = '@00px',
                                           style = "margin:10px 10px"),"Precision Rumen Kinetics Model")),
                  
                  tabsetPanel(
                    
                    tabPanel("Rumen Kinetics Parameters",
                            
                sidebarLayout(
                  sidebarPanel(
                    fileInput(inputId = "dmi_file", label = "Import  DMI file", multiple = F,placeholder = " No file selected", accept = c(".csv")),
                    sliderInput("feeds",
                                "Time Range Select (from data)",
                                min = 0,
                                max = NROW(dmi),
                                value = c(0,20)),
                    sliderInput("interval",
                                "Number of Feedings Per Day",
                                min = 1,
                                max = 24,
                                value = 1),
                    selectInput("cow", "Select unique animal to display",
                                choices=  names(dmi[,-c(1)]),
                                multiple=F),
                    sliderInput("ka",
                                "Ka",
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
                                         h3('Model Output Plot', align='center'),
                                         plotOutput("distPlot"),
                                         h3('DMI Data', align='center'),
                                         DT::dataTableOutput('data')
                                )
                  )
                ),
                tabPanel("Additional Information",
                         h4("The purpose of this app is to show how linear functions, like the kd, ka, and kp rates of 5, 25, and 15%/h can create non-linear behavior when used in a dynamic-mechanistic model. Notice that there is a loss in the model resolution or granularity by how digesta is handled in the model. After a discrete feeding event causes the feed to enter the rumen (S1), a perfect mixing model is applied to the model that remains mathematical and computationally simple enough for the user to draw insights. "),
                         
                         h3("Input Variables"),
                         tags$ul(
                           tags$li("Browse feature: This feature allows the import of daily data from excel which must be saved as CSV,  contain a TIME column and unique animal ID numbers "),
                         
                           tags$li("Feeding Period (Days): This represents the days that an animal is fed over time that will be represented graphically."),
                          
                           tags$li('Number of Feedings Per Day: This represents the number of feeding events at a subdaily level. For example, "2" would mean that the animal is fed twice a day (#).'),
   
                           tags$li("Ka: This respresents the absorption rate of rumen digesta into the rumen wall (%/hour). "),
                           tags$li("Kd: This represents the degradation rate of rumen digesta in the rumen (%/hour). "),
                           tags$li("Kp: This represents the passage rate of rumen digesta out of the rumen (%/hour). ")
                   
                           
                         ),
                         
                         
                         
                         h4("Contributors"),
                         tags$ul(
                           tags$li("Hector Menendez III, Assistant Professor Department of Animal Science, South Dakota State University,", a("hector.menendez@sdstate.edu", href='mailto:hector.menendez@sdstate.edu'),
                                   tags$li("Reid Hensen, Hensen LLC"),
                           tags$li("Jameson Brennan, Assistant Professor Department of Animal Science, South Dakota State University University,", a("jameson.brennan@sdstate.edu", href='mailto:jameson.brennan@sdstate.edu') 
                                   
                                           
                                   ))),
                         
                         h4("Other Great Information"),
                         # p("For additional information on SDSU Extension programs and publications click the button below:"),
                         p("SDSU Extension is an equal opportunity provider and employer in accordance with the nondiscrimination policies of South Dakota State University, the South Dakota Board of Regents and the United States Department of Agriculture. Learn more at extension.sdstate.edu or click on the button below:"),
                         a(h4("Extension Link", class = "btn btn-default action-button" , 
                              style = "fontweight:600"), target = "_blank",
                           href = 'https://extension.sdstate.edu/'),
                         tags$br(),
                         tags$br(),
                         tags$sup("Model adatped from Tedeschi, L. O., & Fox, D. G. (2020). Ruminant Nutrition System. Xanedu Publishing Incorporated. "),
                         tags$br(),
                         tags$br(),
                         tags$li('Figure 1: Example of data format for importing .csv data.'),
                         img(src="input.PNG", height = '400px',
                             width = '1200px',
                             align="left"),
                         tags$br(),
                         tags$br(),
                         
                         tags$br(),
                         tags$br(),
   
                         
                         
                )))
#)

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
  output$data <- renderDataTable(DT::datatable(DMI(),rownames	= F))
}

# Run the application 
shinyApp(ui = ui, server = server)