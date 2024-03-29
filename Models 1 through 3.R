installed.packages("deSolve", "ggplot2", "gridExtra", "scales", "readr", 
                   "data.table")

library(deSolve)
library(ggplot2)
require(gridExtra)
library(scales)
library(readr)
library(data.table)

###MODEL 1######

# Create the start time, finish time, and time step
START<-0; FINISH<-120; STEP<-1

# Create time vector
simtime <- seq(START, FINISH, by=STEP)

# Create stocks vector, with initial values
stocks  <- c(sS1=100)

# Create auxiliaries vector, with values
auxs    <- c(aKa=0.02, aKd=0.25, aKp=0.05)

# Write callback function (model equations)
model <- function(time, stocks, auxs){
  with(as.list(c(stocks, auxs)),{ 
    
    
    KaR <- sS1 * aKa
    KpR <- sS1 * aKp
    KdR <- sS1 * aKd
    
    d_sS1_dt  <- -(KaR+KpR+KdR)
    
    
    
    return (list(c(d_sS1_dt),
                 fKaR=KaR, fKpR=KpR))   
  })
}

# Call Solver, and store results in a data frame
o<-data.frame(ode(y=stocks, times=simtime, func = model, 
                  parms=auxs, method="euler"))

# Plot output
p1<-ggplot()+
  geom_line(data=o,aes(time,o$sS1,color="1. Susceptible"))+
  scale_y_continuous(labels = comma)+
  ylab("System Stocks")+
  xlab("Day") +
  labs(color="")+
  theme(legend.position="bottom")

p1


###Model 3#### 
#This model will have four stocks but only have 100 kg input on day one for DMI

# Create the start time, finish time, and time step
START<-0; FINISH<-120; STEP<-1

# Create time vector
simtime <- seq(START, FINISH, by=STEP)

# Create stocks vector, with initial values
stocks  <- c(sS1=100,sS2=0, sS3=0, sS4=0)

# Create auxiliaries vector, with values
auxs    <- c(aKa=0.02, aKd=0.25, aKp=0.05)

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

# Plot output
p2<-ggplot()+
  geom_line(data=o,aes(time,o$sS1,color="1. S1"))+
  geom_line(data=o,aes(time,o$sS2,color="2. S2"))+
  geom_line(data=o,aes(time,o$sS3,color="3. S3"))+
  geom_line(data=o,aes(time,o$sS4,color="4. S3"))+
  scale_y_continuous(labels = comma)+
  ylab("System Stocks")+
  xlab("Day") +
  labs(color="")+
  theme(legend.position="bottom")

p2

### Model 3 Precision Data Integration ###

#load the data from gitub repo
dmi = read_csv('https://raw.githubusercontent.com/sdsu-cottonwood-precision-ranch/rumen_kinetics_modeling/main/data/dmi_example.csv')

#Option to export data to compture to perform checks
#write.csv(dmi,"C:/Users/hector.menendez/Desktop/SDSU_2020/MODNUT/dmi.csv")

#View data 
head (dmi)

#call the first cow
dmi_2=dmi[1:24,2] #can subset by [row1:20,column]


#get the starting stock for s1. 
dmi_0=as.numeric(dmi_2[1,1])

#this is the function that runs our daily step
full_model= function(stocks){

START<-0; FINISH<-23; STEP<-1

# Create time vector
simtime <- seq(START, FINISH, by=STEP)

# Create stocks vector, with initial values
# Create auxiliaries vector, with values
auxs  <- c(aKa=0.02, aKd=0.25, aKp=0.05)

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
stocks=c(f_0[24,]$sS1,f_0[24,]$sS2,f_0[24,]$sS3,f_0[24,]$sS4)
f_0$day=1
results[[1]]=f_0

#the following loops over the dmi file and adds in those daily intakes to the model in addition to whatever remains in each stock at the end of the day.
for(i in 2:nrow(dmi_2)) {
  dmi_daily=as.numeric(dmi_2[i,1])
  stocks=stocks+(c(dmi_daily,0,0,0))#add dmi to end of the previous day's stocks
  names(stocks)=c('sS1','sS2','sS3','sS4')
  output=full_model(stocks)
  stocks=c(output[24,]$sS1,output[24,]$sS2,output[24,]$sS3,output[24,]$sS4)
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

p3

