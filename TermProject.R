library(regtools)
library(mlbench)
library(geosphere)
library("pracma")
library(polyreg)
library(keras)

#PART A
tripSampleSize <- 50000
hist(rgamma(tripSampleSize, 3, rate=.004),  main = "Gamma Histogram")
porto.data <- read.csv("C://Users/natha/Downloads/train.csv/train.csv") #make sure that this is right for your computer
triptimes <- c()
for(i in sample(1:length(porto.data$POLYLINE), tripSampleSize, replace=TRUE)) {
  triptime <- as.numeric(lapply(strsplit(porto.data$POLYLINE[i], "\\],\\["), length)) 
  #print((triptime * 15)/60) 
  if(triptime * 15 < 4000) {
    triptimes <- append(triptimes, triptime * 15)
  }
}
hist(triptimes, main = "Histogram of Triptimes")
#PART B
#FUNCTIONS:
truncated <- porto.data[sample(nrow(porto.data), 50000),]
for(i in 1:length(truncated$POLYLINE)){
  truncated$TRIPDURATION[i] <- length(unlist(strsplit(truncated$POLYLINE[i], "\\],\\["))) * 15
  truncated$TRIPDISTANCE[i] <- finddistance(truncated$POLYLINE[i])
}
#Now lets clean the data
cleaned <- truncated[complete.cases(truncated$TRIPDISTANCE),] # we don't want na values, or values that are zero
cleaned <- cleaned[cleaned$TRIPDISTANCE!=0,] #it means that a trip didn't happen, or only one pair which is 15 seconds. Most likely a mistake

finddifference<-function(d){ #function finds the time between trips for all trips
  workingtime <- vector()
  for(i in 2:length(d$TRIP_ID)){ # for every pair
    time <- length(unlist(strsplit(d$POLYLINE[i-1], "\\],\\["))) * 15 # take the amount of pairs and multiply by 15
    workingtime[i-1] = d$TIMESTAMP[i] - time - d$TIMESTAMP[i-1] # add to previous time stamp, and subtract from current to find break
  }
  return(workingtime)
}

findworkingtime <-function(d){
  timedriving <- 0
  timeworking <- 0
  index <- 1
  breakoff_time <- 2* 3600
  #If the period between the end of the last trip and the beginning of the next
  #is greater than 2 hours, we can assume that the taxi driver has stopped working
  f<- list()
  for(i in 2:length(d$TIMESTAMP)){
    time <- length(unlist(strsplit(d$POLYLINE[i], "\\],\\["))) * 15 
    timedriving <- timedriving + time 
    # Driver is currently working, he has actively driven for time amount
    # So we add it to timedriving
    if((d$TIMESTAMP[i] - time - d$TIMESTAMP[i-1]) > breakoff_time){ 
      #The driver has not been driving for at least 2 hours.
      #The total time period is the previous timestamp + time just driven
      #subtracted from the beginning of the working period.
      f[length(f) + 1] <- timedriving/(time + d$TIMESTAMP[i] - d$TIMESTAMP[index])# the index is where our period starts
      index <- i #new working period starts at the current timestamp
      timedriving <-0 #reset timedriving
    }
  }
  return (f)
}
# The amount of breaks has few outliers. Lets see if we can explain it.
#In Portugal, it is common for breaks to be 2 hours.

retvec <-vector()
for(i in 1:length(unique(porto.data$TAXI_ID))){ # go through each taxi, and find the amount of working time
  newdata <- porto.data[porto.data$TAXI_ID == unique(porto.data$TAXI_ID)[i],]
  retvec[i] <- sum(finddifference(newdata) > 7200) # find number of breaks for each taxi that are above 7200
}
mean(retvec, na.rm = TRUE)
hist(retvec, main = "Histogram for Number of Breaks Over 2 Hours")
# The mean number of breaks for all taxi drivers over 2 hours
#is 656.9056, nearly twice a day. This makes sense because the drivers will
#take one break, and one longer break once they actually stop working. (go home)
# to find number of trips on average to validate above assumption
numtrips <- vector()
for(i in 1:length(unique(cleaned$TAXI_ID))){
  newid <- unique(cleaned$TAXI_ID)[i]
  newdata <- cleaned[mydata$TAXI_ID == unique(cleaned$TAXI_ID)[i],]
  numtrips[i] <- length(newdata$TRIP_ID)
}
mean(numtrips)
# The mean number of trips was 3818. This is reasonable. Most of the drivers
#are driving around 10 times a day.

times <- vector()
#we are using dataset cleaned which is the cleaned version of mydata. The cleaning was done in p4.r.
for(i in 1:length(unique(cleaned$TAXI_ID))){
  newdata <- cleaned[cleaned$TAXI_ID == unique(cleaned$TAXI_ID)[i],]
  c<-findworkingtime(newdata)
  times[(length(times)+1):(length(times) + length(c))] <- unlist(c)
}
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}
# https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance

betaparams<- estBetaParams(mean(times[times<1], rm.na = TRUE), (var(times[times<1])))
gamma1<- mean(times[times<1], rm.na = TRUE)^2/var(times[times<1])
gamma2<- mean(times[times<1], rm.na = TRUE)/var(times[times<1])
hist(times, main = "Gamma fit of Time Busy",freq = FALSE)
curve(dgamma(x,betaparams$alpha,betaparams$beta),0,1, add= TRUE)
hist(times, main = "Beta fit of Time Busy", freq = FALSE)
curve(dbeta(x,gamma1,gamma2),0,1, add= TRUE)
hist(times, main = "Exponential fit of Time Busy",freq = FALSE)
curve(dexp(x,1/mean(times[times<1], rm.na = TRUE)), 0,1,add = TRUE)
#The exponential distribution seems to be the best fit.
mean(times, rm.na = TRUE)
mean(times[times<1], rm.na = TRUE)

#PART C
expectedValue <- function(values) {
  return(sum(values) / length(values))
}

variance <- function(values) {
  return(expectedValue(values ^ 2) - expectedValue(values) ^ 2)
}

coefficientOfVariation <- function(values) {
  return(sqrt(variance(values)) / expectedValue(values))
}

standardDeviation <- function(values) {
  return(sqrt(variance(values)/length(values)))
}

confidenceInterval95 <- function(values) {
  marginOfError <- 1.96 * standardDeviation(values) / sqrt(length(values))
  conInt <- c(expectedValue(values) - marginOfError, expectedValue(values) + marginOfError)
  return(conInt)
}

printStats <- function(tripTimesA, tripTimesB, tripTimesC) {
  expectedValueA <- expectedValue(tripTimesA)
  expectedValueB <- expectedValue(tripTimesB)
  expectedValueC <- expectedValue(tripTimesC)
  
  confidenceIntervalA <- confidenceInterval95(tripTimesA)
  confidenceIntervalB <- confidenceInterval95(tripTimesB)
  confidenceIntervalC <- confidenceInterval95(tripTimesC)
  
  varianceA <- variance(tripTimesA)
  varianceB <- variance(tripTimesB)
  varianceC <- variance(tripTimesC)
  
  coefficientOfVariationA <- coefficientOfVariation(tripTimesA)
  coefficientOfVariationB <- coefficientOfVariation(tripTimesB)
  coefficientOfVariationC <- coefficientOfVariation(tripTimesC)
  
  print("Trips requested through central dispatch")
  print("Expected value")
  print(expectedValueA)
  print("Confidence interval")
  print(confidenceIntervalA)
  print("Variance")
  print(varianceA)
  print("Cooefficient of variation")
  print(coefficientOfVariationA)
  
  print("Trips requested from a stand")
  print("Expected value")
  print(expectedValueB)
  print("Confidence interval")
  print(confidenceIntervalB)
  print("Variance")
  print(varianceB)
  print("Cooefficient of variation")
  print(coefficientOfVariationB)
  
  print("Trips requested at random, such as on a street corner")
  print("Expected value")
  print(expectedValueC)
  print("Confidence Interval")
  print(confidenceIntervalC)
  print("Variance")
  print(varianceC)
  print("Cooefficient of variation")
  print(coefficientOfVariationC)
}

tripSampleSize <- 50000
tripTimesA <- c() # Dispatch from central
tripTimesB <- c() # Trip requested from stand
tripTimesC <- c() # Trip requested randomly, such as on a street
for(i in sample(1:length(porto.data$POLYLINE), tripSampleSize, replace=TRUE)) {
  beacons <- as.numeric(lapply(strsplit(porto.data$POLYLINE[i], "\\],\\["), length)) 
  realTripTime <- beacons * 15
  callType <- porto.data$CALL_TYPE[i]
  
  if(strcmp("A", callType)) {
    tripTimesA <- append(tripTimesA, realTripTime)
  } else if(strcmp("B", callType)) {
    tripTimesB <- append(tripTimesB, realTripTime)
  } else if(strcmp("C", callType)) {
    tripTimesC <- append(tripTimesC, realTripTime)
  }
}

printStats(tripTimesA, tripTimesB, tripTimesC)

hist(tripTimesA)
hist(tripTimesB)
hist(tripTimesC)

#PART D

finddistance<- function(string){ #this function just needs the polyline, we will then strip away the first two brackets, and last two brackets on the first
  #and last element respectively.
  coordinates <-  unlist(strsplit(string, "\\],\\["))
  firstcoord <- coordinates[1]
  firstcoord <- unlist(strsplit(firstcoord, "\\[\\["))[2]
  pair1 <- vector()
  pair1[1]<- as.numeric(unlist(strsplit(firstcoord,","))[1])
  pair1[2]<- as.numeric(unlist(strsplit(firstcoord,","))[2])
  lastcoord <- coordinates[length(coordinates)]
  lastcoord <- unlist(strsplit(lastcoord, "\\]\\]"))[1]
  pair2 <- vector()
  pair2[1]<- as.numeric(unlist(strsplit(lastcoord,","))[1])
  pair2[2]<- as.numeric(unlist(strsplit(lastcoord,","))[2])
  return (distGeo(pair1,pair2))
  
}

changedays <- function(data, times){
  sorted <- data[order(data$TIMESTAMP),] # need to sort by timestamp for the below logic to work
  holidays <- sort(times)
  timeinterval <- 3600 *24
  for(i in 1:length(sorted$TIMESTAMP)){
    while(sorted$TIMESTAMP[i]> (holidays[1] + timeinterval)){ # we have already sorted the times, so if we are greater than the current holiday, we dont' need to worry about
      if(length(holidays) == 1){ # it anymore
        return (sorted)
      }
      holidays <- holidays[2:length(holidays)]  
    }
    if((sorted$TIMESTAMP[i] >holidays[1]- timeinterval) &(sorted$TIMESTAMP[i] <holidays[1])){ # if it is 24 hours before holiday
      sorted$DAY_TYPE[i] <- "C"
    }
    if(sorted$TIMESTAMP[i] > holidays[1]){ # we know that if it is greater than the holiday, it is not by more than 24 hours
      sorted$DAY_TYPE[i] <-"B"
    }
  }
  return (sorted)
}
# TASK A and B have been taken care of in the second part of the project.
#TASK C
# We need to fix the data as the day types are incorrect
# Use an online unix time to date converter to find holidays and change the days  
#Below are the unix timestamps for Portugal's holidays from 2013 - 2014
holidays <- vector()
holidays[1] <-1401778800
holidays[2] <-1402383600
holidays[3] <- 1402642800
holidays[4] <- 1376550000
holidays[5] <- 1380956400
holidays[6] <- 1383289200
holidays[7] <- 1385884800
holidays[8] <- 1386489600
holidays[9] <- 1387958400
holidays[10] <- 1388563200
holidays[11] <- 1392364800
holidays[12]<-1392451200
holidays[13] <-1396422000
holidays[14] <-1396594800
holidays[15] <-1396681200
holidays[16]<-1398409200
holidays[17] <-1398927600
holidays[18]<-1399014000
holidays[19] <-1399964400

cleaned<- changedays(cleaned, holidays)

#TASK D
# Test linear model
#Try to see Trip Duration from distance, daytype, and call type
callandday <- cleaned[, c(2,7,10,11)]
lin <- qeLin(callandday ,'TRIPDURATION')
print("Mean Absolute Prediction Error from distance, daytype, and call type.")
print(lin$testAcc)
# from distance and daytype
callandday <- cleaned[, c(7,10,11)]
lin <- qeLin(callandday ,'TRIPDURATION')
print("Mean Absolute Prediction Error from distance and daytype.")
print(lin$testAcc)
#from distance and call type
callandday <- cleaned[, c(2,10,11)]
lin <- qeLin(callandday ,'TRIPDURATION')
print("Mean Absolute Prediction Error from distance and call type.")
print(lin$testAcc)
#from calltype and daytype
callandday <- cleaned[, c(2,7,10)]
lin <- qeLin(callandday ,'TRIPDURATION')
print("Mean Absolute Prediction Error from daytype, and call type.")
print(lin$testAcc)
#from distance
callandday <- cleaned[, c(10,11)]
lin <- qeLin(callandday ,'TRIPDURATION')
print("Mean Absolute Prediction Error from distance")
print(lin$testAcc)
#from calltype
callandday <- cleaned[, c(2,10)]
lin <- qeLin(callandday ,'TRIPDURATION')
print("Mean Absolute Prediction Error from call type.")
print(lin$testAcc)
#from daytype
callandday <- cleaned[, c(7,10)]
lin <- qeLin(callandday ,'TRIPDURATION')
print("Mean Absolute Prediction Error from daytype.")
print(lin$testAcc)
#Trying polylin on the distance
callandday <- cleaned[, c(10,11)]
lin <- qePolyLin(callandday ,'TRIPDURATION')
print("Mean Absolute Prediction Error from distance using polylin.")
print(lin$testAcc)
#it seems using all three variables as predictors gives us a reasonable error so lets use that as our model
#set arbitrary u calltype = c, daytype = a, distance = 2000
#having a strict bearing on distance being 2000 is quite useless
#so let's get the values 
#create a new set conditional that will have the above qualities
conditional <- cleaned[cleaned$CALL_TYPE == "C",]
conditional <- conditional[conditional$DAY_TYPE == "A",]
conditional <- conditional[(conditional$TRIPDISTANCE >=(.99*2000)) &(conditional$TRIPDISTANCE <= (1.01*2000)),]
mean(conditional$TRIPDURATION)
#mean trip duration was 649.5056
s2<- var(conditional$TRIPDURATION)
tail <- -qnorm(.05/2)
plusminus <- tail * sqrt(s2/length(conditional$TRIPDURATION))
confidenceinterval <- c(avgtrip - plusminus, avgtrip + plusminus)
print(confidenceinterval)
#95 % chance that values are in between (540, 703)
#if we had a larger sample size we could have much tighter bounds
#but this is simply a tradeoff between time and the tightness of our CI
#Trying machine learning methods
#Let's try K-Nearest Neighbors
#First we need to subfactor DAY_TYPE and CALL_TYPE
callandday <- cleaned[, c(2,5,7,10,11)]
callandday$DAY_TYPE <- toSubFactor(callandday$DAY_TYPE, c("A", "B"))
callandday$CALL_TYPE <- toSubFactor(callandday$CALL_TYPE, c("A", "B"))
lin <- qeKNN(callandday ,'TRIPDURATION')
print(lin$testAcc)
#Testing Error was 320
#Try only using call type and distance
callanddistance <- callandday[, c(1,4,5)]
lin <- qeKNN(callanddistance, "TRIPDURATION")
print(lin$testAcc)
#Error at 275
#Try only using distance
distanceonly <- callandday[, c(4,5)]
lin <- qeKNN(distanceonly, "TRIPDURATION")
print(lin$testAcc)
# Error at 277
#Lets take the most accurate KNN and adjust hyperparameters.
lin <- qeKNN(callanddistance, "TRIPDURATION", k=15)
print("Mean absolute prediction error for k =  15:")
print(lin$testAcc)
#Error- 288
lin <- qeKNN(callanddistance, "TRIPDURATION", k=50)
print("Mean absolute prediction error for k =  50:")
print(lin$testAcc)
#Error- 287
lin <- qeKNN(callanddistance, "TRIPDURATION", k=100)
print("Mean absolute prediction error for k =  100:")
print(lin$testAcc)
#Error- 300.5
lin <- qeKNN(callanddistance, "TRIPDURATION", k=30)
print("Mean absolute prediction error for k =  30:")
print(lin$testAcc)
#Error- 266
lin <- qeKNN(callanddistance, "TRIPDURATION", k=28)
print("Mean absolute prediction error for k =  28:")
print(lin$testAcc)
#Error- 263! Lowest one yet
#Try using Neural Networks

lin <- qeNeural(callandday ,'TRIPDURATION')
print("Mean absolute prediction error using all predictor variables")
print(lin$testAcc)
#When using all predictors, we have an error of 400! Overfitting?
lin <- qeNeural(callanddistance, "TRIPDURATION")
print("Mean absolute prediction error using call type and distance")
print(lin$testAcc)
#Only using calltype and distance, we have an error of 298.
lin <- qeNeural(distanceonly, "TRIPDURATION")
print("Mean absolute prediction error only using distance")
print(lin$testAcc)
#Only using distance we have an error of 287,
qeCompare(cleaned[, c(10,11)], "TRIPDURATION", c('qeLin', 'qePolyLin', 'qeKNN', 'qeNeural'), 25 )