rm(list=ls())
library("ggplot2")
library("gridExtra")
#Author: Eirik Sakariassen

######### 1B

#Integrate (1 + (t/360))^365 * exp(-t) from 0 to infinity
simMonteCarlo <- function(Nsim){
  X <- rexp(Nsim,rate=1)
  e <- (1+ (X/365) )^365
  #sd = (1/sqrt(Nsim))*sqrt(  (1/(Nsim-1))*sum( (e - mean(e))^2 ) )
  sd = ((1/(Nsim-1))*var(e))/sqrt(Nsim)
  #Compute the number of simulations needed for estimate of e=0.01 and confidence of 90%
  num <- ((1.645^2)*(sd^2)*Nsim)/(0.01^2)
  return(list("E" = mean(e),"numsim" = num)) 
}

gcmc <- simMonteCarlo(1000)
firstSim <- gcmc$numsim
cat("The estimated integral with 1000 simulations is: ",gcmc$E)
cat("Run the simulation again with ", firstSim, " simulations")
gcmc <- simMonteCarlo(firstSim)
cat("After  ", firstSim, " simulations, we get an estimate of: ", gcmc$E)


######### 1D
invertedIntensity <- function(t){
  return(sqrt(2*t))
}

#Based on the transformation method
simNHHPbirthdayProblem <- function(Nsim, lambda){
  #Simulate the probability distribution of the birthday problem
  k <- 366
  p <- numeric(k)  #Probabilites for k 1-366
  for (i in 1:k){
    p[i] = 1-(364/365)^( i*(i-1)/2 )
  }
  
  resVecT <- rep(0,Nsim) # Will be used for E(T)
  resVecK <- rep(0,Nsim) # Will be used for E(K)
  for(s in 1:Nsim){
    U <- runif(k)
    #First, simulate HPP with intensity 1 
    arrivalTimes <- cumsum(rexp(5*k*lambda,lambda)) # Simulate the interarrival times T1,T2,..
    arrivalTimes <- arrivalTimes[1:366] #Need maximum 366 arrival times
    #Now calculate s1, s2,...,sN
    arrivalTimes <- invertedIntensity(arrivalTimes)

    #Check if there is two people with birthday the same day after each arrival
    g <- p-U
    resVecT[s] <- arrivalTimes[which(g > 0)[1]] #First time when two have same bday
    resVecK[s] <- which(g > 0)[1] #How many people was in the room
  }

  return(list("T" = resVecT, "K" = resVecK)) 
}

bdp <- simNHHPbirthdayProblem(1000,1) #Maximum 365 arrivals
mean(bdp$K) #E(K)
mean(bdp$T) #E(T)
#As we can see, the estimated value of T and K does not equal each other in the NHPP case

######### 2B

##Generate transition matrices for both cases
Px <- matrix(1/11, nrow=11, ncol=11) #Transition matrix for the complete graph
Py <- matrix(0, nrow=11, ncol=11) #Transition matrix for the discrete circle
Py[1,] <- c(0,1/2,0,0,0,0,0,0,0,0,1/2)

for(i in 2:11){
  newRow <- rep(0,11)
  newRow[(which(Py[i-1,1:11] > 0)+1)%%11] = 1/2
  Py[i,] <- newRow
}

##Generate the given zipfs function
s <- 3
n <- 11
calc <-seq(11)
for (i in 1:11){
  calc[i] <- (1/(i^s))
}
divider <-sum(calc)
pmfZipf <-calc/divider


#This function is insipred by markovLLN in Markov.R
mha_zipfs  <- function(start, Nsim, P, proposal){ #(starting state, #simulations, transition matrix, zipfs pdf)
  nRow <- dim(P)[[1]]
  statesVec <- seq(Nsim+1) #Store states
  statesVec[1] <- start
  for(s in 2:(Nsim+1)){
    i <- statesVec[s-1] #current step
    j <- sample(1:nRow,1,prob=P[statesVec[s-1],]) #proposed next step sampled from the probability matrix based on the last state
    #Based on this sample j in state space S, compute the acceptance function
    acceptance <- (proposal[j]*P[j,i])/(P[i,j]*proposal[i])
    statesVec[s] <- ifelse(runif(1) <= acceptance, j, i) #accept j if the acceptance function is larger than U[0,1]
  }
  sim <-seq(nRow)
  for (k in 1:nRow){
    sim[k]<-sum(statesVec==k)/Nsim
  }
  return(list("sim" = sim, "dtv" = sum((proposal-sim)^2)/2 )) 
}

simZipfX <- mha_zipfs(1,10000,Px,pmfZipf)
simZipfY <- mha_zipfs(1,10000,Py,pmfZipf)


##### Plot 2B #####
dataX <-data.frame(x=c(1:11), probability = simZipfX$sim)
dataY <-data.frame(x=c(1:11), probability = simZipfY$sim)

titleX <- sprintf("Simulation by X with dtv %f", simZipfX$dtv)
titleY <- sprintf("Simulation by Y with dtv %f", simZipfY$dtv)

p1 <-ggplot(dataX, aes(x,probability)) + labs(title=titleX) + geom_bar(stat = "identity") 
p2 <-ggplot(dataY, aes(x,probability)) + labs(title=titleY) + geom_bar(stat = "identity") 

grid.arrange(p1, p2, ncol = 1)
##### end plot #####


######### 2C

Psurfer <- matrix(c(0,1/3,1/3,0,1/3,0,1/2,0,1/2,0,0,0,1/3,1/3,0,0,0,1/3,1/3,1/3,1/3,0,0,0,1/2,0,0,1/2,0,0,1,0,0,0,0,0),
                  nrow=6,ncol=6,byrow=TRUE)


simRandomSurfer <- function(Nsim,P){
  start <- round(runif(1,min=1,max=6))
  nRow <- dim(P)[[1]]
  statesVec <- seq(Nsim+1)
  statesVec[1] <- start
  for(i in 2:(Nsim+1)){
    statesVec[i] <- sample(1:nRow,1,prob=P[statesVec[i-1],])
  }
  #statesVec <- statesVec[-c(1:1000)] We may remove the first 10% for more accurate result
  sim <-seq(nRow)
  for (k in 1:nRow){
    sim[k]<-sum(statesVec==k)/Nsim
  }
  return(sim)
}

simSurf <- simRandomSurfer(10000,Psurfer)
#Rank of the webpages based on the simulation
x <- order(simSurf,decreasing = TRUE)
y <- sort(simSurf, decreasing = TRUE)

barplot(y, names.arg=x, main="Ranking of websites", xlab="Website number", ylab="Probability")

######### 2D


simModifiedRandomSurfer <- function(Nsim,P,prob){
  start <- round(runif(1,min=1,max=6))
  nRow <- dim(P)[[1]]
  statesVec <- seq(Nsim+1)
  statesVec[1] <- start
  for(i in 2:(Nsim+1)){
    if(runif(1) <= prob){ #if heads, proceed as in c)
      statesVec[i] <- sample(1:nRow,1,prob=P[statesVec[i-1],])
    }else{
      statesVec[i] <- round(runif(1,min=1,max=6))
    }
  }
  #statesVec <- statesVec[-c(1:1000)] We may remove the first 10% for more accurate result
  sim <-seq(nRow)
  for (k in 1:nRow){
    sim[k]<-sum(statesVec==k)/Nsim
  }
  return(sim)
}
simModSurf1 <- simModifiedRandomSurfer(10000,Psurfer,0.2)
simModSurf2 <- simModifiedRandomSurfer(10000,Psurfer,0.85)

#Rank both simulations
x1 <- order(simModSurf1,decreasing = TRUE)
y1 <- sort(simModSurf1, decreasing = TRUE)
x2 <- order(simModSurf2,decreasing = TRUE)
y2 <- sort(simModSurf2, decreasing = TRUE)

#Expand the plot in order to see the x-labels and y-labels better
par(mfrow = c(1,2))
barplot(y1, names.arg=x1, main="Ranking of websites p=0.2", xlab="Website number", ylab="Probability")
barplot(y2, names.arg=x2, main="Ranking of websites p=0.85", xlab="Website number", ylab="Probability")

#As we can see, as the probability of choosing head increases, the distribution approaches the case in c).
#This makes a lot of sense since we sample the same way as in c) whenever head is tossed.
#As p approaches 0, visits becomes more random and the probabilities for all 
#websites will approach 1/(number of websites).

######### 3B

simCaseB <- function(Nsim,lambda){
  resVec <- vector(length=10)
  for(i in 1:10){
    theta <- i/10
    profitVec <- rep(0,Nsim)
    for(sim in 1:Nsim){ #simulate for theta 
      arrivalTimes <- cumsum(rexp(5,lambda))
      arrivalTimes <- arrivalTimes[arrivalTimes <= 1] #arrivaltimes in the window [0,1]
      if(length(arrivalTimes) > 0){
        for(t in 1:length(arrivalTimes)){
          bid <- runif(1)
          thresh <- theta*(1- ( (arrivalTimes[t]^2)/2 ) )
          if(bid > thresh){
            profitVec[sim] <- bid #profit in this simulation
          }
        }
      }
    }
    resVec[i] <- mean(profitVec)
  }
  return(resVec)
}


#Estimate the integral in d) based on hit or miss method
simHitorMiss <- function(Nsim,a,b){
  c <- exp(0.5) #Maximum y value of the function within [a,b]
  X <- runif(Nsim,a,b)
  Y <- runif(Nsim,0,c)
  Z <- Y <= exp((X^2)/2)
  integral <- (b*c)*mean(Z)
  return(integral)
}


#Simulate case C by code
simCaseC <- function(Nsim,lambda){
  resVec <- rep(0,Nsim)
  for(i in 1:Nsim){
    arrivalTimes <- cumsum(rexp(5,lambda))
    arrivalTimes <- arrivalTimes[arrivalTimes <= 1]
    if(length(arrivalTimes) > 0){
      for(t in 1:length(arrivalTimes)){
        bid <- runif(1)
        if(bid > (1-arrivalTimes[t])){
          resVec[i] <- bid
        }
      }
    }
  }
  return(mean(resVec))
}
#Compare the hit or miss integral with the simulated function of case C
caseCsim <- simCaseC(100000,1)
caseCintegral <- 1-exp(-0.5)*simHitorMiss(100000,0,1)
cat("Simulated function gives: ", caseCsim, ", and hit or miss estimation gives: ", caseCintegral)

#Estimated profit in case A for theta = i/10, i=1,..,10
simCaseA <- function(){
  resVec <- rep(0,10)
  for(i in 1:10){
    theta <- i/10
    resVec[i] <- ((1+theta)/2)*(1-exp(theta-1))
  }
  return(resVec)
}

caseA <- simCaseA()
caseB <- simCaseB(10000,1)
caseC <- 1-exp(-0.5)*simHitorMiss(100000,0,1)

##### Plot 3E #####
dataA <-data.frame(x=c(1:10)/10, probability = caseA)
dataB <-data.frame(x=c(1:10)/10, probability = caseB)
dataC <-data.frame(x=c(1:10)/10, probability = rep(caseC,10))

p1 <-ggplot(dataA, aes(x,probability)) + labs(title="Case A") + geom_bar(stat = "identity") + labs(x = "price theta")
p2 <-ggplot(dataB, aes(x,probability)) + labs(title="Case B") + geom_bar(stat = "identity") + labs(x = "price theta")
p3 <-ggplot(dataC, aes(x,probability)) + labs(title="Case C") + geom_bar(stat = "identity") + labs(x = "price theta")
grid.arrange(p1, p2, p3, ncol = 1)
##### end plot #####


#Estimates for all the three strategiesin a table.
estimationTable <- merge(dataA,dataB,by="x")
estimationTable <- merge(estimationTable,dataC,by="x")
names(estimationTable) <- c("Theta","Case A", "Case B", "Case C")
estimationTable

###Which strategy gives the highest expected profit?

##If we consider the case where theta varies
strategyMeans <- colMeans(estimationTable[2:4]) #The mean for 10 different theta prices
strategyMeans #Table of mean per strategy
which.max(strategyMeans) #The best strategy (C)
#Strategy C is the best pick if the theta price varies a lot. This is because strategy C
#does not depend on the price (constant profit)


##If we fix theta such that we maximize profit for a strategy
maxValues <- apply(estimationTable[2:4],2,max) #The maximum profit for each strategy
maxValues
which.max(maxValues) #The best one in this case (A)
#Strategy A is the best pick if we fix a price theta. For A, theta=0.2 maximized the profit











