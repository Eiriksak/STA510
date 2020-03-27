rm(list=ls())
library("ggplot2")
library("gridExtra")
##Author: Eirik Sakariassen
######### 1C

#N = num of cards, k=num of success states among the cards, n=num of cards drawn
simCardDraws <- function(Nsim, N, k, n){
  #Generate a deck of cards
  suits <- c('Red','Failure') # Red=success, black+yellow=failure
  ResVec <- vector(length=Nsim)  # Vector to store the simulated draw results
  for(i in 1:Nsim){
    deck <- rep(suits, times = c(k,N-k)) #creates a deck of k (18) red cards, then the remainding ones becomes failure cards
    drawn_cards <- c() #vector of drawn cards
    #Now pick n random cards from the deck and append to drawn_cards
    for(j in 1:n){
      index <- round(runif(1, 1, length(deck))) #random number between 1 and length of deck
      drawn_cards <- c(drawn_cards, deck[index]) #append to drawn_cards vector
      deck <- deck[-index] #remove the card from the deck
    }
    ResVec[i] = sum(drawn_cards == 'Red') #add the number of red cards drawn into the result vector
  }
  return(ResVec)
}

cardSim <- simCardDraws(Nsim=10000, N=36, k=18, n=6)
hist(cardSim)
sum(cardSim >= 2)/length(cardSim) #P(Y6 >= 2), 0.9113 from 1a
mean(cardSim) #expected value, 3 from 1a


######### 1E

simCardGame <- function(Nsim, condition1){
  resVec <- vector(length = 11)
  for(i in 0:10){ #simulate 10 probabilities for p
    p <- i/10
    simWinner <- vector(length = Nsim)
    for(j in 1:Nsim){
      deck <- rep(c('Red','Black', 'Yellow'), times = c(18,18,1)) #creates a deck of 18 red, 18 black and 1 yellow
      
      player1 = c()
      player2 = c()
      
      for(k in 1:8) { #pick 8 cards
        index <- round(runif(1, 1, length(deck))) #random number between 1 and length of deck
        if(k%%2 == 0){ #the players pick a new card every other time
          player2 <- c(player2, deck[index]) 
        }else{
          if(condition1){ #First three cards are red, the last one random
            if(k==7){
              player1 <- c(player1, deck[index]) #last random
            }else{
              player1 <- c(player1, deck[1]) #First value allways red in the deck since there is no shuffle and k<18
              deck <- deck[-1] 
              next
            }
          }else{
            player1 <- c(player1, deck[index]) 
          }
        }
        deck <- deck[-index] 
      }
      winner = ""
      
      #Coin toss  (ensures a winner if rest of the conditions does not appair)
      toss <- ifelse(runif(1) < p, "Heads", "Tails")
      if(toss == 'Heads'){
        winner = "Player1"
      }else{
        winner = "Player2"
      }
      
      #Player with most red cards 
      if (sum(player1 == 'Red') > sum(player2 == 'Red')){
        winner = "Player1"
      } else if (sum(player2 == 'Red') > sum(player1 == 'Red')){
        winner = "Player2"
      }
      
      #Check if anyone has a yellow card
      if('Yellow' %in% player1){
        winner = "Player1"
      } else if('Yellow' %in% player2){
        winner = "Player2"
      }
      simWinner[j] = winner #Store all winners for this simulation of p
    }
    resVec[i+1] <- (sum(simWinner == 'Player2'))/Nsim #Add the probability that p2 won based on one simulation for p
  }
  return(resVec)
}

CardGameC1 <- simCardGame(10000, TRUE) #with condition 1
CardGameC2 <- simCardGame(10000, FALSE) #with condition 2

##### Plot 1E #####
pValues <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
mydata1 <-data.frame(pValues, CardGameC1)
mydata2 <-data.frame(pValues, CardGameC2)

p1 <-ggplot(mydata1, aes(pValues, CardGameC1)) + 
  labs(x = "p values", y="Probability that Player 2 wins", title="Condition 1") + geom_bar(stat = "identity")
p2 <-ggplot(mydata2, aes(pValues, CardGameC2)) + 
  labs(x = "p values", y="Probability that Player 2 wins", title="Condition 2") + geom_bar(stat = "identity")
grid.arrange(p1, p2, ncol = 2)
##### end plot #####


######### 2C
simTimeToHospital <- function(Nsim){
  resVec <- vector(length = Nsim)
  for(i in 1:Nsim){
    #ACCIDENT JUST HAPPENED
    peopleInAccident <- sample(c(1:4), size = 1, prob=c(0.4,0.3,0.2,0.1), replace = TRUE) 
    timeToLastPerson <- sum(rexp(peopleInAccident,1/15))
    resVec[i] <- timeToLastPerson
  }
  return(resVec)
}

T <- simTimeToHospital(10000)
sum(T > 30)/length(T) #P(T>30)
median(T)
max(T)


######### 2E
simBusyAmbulance <- function(Nsim, stoptime, lambda){
  resVec <- vector(length=Nsim)
  for(i in 1:Nsim){
    #Simulate a poisson process with intensity lambda over a time period
    #The inter arrival times can be simulated from an exponential distribution with parameter lambda
    #This part is directly taken from stochastic_processes_examples.R
    expectednumber <- stoptime*lambda  # Expected number of event until stoptime
    sim <- 3*expectednumber  # Simulate more than the expected number to be certain to exceed stoptime
    timesBetween <- rexp(sim,lambda) # Simulate the interarrival times T1,T2,..
    arrivalTimes <- cumsum(timesBetween)   # Calculate arrival times 
    arrivalTimes <- arrivalTimes[arrivalTimes<stoptime] # Discard the times larger than stoptime
    if(length(arrivalTimes) < 2){ #go to next simulation if it was only one accident this time period
      resVec[i] <- 0
      next
    }
    #Now simulate a queue system for the ambulance based on the arrival times
    isBusy <- FALSE
    for(accident in 1:(length(arrivalTimes)-1)){ #Last arrival has no accident after that may crash
      peopleInAccident <- sample(c(1:4), size = 1, prob=c(0.4,0.3,0.2,0.1), replace = TRUE) 
      timeToLastPerson <- sum(rexp(peopleInAccident,1/15))/60 #represent accident times as hours
      
      if(arrivalTimes[accident+1]-arrivalTimes[accident] < timeToLastPerson){ #Check if ambulance is done before next arrival
        isBusy <- TRUE
        break
      }
    }
    if(isBusy){
      resVec[i] <- 1
    }else{resVec[i] <- 0}
  }
  return(sum(resVec)/Nsim)
}
T24 <- simBusyAmbulance(10000,24,0.2)
T24 #Probability that the ambulance is being out on a mission at the time when an accident occurs within 24h 
TWEEK <- simBusyAmbulance(10000,168,0.2)
TWEEK 
TMONTH <- simBusyAmbulance(10000,720,0.2)
TMONTH 

######### 3A
simRenyiRandomGraph <- function(n, p, Nsim){
  degreeDis <- rep(0,n) #Add distribution after each simulation
  for(s in 1:Nsim){
    degreeVec <- rep(0,n) #store the degree of each vertex
    
    #Represent the graph as a matrix
    adjacencyMatrix <- matrix(1, nrow=n, ncol=n) #Create a nxn adjacency matrix 
    diag(adjacencyMatrix) <- 0 #Set all nodes to be connected to each other, except connected to itself (0 at the diagonal)
    
    #Since this is an undirected graph, we dont need to look at more than the upper half of the diagonal in the matrix
    for(i in 2:n){ #n columns (skip the first since thats 0 anyway)
      for(j in 1:i-1){ #stop at the upper half of the diagonal
        #Now evaluate each edge with a coin toss
        toss <- ifelse(runif(1) < p, "Head", "Tail")
        if(toss =="Tail"){
          adjacencyMatrix[j,i] <- 0 #remove the edge
        }else{
          #Increment both the i and j nodes degree (since its undirected)
          degreeVec[i] <- degreeVec[i] +1
          degreeVec[j] <- degreeVec[j] +1
        }
      }
    }
    #Create a degree distribution for this simulation
    simDegreeDis <- rep(0,n)
    for(i in 1:n){
      simDegreeDis[degreeVec[i]] <- simDegreeDis[degreeVec[i]] +1
    }
    degreeDis <- degreeDis + simDegreeDis/n #Add to the total distribution
  }
  return(list("distribution" = degreeDis/n, "mean" = mean(degreeVec), "max" = max(degreeVec), "var" = var(degreeVec))) 
}

rg1 = simRenyiRandomGraph(100,0.5,100)
rg2 = simRenyiRandomGraph(100,0.1,100)
rg3 = simRenyiRandomGraph(100,0.01,100)

##### Plot 3A #####
mydata1 <-data.frame(x=c(1:100), probability = rg1$distribution)
mydata2 <-data.frame(x=c(1:100), probability = rg2$distribution)
mydata3 <-data.frame(x=c(1:100), probability = rg3$distribution)

title1 <- sprintf("p=0.5, mean=%2.1f, var=%2.1f, max=%2.1f", rg1$mean, rg1$var, rg1$max)
title2 <- sprintf("p=0.1, mean=%2.1f, var=%2.1f, max=%2.1f", rg2$mean, rg2$var, rg2$max)
title3 <- sprintf("p=0.01, mean=%2.1f, var=%2.1f, max=%2.1f", rg3$mean, rg3$var, rg3$max)

p1 <-ggplot(mydata1, aes(x,probability)) + labs(title=title1) + geom_bar(stat = "identity") 
p2 <-ggplot(mydata2, aes(x,probability)) + labs(title=title2) + geom_bar(stat = "identity") 
p3 <-ggplot(mydata3, aes(x,probability)) + labs(title=title3) + geom_bar(stat = "identity") 
grid.arrange(p1, p2, p3, ncol = 1)
##### end plot #####


######### 3B

#the adjacencyMatrix is commented out due to the compute time. Could be used as a representation of the graph if needed
simPreferentialAttachment <- function(n, p, Nsim){
  degreeDis <- rep(0,n) #Add distribution after each simulation
  for(s in 1:Nsim){
    adjacencyMatrix <- matrix(c(0,1,0,0), nrow=2, ncol=2) #Start with 2 vertices, where vertex2 is connected to vertex1
    V <- seq(2) #Collection of vertices
    degreeVec <- c(1,0) #store the degree of each vertex
    for(i in 3:n){ #down the rows
      toss <- ifelse(runif(1) < p, "Head", "Tail")
      #Existing rows in the matrix only needs a new zero column to be inserted
      #adjacencyMatrix <- cbind(adjacencyMatrix, 0) 
      #The new row (representing newly inserted vertex) need to fill in the entire row with i columns
      #adjacencyMatrix <- rbind(adjacencyMatrix, rep(0,i)) 
      if(toss == "Head"){
        #Connect to an existing vertex drawn uniformly at random
        vertex <- sample(V, 1, replace = TRUE) 
        #adjacencyMatrix[i,vertex] <- 1 #At the newly inserted row, add which column (vertex) it is connected to
        degreeVec[vertex] <- degreeVec[vertex] + 1 #New edge at the drawn vertex
      } else{
        #Connect to an existing vertex (the number of edges with end-point in an existing vertex determines its probability)
        probVec <- degreeVec/sum(degreeVec) #calculates the probabilty for each vertex
        vertex <- sample(V, size = 1, prob=probVec, replace = TRUE) #pick a vertex based on each probability
        #adjacencyMatrix[i,vertex] <- 1
        degreeVec[vertex] <- degreeVec[vertex] + 1 #New edge at the drawn vertex
      }
      V <- seq(i) #new vertex in the collection of vertices
      degreeVec <- c(degreeVec,0) #new vertex added to the degree vector (0 value since no one is connected atm)
    }
    simDegreeDis <- rep(0,n)
    for(l in 1:n){
      simDegreeDis[degreeVec[l]+1] <- simDegreeDis[degreeVec[l]+1] + 1
    }
    degreeDis <- degreeDis + simDegreeDis/n #Add to the total distribution
  }
  return(list("distribution" = degreeDis/n, "degrees" = degreeVec, "mean" = mean(degreeVec), "max" = max(degreeVec), "var" = var(degreeVec))) 
}

pam1 <- simPreferentialAttachment(10000, 0.5, 1)
pam2 <- simPreferentialAttachment(10000,0.1,1)
pam3 <- simPreferentialAttachment(10000,0.01,1)

##### Plot 3B part 1##### (if just 1 of 6 plots shows up, click on the zoom button or run plot again)
mydata1 <-data.frame(age=c(1:10000), degree = pam1$degrees)
mydata2 <-data.frame(age=c(1:10000), degree = pam2$degrees)
mydata3 <-data.frame(age=c(1:10000), degree = pam3$degrees)

p1 <- ggplot(mydata1,aes(x=age,y=degree)) + geom_point(alpha = 1) + labs(title="p=0.5")

logp1 <- ggplot(mydata1,aes(x=log(age),y=log(degree))) + 
  labs(title="Log-log plot with regression") + geom_point(alpha = 1) +
  geom_smooth(method = "lm", se = FALSE)

p2 <- ggplot(mydata2,aes(x=age,y=degree)) + geom_point(alpha = 1) + labs(title="p=0.1")

logp2 <- ggplot(mydata2,aes(x=log(age),y=log(degree))) +
  labs(title="Log-log plot with regression") + geom_point(alpha = 1) +
  geom_smooth(method = "lm", se = FALSE)

p3 <- ggplot(mydata3,aes(x=age,y=degree)) + geom_point(alpha = 1) + labs(title="p=0.01")

logp3 <- ggplot(mydata3,aes(x=log(age),y=log(degree))) +
  labs(title="Log-log plot with regression") + geom_point(alpha = 1) +
  geom_smooth(method = "lm", se = FALSE)

grid.arrange(p1, logp1, p2, logp2, p3, logp3, ncol = 2)
##### end plot part 1#####

#3B part 2
#simulate 100 times with N=100
pam4 <- simPreferentialAttachment(100,0.5,100)
pam5 <- simPreferentialAttachment(100,0.1,100)
pam6 <- simPreferentialAttachment(100,0.01,100)
##### Plot 3B part 2 #####
mydata4 <-data.frame(x=c(1:100), probability = pam4$distribution)
mydata5 <-data.frame(x=c(1:100), probability = pam5$distribution)
mydata6 <-data.frame(x=c(1:100), probability = pam6$distribution)

title4 <- sprintf("p=0.5, mean=%2.1f, var=%2.1f, max=%2.1f", pam4$mean, pam4$var, pam4$max)
title5 <- sprintf("p=0.1, mean=%2.1f, var=%2.1f, max=%2.1f", pam5$mean, pam5$var, pam5$max)
title6 <- sprintf("p=0.01, mean=%2.1f, var=%2.1f, max=%2.1f", pam6$mean, pam6$var, pam6$max)

p4 <-ggplot(mydata4, aes(x,probability)) + labs(title=title4) + geom_bar(stat = "identity") 
p5 <-ggplot(mydata5, aes(x,probability)) + labs(title=title5) + geom_bar(stat = "identity") 
p6 <-ggplot(mydata6, aes(x,probability)) + labs(title=title6) + geom_bar(stat = "identity") 
grid.arrange(p4, p5, p6, ncol = 1)
##### end plot 3B part 2 #####


######### 4B
#Implemented the discrete inverse transfer algorithm from 4a

simDiscInverseTransf = function(Nsim, n, s) {
  pmfVec <- vector(length=n) #Store n values calculated from the pmf 
  #Generate the constant divider
  divider <- 0
  for(l in 1:n){
    divider <- divider + 1/l^s
  }
  #Calcuate the k values from pmf
  for(v in 1:n){
    pmfVec[v] <- (1/v^s)/divider
  }
  
  U <- runif(Nsim)
  resVec <- rep(0,Nsim)
  
  for (i in 1:Nsim) {
    if(U[i] <= pmfVec[1]){
      resVec[i] = 1
    }
    #Check if the condition applies for a constant k
    for(k in 2:length(pmfVec)) {
      if(sum(pmfVec[1:(k-1)]) < U[i] && U[i] <= sum(pmfVec[1:k]) ) {
        resVec[i] <- k
      }
    }
  }
  
  return(list("result" = resVec, "theoretical" = pmfVec))
}

ditf <- simDiscInverseTransf(10000,10,3)
##### Plot 4B #####
theory <- data.frame(x=c(1:10), probability = ditf$theoretical)
res <- as.data.frame(table(ditf$result))
colnames(res) <- c("x", "probability")

p1 <- ggplot(theory, aes(x,probability)) + labs(title="Theory") + geom_bar(stat = "identity", fill="steelblue") 
p2 <- ggplot(res, aes(x,probability/10000)) + labs(y="probability", title="Simulation") + geom_bar(stat = "identity", fill="red")
 
grid.arrange(p1,p2, ncol = 1)
##### Plot 4B #####


######### 4C
#Simulate K by inverse transfer method
#For a geometric distribution, P(X=k) = p(1-p)^k-1, k=1,2,..
simK <- function(Nsim){
  pmfVec <- vector(length=100)
  for(v in 1:100){
    W <- rexp(1, 1/2)
    p <- 1-exp(-W)
    pmfVec[v] <- p*(1-p)^(v-1)
  }
  U <- runif(Nsim)
  resVec <- rep(0,Nsim)
  
  for (i in 1:Nsim) {
    #Check if the condition applies for a constant k
    if(U[i] <= pmfVec[1]){
      resVec[i] = 1
    }
    for(k in 2:length(pmfVec)) {
      if(sum(pmfVec[1:(k-1)]) < U[i] && U[i] <= sum(pmfVec[1:k]) ) {
        resVec[i] <- k
      }
    }
  }
  return(list("result" = resVec, "theoretical" = pmfVec))
}

sk <- simK(10000)
hist(sk$result)
#The plot of the simulated rv K proves that the
#preferential  attachment  model  from  Problem  3 converges K as N -> inf. (and a low p value)


##Compare simulation with real generated if needed
#W <- rexp(10000, 1/2)
#Z = 1-exp(-W)
#K = rgeom(10000,Z)
#hist(K)
