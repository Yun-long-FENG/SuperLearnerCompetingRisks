# Simulation: Super Learner for Competing Risks


# Continuous-time Super Learner
library(survival)
library(pch) 
library(mgcv) 
library(randomForestSRC) 
library(cmprsk) # to generate theoretical cif graph

source("functionC.R")
source("data_graph.R")

result1 <- list()
libraryC <- c("coxph","km","gam","pch","loglog","weibull")
#libraryC <- c("coxph","km","gam","pch")  #Choose the candidate learners as we described in the paper and supplementary materials.

#Margical CIFs, RR and RD:
for (i in 1:100) {
   print(i)
   set.seed(123+i)
   data <- W2data2(1000) # For data without unmeasured common risk factor
   #data <- W2dataU_D(1000) #For supplementary materials: data with unmeasured common risk factor
   X = data.frame(Z=data$Z, A=data$A)
   NewX = X   
   resultC  <- SuperLearnerC(Time = data$Time,
                         X = X,
                         Event = data$Event,
                         NewX = NewX,
                         NewX1 = data.frame(Z=data$Z,A=1),
                         NewX0= data.frame(Z=data$Z,A=0),
                         # V = 3,
                         library = libraryC,
                         #k = 1,
                         #w_breaks = NULL,
                         #pch_breaks=20, ObsWeights=NULL, rsf_nsplit=3, rsf_ntree=100
                         )                    
   mean_cif <- apply(resultC$cif, c(2, 3), mean) 
   mean_cif1 <- apply(resultC$cif1, c(2, 3), mean) 
   mean_cif0 <- apply(resultC$cif0, c(2, 3), mean) 
   result1[[i]] <- list(cif = mean_cif, cif1 = mean_cif1, cif0 = mean_cif0,
                        time = resultC$time, weights = resultC$weights)
}

drawSL(W2data2, result1)
drawrrrd(W2data2, result1)
#drawSL(W2dataU_D, result1)
#drawrrrd(W2dataU_D, result1)


## Conditional CIFs:
for (i in 1:100) {
   print(i)
   set.seed(123+i)
   data <- W2data2(1000)  # For data without unmeasured common risk factor
   #data <- W2dataU_D(1000) #For supplementary materials: data with unmeasured common risk factor
   #data <- W2dataU_D(100) #For supplementary materials: data with unmeasured common risk factor and small sample size
   #data <- W2dataU_D_H(1000) #For supplementary materials: data with unmeasured common risk factor and high-dimensional sparse covariates
   X = data.frame(Z=data$Z, A=data$A)
   NewX = data.frame(Z=0.5, A=1)  # Conditional CIFs
   resultC  <- SuperLearnerC(Time =data$Time,
                         X = X,
                         Event = data$Event,
                         NewX = NewX,
                         # V = 3,
                         library = libraryC,
                         #k = 1,
                         #w_breaks=NULL,
                         #pch_breaks=20, ObsWeights=NULL, rsf_nsplit=3, rsf_ntree=100
                         )                    
   mean_cif <- apply(resultC$cif, c(2, 3), mean) 
   result1[[i]] <- list(cif= mean_cif, time=resultC$time, weights =resultC$weights)
}

drawSL(W2data2_fixed, result1) # default value of W2data2_fixed: Z=0.5,A=1   
#drawSL(W2dataU_fixed, result1) # For supplementary materials: data with unmeasured common risk factor, default value of W2dataU_fixed: Z=0.5,A=1
  


#########################################################

# Super Learner for Discrete Time

library(VGAM) 
library(xgboost)  
library(caret) 
library(cmprsk) # to generate theoretical cif graph

source("functionD.R")
source("data_graph.R")

result2 <- list()
time_intervals <- seq(0, 1, length.out = 40 + 1) 
libraryD <- c("VGAM.mn","knn","xgboost","vgam_s")
#libraryD <- c("VGAM.mn")  # Only choose the multinomial regression as the candidate learner, which is our comparison method in the paper and supplementary materials.

for (i in 1:100) {
  print(i)
  set.seed(123+i)
  data<- W2data2(1000)
  resultD <- SuperLearnerD( Time = data$Time,
                             X = data.frame(Z=data$Z, A=data$A),
                             Event = data$Event, #Event_indicator = NULL,
                             time_intervals = time_intervals,
                             NewX = data.frame(Z=data$Z, A=data$A),
                             NewX0= data.frame(Z=data$Z,A=0),
                             NewX1= data.frame(Z=data$Z,A=1),
                             #V = 3, 
                             library = libraryD,
                             #k = 1, w_breaks = NULL,
                             #VGAM.mn_maxit=200, xgb_nrounds=200, xgb_eta=0.1, 
                             #xgb_maxdep=8, xgb_verbose=0, knn_tuneLength= 10, rf_tuneLength =5
                           )

  mean_cif <- apply(resultD$cif, c(2, 3), mean) 
  mean_cif1 <- apply(resultD$cif1, c(2, 3), mean) 
  mean_cif0 <- apply(resultD$cif0, c(2, 3), mean) 
  result2[[i]] <- list(cif= mean_cif, cif1 = mean_cif1, cif0 = mean_cif0, 
                       time = resultD$time*0.025, weights = resultD$weights)
}

drawSL(W2data2, result2) 
drawrrrd(W2data2,result2)
#drawSL(W2dataU_D, result2)
#drawrrrd(W2dataU_D, result2)





########################################################################
######## State learner
# https://github.com/amnudn/joint-survival-super-learner
library(here)
library(targets)
tar_source(here("statelearner/R-code/functions"))

source("functionC.R")
source("data_graph.R")

library(tarchetypes)
library(parallel)
library(data.table)
library(riskRegression)
library(survival)
library(randomForestSRC)
library(prodlim)
library(MASS)


learners <- list(
  cox = list(model = "cox", x_form = ~Z+A),
  cox_penalty = list(model = "GLMnet", x_form = ~Z+A),
  N_Aa = list(model = "cox", x_form = ~1),
  #N_Aa_strat = list(model = "cox", x_form = ~strata(A)),
  rf = list(model = "rfsrc", x_form = ~Z+A, ntree = 50)
)

result_list2 <- list()
times <- seq(0, 1, length.out = 100) #Predict times, but we will only use and draw the min of the max of the time in the 100 simulations
for (i in 1:100) {
  time0<- Sys.time()
  print(i)
  set.seed(123+i)
  data <- W2data2(1000) # Choose different data generating processes as in the Super Learner
  #data <- W2dataU_D(1000)
  data <- as.data.table(data)
  sl = statelearner(learners = list(cause1 = learners,
                                    cause2 = learners,
                                    censor = learners),
                    data = data,
                    time = 1,
                    time_name = "Time",
                    status_name = "Event")
  #Mariginal CIFs
  CHF1 <- predictCHF(sl$fitted_winners$cause1,newdata = data.frame(Z=data$Z,A=data$A), times = times)
  CHF2 <- predictCHF(sl$fitted_winners$cause2,newdata = data.frame(Z=data$Z,A=data$A), times = times)
  #Conditional CIFs
  #CHF1 <- predictCHF(sl$fitted_winners$cause1,newdata = data.frame(Z=rep(0.5,1000),A=rep(1,1000)), times = times)
  #CHF2 <- predictCHF(sl$fitted_winners$cause2,newdata = data.frame(Z=rep(0.5,1000),A=rep(1,1000)), times = times)

  H1 <- compute_hazard_from_CHF(CHF1, times) #This function is in "data_graph.R"
  H2 <- compute_hazard_from_CHF(CHF2, times)

  h_matrix <- array(c(H1, H2), dim = c(nrow(H1), ncol(H1), 2))
  cif <- compute_cif(h_matrix,times) #This function is in "functionC.R"

  mean_hazard <- apply(h_matrix, c(2, 3), mean)     
  mean_cif <- apply(cif, c(2, 3), mean) 
  result_list2[[i]] <- list(
  hazard = mean_hazard,
  cif = mean_cif,
  time = times) 
  print(Sys.time()-time0)
}

drawState(W2data2, result_list2)
#drawState(W2dataU_D, result_list2)




#############################################################
############## Random Survival Forest
library(randomForestSRC)
library(cmprsk)

result_rf1 <- list()
for (i in 1:100) {
  set.seed(123+i)
  data<- W2data2(1000)
  #data <- W2dataU_D(1000) #For supplementary materials: data with unmeasured common risk factor
  #data <- W2dataU_D(100) #For supplementary materials: data with unmeasured common risk factor and small sample size
  #data <- W2dataU_D_H(1000) #For supplementary materials: data with unmeasured common risk factor and high-dimensional sparse covariates  
  rf <- rfsrc(Surv(Time, Event) ~ Z + A, data, nsplit = 3, ntree = 100)
  
  #Marginal CIFs
  cif <- rf$cif.oob
  time <-rf$time.interest
  cif1<- apply(cif[,,1], 2, mean)
  cif2<- apply(cif[,,2], 2, mean)
  cifc<- cbind(cif1,cif2)
  result_rf1[[i]]<- list(cif = cifc, times = time) 
}


# The predict time of RSF is fixed 150 points, we need to calculate the minimum of the maximum of time in 100 simulations
time <- rep(0,100)
for (i in 1:100){
  set.seed(123+i)
  data <- W2data2(1000)  #Choose different data generating functions
  time[i] <- max(data$Time)
}
mintime <- min(time)  
drawRSF(W2data2, result_rf1, mintime) 


##Condtiional CIF with high-dimensional sparse data 
result_rf1<- list()
for (i in 1:100){
   set.seed(123+i)
   data<- W2dataU_D_H(1000)
   rf <- rfsrc(Surv(Time, Event) ~ paste0("X", 1:100) +A, data, nsplit = 10, ntree = 300)
   #conditional cif
   newdata <- data.frame(matrix(0.5, nrow = 1, ncol = 100))
   colnames(newdata) <- paste0("X", 1:100)
   newdata$A <- 1 
   rf_pred <- predict(rf, newdata = newdata)

   cif_pred <- rf_pred$cif       
   time_pred <- rf_pred$time.interest 

   cif1 <- as.vector(cif_pred[1, , 1])  # Event 1 CIF
   cif2 <- as.vector(cif_pred[1, , 2])  # Event 2 CIF
   cifc<- cbind(cif1,cif2)

   result_rf1[[i]]<- list(cif = cifc, time = time_pred) 
}

time <- rep(0,100)
for (i in 1:100){
  set.seed(123+i)
  data <- W2dataU_D_H(1000)  
  time[i] <- max(data$Time)
}
mintime <- min(time)  
drawRSF(W2dataU_fixed, result_rf1, mintime) 

