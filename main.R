
# Super Learner for Continuous Time 


###########Input:

#### The observed data should be input as:
# Time: Continuous observed time, The time origin should be 0 or has adjusted to 0
# X: Covariates, must be data.frame even only has 1 column  
# Event: (num or factor): 0,1,2...  0 should be censoring.

#### Super Learner will give a prediction on each row a NewX on the time period from 0 to max(Observed time).  

#NewX: new covariates for prediction, must be data.frame and the same covariates names as X
#Super Learner can compute at most 3 groups of new X, so we just need fit and optimize once:
#(1) only 1 NewX: 
#     Conditional CIF:  NewX=data.frame(Z=0.5, A=1)  (output name: CIF)
#   Or marginal CIF:    NewX=data.frame(Z=data$Z,A=data$A)  (output name: CIF)
#(2) Compute two NewX:  A=1:  NewX=data.frame(Z=data$Z, A=1)  (output name:CIF1)
#                       A=0:  NewX0=data.frame(Z=data$Z, A=0)  (output name:CIF0)
#(3) Compute three NewX:  Conditional/marginal CIF: NewX=data.frame(Z=0.5, A=1) (output name: CIF)
#                          A=1:  NewX1=data.frame(Z=data$Z, A=1)  (output name:CIF1)
#                          A=0:  NewX0=data.frame(Z=data$Z, A=0)  (output name:CIF0)


#### V: V folds cross-validation, default 3.


#### Candidate learners
# library: "coxph","km","gam","pch","loglog","weibull","exp", 
          #"gam_s"(gam_spline in mgcv), "rsf"(random survival forest, but it will be used to compute single event)
          #"gam_spline" (another gam_spline in survPen, unstable, slow, and may report errors)
          #"aalen" (Aalen additive hazard, may output negative hazard)
#The number of candidate learners may influence the running speed


#### Time-dependent Weights: 
#k and w_breaks: provied one (or none) of the two 
                 # k is the intervals of time for the combination of weights, the interval will
                 # be average splitted from 0 to max time, default k=1;
                 # w_breaks: a sequence from 0 to a time point, 
                 # eg: w_breaks=c(0,0.2,0.4,0.6,0.8,1.0), different intervals will use different weights
                 # The performance will be slow for lots of intervals and the results might
                 # not improve a lot, we suggest begin with default k=1.

# ObsWeights: the weights of observation, defualt the weights for all observations
# Parameters for random forest: default rsf_nsplit=3, rsf_ntree=100
# Parameters for pch: default pch_breaks=10
           
      

#######Output: 
# Output is a list include:
#library: the names of candidate learners
#weights: a list of the weights for each event, 
          # each list is a 2 dimensional matrix: time intervals for weights * candidate learners
          # the order is consistent with the learners in the library.

#hazard: cause-specific hazard, a 3 dimensional matrix: predicted observations (in NewX) * Time points * Events
#cif: a 3 dimensional matrix: predicted observations (in NewX) * Time points * Events
#time: 101 predicted time points.


library(survival)
library(pch) #if use pch
library(mgcv) #if use GAM or GAM_Spline
library(timereg) #if use Aalen Additive Hazard
#library(survPen) #if use GAM_Spline in SurvPen, it is unstable
library(randomForestSRC) #if use Random Forrest


library(cmprsk) # to generate theoretical cif graph


source("functionC.R")
source("data_graph.R")


#simulation example
result1 <- list()
#libraryC <- c("coxph","km","gam","pch","loglog","weibull","exp")
libraryC <- c("coxph","pch","km","gam","weibull","loglog")
libraryC <- c("weibull")
libraryC <- c("pch","weibull")
libraryC<- c("true2","coxph","km")
 


#High Dimension Sparsity
#X=data.frame(data[, paste0("X", 1:100)], A = data$A)
#NewX <- data.frame(matrix(0.5, nrow = 1, ncol = 100))
#colnames(NewX) <- paste0("X", 1:100)
#NewX$A <- 1


 for (i in 1:68) {
   print(i)
set.seed(123+i)
data <- W2data2(1000) 

# data <- W2dataU_I(1000) #independent frailty 
# data <- W2dataU_D(1000) #dependent frailty
# data <- W2dataU_D(100)  #small sample size
# data <- W2dataU_D_H(1000)   # high-dimensional sparsity


X = data.frame(Z=data$Z, A=data$A)
NewX = X    # Marginal CIFs
# NewX = data.frame(Z=0.5, A=1)  # Conditional CIFs

resultC  <- SuperLearnerC(Time=data$Time,
                         X=X,
                         Event=data$Event,
                         NewX=NewX,
                         # V=3,
                         library=libraryC,
                         #k=1,
                         # w_breaks=NULL,
                         #pch_breaks=20, ObsWeights=NULL, rsf_nsplit=3, rsf_ntree=100
                         ) 


# To only compute CIF(t;Z,A=1) and  CIF(t;Z,A=0), use:
#               SuperLearnerC(...,NewX = data.frame(Z=data$Z,A=1),
#                                 NewX0= data.frame(Z=data$Z,A=0),...)
#The output is result$CIF1 (for NewX) and result$CIF0 (for NewX0)


# To compute Conditional (or marginal) CIF,  CIF(t;Z,A=1)  and  CIF(t;Z,A=0), use:
#               SuperLearnerC(...,NewX = data.frame(Z=0.5,A=1),
#                                 NewX1= data.frame(Z=data$Z,A=1),...
#                                 NewX0= data.frame(Z=data$Z,A=0),...)
#The output is result$CIF (for NewX), result$CIF1 (for NewX1), result$CIF0 (for NewX0)

mean_cif <- apply(resultC$cif, c(2, 3), mean) 

result1[[i]] <- list(cif= mean_cif, time=resultC$time, weights =resultC$weights)

}

drawSL(W2data2, result1)
# drawSL(W2dataU_D, result1)
## Conditional 
# drawSL(W2data2_fixed, result1) # default value of W2data2_fixed: Z=0.5,A=1   

# drawSL_NULL(result1) without theoretical results, xlim=ylim=c(0,1)  
  


##For RD and RR:
libraryC <- c("coxph","gam","km")

result1 <- list()
result2 <- list()
for (i in 1:100) {
  set.seed(123+i)
  print(i)
  data <- W2dataU_D(1000) 
  result  <- SuperLearnerC(Time=data$Time,
                           X=data.frame(Z = data$Z,
                                        A = data$A),
                           Event=data$Event,
                           NewX=data.frame(Z = data$Z,
                                           A = 1),
                           NewX0=data.frame(Z = data$Z,
                                            A = 0),
                           # V=3,
                           library=libraryC,
                           # k=1,
                           # w_breaks=NULL,
                           #pch_breaks=10, ObsWeights=NULL, rsf_nsplit=3, rsf_ntree=100
  ) 
  
  
  
  # Draw graph
  mean_cif1 <- apply(result$cif1, c(2, 3), mean) 
  mean_cif0 <- apply(result$cif0, c(2, 3), mean) 
  result1[[i]] <- list(cif1= mean_cif1,cif0=mean_cif0, time=result$time)
  result2[[i]] <- result
  
}
#save(result1,result2, file = "rrrdC_U.Rdata")


drawrrrd(W2dataU_D,result1)






#########################################################

# Super Learner for Discrete Time

######Input:
#### The observed data should be input as:
# Time: Observed time, The time origin should be 0 or has adjusted to 0
# X: Covariates, must be data.frame even only has 1 column  
# Event or Event_indicator: One of the two must be specified.
#     Event: (num): 0,1,2...  0 should be censoring.
#     Event matrix: a matrix or data.frame, the columns of the events with such order: 
#                   Event 0 (Censoring), Event 1, Event 2... the values are 1 or 0.\

# time_intervals: e.g. seq(0, 1, length.out = 40 + 1) , the sequence to discreet time

#### Super Learner will give a prediction on each row a NewX on the time period 
####    from 0 to max(Observed time sequence).  
# NewX: new covariates for prediction, must be data.frame and the same covariates names as X
#Super Learner can compute at most 3 groups of new X, so we just need fit and optimize once:
#(1) only 1 NewX: 
#     Conditional CIF:  NewX=data.frame(Z=0.5, A=1)  (output name: CIF)
#   Or marginal CIF:    NewX=data.frame(Z=data$Z,A=data$A)  (output name: CIF)
#(2) Compute two NewX:  A=1:  NewX=data.frame(Z=data$Z, A=1)  (output name:CIF1)
#                       A=0:  NewX0=data.frame(Z=data$Z, A=0)  (output name:CIF0)
#(3) Compute three NewX:  Conditional/marginal CIF: NewX=data.frame(Z=0.5, A=1) (output name: CIF)
#                          A=1:  NewX1=data.frame(Z=data$Z, A=1)  (output name:CIF1)
#                          A=0:  NewX0=data.frame(Z=data$Z, A=0)  (output name:CIF0)


#### V: V folds cross-validation, default 3.


#### Candidate learners
# library: "VGAM.mn": multinomial logit model using VGAM package
#          "vgam_s":  multinomial logit model with smooth using VGAM package  
#          "nnet.mn": multinomial logit model using nnet package
#          "cv.glmnet.mn": multinomial logit model with elastic net
#          "lda", "qda", "naiveBayes", "svm", "xgboost", "knn", "rf" (random forest)

#The number of candidate learners may influence the running speed


#### Time-dependent Weights: 
#k and w_breaks: only specify 1 (or none) of the two:
# k is the intervals of time for the combination of weights, the interval will
# be average splitted from 0 to max time, default K=1;
# w_breaks: a sequence from 0 to a time point, 
# eg: w_breaks=c(0,0.2,0.4,0.6,0.8,1.0), different intervals will use different weights
# The performance will be slow for lots of intervals and the results might
# not improve a lot, we suggest begin with default k=1.

#### Others are the tuning parameters named after the model names, 
#### and default values are in the example.



######Output:
# Output is a list include:
#library: the names of candidate learners
#weights: a 2 dimensional matrix: time intervals (for weights) * candidate learners
#            the order is consistent with the learners in the library.

#hazard: cause-specific hazard, a 3 dimensional matrix: predicted observations (in NewX) * Time intervals * Events
#               (including Censoring as Event 0: hazards for all events each row == 1) 

#cif: a 3 dimensional matrix: predicted observations (in NewX) * Time intervals * Events 


###
#time: right boundaries of each time interval. e.g., hazard[,i,] corresponds the time interval [i-1,i).

library(VGAM) # if use multinomial logit regression in VGAM or vgam_s
library(nnet) # if use multi-logit regression in "nnet"
library(glmnet) # if use penalised multi-logti regression in "glmnet"
library(MASS) # if use lda or qda
library(e1071) # if use naive Bayes, or SVM
library(xgboost)  # if use xgboost
library(caret) # if use knn or random forest


library(cmprsk) # to generate theoretical cif graph

source("functionD.R")
source("data_graph.R")

#simulation  example for CIF

result2 <- list()
time_intervals <- seq(0, 1, length.out = 20 + 1) 
libraryD <- c("VGAM.mn","knn","xgboost","vgam_s")

libraryD <- c("rf")
libraryD <- c("VGAM.mn")

libraryD <- c("cv.glmnet.mn")
libraryD <- c("svm")

for (i in 1:100) {
  print(i)
set.seed(123+i)
data<- W2dataU_D(1000)

resultD <- SuperLearnerD( Time = data$Time,
                          X = data.frame(Z=data$Z, A=data$A),
                          Event = data$Event, #Event_indicator = NULL,
                          time_intervals = time_intervals,
                          NewX = data.frame(Z=data$Z, A=data$A),
                          #NewX = data.frame(Z=0.5, A=1),
                          #NewX0= data.frame(Z=data$Z,A=0),
                          #NewX1= data.frame(Z=data$Z,A=1),
                          #V=3, 
                          library = libraryD,
                          #k=1, w_breaks=NULL,
                          #VGAM.mn_maxit=200, xgb_nrounds=200, xgb_eta=0.1, 
                          #xgb_maxdep=8, xgb_verbose=0, knn_tuneLength= 10, rf_tuneLength =5
                          )


mean_cif <- apply(resultD$cif, c(2, 3), mean) 

# Cumulative hazard
#hazard <- resultD$hazard
#cumhazard <- array(NA, dim = c(dim(hazard)[1], dim(hazard)[2], dim(hazard)[3] - 1))

#for (j in 1:2) {
  # 对第 j 个事件，计算每个个体的累积 hazard（对时间维度进行 cumsum）
#  cumhazard[1,,j] <- cumsum(hazard[1, , j+1])
# }

result2[[i]] <- list(cif= mean_cif, 
                     #Hazard = cumhazard,
                     time=resultD$time*0.05,
                     
                     weights = resultD$weights)

}

for (i in 1:100){ result2[[i]]$time = result2[[i]]$time *2 }
drawSL(W2data2, result2) # CIF


#drawh(result2)  # Cumulative Hazard



##RD AND RR
time_intervals <- seq(0, 1, length.out = 40 + 1) 
# libraryD <- c("VGAM.mn","lda","naiveBayes","svm","xgboost","vgam_s")
libraryD <- c("VGAM.mn","lda","knn","xgboost","vgam_s")


result3<- list()
#result4<- list()

for (i in 1:100) {
  print(i)
  set.seed(123+i)
  data<- W2dataU_D(1000)
  
  resultD <- SuperLearnerD( Time = data$Time,
                            X = data.frame(Z=data$Z,A=data$A),
                            Event = data$Event, #Event_indicator = NULL,
                            time_intervals = time_intervals,
                            NewX = data.frame(Z=data$Z, A=data$A),
                            NewX0= data.frame(Z=data$Z,A=0),
                            NewX1= data.frame(Z=data$Z,A=1),
                            #V=3, 
                            library = libraryD,
                            #k=1, w_breaks=NULL,
                            #VGAM.mn_maxit=200, xgb_nrounds=200, xgb_eta=0.1, 
                            #xgb_maxdep=8, xgb_verbose=0, knn_tuneLength= 10, rf_tuneLength =5
  )
  
  mean_cif <- apply(resultD$cif, c(2, 3), mean) 
  mean_cif1 <- apply(resultD$cif1, c(2, 3), mean) 
  mean_cif0 <- apply(resultD$cif0, c(2, 3), mean) 
  result3[[i]] <- list(cif = mean_cif, cif1= mean_cif1,cif0=mean_cif0, time=resultD$time*0.025,
                       weights = resultD$weights)
  #result4[[i]] <- resultD
  
}

drawrrrd(W2dataU_D,result3)



drawrrrd(W2dataU_D,result3)


