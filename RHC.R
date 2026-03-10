# Follow https://htmlpreview.github.io/?https://github.com/Bella2001/causalCmprsk/blob/dev/index.html to download and clean the data
library(Hmisc)
library(naniar)
library(kableExtra)
library(dplyr)

Hmisc::getHdata(rhc) # loading data using the Hmisc package
rhc_raw <- rhc

miss.report <- rhc_raw %>% miss_var_summary()
kable(miss.report[1:10,]) 
kable(miss.report[1:10,])  %>% kable_styling(bootstrap_options = "striped", full_width = F)


rhc_cleaning1 <- rhc_raw %>%
  mutate(RHC = as.numeric(swang1 == "RHC"), 
         trt=swang1,
         E = ifelse(is.na(dthdte), 1, ifelse(dschdte==dthdte, 2, 1)), 
         Time = dschdte - sadmdte, #=time-to-"Discharge" (either alive or due to death in a hospital)
         T.death = ifelse(is.na(dthdte), lstctdte - sadmdte, dthdte - sadmdte), #=min(Time to death before or after discharge, Time to a last follow-up visit)
         D = ifelse(death =="No", 0, 1), # death indicator
         sex_Male = as.numeric(sex == "Male"),
         race_black = as.numeric(race == "black"),
         race_other = as.numeric(race == "other"),
         income_11_25K = as.numeric(income == "$11-$25k"),
         income_25_50K = as.numeric(income == "$25-$50k"),
         income_50K = as.numeric(income == "> $50k"),
         ninsclas_Private_Medicare = as.numeric(ninsclas == "Private & Medicare"),
         ninsclas_Medicare = as.numeric(ninsclas == "Medicare"),
         ninsclas_Medicare_Medicaid = as.numeric(ninsclas == "Medicare & Medicaid"),
         ninsclas_Medicaid = as.numeric(ninsclas == "Medicaid"),
         ninsclas_No_Insurance = as.numeric(ninsclas == "No insurance"),
         # combine cat1 with cat2, i.e the primary disease category with the secondary disease category:         
         cat_CHF = as.numeric(cat1 == "CHF" | (!is.na(cat2))&(cat2 == "CHF") ),
         cat_Cirrhosis = as.numeric(cat1 == "Cirrhosis" | (!is.na(cat2))&(cat2 == "Cirrhosis")),
         cat_Colon_Cancer = as.numeric(cat1 == "Colon Cancer" | (!is.na(cat2))&(cat2 == "Colon Cancer")),
         cat_Coma = as.numeric(cat1 == "Coma" | (!is.na(cat2))&(cat2 == "Coma")),
         cat_COPD = as.numeric(cat1 == "COPD" | (!is.na(cat2))&(cat2 == "COPD")),
         cat_Lung_Cancer = as.numeric(cat1 == "Lung Cancer" | (!is.na(cat2))&(cat2 == "Lung Cancer")),
         cat_MOSF_Malignancy = as.numeric(cat1 == "MOSF w/Malignancy" | (!is.na(cat2))&(cat2 == "MOSF w/Malignancy")),
         cat_MOSF_Sepsis = as.numeric(cat1 == "MOSF w/Sepsis" | (!is.na(cat2))&(cat2 == "MOSF w/Sepsis")),
         dnr1_Yes = as.numeric(dnr1 == "Yes"),
         card_Yes = as.numeric(card == "Yes"),
         gastr_Yes = as.numeric(gastr == "Yes"),
         hema_Yes = as.numeric(hema == "Yes"),
         meta_Yes = as.numeric(meta == "Yes"),
         neuro_Yes = as.numeric(neuro == "Yes"),
         ortho_Yes = as.numeric(ortho == "Yes"),
         renal_Yes = as.numeric(renal == "Yes"),
         resp_Yes = as.numeric(resp == "Yes"),
         seps_Yes = as.numeric(seps == "Yes"),
         trauma_Yes = as.numeric(trauma == "Yes"),
         ca_Yes = as.numeric(ca == "Yes"),
         ca_Metastatic = as.numeric(ca == "Metastatic")
  )

# variables selection and data reordering:
rhc_full <- rhc_cleaning1 %>% 
  select(ptid, RHC, trt, Time, T.death, E, D, sex_Male, 
         age, edu, race_black, race_other, income_11_25K, income_25_50K, income_50K,
         ninsclas_Private_Medicare, ninsclas_Medicare, ninsclas_Medicare_Medicaid,
         ninsclas_Medicaid, ninsclas_No_Insurance, 
         cat_CHF, cat_Cirrhosis, cat_Colon_Cancer, cat_Coma, cat_COPD, cat_Lung_Cancer, 
         cat_MOSF_Malignancy, cat_MOSF_Sepsis,
         # diagnoses:
         dnr1_Yes, card_Yes, gastr_Yes, hema_Yes, meta_Yes, neuro_Yes, ortho_Yes, renal_Yes, 
         resp_Yes, seps_Yes, trauma_Yes,
         ca_Yes, ca_Metastatic,
         # lab tests:
         wtkilo1, hrt1, meanbp1, resp1, temp1,
         aps1, das2d3pc, scoma1, 
         surv2md1, alb1, bili1, crea1, hema1, paco21, 
         pafi1, ph1, pot1, sod1, wblc1,
         # all variables with "hx" are preexisting conditions: 
         amihx, cardiohx, chfhx, chrpulhx, dementhx, 
         gibledhx, immunhx, liverhx, malighx, psychhx, 
         renalhx, transhx, 
         death, sadmdte, dschdte, dthdte, lstctdte)


# omit 1 obs with missing discharge date for the length-of-stay analysis:
rhc_full <- rhc_full[!is.na(rhc_full$Time), ]

# We only use the first 30 days for the analysis, the observed time is censored at 30 days
rhc_full$E[rhc_full$Time > 30] <- 0
rhc_full$Time[rhc_full$Time > 30] <- 30

covs.names <- c("age", "sex_Male", "edu", "race_black", "race_other",
                "income_11_25K", "income_25_50K", "income_50K", 
                "ninsclas_Private_Medicare", "ninsclas_Medicare", "ninsclas_Medicare_Medicaid",
                "ninsclas_Medicaid", "ninsclas_No_Insurance",
                "cat_CHF", "cat_Cirrhosis", "cat_Colon_Cancer", "cat_Coma", "cat_COPD", "cat_Lung_Cancer", 
                "cat_MOSF_Malignancy", "cat_MOSF_Sepsis", 
                "dnr1_Yes", "wtkilo1", "hrt1", "meanbp1", 
                "resp1", "temp1", 
                "card_Yes", "gastr_Yes", "hema_Yes", "meta_Yes", "neuro_Yes", "ortho_Yes", 
                "renal_Yes", "resp_Yes", "seps_Yes", "trauma_Yes", 
                "ca_Yes", "ca_Metastatic",
                "amihx", "cardiohx", "chfhx", "chrpulhx",
                "dementhx", "gibledhx", "immunhx", "liverhx", 
                "malighx", "psychhx", "renalhx", "transhx",
                "aps1", "das2d3pc", "scoma1", "surv2md1",
                "alb1", "bili1", "crea1", "hema1", "paco21", "pafi1", 
                "ph1", "pot1", "sod1", "wblc1")

X=data.frame(rhc_full[, c(covs.names,"RHC")])
NewX = X
NewX0 = data.frame(rhc_full[, covs.names], "RHC" = 0)
NewX1 = data.frame(rhc_full[, covs.names], "RHC" = 1)



## Continuous-time Super Learner
library(survival)
#library(pch) #if use pch
library(mgcv) #if use GAM or GAM_Spline
#library(timereg) #if use Aalen Additive Hazard
#library(survPen) #if use GAM_Spline in SurvPen, it is unstable
#library(randomForestSRC) #if use Random Forrest

source("functionC.R")

libraryC <- c("coxph","km","gam","weibull","gam_s")

resultC  <- SuperLearnerC(Time=rhc_full$Time,
                            X=X,
                            Event=rhc_full$E,
                            NewX=NewX,
                            NewX0= NewX0,
                            NewX1= NewX1,
                            # V=3,
                            library=libraryC,
                            # w_breaks=NULL,
                            pch_breaks=20, ObsWeights=NULL, rsf_nsplit=3, rsf_ntree=100
                          ) 
  
mean_cif <- apply(resultC$cif, c(2, 3), mean) 
mean_cif1 <- apply(resultC$cif1, c(2, 3), mean) 
mean_cif0 <- apply(resultC$cif0, c(2, 3), mean) 

result1 <- list(cif = mean_cif, cif1 = mean_cif1, cif0 = mean_cif0, 
                   time = resultC$time, weights = resultC$weights)



## Discrete-time Super Learner
library(VGAM) # if use multinomial logit regression in VGAM or vgam_s
library(xgboost)  # if use xgboost
library(caret) # if use knn or random forest

source("functionD.R")

time_intervals <- seq(0, 30, length.out = 30 + 1) 
libraryD <- c("VGAM.mn","xgboost","knn","vgam_s")
resultD <- SuperLearnerD( Time = rhc_full$Time,
                          X = X,
                          Event = rhc_full$E, #Event_indicator = NULL,
                          time_intervals = time_intervals,
                          NewX = X,
                          NewX0 = NewX0,
                          NewX1 = NewX1,
                          #V=3, 
                          library = libraryD,
                          #k=1, w_breaks=NULL,
                          #VGAM.mn_maxit=200, xgb_nrounds=200, xgb_eta=0.1, 
                          #xgb_maxdep=8, xgb_verbose=0, knn_tuneLength= 10, rf_tuneLength =5
                       )

mean_cif <- apply(resultD$cif, c(2, 3), mean) 
mean_cif1 <- apply(resultD$cif1, c(2, 3), mean) 
mean_cif0 <- apply(resultD$cif0, c(2, 3), mean) 

result2 <- list(cif= mean_cif, cif1 = mean_cif1, cif0 = mean_cif0, time=resultD$time, weights =resultC$weights)



########Compare with causalCmprsk: code from https://htmlpreview.github.io/?https://github.com/Bella2001/causalCmprsk/blob/dev/index.html 
library("causalCmprsk")  
form.txt <- paste0("RHC", " ~ ", paste0(covs.names, collapse = "+"))
trt.formula <- as.formula(form.txt)
res.stab.ATE <- fit.nonpar(df=rhc_full, X="Time", E="E", trt.formula=trt.formula, 
                             wtype="stab.ATE", cens=0, conf.level=0.95, bs=TRUE, nbs.rep=50,
                             seed=17, parallel = FALSE)
  
df <- rbind(summary(res.stab.ATE, event=1, estimand="CIF"),
             summary(res.stab.ATE, event=2, estimand="CIF"))
df$Event_TRT <- factor(2*(df$Event==2) + 1*(df$TRT==0))
levels(df$Event_TRT) <- c("Discharge-RHC", "Discharge-No RHC",
                            "In-hospital death-RHC", "In-hospital death-No RHC")
  
df2 <- rbind(summary(res.stab.ATE, event=1, estimand="RD"),
            summary(res.stab.ATE, event=2, estimand="RD"))
df2$Event <- as.factor(df2$Event)
levels(df2$Event) <- c("Discharge", "In-hospital death")

df3 <- rbind(summary(res.stab.ATE, event=1, estimand="RR"),
            summary(res.stab.ATE, event=2, estimand="RR"))
df3$Event <- as.factor(df3$Event)
levels(df3$Event) <- c("Discharge", "In-hospital death")
  
  
  
### Plotting: RHC and No RHC:  
plot_compare <- function(result1, df) {
    time_grid <- result1$time
    
    plot(time_grid, result1$cif1[,1], type="l", col="blue", lwd=2.5, lty=1,
         xlab="Time (days)", ylab="CIF", xlim=c(0,30), ylim=c(0,0.55),
         main="CIFs - RHC")
    lines(time_grid, result1$cif1[,2], col="red", lwd=2.5, lty=2)
    
    lines(df$time[df$Event_TRT=="Discharge-RHC"], 
          df$CIF[df$Event_TRT=="Discharge-RHC"], col="green", lty=3, lwd=2)
    lines(df$time[df$Event_TRT=="In-hospital death-RHC"], 
          df$CIF[df$Event_TRT=="In-hospital death-RHC"], col="orange", lty=4, lwd=2)
    
    s <- subset(df, Event_TRT=="Discharge-RHC")
    polygon(c(s$time, rev(s$time)),
            c(s$CIU.CIF, rev(s$CIL.CIF)),
            col=adjustcolor("green", alpha.f=0.3), border=NA)
    
    s <- subset(df, Event_TRT=="In-hospital death-RHC")
    polygon(c(s$time, rev(s$time)),
            c(s$CIU.CIF, rev(s$CIL.CIF)),
            col=adjustcolor("orange", alpha.f=0.3), border=NA)
    
    legend("bottomright",
           legend=c("Discharge alive-Super Learner","In-hospital death-Super Learner","Discharge-alive-causalCmprsk","In-hospital death-causalCmprsk"),
           col=c("blue","red","green","orange"),
           lty=c(1,2,3,4), lwd=2, bty="o")
    
    plot(time_grid, result1$cif0[,1], type="l", col="blue", lwd=2.5, lty=1,
         xlab="Time (days)", ylab="CIF", xlim=c(0,30), ylim=c(0,0.55),
         main="CIFs - No RHC")
    lines(time_grid, result1$cif0[,2], col="red", lwd=2.5, lty=2)
    
    lines(df$time[df$Event_TRT=="Discharge-No RHC"], 
          df$CIF[df$Event_TRT=="Discharge-No RHC"], col="green", lty=3, lwd=2)
    lines(df$time[df$Event_TRT=="In-hospital death-No RHC"], 
          df$CIF[df$Event_TRT=="In-hospital death-No RHC"], col="orange", lty=4, lwd=2)
    
    s <- subset(df, Event_TRT=="Discharge-No RHC")
    polygon(c(s$time, rev(s$time)),
            c(s$CIU.CIF, rev(s$CIL.CIF)),
            col=adjustcolor("green", alpha.f=0.3), border=NA)
    
    s <- subset(df, Event_TRT=="In-hospital death-No RHC")
    polygon(c(s$time, rev(s$time)),
            c(s$CIU.CIF, rev(s$CIL.CIF)),
            col=adjustcolor("orange", alpha.f=0.3), border=NA)
    
    legend("bottomright",
           legend=c("Discharge alive-Super Learner","In-hospital death-Super Learner","Discharge-alive-causalCmprsk","In-hospital death-causalCmprsk"),
           col=c("blue","red","green","orange"),
           lty=c(1,2,3,4), lwd=2, bty="o")
  }
  
plot_compare(result1, df)  
plot_compare(result2, df)  

# Plot RD and RR:
draw_compare_rd_rr <- function(result1, df_rd, df_rr) {
  time <- result1$time
  cif1 <- result1$cif1
  cif0 <- result1$cif0
  
  rd_my <- cif1 - cif0
  rr_my <- cif1 / cif0
  
  # -------------------- RD --------------------

  plot(time, rd_my[,1], type="n", xlim=c(0,30), ylim=range(c(rd_my, df_rd$CIL, df_rd$CIU), na.rm=TRUE),
       xlab="Time (days)", ylab="Risk Difference", main="Risk Differences (RD)")

  s <- subset(df_rd, Event=="Discharge")
  polygon(c(s$time, rev(s$time)), c(s$CIU, rev(s$CIL)), 
            col=adjustcolor("green", alpha.f=0.3), border=NA)
  
  s <- subset(df_rd, Event=="In-hospital death")
  polygon(c(s$time, rev(s$time)), c(s$CIU, rev(s$CIL)), 
            col=adjustcolor("orange", alpha.f=0.3), border=NA)

  lines(time, rd_my[,1], col="blue", lwd=2.5, lty=1)
  lines(time, rd_my[,2], col="red",  lwd=2.5, lty=2)
  
  rd <- df_rd$RD[df_rd$Event=="Discharge"]
  time1 <- df_rd$time[df_rd$Event=="Discharge"]
  lines(time1, rd, col = "green", lwd=2, lty=3)
  
  rd2 <- df_rd$RD[df_rd$Event=="In-hospital death"]
  time2 <- df_rd$time[df_rd$Event=="In-hospital death"]
  lines(time2, rd2, col = "orange", lwd=2, lty=4)
  
  legend("right", legend=c("Discharge alive-Super Learner","In-hospital death-Super Learner","Discharge alive-causalCmprsk","In-hospital death-causalCmprsk"),
         col=c("blue","red","green","orange"), lty=c(1,2,3,4), lwd=2, bty="n",cex=0.8,pt.cex=0.8)
  
  for (event in c("Discharge","In-hospital death")) {
    s <- subset(df_rd, Event==event & abs(time-30)==min(abs(s$time-30)))
    cat("Other RD at 30 days (", event, "): Estimate=", s$RD, 
        ", CI=(", s$CIL, ", ", s$CIU, ")\n", sep="")
  }
  
  # ------------------ RR --------------------
  if (!is.null(df_rr)) {
    rr_my <- cif1 / cif0
    plot(time, rr_my[,1], type="n", xlim=c(0,30), ylim=range(c(rr_my, df_rr$CIL, df_rr$CIU), na.rm=TRUE),
         xlab="Time (days)", ylab="Risk Ratio", main="Risk Ratios (RR)")

      
    s <- subset(df_rr, Event=="Discharge")
    polygon(c(s$time, rev(s$time)), c(s$CIU, rev(s$CIL)), 
              col=adjustcolor("green", alpha.f=0.3), border=NA)
    
    s <- subset(df_rr, Event=="In-hospital death")
    polygon(c(s$time, rev(s$time)), c(s$CIU, rev(s$CIL)), 
              col=adjustcolor("orange", alpha.f=0.3), border=NA)
    
    rr <- df_rr$RR[df_rr$Event=="Discharge"]
    time1 <- df_rr$time[df_rr$Event=="Discharge"]
    lines(time1, rr, col = "green", lwd=2, lty=3)
    
    rr2 <- df_rr$RR[df_rr$Event=="In-hospital death"]
    time2 <- df_rr$time[df_rr$Event=="In-hospital death"]
    lines(time2, rr2, col = "orange", lwd=2, lty=4)
    
    lines(time, rr_my[,1], col="blue", lwd=2, lty=1)
    lines(time, rr_my[,2], col="red",  lwd=2, lty=2)
    
    legend("bottomright", legend=c("Discharge alive-Super Learner","In-hospital death-Super Learner","Discharge alive-causalCmprsk","In-hospital death-causalCmprsk"),
           col=c("blue","red","green","orange"), lty=c(1,2,3,4), lwd=2, bty="o")
    
    for (event in c("Discharge","In-hospital death")) {
      s <- subset(df_rr, Event==event & abs(time-30)==min(abs(s$time-30)))
      cat("Other RR at 30 days (", event, "): Estimate=", s$RR, 
          ", CI=(", s$CIL, ", ", s$CIU, ")\n", sep="")
    }
  }
}

draw_compare_rd_rr(result1, df2, df3)
draw_compare_rd_rr(result2, df2, df3)
