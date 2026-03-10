
create_event_data <- function(event) {
  event_levels <- sort(unique(event))
  event_data_list <- list()
  
  for (n in event_levels) {
    if (n == 0) next  
    En <- ifelse(event == n, 1, 0)
    event_data_list[[paste0("E", n)]] <- En
  }
  
  event_data_list <- event_data_list[sort(names(event_data_list))]
  return(event_data_list)
}


softmax <- function(params) {
  exp_params <- exp(params)
  denominator <- sum(exp_params) + 1

  weights <- exp_params / denominator
  additional_weight <- 1 / denominator
  
  result <- c(weights, additional_weight)
  
  return(result)
}

expand <- function(w) {
  # w: [0, 1]  project to [-0.001, 1.001]
  w1 <- w * (1.001 + 0.001) - 0.001
  w1 <- pmax(w1,0)
  w2 <- w1 / sum(w1)
  return(w2)
}


hazard_function <- function(survival, time, min_surv = 1e-6) {
  survival <- pmax(survival, min_surv)
  logS <- log(survival)
  n <- length(survival)
  
  hazard <- numeric(n)
  
  for (i in 2:(n - 1)) {
    hazard[i] <- - (logS[i + 1] - logS[i - 1]) / (time[i + 1] - time[i - 1])
  }
  
  hazard[1]     <- - (logS[2] - logS[1]) / (time[2] - time[1])
  hazard[n]     <- - (logS[n] - logS[n - 1]) / (time[n] - time[n - 1])
  
  return(hazard)
}



SL_fit <- function(Time, E, X, library,ObsWeights=NULL,pch_breaks=10,rsf_nsplit=3,rsf_ntree=100){   
  
  fits <- list()
  data <- data.frame(Time =Time, E=E, X)
  
  covariates <- colnames(X)
  formula <- as.formula(paste("Surv(Time, E) ~", 
                              paste(covariates, collapse = " + ")))
  

  for (learner in library) {
    if (learner == "coxph") {
      fits[[learner]] <- coxph(formula, data = data, 
                               control = coxph.control(iter.max = 1000), weights= ObsWeights,x = TRUE)
    } else if (learner == "km") {
      fits[[learner]] <- survfit(Surv(Time, E) ~ 1, data = data, weights = ObsWeights)
    } else if (learner == "pch") {
    
      fits[[learner]] <- pch::pchreg(formula, breaks = pch_breaks, data = data, weights =  ObsWeights)
    } else if (learner == "gam") {
      gam_model <- as.formula(paste("Time ~", paste(covariates, collapse = " + ")))
      fits[[learner]] <- mgcv::gam(gam_model, family = mgcv::cox.ph(), 
                                   data = data, weights = E) 
    } else if (learner == "gam_s") {
      unique_counts <- sapply(X, function(col) length(unique(col)))
      smooth_vars <- names(unique_counts[unique_counts > 6])
      linear_vars <- names(unique_counts[unique_counts <= 6])
      rhs_terms <- c(
        paste0("s(", smooth_vars, ")"),
        linear_vars
      )
      gam_model2 <- as.formula(
        paste("Time ~", paste(rhs_terms, collapse = " + "))
      )
      fits[[learner]] <- mgcv::gam(gam_model2, family = mgcv::cox.ph(), 
                                   data = data, weights = E)
    } else if (learner =="exp") { 
      if (is.null(fits[["fit.pos"]])) {
        if(any(Time == 0 & E == 1)) {
          timepos <- as.numeric(Time > 0 & E == 1)
          fits[["fit.pos"]] <- stats::glm(timepos ~ ., data=cbind(timepos, X)[E == 1,], family='binomial', weights = ObsWeights[E == 1])
        } else {
          fits[["fit.pos"]] <- 1
        }  #Parametric regression learners we use the code from Westling et el: survSuperLearner https://github.com/tedwestling/survSuperLearner/blob/master/R/SL_wrappers.R
      }
      fits[[learner]] <- survival::survreg(survival::Surv(Time[Time > 0], E[Time > 0]) ~ .,
                                           data = X[Time > 0, , drop = FALSE],
                                           weights = ObsWeights[Time > 0], dist = 'exponential')
    } else if (learner == "weibull"){
      if (is.null(fits[["fit.pos"]])) {
        if(any(Time == 0 & E == 1)) {
          timepos <- as.numeric(Time > 0 & E == 1)
          fits[["fit.pos"]] <- stats::glm(timepos ~ ., data=cbind(timepos, X)[E == 1,], family='binomial', weights = ObsWeights[E == 1])
        } else {
          fits[["fit.pos"]] <- 1
        }
      }
      fits[[learner]] <- survival::survreg(survival::Surv(Time[Time > 0], E[Time > 0]) ~ .,
                                           data = X[Time > 0,,drop=FALSE],
                                           weights = ObsWeights[Time > 0], dist = 'weibull')
    }  else if (learner == "loglog") {
      if (is.null(fits[["fit.pos"]])) {
        if(any(Time == 0 & E == 1)) {
          timepos <- as.numeric(Time > 0 & E == 1)
          fits[["fit.pos"]] <- stats::glm(timepos ~ ., data=cbind(timepos, X)[E == 1,], family='binomial', weights = ObsWeights[E == 1])
        } else {
          fits[["fit.pos"]] <- 1
        }
      }
      fits[[learner]] <- survival::survreg(survival::Surv(Time[Time > 0], E[Time > 0]) ~ .,
                                           data = X[Time > 0,,drop=FALSE],
                                           weights = ObsWeights[Time > 0], dist = 'loglogistic')
    }  else if (learner == "aalen"){
      fits[[learner]] <- aalen(formula, data = data, resample.iid = TRUE, weights = ObsWeights)
    } else if (learner == "gam_spline") {
      #library(survPen)
      spline_formula <- as.formula(paste("~ tensor(Time, ", paste(covariates, collapse = ", "), ")"))
      fits[[learner]] <- survPen(spline_formula, data = data, 
                                 t1 = Time, event = E,
                                 max.it.beta = 3000,
                                 tol.beta = 1e-2)
    } else if (learner == "rsf"){
      fits[[learner]]  <- rfsrc(formula, data, nsplit = rsf_nsplit, ntree = rsf_ntree)
    } 
  }
  return(fits)
}


SL_pred  <- function (newTime, newX, fits,  library) {

  compute_hazard <- function(X_row) {
    
    rownames(X_row) <- NULL  
    X_expanded <- X_row[rep(1, length(newTime)),,  drop = FALSE ]
    newdata <- data.frame(Time = newTime, X_expanded)

    hazards <- vector("list", length(library))
    for (k in seq_along(library)) {
      learner <- library[k]
      if (learner == "coxph") {
        survfit_cox <- survfit(fits[["coxph"]], newdata = X_row)
        time_cox <- survfit_cox$time
        surv_cox <- survfit_cox$surv
        stepfun_cox <- stepfun(time_cox, c(1, surv_cox))
        surv_cox <- stepfun_cox(newTime)
        hazards[[k]] <- hazard_function(surv_cox, newTime)
      } else if (learner == "km") {
        summary_km <- summary(fits[["km"]], times = newTime)
        stepfun_h2 <- stepfun(summary_km$time, c(1, summary_km$surv))
        surv_km <- stepfun_h2(newTime)
        hazards[[k]] <- hazard_function(surv_km, newTime)
      } else if (learner == "pch") {
        pre_pch <- predict(fits[["pch"]], type = "distr", 
                           newdata = newdata)
        hazards[[k]] <- pre_pch$haz 
      } else if (learner == "gam") {
        pred <- predict(fits[["gam"]], newdata = newdata, type = "response", se = FALSE)
        hazards[[k]] <- hazard_function(as.vector(pred), newTime)
      } else if (learner == "gam_s") {
        pred <- predict(fits[["gam_s"]], newdata = newdata, type = "response", se = FALSE)
        hazards[[k]] <- hazard_function(as.vector(pred), newTime)
      } else if (learner == "exp") {
        if(!exists("pos.pred")) {
          if(fits[["fit.pos"]] == 1) {
            pos.pred <- 1
          } else {
            pos.pred <- predict(fits[["fit.pos"]], newdata = X_row, type = 'response')
          } 
        }
        surv <- predict(fits[["exp"]], newdata = X_row, type = 'quantile', p = seq(0, .999, by=.001))
        surv <- pos.pred * (1-stats::approx(surv, seq(0, .999, by=.001), xout = newTime, method = 'linear', rule = 2)$y)  
        hazards[[k]] <- hazard_function(surv, newTime)
      } else if (learner == "weibull")  {
        if(!exists("pos.pred")) {
          if(fits[["fit.pos"]] == 1) {
            pos.pred <- 1
          } else {
            pos.pred <- predict(fits[["fit.pos"]], newdata = X_row, type = 'response')
          } }
        surv <- predict(fits[["weibull"]], newdata = X_row, type = 'quantile', p = seq(0, .999, by=.001))
        surv <- pos.pred * (1-stats::approx(surv, seq(0, .999, by=.001), xout = newTime, method = 'linear', rule = 2)$y)
        hazards[[k]] <- hazard_function(surv, newTime)
      } else if (learner == "loglog") {
        if(!exists("pos.pred")) {
          if(fits[["fit.pos"]] == 1) {
            pos.pred <- 1
          } else {
            pos.pred <- predict(fits[["fit.pos"]], newdata = X_row, type = 'response')
          } }
        surv <- predict(fits[["loglog"]], newdata = X_row, type = 'quantile', p = seq(0, .999, by=.001))
        surv <- pos.pred * (1-stats::approx(surv, seq(0, .999, by=.001), xout = newTime, method = 'linear', rule = 2)$y)
        hazards[[k]] <- hazard_function(surv, newTime)
      }  else if (learner == "aalen") {
        pred <- predict(fits[["aalen"]],  newdata=X_row, times=newTime)
        hazards[[k]] <- hazard_function(cummin(as.vector(pred$S0)), newTime)
      } else if (learner == "gam_splines") {
        pred <- predict(fits[["gam_splines"]],data.frame(Time= newTime, X_row) )
        hazards[[k]] <- hazard_function(pred$surv, newTime)
      } else if (learner =="rsf"){
        pred <- predict(fits[["rsf"]], newdata = X_row)
        surv_rsf <- approx(x = pred$time.interest, 
                           y = pred$survival, 
                           xout = newTime, 
                           method = "linear", 
                           rule = 2)$y
        hazards[[k]] <- hazard_function(surv_rsf, newTime)
      }  
    }
    return (hazards)
  }
  hazard_list <- lapply(seq_len(nrow(newX)), function(i) compute_hazard( newX[i, , drop = FALSE]))
  return(hazard_list)
}

# A temporary solution to fix the problem that some learners themselves cannot converge: Use 95% of the data to refit:
fix_fits <- function(fits_all, Time, E, X, ObsWeights, library, pch_breaks, k = 20) {
  check_learners <- setdiff(names(fits_all), c("fit.pos", "aalen", "km", "pch","true1","true2"))
  if (length(check_learners) == 0) {
    #cat("No learners to check.\n")
    return(fits_all)
  }

  problematic_learners <- check_learners[sapply(check_learners, function(learner) {
    fit <- fits_all[[learner]]
    is.null(fit) || any(is.na(fit$coefficients)) 
  })]
  
  if (length(problematic_learners) == 0) {
    #cat("All learners are fine.\n")
    return(fits_all)
  }
  
  #cat("Problematic learners:", paste(problematic_learners, collapse = ", "), "\n")
  
  n <- nrow(X)
  folds <- split(sample(1:n), rep(1:k, length.out = n))
  
  for (learner in problematic_learners) {
    for (i in 1:k) {
      training_indices <- unlist(folds[-i])

      train_Time <- Time[training_indices]
      train_E <- E[training_indices]
      train_X <- X[training_indices, , drop = FALSE]
      train_ObsWeights <- ObsWeights[training_indices]
      
      cat("Refitting learner:", learner, "using fold", i, "\n")
      refit <- try(SL_fit(
        Time = train_Time, 
        E = train_E, 
        X = train_X, 
        library = c(learner), 
        ObsWeights = train_ObsWeights, 
        pch_breaks = pch_breaks
      )[[learner]])
      
      if (!inherits(refit, "try-error") && !is.null(refit) && 
          all(!is.na(unlist(refit)))) {
        fits_all[[learner]] <- refit
        #cat("Successfully refitted learner:", learner, "\n")
        break
      }
      
      if (i == k) {
        warning("Failed to refit learner:", learner)
      }
    }
  }
  
  return(fits_all)
}



compute_hazard <- function( Time, E, X, folds, library, ObsWeights=NULL, pch_breaks =10 ,rsf_nsplit=3,rsf_ntree=100 ) {
  
  hazard_list <- vector("list", length = nrow(X))
  
  event_times <- sort(unique(Time)) #NewTime
  
  for (fold in seq_along(folds)) {
    test_idx <- folds[[fold]]
    train_idx <- setdiff(1:length(Time), test_idx)
    
    T_train <- Time[train_idx]
    E_train <- E[train_idx]
    X_train <- X[train_idx, , drop = FALSE]
    ObsWeights_train <- ObsWeights[train_idx]
    
    X_test <- X[test_idx, , drop = FALSE]
    T_test <- Time[test_idx]
    E_test <-E[test_idx]
    
    fits <- SL_fit(T_train,E_train,X_train,library,ObsWeights_train,pch_breaks,rsf_nsplit = rsf_nsplit,rsf_ntree = rsf_ntree)
    fits <- fix_fits(fits, T_train, E_train, X_train, ObsWeights, library, pch_breaks)
    hazard_list[test_idx] <- SL_pred(event_times,X_test,fits,library)
  }
  
  return(hazard_list)
}



loss <- function(params, hazard, E  ,Time, w_breaks ,p ) {
  k = length(w_breaks)-1
  event_times <- sort(unique(Time))

  w_vector <- vector("numeric", length = k*p)
  for (m in 1:k){
    w <- params[((m-1)*(p-1)+1) : (m*(p-1))]
    w <- softmax(w) 
    w <- expand(w)
    w_vector[((m-1)*p+1) : (m*p)] <- w
  }
  
  w_matrix <- matrix(w_vector,  ncol = k)
  
  n_samples <- length(E)
  n_times <- length(event_times)  
  

  LOSS <- sum(sapply(1:length(Time), function(i){
    
    hazards <- hazard[[i]]  
    hazards_matrix <- do.call(rbind, hazards)  
    hazard_mix <- numeric(ncol(hazards_matrix))
    for (j in seq_len(k)) {
      time_mask <- event_times >= w_breaks[j] & event_times < w_breaks[j + 1]
      idx <- which(time_mask)  
      if (length(idx) > 0) {
        hazard_mix[idx] <- colSums(w_matrix[, j] * hazards_matrix[, idx, drop = FALSE])
      }
    }
    Ei <- E[i]
    index <- which.min(abs(event_times - Time[i]))
    
    log_likelihood <- 0
    if (Ei == 1) {
      if (hazard_mix[index] == 0) {
        hazard_mix[index] <- 0.00001
      }
      log_likelihood <-   log(hazard_mix[index])
    }
    cumsum_hazard <- cumsum(hazard_mix* c(diff(event_times),0))
    integral_term <- cumsum_hazard[index]
    
    - (log_likelihood - integral_term)  
  }))
  
  return(LOSS)
}


# Include loss() k: breaks of w， p: num of learners
optimize <- function(hazard,E,Time ,w_breaks,p) {
  k = length(w_breaks)-1
  optim(
    par = rep(1, (p - 1)*k), 
    fn = loss,
    hazard =hazard,
    E =E,
    Time = Time,
    w_breaks = w_breaks,
    p = p,
    method = "L-BFGS-B",
    lower = rep(-7, (p - 1)*k),  
    upper = rep(7, (p - 1)*k),  
    #lower = 0,
    #upper = 1,
    control = list(fnscale = 1)
  )$par
}



compute_cif <- function(hazard_results, event_times) {
  n <- dim(hazard_results)[1]          
  m <- dim(hazard_results)[2]          
  J <- dim(hazard_results)[3]          
  
  # time grid: t0 = 0, t1,...,tm
  event_times0 <- c(0, event_times)     # 
  dt <- diff(event_times0)              #  t1 - t0, ..., tm - t_{m-1}
  
  cif_results <- array(0, dim = c(n, m, J))  #  n x m x J
  
  for (j in 1:n) {
    hazard_sample <- hazard_results[j, , ]   # m × J
    hazard_sample <- matrix(hazard_sample, nrow = m, ncol = J)
    
    # total hazard: h^1 + h^2 + ...
    total_hazard <- rowSums(hazard_sample)
    
    cumhaz_integrand <- total_hazard * dt
    cumhaz <- c(0, cumsum(cumhaz_integrand))  # shift right by one
    S <- exp(-cumhaz)   # S[1] = 1 (t0)
    
    for (k in 1:J) {
      h_k <-c(0, hazard_sample[, k])  # hazard at t1,...,tm，length m
      
      integrand <- 0.5 * (S[-(m+1)] * h_k[-m-1] + S[-1] * h_k[-1]) * dt  #length m
      
      cif_k <- cumsum(integrand)  
      
      cif_results[j, , k] <- cif_k
    }
  }
  
  return(cif_results)
}


predict_final <- function(Time,Ev, X, weights, newX,newX0=NULL,newX1=NULL,newTime,library, w_breaks, ObsWeights=NULL,  pch_breaks =10,rsf_nsplit=3, rsf_ntree=100) {
  num_events <- length(Ev)  
  n_X <- nrow(newX)           
  k <- length(w_breaks)-1
  
  if (is.null(newTime)) {
    event_times <- sort(unique(Time[do.call(pmax, Ev) > 0])) 
  } else {
    event_times <- newTime
  }
  
  hazard_list <- vector("list", length = n_X)   
  hazard_results <- array(0, dim = c(n_X, length(event_times), num_events))  
  
  
  if (is.null(newX0)){
  for (i in 1:num_events) {
    cat("Processing Event", i, "\n")
    weights_ev <- weights[[i]]  
    w_matrix <- matrix(weights_ev,  ncol = k)
    Ei <- Ev[[i]]                  
    
    fits_all <- SL_fit(Time = Time, E=Ei, X=X,library=library,ObsWeights=ObsWeights,
                       pch_breaks=pch_breaks, rsf_nsplit = rsf_nsplit, rsf_ntree = rsf_ntree)
    fits_all <- fix_fits(fits_all, Time, Ei, X, ObsWeights, library, pch_breaks)
   
     hazard_list  <-   SL_pred(newTime = event_times, newX = newX, 
                              fits = fits_all, library = library)
    for (j in seq_along(hazard_list)) {
      hazards <- hazard_list[[j]] 
      
      hazards_matrix <- do.call(rbind, hazards)  # n_learners x n_events matrix
      
      for (q in 1:k) {
        if (q == k) {
          time_mask <- event_times >= w_breaks[q]  
        } else {
          time_mask <- event_times >= w_breaks[q] & event_times < w_breaks[ q+1 ]
        }
        idx <- which(time_mask)  
        if (length(idx) > 0) {
          hazard_results[j,idx,i] <- colSums(w_matrix[, q] * hazards_matrix[, idx, drop = FALSE])
        }
      }
    }
  } 
  cif_results <- compute_cif(hazard_results, event_times)
  return(list(hazard = hazard_results, cif = cif_results, time = event_times))
  } else {
    if (is.null(newX1)) {
    n_X0 <- nrow(newX0)
    hazard_list0 <- vector("list", length = n_X0)   
    hazard_results0 <- array(0, dim = c(n_X0, length(event_times), num_events))  
    for (i in 1:num_events) {
      cat("Processing Event", i, "\n")
      weights_ev <- weights[[i]]  
      w_matrix <- matrix(weights_ev,  ncol = k)
      Ei <- Ev[[i]]                 
      fits_all <- SL_fit(Time = Time, E=Ei, X=X,library=library,ObsWeights=ObsWeights,
                         pch_breaks=pch_breaks, rsf_nsplit = rsf_nsplit, rsf_ntree = rsf_ntree)
      fits_all <- fix_fits(fits_all, Time, Ei, X, ObsWeights, library, pch_breaks)
      hazard_list  <-   SL_pred(newTime = event_times, newX = newX, 
                                fits = fits_all, library = library)
      hazard_list0  <-   SL_pred(newTime = event_times, newX = newX0, 
                                fits = fits_all, library = library)
      for (j in seq_along(hazard_list)) {
        hazards <- hazard_list[[j]]  
        hazards_matrix <- do.call(rbind, hazards) 
        for (q in 1:k) {
          if (q == k) {
            time_mask <- event_times >= w_breaks[q] 
          } else {
            time_mask <- event_times >= w_breaks[q] & event_times < w_breaks[ q+1 ]
          }
          idx <- which(time_mask)  
          if (length(idx) > 0) {
            hazard_results[j,idx,i] <- colSums(w_matrix[, q] * hazards_matrix[, idx, drop = FALSE])
          }
        }
      }
      for (j in seq_along(hazard_list0)) {
        hazards0 <- hazard_list0[[j]]  
        hazards_matrix0 <- do.call(rbind, hazards0)  
        
        for (q in 1:k) {
          if (q == k) {
            time_mask <- event_times >= w_breaks[q]  
          } else {
            time_mask <- event_times >= w_breaks[q] & event_times < w_breaks[ q+1 ]
          }
          idx <- which(time_mask)  
          if (length(idx) > 0) {
            hazard_results0[j,idx,i] <- colSums(w_matrix[, q] * hazards_matrix0[, idx, drop = FALSE])
          }
        }
      }
    } 
    
    cif_results <- compute_cif(hazard_results, event_times)
    cif_results0 <- compute_cif(hazard_results0, event_times)
    return(list(hazard1 = hazard_results, 
                hazard0 = hazard_results0, 
                cif1 = cif_results, 
                cif0 = cif_results0, 
                time = event_times))
  } else {
    n_X0 <- nrow(newX0)
    hazard_list0 <- vector("list", length = n_X0)   
    hazard_results0 <- array(0, dim = c(n_X0, length(event_times), num_events))  
    
    n_X1 <- nrow(newX1)
    hazard_list1 <- vector("list", length = n_X1)   
    hazard_results1 <- array(0, dim = c(n_X1, length(event_times), num_events)) 
    for (i in 1:num_events) {
      cat("Processing Event", i, "\n")
      weights_ev <- weights[[i]]  
      w_matrix <- matrix(weights_ev,  ncol = k)
      Ei <- Ev[[i]]                
      fits_all <- SL_fit(Time = Time, E=Ei, X=X,library=library,ObsWeights=ObsWeights,
                         pch_breaks=pch_breaks, rsf_nsplit = rsf_nsplit, rsf_ntree = rsf_ntree)
      fits_all <- fix_fits(fits_all, Time, Ei, X, ObsWeights, library, pch_breaks)
      
      hazard_list  <-   SL_pred(newTime = event_times, newX = newX, 
                                fits = fits_all, library = library)
      hazard_list0  <-   SL_pred(newTime = event_times, newX = newX0, 
                                fits = fits_all, library = library)
      hazard_list1  <-   SL_pred(newTime = event_times, newX = newX1, 
                                 fits = fits_all, library = library)
      for (j in seq_along(hazard_list)) {
        hazards <- hazard_list[[j]]  
        hazards_matrix <- do.call(rbind, hazards) 
        for (q in 1:k) {
          if (q == k) {
            time_mask <- event_times >= w_breaks[q] 
          } else {
            time_mask <- event_times >= w_breaks[q] & event_times < w_breaks[ q+1 ]
          }
          idx <- which(time_mask)  
          if (length(idx) > 0) {
            hazard_results[j,idx,i] <- colSums(w_matrix[, q] * hazards_matrix[, idx, drop = FALSE])
          }
        }
      }
      for (j in seq_along(hazard_list0)) {
        hazards0 <- hazard_list0[[j]]  
        hazards_matrix0 <- do.call(rbind, hazards0) 
        for (q in 1:k) {
          if (q == k) {
            time_mask <- event_times >= w_breaks[q]  
          } else {
            time_mask <- event_times >= w_breaks[q] & event_times < w_breaks[ q+1 ]
          }
          idx <- which(time_mask)  
          if (length(idx) > 0) {
            hazard_results0[j,idx,i] <- colSums(w_matrix[, q] * hazards_matrix0[, idx, drop = FALSE])
          }
        }
      }
      for (j in seq_along(hazard_list1)) {
        hazards1 <- hazard_list1[[j]]  
        hazards_matrix1 <- do.call(rbind, hazards1)  
        for (q in 1:k) {
          if (q == k) {
            time_mask <- event_times >= w_breaks[q]  
          } else {
            time_mask <- event_times >= w_breaks[q] & event_times < w_breaks[ q+1 ]
          }
          idx <- which(time_mask)  
          if (length(idx) > 0) {
            hazard_results1[j,idx,i] <- colSums(w_matrix[, q] * hazards_matrix1[, idx, drop = FALSE])
          }
        }
      }
    } 
    
    cif_results <- compute_cif(hazard_results, event_times)
    cif_results0 <- compute_cif(hazard_results0, event_times)
    cif_results1 <- compute_cif(hazard_results1, event_times)
    return(list(hazard = hazard_results, 
                hazard0 = hazard_results0,
                hazard1 = hazard_results1,
                cif = cif_results, 
                cif0 = cif_results0, 
                cif1 =cif_results1,
                time = event_times))
    } 
  }
}




SuperLearnerC <- function(Time, X, Event, NewX, NewX0=NULL, NewX1=NULL, V=3, k=1, w_breaks=NULL,library, pch_breaks=10, ObsWeights=NULL, rsf_nsplit=3, rsf_ntree=100) {
  
  n <- length(Time)
  folds <- split(sample(1:n), rep(1:V, length.out = n))
  Ev <-  create_event_data(Event)
  
  if (is.null(w_breaks)) {
    w_breaks <- seq(0, max(Time), length.out = k + 1)
  } else {
    k <- length(w_breaks)-1  
  }
  
  p <- length(library)
  weights <- list()
  if(p>1){
    for (j in 1: length(Ev)) {
      cat("Optimizing Event", j, "\n")
      hazard <- compute_hazard( Time=Time, 
                                E= Ev[[j]], 
                                X=X,
                                folds=folds,
                                library=library,
                                ObsWeights=ObsWeights, 
                                pch_breaks=pch_breaks,
                                rsf_nsplit= rsf_nsplit,
                                rsf_ntree = rsf_ntree)
      
      
      weight <- optimize(hazard= hazard,
                         E = Ev[[j]],
                         Time = Time, 
                         w_breaks = w_breaks,
                         p=p)
      
      w_vector <- vector("numeric", length = k*p)
      for (m in 1:k){
        w <- weight[((m-1)*(p-1)+1) : (m*(p-1))]
        
        if (length(w) == 1) {
          w <- expand(softmax(c(w)))  
        } else {
          w <- expand(softmax(w))    
        }
        
        w_vector[((m-1)*p+1) : (m*p)] <- w
      }
    
      weights[[j]] <- w_vector
    }
  }
  else
  {       w_vector <- 1 
  for (j in 1: length(Ev)) {
    weights[[j]] <- w_vector}}
  names(weights) <- paste("Event", 1:length(Ev))
  
  newTime <- sort(unique(Time))
  
  
  result <- predict_final(Time = Time, 
                          Ev,
                          X,
                          weights = weights,
                          newX = NewX,
                          newX0= NewX0,
                          newX1 = NewX1,
                          newTime =  newTime,
                          library,
                          w_breaks = w_breaks,
                          ObsWeights= ObsWeights,
                          pch_breaks=pch_breaks,
                          rsf_nsplit= rsf_nsplit,
                          rsf_ntree = rsf_ntree)
  
    final_result <- c(
      list(library = library),
      list(weights = weights),
      result)
  
  return(final_result)
}


