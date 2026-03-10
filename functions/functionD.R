


softmax <- function(params) {
  exp_params <- exp(params)
  denominator <- sum(exp_params) + 1
  weights <- exp_params / denominator
  additional_weight <- 1 / denominator
  result <- c(weights, additional_weight)
  return(result)
}

expand <- function(w) {
  w1 <- w * (1.001 + 0.001) - 0.001
  w1 <- pmax(w1,0)
  w2 <- w1 / sum(w1)
  return(w2)
}

stack_data <- function(time, event = NULL, indicators = NULL, covariates, time_intervals) {
  if (is.vector(covariates)) {
    covariates <- data.frame(covariates)
  }
  
  if (!is.null(event)) {
    stacked_data_list <- lapply(seq_along(time), function(i) {
      t_i <- time[i]
      event_i <- event[i]
      covariate_i <- covariates[i, , drop = FALSE]
      valid_intervals <- which(time_intervals[-length(time_intervals)] <= t_i)
      if (length(valid_intervals) == 0) {
        return(NULL)
      }
      rows <- length(valid_intervals)
      time_start <- time_intervals[valid_intervals]
      time_end <- time_intervals[valid_intervals + 1]
      event_status <- ifelse(time_end > t_i, event_i, 0)
      cbind(
        ID = i,  
        Time_interval = valid_intervals,
        Time_start = time_start,
        Time_end = time_end,
        Event = event_status,
        covariate_i[rep(1, rows), , drop = FALSE]
      )
    })
    stacked_data <- do.call(rbind, stacked_data_list)
    unique_events <- sort(unique(c(0, event)))
    for (e in unique_events) {
      stacked_data[[paste0("e", e)]] <- ifelse(stacked_data$Event == e, 1, 0)
    }
  } else if (!is.null(indicators)) {
    event_matrix <- as.matrix(indicators)
    e0 <- rowSums(event_matrix) == 0
    event_matrix <- cbind(e0 = as.integer(e0), event_matrix)
    colnames(event_matrix) <- paste0("e", 0:(ncol(event_matrix) - 1))
    stacked_data_list <- lapply(seq_along(time), function(i) {
      t_i <- time[i]
      covariate_i <- covariates[i, , drop = FALSE]
      indicator_row <- event_matrix[i, ]
      event_i <- if (all(indicator_row == 0)) { 0 } else { which(indicator_row == 1)[1] - 1 }
      valid_intervals <- which(time_intervals[-length(time_intervals)] <= t_i)
      rows <- length(valid_intervals)
      time_start <- time_intervals[valid_intervals]
      time_end <- time_intervals[valid_intervals + 1]
      event_status <- ifelse(time_end > t_i, event_i, 0)
      cbind(
        ID = i,  
        Time_interval = valid_intervals,
        Time_start = time_start,
        Time_end = time_end,
        Event = event_status,
        covariate_i[rep(1, rows), , drop = FALSE]
      )
    })
    stacked_data <- do.call(rbind, stacked_data_list)
    for (k in seq_len(ncol(event_matrix))) {
      stacked_data[[colnames(event_matrix)[k]]] <- ifelse(stacked_data$Event == (k - 1), 1, 0)
    }
  }
  stacked_data$Time_factor = factor(stacked_data$Time_interval)
  return(stacked_data)
}


survival <- function(h) {
  m <- nrow(h)
  S <- numeric(m)
  S[1] <- h[1,1]
  if (m >1)   {
    for (s in 2:m) {
      S[s] <- S[s - 1] * h[s, 1] 
  }
  return(S)
  }
}


cif_disc <- function(h) {
  S <- survival(h)

  m <- nrow(h)
  j <- ncol(h)

  cif <- matrix(0, nrow = m, ncol = j - 1)
  
  for (e in 2:j) { 
    for (t in 1:m) {
      cif[t, e - 1] <- sum(h[1:t, e] * c(1, S[1:t-1]))
    }
  }
  return(cif)
}


SL_disc_fit <- function(Time, levels, Ev, Indi, X, library,
                        VGAM.mn_maxit=200, vgam_maxit=200,
                        xgb_nrounds=200,  xgb_eta=0.1, xgb_maxdep=8,xgb_verbose=0,
                        knn_tuneLength=10, rf_tuneLength =5) {   
  fits <- list()
  
  covariates <- colnames(X)
  
  data <- data.frame(
    Time = factor(Time, levels=levels),
    Event = Ev,
    Indi,
    X
  )
  
  data$Event <- as.factor(data$Event)
  formula <- as.formula(paste("Event ~ Time +", paste(covariates, collapse = " + ")))
  for (learner in library) {
    if (learner == "VGAM.mn") {
      fits[["VGAM.mn"]] <- VGAM::vglm(
        formula = as.formula(paste("cbind(", paste(colnames(Indi), collapse = ", "), ") ~ -1 + Time +", paste(covariates, collapse = " + "))),
        data = data,
        family = VGAM::multinomial(refLevel = "e0"),
        maxit = VGAM.mn_maxit
      )
    } else if (learner == "vgam_s") {
      all_vars <- c("Time", covariates)
      smoothable_vars <- setdiff(
        all_vars[sapply(all_vars, function(v) length(unique(data[[v]])) >= 7)], "Time" 
      )
      linear_vars <- setdiff(all_vars, smoothable_vars)
      
      rhs_terms <- c(
        paste0("s(", smoothable_vars, ")"),  
        linear_vars          
      )
      vgam_formula <- as.formula(paste(
        "cbind(", paste(colnames(Indi), collapse = ", "), ") ~ ",
        paste(rhs_terms, collapse = " + ")
      ))
      fits[["vgam_s"]] <- VGAM::vgam(
        formula = vgam_formula,
        data = data,
        family = VGAM::multinomial(refLevel = "e0"),
        maxit = vgam_maxit
      )
    }  else if (learner == "nnet.mn") {
      fits[["nnet.mn"]] <- nnet::multinom(
        formula = formula,
        data = data,
        trace=FALSE
      )
    }  else if (learner == "lda") { 
      fits[["lda"]] <- lda(
        formula = formula,
        data = data
      )
    } else if (learner == "naiveBayes") {
      fits[["naiveBayes"]] <- naiveBayes(
        formula = formula,
        data = data
      )
    } else if (learner =="svm") {
      fits[["svm"]] <- svm(formula, data =data ,
                           probability = TRUE)
    } else if (learner == "xgboost"){
      formula2 = as.formula(paste("as.factor(Event) ~ Time +", 
                                  paste(covariates, collapse = " + "), "- 1"))
      model_matrix <- model.matrix(formula2, data = data,)
      fits[["xgboost"]] <- xgboost(data = model_matrix, label = as.numeric(data$Event)-1,
                                   objective = "multi:softprob",
                                   num_class = length(unique(Ev)),                
                                   nrounds = xgb_nrounds,
                                   eta = xgb_eta,
                                   max_depth = xgb_maxdep,
                                   verbose = xgb_verbose)
    } else if (learner =="qda" ) {
      fits[["qda"]] <- qda(
        formula = formula,
        data = data
      )
    } else if (learner =="knn" ) {
      train_control <- trainControl(method = "cv", number = 10)  # 10 folds

      fits[["knn"]] <- train(formula, data = data,
                             method = "knn",
                             tuneLength = knn_tuneLength,  
                             trControl = train_control)
    } else if (learner =="rf" ) {
      train_control <- trainControl(method = "cv", number = 10)  
      fits[["rf"]] <- train(formula, data = data,
                            method = "rf",
                            tuneLength = rf_tuneLength,  
                            trControl = train_control)
    } 
  }
  return(fits)
}


SL_disc_pred  <- function (newTime, levels, newX, fits,  library) {
  
  compute_hazard <- function(Time, X_row) {
    new_data <- data.frame(Time = factor(Time, levels= levels), X_row)
    hazards <- vector("list", length(library))
    for (k in seq_along(library)) {
      learner <- library[k]
      if (learner == "VGAM.mn") {
        hazards[[k]] <- VGAM::predictvglm(fits[["VGAM.mn"]],
                                          newdata = new_data,
                                          type = "response")
      }  else if (learner == "vgam_s") {
        hazards[[k]] <- VGAM::predict(fits[["vgam_s"]],
                                      newdata = new_data,
                                      type = "response")
      } else if (learner == "nnet.mn") {
        hazards[[k]] <- t(as.matrix(predict(fits[["nnet.mn"]], newdata= new_data,type="prob")))
      } else if (learner == "cv.glmnet.mn") {
        new_x <- model.matrix(~ Time + . -1, data = new_data)  
        hazards[[k]] <- t(as.matrix(predict(fits[["cv.glmnet.mn"]], newx = new_x, type="response")))
      } else if (learner == "lda") {
        predictions <- predict(fits[["lda"]], newdata = new_data)
        hazards[[k]]  <- predictions$posterior
      } else if (learner == "naiveBayes") {
        hazards[[k]] <- predict(fits[["naiveBayes"]], newdata = new_data, type = "raw")
      } else if (learner == "svm") {
        h_svm <- predict(fits[["svm"]], newdata = new_data,  probability = TRUE)
        hazards[[k]] <- attr(h_svm, "probabilities")
      } else if (learner == "xgboost") {
        new_matrix <- model.matrix(~ . - 1, data = new_data)
        h_xgb <- predict(fits[["xgboost"]], newdata = new_matrix)
        hazards[[k]] <- matrix(h_xgb, ncol = fits[["xgboost"]]$params$num_class, byrow = TRUE)
      } else if (learner == "qda" ) {
        predictions <- predict(fits[["qda"]], newdata = new_data)
        hazards[[k]]  <- predictions$posterior
      } else if (learner == "knn" ) {
        hazards[[k]]  <- as.matrix(predict(fits[["knn"]], newdata = new_data,type = "prob"))
      } else if (learner == "rf" ) {
        hazards[[k]]  <- as.matrix(predict(fits[["rf"]], newdata = new_data,type = "prob"))
      } 
    }
    
    return(hazards)
  }
  
  hazard_list <- lapply(seq_len(nrow(newX)), function(i)
    compute_hazard( factor(newTime[i],levels=levels), newX[i, , drop = FALSE]))
  
  return(hazard_list)
  
}



hazard_disc <- function(stacked_data, folds, library,
                        VGAM.mn_maxit=200,
                        xgb_nrounds=200,  xgb_eta=0.1, xgb_maxdep=8,xgb_verbose=0,
                        knn_tuneLength=10, rf_tuneLength =5) {
  
  hazard_list_total <- vector("list", length = nrow(stacked_data))
  
  indicator_columns <- grep("^e[0-9]+$", colnames(stacked_data), value = TRUE)
  covariate_columns <- colnames(stacked_data)[!colnames(stacked_data) %in% c("ID", "Time_interval", "Time_start",  "Time_end", "Event","Time_factor", indicator_columns)]
  levels <- levels(stacked_data$Time_factor)

  for (fold in seq_along(folds)) {
    test_idx <- folds[[fold]]
    train_idx <- setdiff(unique(stacked_data$ID), test_idx)
    
    # Saparate data according to ID
    train_data <- stacked_data[stacked_data$ID %in% train_idx, ]
    test_data <- stacked_data[stacked_data$ID %in% test_idx, ]
    
    max_train_time <- max(train_data$Time_interval)
    
    test_data <- test_data[test_data$Time_interval <= max_train_time, ]
    
    T_train <- factor(train_data$Time_interval, levels = levels)
    E_train <- train_data$Event
    Ind_train <- train_data[, indicator_columns, drop = FALSE]
    X_train <- train_data[, covariate_columns, drop = FALSE]
    
    T_test <- factor(test_data$Time_interval, levels = levels)
    Ind_test <- test_data[, indicator_columns, drop = FALSE]
    X_test <- test_data[, covariate_columns, drop = FALSE]
    
    fits <- SL_disc_fit(T_train,levels, E_train, Ind_train, X_train, library,VGAM.mn_maxit=VGAM.mn_maxit, xgb_nrounds=xgb_nrounds,  xgb_eta=xgb_eta, xgb_maxdep=xgb_maxdep,xgb_verbose=xgb_verbose, knn_tuneLength=knn_tuneLength, rf_tuneLength =rf_tuneLength)
    
    hazards <- SL_disc_pred(T_test,levels, X_test, fits, library)
    hazard_list_total[which(stacked_data$ID %in% test_idx)] <- hazards
    
  }
  
  return(hazard_list_total)
}


loss_disc <- function(params, hazard_list_total, Ind,Time,w_breaks,p) {
  k = length(w_breaks)-1
  w_vector <- vector("numeric", length = k*p)
  for (m in 1:k){
    w <- params[((m-1)*(p-1)+1) : (m*(p-1))]
    w <- softmax(w) 
    w <- expand(w)
    w_vector[((m-1)*p+1) : (m*p)] <- w
  }
  w_matrix <- matrix(w_vector,  ncol = k)
  
  n_samples <- length(hazard_list_total)
  n_events <- ncol(Ind)  
  
  hazard_mix <- matrix(0, nrow = n_samples, ncol = n_events)
  
  for (i in seq_along(hazard_list_total)) {
    hazards <- hazard_list_total[[i]]  
    hazards_matrix <- do.call(rbind, hazards)  # n_learners x n_events 
    
    for (j in seq_len(k)) {
      time_mask <- Time[i] >= w_breaks[j] & Time[i] < w_breaks[j + 1] 
      if (time_mask) {  
        hazard_mix[i,] <- colSums(w_matrix[, j] * hazards_matrix)
      }
    }
  }
  
  hazard_mix[hazard_mix == 0] <- 1e-6
  cross_entropy_loss <- -mean(rowSums(Ind * log(hazard_mix)))
  return(cross_entropy_loss)
}


optimize_disc <- function(hazard_list_total,  Ind,Time,w_breaks, p) {
  k = length(w_breaks)-1
  optim(
    par = rep(1, (p - 1)*k),  
    fn = loss_disc,
    hazard_list_total= hazard_list_total,
    Ind = Ind,
    #ID = ID,
    Time = Time,
    w_breaks = w_breaks,
    p = p,
    method = "L-BFGS-B",
    lower = rep(-7, p - 1),  
    upper = rep(7, p - 1),  
    control = list(fnscale = 1)
  )$par
}




predict_disc <- function(stacked_data, library,
                         weights, newTime, newX,newX0=NULL, newX1 = NULL,
                         w_breaks,
                         VGAM.mn_maxit = 200,
                         xgb_nrounds = 200, xgb_eta = 0.1, xgb_maxdep = 8, xgb_verbose = 0,
                         knn_tuneLength=10, rf_tuneLength =5) {
  
  n_X <- nrow(newX)           
  n_time <- length(newTime)
  n_event <- length(unique(stacked_data$Event))
  
  k <- length(w_breaks) - 1
  w_matrix <- matrix(weights, ncol = k)
  
  hazard_results <- array(0, dim = c(n_X, n_time, n_event))  
  
  indicator_columns <- grep("^e[0-9]+$", colnames(stacked_data), value = TRUE)
  Ind <- as.matrix(stacked_data[, indicator_columns, drop = FALSE])
  covariate_columns <- colnames(stacked_data)[!colnames(stacked_data) %in% c("ID", "Time_interval", "Time_start",  "Time_end", "Event","Time_factor" ,indicator_columns)]
  levels <- levels(stacked_data$Time_factor)

  fits <- SL_disc_fit(Time=stacked_data$Time_factor, 
                      levels,
                      Ev=stacked_data$Event, 
                      Ind = Ind, 
                      X=stacked_data[, covariate_columns, drop = FALSE] , 
                      library = library, 
                      VGAM.mn_maxit = VGAM.mn_maxit,
                      xgb_nrounds = xgb_nrounds, xgb_eta = xgb_eta, 
                      xgb_maxdep = xgb_maxdep, xgb_verbose = xgb_verbose,
                      knn_tuneLength=knn_tuneLength, rf_tuneLength =rf_tuneLength)
  

  compute_hazard <- function(X_row) {
    rownames(X_row) <- NULL  
    X_expanded <- X_row[rep(1, length(newTime)), , drop = FALSE]
    newTime_factor <- factor(newTime, levels = levels)
    hazard_list <- SL_disc_pred(newTime_factor,levels, X_expanded, fits, library)
    hazard_mix <- array(0, dim = c(n_time, n_event))
    for (l in 1:n_time) {
      hazards <- hazard_list[[l]]
      hazards_matrix <- do.call(rbind, hazards)  
   
      for (j in 1:k) {
        if (j == k) {
          time_mask <- newTime[l] >= w_breaks[j]
        } else {
          time_mask <- newTime[l] >= w_breaks[j] & newTime[l] < w_breaks[j + 1]
        }
        
        if (time_mask) {  # 如果当前时间点落在当前区间内
          # 直接加权
          hazard_mix[l,] <- colSums(w_matrix[, j] * hazards_matrix)
        }
      }
    }
    return(hazard_mix)
  }
  
  if(is.null(newX0)) {
  for (a in 1:n_X) {
    hazard_results[a, , ] <- compute_hazard(newX[a, , drop = FALSE])
  }
  
  cif <- array(0, dim = c(n_X, n_time, n_event-1)) 
  for (c in 1:n_X){
    hazard <- hazard_results[c,,]
    cif[c,,] <- cif_disc(hazard)
  }

    
  return (list(hazard = hazard_results, cif = cif))
  
} else{
  if(is.null(newX1)) {
  n_X0 <- nrow(newX0) 
  hazard_results0 <- array(0, dim = c(n_X0, n_time, n_event))
  for (a in 1:n_X) {
    hazard_results[a, , ] <- compute_hazard(newX[a, , drop = FALSE])
  }
  for (a in 1:n_X0) {
    hazard_results0[a, , ] <- compute_hazard(newX0[a, , drop = FALSE])
  }
  
  cif1 <- array(0, dim = c(n_X, n_time, n_event-1)) 
  cif0 <- array(0, dim = c(n_X0, n_time, n_event-1)) 
  for (c in 1:n_X){
    hazard1 <- hazard_results[c,,]
    cif1[c,,] <- cif_disc(hazard1)
  }
  for (c in 1:n_X0){
    hazard0 <- hazard_results0[c,,]
    cif0[c,,] <- cif_disc(hazard0)
  }
  
  return (list(hazard1 = hazard_results,
               hazard0 = hazard_results0,
               cif1 = cif1,
               cif0 = cif0))
  }
  else {
    n_X0 <- nrow(newX0) 
    n_X1 <- nrow(newX1) 
    hazard_results0 <- array(0, dim = c(n_X0, n_time, n_event))
    hazard_results1 <- array(0, dim = c(n_X1, n_time, n_event))
    for (a in 1:n_X) {
      hazard_results[a, , ] <- compute_hazard(newX[a, , drop = FALSE])
    }
    for (a in 1:n_X0) {
      hazard_results0[a, , ] <- compute_hazard(newX0[a, , drop = FALSE])
    }
    for (a in 1:n_X1) {
      hazard_results1[a, , ] <- compute_hazard(newX1[a, , drop = FALSE])
    }
    cif <- array(0, dim = c(n_X, n_time, n_event-1)) 
    cif0 <- array(0, dim = c(n_X0, n_time, n_event-1)) 
    cif1 <- array(0, dim = c(n_X1, n_time, n_event-1)) 
    
    for (c in 1:n_X){
      hazard <- hazard_results[c,,]
      cif[c,,] <- cif_disc(hazard)
    }
    for (c in 1:n_X0){
      hazard0 <- hazard_results0[c,,]
      cif0[c,,] <- cif_disc(hazard0)
    }
    for (c in 1:n_X1){
      hazard1 <- hazard_results1[c,,]
      cif1[c,,] <- cif_disc(hazard1)
    }
    return (list(hazard = hazard_results,
                 hazard0 = hazard_results0,
                 hazard1 = hazard_results1,
                 cif = cif,
                 cif0 = cif0,
                 cif1 = cif1 ))
    }
  }
}


SuperLearnerD <- function( Time,X,Event=NULL,Event_indicator=NULL,time_intervals,
                           NewX, NewX0=NULL, NewX1 = NULL,
                           V=3, library,k=1,w_breaks=NULL,
                           VGAM.mn_maxit=200,
                           xgb_nrounds=200, xgb_eta=0.1, xgb_maxdep=8, xgb_verbose=0,
                           knn_tuneLength= 10, rf_tuneLength =5) {
  
  
  cat("Preparing...", "\n")
  
  n <- length(Time)
  folds <- split(sample(1:n), rep(1:V, length.out = n))
  p <- length(library)
  
  stacked_data <- stack_data(time = Time,
                             event = Event, 
                             covariates= X, 
                             time_intervals = time_intervals) 

  
  if(is.null(w_breaks) ) {
    w_breaks <- seq(0, max(stacked_data$Time_interval), length.out = k + 1)
  } else {
    k = length(w_breaks)-1
  }
  q <- max(stacked_data$Time_interval)
  
  if (p > 1){
    indicator_columns <- grep("^e[0-9]+$", colnames(stacked_data), value = TRUE)
    Ind <- as.matrix(stacked_data[, indicator_columns, drop = FALSE])
  
    cat("CV fitting and computing hazards...", "\n")
    hazards <- hazard_disc(stacked_data= stacked_data,
                           folds =folds,
                           library= library,
                           VGAM.mn_maxit=VGAM.mn_maxit,
                           xgb_nrounds=xgb_nrounds,  xgb_eta=xgb_eta, 
                           xgb_maxdep=xgb_maxdep, xgb_verbose=xgb_verbose,
                           knn_tuneLength= knn_tuneLength, rf_tuneLength= rf_tuneLength)
    
    cat("Optimizing...", "\n")
    
    weight <- optimize_disc(hazards, Ind, 
                            Time=stacked_data$Time_interval ,
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
  } else {
    w_vector <- 1
  }
  
  
  newTime <- seq(1:q)
  cat("Fitting with optimal weights...", "\n")

  result <- predict_disc(stacked_data=stacked_data, 
                         library=library,
                         weights= w_vector, 
                         newTime= newTime, 
                         newX=NewX ,
                         newX0 =NewX0,
                         newX1 = NewX1,
                         w_breaks=w_breaks,
                         VGAM.mn_maxit=VGAM.mn_maxit,
                         xgb_nrounds=xgb_nrounds,  xgb_eta=xgb_eta, 
                         xgb_maxdep=xgb_maxdep, xgb_verbose=xgb_verbose)
  
  final_result <- c(list(library = library), list(weights = w_vector), result,
                    list(time = newTime)) 
  
  return(final_result)
  
}

