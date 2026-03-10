#Some functions for data generation and graphing in the simulation study and real data analysis

# Date without frailty
W2data2 <- function(n) {
  X <- rnorm(n, mean = 0.5, sd = 0.2)
  
  p_treat <- 0.7 
  Treat <- rbinom(n, size = 1, prob = p_treat)
  
  shape <- 5
  scale <- 1
  beta1 <- 3
  beta2 <- 2
  beta_t1 <- 1.5
  beta_t2 <- -1.0
  
  scale1 <- scale * exp(-(beta1 * X + beta_t1 * Treat) / shape)
  scale2 <- scale * exp(-(beta2 * X + beta_t2 * Treat) / shape)
  
  time_state1 <- rweibull(n, shape = shape, scale = scale1)
  time_state2 <- rweibull(n, shape = shape, scale = scale2)
  
  censoring_time <- rweibull(n, shape = 5, scale = 0.75)
  
  time <- pmin(time_state1, time_state2)
  status <- ifelse(time_state1 < time_state2, 1, 2)
  observed_time <- pmin(time, censoring_time)
  status <- ifelse(time > censoring_time, 0, status)
  
  data <- data.frame(Time = observed_time, Event = status, Z = X, A = Treat)
  return(data)
}


## Dependent Frailty
W2dataU_D <- function(n) {
  X <- rnorm(n, mean = 0.5, sd = 0.2)
  e <- rnorm(n, mean = 0, sd =0.1)
  Y <-  X^3 + X^2 + e
  
  p_treat <- 0.7 
  Treat <- rbinom(n, size = 1, prob = p_treat)

  shape <- 5
  scale <- 1
  beta1 <- 3
  beta2 <- 2
  beta_t1 <- 1.5
  beta_t2 <- -1.0
  betay1 <- -1
  betay2 <- 2

  scale1 <- scale * exp(-(beta1 * X + beta_t1 * Treat + betay1 * Y ) / shape)
  scale2 <- scale * exp(-(beta2 * X + beta_t2 * Treat + betay2 * Y ) / shape)

  time_state1 <- rweibull(n, shape = shape, scale = scale1)
  time_state2 <- rweibull(n, shape = shape, scale = scale2)

  censoring_time <- rweibull(n, shape = 5, scale = 0.75)

  time <- pmin(time_state1, time_state2)
  status <- ifelse(time_state1 < time_state2, 1, 2)
  observed_time <- pmin(time, censoring_time)
  status <- ifelse(time > censoring_time, 0, status)

  data <- data.frame(Time=observed_time,Event= status,Z= X,A= Treat)
  return(data)
}

# High-dimensional sparse data
W2dataU_D_H <- function(n) {
  p=100
  Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  X <- rnorm(n, mean = 0.5, sd = 0.2)
  e <- rnorm(n, mean = 0, sd =0.1)
  Y <-  X^3 + X^2 + e
  
  Z[,1] <- X
  
  p_treat <- 0.7 
  Treat <- rbinom(n, size = 1, prob = p_treat)

  shape <- 5
  scale <- 1
  beta1 <- 3
  beta2 <- 2
  beta_t1 <- 1.5
  beta_t2 <- -1.0
  betay1 <- -1
  betay2 <- 2

  scale1 <- scale * exp(-(beta1 * X + beta_t1 * Treat + betay1 * Y ) / shape)
  scale2 <- scale * exp(-(beta2 * X + beta_t2 * Treat + betay2 * Y ) / shape)

  time_state1 <- rweibull(n, shape = shape, scale = scale1)
  time_state2 <- rweibull(n, shape = shape, scale = scale2)
  
  censoring_time <- rweibull(n, shape = 5, scale = 0.75)

  time <- pmin(time_state1, time_state2)
  status <- ifelse(time_state1 < time_state2, 1, 2)
  observed_time <- pmin(time, censoring_time)
  status <- ifelse(time > censoring_time, 0, status)
  
  data <- data.frame(Time=observed_time,Event= status,Z,A= Treat)
  return(data)
}

cif_fit_n <- function(func) {
  cif_data <- func(300000)
  cif_fit <- cuminc(ftime = cif_data$Time, fstatus = cif_data$Event,
                    cencode=0)
  return (cif_fit)
}

compute_hazard_from_CHF <- function(chf_matrix, times) {
  dt <- diff(times)
  hazard_matrix <- t(apply(chf_matrix, 1, function(chf_row) {
    diff_chf <- diff(chf_row)
    hazard <- diff_chf / dt
    c(hazard, 0)
  }))
  return(hazard_matrix)
}

#CIF
drawSL <- function(func, result) {
  
  cif_fit <- cif_fit_n(func)
  ###min max time
  max_times <- sapply(result, function(x) max(x$time))
  time_max_common <- min(max_times)
  common_time <- seq(0, time_max_common, by = 0.01)
  n_time <- length(common_time)
  B <- length(result)
  n_event <- dim(result[[1]]$cif)[2]
  step_interp <- function(orig_time, orig_val, new_time) {
    approx(orig_time, orig_val, xout = new_time,
           method = "constant", f = 1, rule = 2)$y
  }
  
  all_CIF      <- array(NA, dim = c(B, n_time, n_event))
  for (i in 1:B) {
    for (k in 1:n_event) {
      all_CIF[i, , k]       <- step_interp(result[[i]]$time, result[[i]]$cif[ , k],       common_time)
    }
  }
  
 
  time_seq <- common_time

  avg_CIF <- apply(all_CIF, c(2, 3), mean)
  ci_CIF  <- apply(all_CIF, c(2, 3), function(x) quantile(x, probs = c(0.025, 0.975)))

  avg_CIF_1_SL <- avg_CIF[, 1]
  avg_CIF_2_SL <- avg_CIF[, 2]
  
  ci_CIF_1_SL <- ci_CIF[, , 1]  # (2 x time_points)
  ci_CIF_2_SL <- ci_CIF[, , 2]  # (2 x time_points)
  
  plot(cif_fit, xlab = "t", ylab = "CIF",
       main = "Comparison", 
       col = c("darkblue", "red"),
       lwd = c(2, 2),
       lty = c(1, 1)
  )
  
  polygon(c(time_seq, rev(time_seq)), 
          c(ci_CIF_1_SL[1, ], rev(ci_CIF_1_SL[2, ])),
          col = rgb(0, 1, 0, 0.2), border = NA)
  
  polygon(c(time_seq, rev(time_seq)), 
          c(ci_CIF_2_SL[1, ], rev(ci_CIF_2_SL[2, ])),
          col = rgb(1, 0, 1, 0.2), border = NA)

  lines(time_seq, avg_CIF_1_SL, col = "green", lwd = 2, lty = 2)
  lines(time_seq, avg_CIF_2_SL, col = "purple", lwd = 2, lty = 2)
  
  lines(time_seq, ci_CIF_1_SL[1, ], col = "darkgreen", lty = 3, lwd = 1.5)
  lines(time_seq, ci_CIF_1_SL[2, ], col = "darkgreen", lty = 3, lwd = 1.5)
  
  lines(time_seq, ci_CIF_2_SL[1, ], col = "darkblue", lty = 3, lwd = 1.5)
  lines(time_seq, ci_CIF_2_SL[2, ], col = "darkblue", lty = 3, lwd = 1.5)

  legend("topleft", legend = c("CIF1 (Theoretical)", "CIF2 (Theoretical)", 
                               "CIF1 (SL)", "CIF2 (SL)",
                               "SL CIF1 95% interval", "SL CIF2 95% interval"),
         col = c("darkblue", "red", "green", "purple",
                 rgb(0, 1, 0, 0.2), rgb(1, 0, 1, 0.2)),
         lwd = c(2, 2, 2, 2, NA, NA),
         lty = c(1, 1, 2, 2, NA, NA),
         fill = c(NA, NA, NA, NA, 
                  rgb(0, 1, 0, 0.2), rgb(1, 0, 1, 0.2)))
}

#result without theoretical
drawSL_NULL <- function(result){
  
  ###min max time
  max_times <- sapply(result, function(x) max(x$time))
  time_max_common <- min(max_times)
  common_time <- seq(0, time_max_common, by = 0.01)
  n_time <- length(common_time)
  B <- length(result)
  n_event <- dim(result[[1]]$cif)[2]
  step_interp <- function(orig_time, orig_val, new_time) {
    approx(orig_time, orig_val, xout = new_time,
           method = "constant", f = 1, rule = 2)$y
  }
  
  all_CIF <- array(NA, dim = c(B, n_time, n_event))
  
  for (i in 1:B) {
    for (k in 1:n_event) {
      all_CIF[i, , k] <- step_interp(result[[i]]$time, result[[i]]$cif[ , k], common_time)
      
    }
  }
  
  
  time_seq <- common_time
  # 计算平均和置信区间
  avg_CIF <- apply(all_CIF, c(2, 3), mean)
  ci_CIF  <- apply(all_CIF, c(2, 3), function(x) quantile(x, probs = c(0.025, 0.975)))
  
  # 分别提取每个事件的平均值和置信区间
  avg_CIF_1_SL <- avg_CIF[, 1]
  avg_CIF_2_SL <- avg_CIF[, 2]
  
  ci_CIF_1_SL <- ci_CIF[, , 1]  # (2 x time_points)
  ci_CIF_2_SL <- ci_CIF[, , 2]  # (2 x time_points)
  
  # 开始绘图
  plot(NULL, xlab = "t", ylab = "CIF",
       xlim = c(0,1),
       ylim = c(0,1),
       main = "SL", 
       
  )
  
  # 添加置信区间的阴影区域
  polygon(c(time_seq, rev(time_seq)), 
          c(ci_CIF_1_SL[1, ], rev(ci_CIF_1_SL[2, ])),
          col = rgb(0, 1, 0, 0.2), border = NA)
  
  polygon(c(time_seq, rev(time_seq)), 
          c(ci_CIF_2_SL[1, ], rev(ci_CIF_2_SL[2, ])),
          col = rgb(1, 0, 1, 0.2), border = NA)
  
  # 添加 Super Learner 的 CIF 平均曲线
  lines(time_seq, avg_CIF_1_SL, col = "green", lwd = 2, lty = 2)
  lines(time_seq, avg_CIF_2_SL, col = "purple", lwd = 2, lty = 2)
  
  # 添加置信区间边界
  lines(time_seq, ci_CIF_1_SL[1, ], col = "darkgreen", lty = 3, lwd = 1.5)
  lines(time_seq, ci_CIF_1_SL[2, ], col = "darkgreen", lty = 3, lwd = 1.5)
  
  lines(time_seq, ci_CIF_2_SL[1, ], col = "darkblue", lty = 3, lwd = 1.5)
  lines(time_seq, ci_CIF_2_SL[2, ], col = "darkblue", lty = 3, lwd = 1.5)
  
  # 图例
  legend("topleft", legend = c(
    "CIF1 (SL)", "CIF2 (SL)",
    "SL CIF1 95% interval", "SL CIF2 95% interval"),
    col = c( "green", "purple",
             rgb(0, 1, 0, 0.2), rgb(1, 0, 1, 0.2)),
    lwd = c( 2, 2, NA, NA),
    lty = c( 2, 2, NA, NA),
    fill = c( NA, NA, 
              rgb(0, 1, 0, 0.2), rgb(1, 0, 1, 0.2)))
  
}


drawMLM <- function(func, result){
  
  cif_fit <- cif_fit_n(func)
  ###min max time
  max_times <- sapply(result, function(x) max(x$time))
  time_max_common <- min(max_times)
  common_time <- seq(0, time_max_common, by = 0.01)
  n_time <- length(common_time)
  B <- length(result)
  n_event <- dim(result[[1]]$cif)[2]
  step_interp <- function(orig_time, orig_val, new_time) {
    approx(orig_time, orig_val, xout = new_time,
           method = "constant", f = 1, rule = 2)$y
  }
  
  all_CIF      <- array(NA, dim = c(B, n_time, n_event))
  
  for (i in 1:B) {
    for (k in 1:n_event) {
      all_CIF[i, , k]       <- step_interp(result[[i]]$time, result[[i]]$cif[ , k],       common_time)
      
    }
  }
  
  time_seq <- common_time
  avg_CIF <- apply(all_CIF, c(2, 3), mean)
  ci_CIF  <- apply(all_CIF, c(2, 3), function(x) quantile(x, probs = c(0.025, 0.975)))
  
  avg_CIF_1_SL <- avg_CIF[, 1]
  avg_CIF_2_SL <- avg_CIF[, 2]
  
  ci_CIF_1_SL <- ci_CIF[, , 1]  # (2 x time_points)
  ci_CIF_2_SL <- ci_CIF[, , 2]  # (2 x time_points)
  
  plot(cif_fit, xlab = "t", ylab = "CIF",
       main = "Comparison", 
       col = c("darkblue", "red"),
       lwd = c(2, 2),
       lty = c(1, 1)
  )
  
  polygon(c(time_seq, rev(time_seq)), 
          c(ci_CIF_1_SL[1, ], rev(ci_CIF_1_SL[2, ])),
          col = rgb(0, 1, 0, 0.2), border = NA)
  
  polygon(c(time_seq, rev(time_seq)), 
          c(ci_CIF_2_SL[1, ], rev(ci_CIF_2_SL[2, ])),
          col = rgb(1, 0, 1, 0.2), border = NA)
  
  lines(time_seq, avg_CIF_1_SL, col = "green", lwd = 2, lty = 2)
  lines(time_seq, avg_CIF_2_SL, col = "purple", lwd = 2, lty = 2)
  
  lines(time_seq, ci_CIF_1_SL[1, ], col = "darkgreen", lty = 3, lwd = 1.5)
  lines(time_seq, ci_CIF_1_SL[2, ], col = "darkgreen", lty = 3, lwd = 1.5)
  
  lines(time_seq, ci_CIF_2_SL[1, ], col = "darkblue", lty = 3, lwd = 1.5)
  lines(time_seq, ci_CIF_2_SL[2, ], col = "darkblue", lty = 3, lwd = 1.5)
  
  legend("topleft", legend = c("CIF1 (Theoretical)", "CIF2 (Theoretical)", 
                               "CIF1 (MLM)", "CIF2 (MLM)",
                               "CIF1 95% interval", "CIF2 95% interval"),
         col = c("darkblue", "red", "green", "purple",
                 rgb(0, 1, 0, 0.2), rgb(1, 0, 1, 0.2)),
         lwd = c(2, 2, 2, 2, NA, NA),
         lty = c(1, 1, 2, 2, NA, NA),
         fill = c(NA, NA, NA, NA, 
                  rgb(0, 1, 0, 0.2), rgb(1, 0, 1, 0.2)))
}

drawState <- function(func, result_list) {
  cif_fit <- cif_fit_n(func = func) 
  num_simulations <- length(result_list)  
  max_times <- sapply(result_list, function(x) max(x$time))
  time_max_common <- min(max_times)
  time_grid <- seq(0, time_max_common, by = 0.01)
  
  num_events <- dim(result_list[[1]]$cif)[2]  
  time_points <- length(time_grid)           

  interpolated_CIF <- array(NA, dim = c(num_simulations, time_points, num_events))
  
  for (i in 1:num_simulations) {
    interpolated_CIF[i, , ] <- sapply(1:num_events, function(event) {
      step_fun <- stepfun(result_list[[i]]$times, c(0, result_list[[i]]$cif[, event]))
      step_fun(time_grid)  
    })
  }

  avg_CIF <- apply(interpolated_CIF, c(2, 3), mean, na.rm = TRUE)  #  (time_points x num_events)
  ci_CIF <- apply(interpolated_CIF, c(2, 3), function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))  #  (2 x time_points x num_events)
  

  avg_CIF_1_SL <- avg_CIF[, 1]
  avg_CIF_2_SL <- avg_CIF[, 2]
  
  ci_CIF_1_SL <- ci_CIF[, , 1]  # (2 x time_points)
  ci_CIF_2_SL <- ci_CIF[, , 2]  # (2 x time_points)
  
  time_seq <- time_grid
  
  plot(cif_fit, xlab = "t", ylab = "CIF",
       main = "Comparison", 
       col = c("darkblue", "red") ,
       lwd =c(2,2),
       lty =c(1,1)
  )  

  polygon(c(time_seq, rev(time_seq)), 
          c(ci_CIF_1_SL[1, ], rev(ci_CIF_1_SL[2, ])),
          col = rgb(0, 1, 0, 0.2), border = NA)  
  polygon(c(time_seq, rev(time_seq)), 
          c(ci_CIF_2_SL[1, ], rev(ci_CIF_2_SL[2, ])),
          col = rgb(1, 0, 1, 0.2), border = NA) 
  
  lines(time_seq, avg_CIF_1_SL, col = "green", lwd = 2, lty = 2)
  lines(time_seq, avg_CIF_2_SL, col = "purple", lwd = 2, lty = 2)
  
  lines(time_seq, ci_CIF_1_SL[1, ], col = "darkgreen", lty = 3, lwd = 1.5)  # 上边界
  lines(time_seq, ci_CIF_1_SL[2, ], col = "darkgreen", lty = 3, lwd = 1.5)  # 下边界
  
  lines(time_seq, ci_CIF_2_SL[1, ], col = "darkblue", lty = 3, lwd = 1.5)  # 上边界
  lines(time_seq, ci_CIF_2_SL[2, ], col = "darkblue", lty = 3, lwd = 1.5)  # 下边界
  
  legend("topleft", legend = c("CIF1 (Theoretical)", "CIF2 (Theoretical)", 
                               "CIF1 (State Learner)", "CIF2 (State Learner)",
                               "StL CIF1 95% interval", "StL CIF2 95% interval"),
         col = c("blue", "red", "green", "purple",
                 rgb(0, 1, 0, 0.2), rgb(1, 0, 1, 0.2)
         ), 
         lwd = c(2, 2, 2, 2, 2, 2, NA, NA, NA, NA),
         lty = c(1, 1, 2, 2, NA, NA),
         fill = c(NA, NA, NA, NA, 
                  rgb(0, 1, 0, 0.2), rgb(1, 0, 1, 0.2) 
         )) 
}


drawRSF <- function(func, result, mintime) {
  
  cif_fit <- cif_fit_n(func)
  common_time <- seq(0, mintime, by = 0.01)
  n_time <- length(common_time)
  B <- length(result)
  n_event <- dim(result[[1]]$cif)[2]
  step_interp <- function(orig_time, orig_val, new_time) {
    approx(orig_time, orig_val, xout = new_time,
           method = "constant", f = 1, rule = 2)$y
  }
  
  all_CIF      <- array(NA, dim = c(B, n_time, n_event))
  
  for (i in 1:B) {
    for (k in 1:n_event) {
      all_CIF[i, , k]       <- step_interp(result[[i]]$time, result[[i]]$cif[ , k],       common_time)
      
    }
  }
  
  time_seq <- common_time

  avg_CIF <- apply(all_CIF, c(2, 3), mean)
  ci_CIF  <- apply(all_CIF, c(2, 3), function(x) quantile(x, probs = c(0.025, 0.975)))
  
  avg_CIF_1_SL <- avg_CIF[, 1]
  avg_CIF_2_SL <- avg_CIF[, 2]
  
  ci_CIF_1_SL <- ci_CIF[, , 1]  # (2 x time_points)
  ci_CIF_2_SL <- ci_CIF[, , 2]  # (2 x time_points)
  

  plot(cif_fit, xlab = "t", ylab = "CIF",
       main = "Comparison", 
       col = c("darkblue", "red"),
       lwd = c(2, 2),
       lty = c(1, 1)
  )
  
  polygon(c(time_seq, rev(time_seq)), 
          c(ci_CIF_1_SL[1, ], rev(ci_CIF_1_SL[2, ])),
          col = rgb(0, 1, 0, 0.2), border = NA)
  
  polygon(c(time_seq, rev(time_seq)), 
          c(ci_CIF_2_SL[1, ], rev(ci_CIF_2_SL[2, ])),
          col = rgb(1, 0, 1, 0.2), border = NA)
  

  lines(time_seq, avg_CIF_1_SL, col = "green", lwd = 2, lty = 2)
  lines(time_seq, avg_CIF_2_SL, col = "purple", lwd = 2, lty = 2)

  lines(time_seq, ci_CIF_1_SL[1, ], col = "darkgreen", lty = 3, lwd = 1.5)
  lines(time_seq, ci_CIF_1_SL[2, ], col = "darkgreen", lty = 3, lwd = 1.5)
  
  lines(time_seq, ci_CIF_2_SL[1, ], col = "darkblue", lty = 3, lwd = 1.5)
  lines(time_seq, ci_CIF_2_SL[2, ], col = "darkblue", lty = 3, lwd = 1.5)
  
  legend("topleft", legend = c("CIF1 (Theoretical)", "CIF2 (Theoretical)", 
                               "CIF1 (RSF)", "CIF2 (RSF)",
                               "RSF CIF1 95% interval", "RSF CIF2 95% interval"),
         col = c("darkblue", "red", "green", "purple",
                 rgb(0, 1, 0, 0.2), rgb(1, 0, 1, 0.2)),
         lwd = c(2, 2, 2, 2, NA, NA),
         lty = c(1, 1, 2, 2, NA, NA),
         fill = c(NA, NA, NA, NA, 
                  rgb(0, 1, 0, 0.2), rgb(1, 0, 1, 0.2)))
  
}


#Theoretical RD and RR
rrrd <- function(func) {
  dat <- func(300000)
  
  cif_fit <- cuminc(ftime = dat$Time, fstatus = dat$Event,
                    group = dat$A, cencode = 0)
  
  extract_cif <- function(obj) data.frame(time = obj$time, est = obj$est)
  
  remove_duplicate_time <- function(df) {
    df[!duplicated(df$time, fromLast = TRUE), ]
  }
  
  cif_0_1 <- remove_duplicate_time(extract_cif(cif_fit[["0 1"]]))
  cif_1_1 <- remove_duplicate_time(extract_cif(cif_fit[["1 1"]]))
  cif_0_2 <- remove_duplicate_time(extract_cif(cif_fit[["0 2"]]))
  cif_1_2 <- remove_duplicate_time(extract_cif(cif_fit[["1 2"]]))
  
  t_max <- min(
    max(cif_0_1$time),
    max(cif_1_1$time),
    max(cif_0_2$time),
    max(cif_1_2$time)
  )
  t_all <- seq(0, t_max, length.out = 101)
  

  step_interp <- function(df) {
    approx(df$time, df$est, xout = t_all,
           method = "constant", f = 1, rule = 2)$y
  }
  
  cif0_e1 <- step_interp(cif_0_1)
  cif1_e1 <- step_interp(cif_1_1)
  cif0_e2 <- step_interp(cif_0_2)
  cif1_e2 <- step_interp(cif_1_2)
  
  rd_e1 <- cif1_e1 - cif0_e1
  rr_e1 <- cif1_e1 / cif0_e1
  rr_e1[!is.finite(rr_e1)] <- NA
  
  rd_e2 <- cif1_e2 - cif0_e2
  rr_e2 <- cif1_e2 / cif0_e2
  rr_e2[!is.finite(rr_e2)] <- NA
  
  return(list(
    time = t_all,
    rd_e1 = rd_e1,
    rr_e1 = rr_e1,
    rd_e2 = rd_e2,
    rr_e2 = rr_e2
  ))
}


drawrrrd <- function(func,result1){
  
  theory_rr_rd  <- rrrd(func)
  max_times <- sapply(result1, function(x) max(x$time))
  common_max_time <- min(max_times)
  common_time <- seq(0, common_max_time, length.out = 101)
  
  step_interp <- function(time_vec, cif_mat, new_time) {
    apply(cif_mat, 2, function(est) {
      approx(time_vec, est, xout = new_time,
             method = "constant", f = 1, rule = 2)$y
    })
  }  
  n_sim <- length(result1)
  n_time <- length(common_time)
  n_event <- ncol(result1[[1]]$cif1)
  
  cif1 <- array(NA, dim = c(n_sim, n_time, n_event))
  cif0 <- array(NA, dim = c(n_sim, n_time, n_event))
  
  for (i in 1:n_sim) {
    cif1[i,,] <- step_interp(result1[[i]]$time, result1[[i]]$cif1, common_time)
    cif0[i,,] <- step_interp(result1[[i]]$time, result1[[i]]$cif0, common_time)
  }
  
  rd_mean <- matrix(NA, nrow = n_time, ncol = n_event)
  rr_mean <- matrix(NA, nrow = n_time, ncol = n_event)
  rd_lower <- matrix(NA, nrow = n_time, ncol = n_event)
  rd_upper <- matrix(NA, nrow = n_time, ncol = n_event)
  rr_lower <- matrix(NA, nrow = n_time, ncol = n_event)
  rr_upper <- matrix(NA, nrow = n_time, ncol = n_event)
  
  for (k in 1:n_event) {
    rd_mat <- cif1[,,k,drop=FALSE] - cif0[,,k,drop=FALSE]
    denom <- cif0[,,k,drop=FALSE]
    
    rr_mat <- ifelse(denom == 0, NA, cif1[,,k,drop=FALSE] / denom)
    
    rd_mat2 <- array(rd_mat, dim = c(dim(rd_mat)[1], dim(rd_mat)[2]))
    rr_mat2 <- array(rr_mat, dim = c(dim(rr_mat)[1], dim(rr_mat)[2]))
    
    rd_mean[,k] <- colMeans(rd_mat, na.rm = TRUE)
    rr_mean[,k] <- colMeans(rr_mat, na.rm = TRUE)
    
    rd_lower[,k] <- apply(rd_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
    rd_upper[,k] <- apply(rd_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
    rr_lower[,k] <- apply(rr_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
    rr_upper[,k] <- apply(rr_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
  }
  
  sim_res<- list(
    rd = rd_mean,
    rr = rr_mean,
    rd_lower = rd_lower,
    rd_upper = rd_upper,
    rr_lower = rr_lower,
    rr_upper = rr_upper
  )

  time_sim <- common_time

  time_theory <- theory_rr_rd$time
  
  rd_sim_e1 <- sim_res$rd[,1]; rd_sim_e1_lower <- sim_res$rd_lower[,1]; rd_sim_e1_upper <- sim_res$rd_upper[,1]
  rd_sim_e2 <- sim_res$rd[,2]; rd_sim_e2_lower <- sim_res$rd_lower[,2]; rd_sim_e2_upper <- sim_res$rd_upper[,2]

  rd_theory_e1 <- theory_rr_rd$rd_e1
  rd_theory_e2 <- theory_rr_rd$rd_e2

  ymin <- min(c(rd_sim_e1_lower, rd_sim_e2_lower, rd_theory_e1, rd_theory_e2), na.rm = TRUE)
  ymax <- max(c(rd_sim_e1_upper, rd_sim_e2_upper, rd_theory_e1, rd_theory_e2), na.rm = TRUE)
  
  plot(time_sim, rd_sim_e1, type = "n",
       ylim = c(ymin, ymax),
       xlab = "Time", ylab = "Risk Difference",
       main = "Risk Difference (RD) over Time")
  
  polygon(c(time_sim, rev(time_sim)), c(rd_sim_e1_lower, rev(rd_sim_e1_upper)),
          col = "lightblue", border = "darkgreen", lwd = 1)
  polygon(c(time_sim, rev(time_sim)), c(rd_sim_e2_lower, rev(rd_sim_e2_upper)),
          col = "yellow", border = "orange", lwd = 1)
  
  lines(time_sim, rd_sim_e1, col = "darkblue", lwd = 1, lty = 2)      # lty = 2: dashed
  lines(time_sim, rd_sim_e2, col = "brown", lwd = 1, lty = 2)
  
  lines(time_theory, rd_theory_e1, col = "blue", lwd = 2, lty = 1)    # lty = 1: solid
  lines(time_theory, rd_theory_e2, col = "red", lwd = 2, lty = 1)
  
  legend("topleft",
         legend = c(
           "Event1 Theory", 
           "Event2 Theory",
           "Event1 SL Simulation", 
           "Event2 SL Simulation", 
           "95% Empirical Interval",
           "95% Empirical Interval"
         ),
         col = c("blue", "red", "darkblue", "brown", "darkgreen", "orange"),
         lty = c(1, 1, 2, 2, NA, NA),
         lwd = c(2, 2, 1, 1, NA, NA),
         fill = c(NA, NA, NA, NA, "lightblue", "yellow"),
         border = c(NA, NA, NA, NA, "darkgreen", "orange"),
         bg = "white",
         cex = 0.8)

  rr_sim_e1 <- sim_res$rr[,1]; rr_sim_e1_lower <- sim_res$rr_lower[,1]; rr_sim_e1_upper <- sim_res$rr_upper[,1]
  rr_sim_e2 <- sim_res$rr[,2]; rr_sim_e2_lower <- sim_res$rr_lower[,2]; rr_sim_e2_upper <- sim_res$rr_upper[,2]
  
  rr_theory_e1 <- theory_rr_rd$rr_e1
  rr_theory_e2 <- theory_rr_rd$rr_e2
  
  ymin <- min(c(rr_sim_e1_lower, rr_sim_e2_lower, rr_theory_e1, rr_theory_e2), na.rm = TRUE)
  ymax <- max(c(rr_sim_e1_upper, rr_sim_e2_upper, rr_theory_e1, rr_theory_e2), na.rm = TRUE)
  
  plot(time_sim, rr_sim_e1, type = "n",
       xlim =c(0.4,0.85),
       ylim = c(ymin, 8),
       xlab = "Time", ylab = "Risk Ratio",
       main = "Risk Ratio (RR) over Time")

  polygon(c(time_sim, rev(time_sim)), c(rr_sim_e1_lower, rev(rr_sim_e1_upper)),
          col = "lightblue", border = "darkgreen", lwd = 1)
  polygon(c(time_sim, rev(time_sim)), c(rr_sim_e2_lower, rev(rr_sim_e2_upper)),
          col = "yellow", border = "orange", lwd = 1)
  
  lines(time_sim, rr_sim_e1, col = "darkblue", lwd = 1, lty = 2)
  lines(time_sim, rr_sim_e2, col = "brown", lwd = 1, lty = 2)
  
  lines(time_theory, rr_theory_e1, col = "blue", lwd = 2, lty = 1)
  lines(time_theory, rr_theory_e2, col = "red", lwd = 2, lty = 1)
  
  legend("topright",
         legend = c(
           "Event1 Theory", 
           "Event2 Theory",
           "Event1 SL Simulation", 
           "Event2 SL Simulation", 
           "95% Empirical Interval",
           "95% Empirical Interval"
         ),
         col = c("blue", "red", "darkblue", "brown", "darkgreen", "orange"),
         lty = c(1, 1, 2, 2, NA, NA),
         lwd = c(2, 2, 1, 1, NA, NA),
         fill = c(NA, NA, NA, NA, "lightblue", "yellow"),
         border = c(NA, NA, NA, NA, "darkgreen", "orange"),
         bg = "white",
         cex = 0.8)  
  
}


W2dataU_fixed <- function(n, Z_value = 0.5, A_value = 1) {
  X <- rep(Z_value, n)
  Treat <- rep(A_value, n)
  e <- rnorm(n, mean = 0, sd = 0.1)
  Y <-  X^3 + X^2 + e
  
  shape <- 5
  scale <- 1
  beta1 <- 3; beta2 <- 2
  beta_t1 <- 1.5; beta_t2 <- -1.0
  betay1 <- -1; betay2 <- 2
  
  scale1 <- scale * exp(-(beta1 * X + beta_t1 * Treat + betay1 * Y) / shape)
  scale2 <- scale * exp(-(beta2 * X + beta_t2 * Treat + betay2 * Y) / shape)
  
  time_state1 <- rweibull(n, shape = shape, scale = scale1)
  time_state2 <- rweibull(n, shape = shape, scale = scale2)
  
  censoring_time <- rweibull(n, shape = 5, scale = 0.75)
  time <- pmin(time_state1, time_state2)
  status <- ifelse(time_state1 < time_state2, 1, 2)
  observed_time <- pmin(time, censoring_time)
  status <- ifelse(time > censoring_time, 0, status)
  
  data.frame(Time = observed_time, Event = status, Z = X, A = Treat, U=Y)
}

W2data2_fixed <- function(n,Z=0.5,A=1) {
 
  X <- rep(Z, n)
  Treat <- rep(A, n)
  

  shape <- 5
  scale <- 1
  beta1 <- 3
  beta2 <- 2
  beta_t1 <- 1.5
  beta_t2 <- -1.0
  
 
  scale1 <- scale * exp(-(beta1 * X + beta_t1 * Treat) / shape)
  scale2 <- scale * exp(-(beta2 * X + beta_t2 * Treat) / shape)
  

  time_state1 <- rweibull(n, shape = shape, scale = scale1)
  time_state2 <- rweibull(n, shape = shape, scale = scale2)
  

  censoring_time <- rweibull(n, shape = 5, scale = 0.75)
  

  time <- pmin(time_state1, time_state2)
  status <- ifelse(time_state1 < time_state2, 1, 2)
  observed_time <- pmin(time, censoring_time)
  status <- ifelse(time > censoring_time, 0, status)
  

  data <- data.frame(Time = observed_time, Event = status, Z = X, A = Treat)
  return(data)
}




