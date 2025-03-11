## Calibriation

library(ggplot2)
library(patchwork)
library(randomForest)
library(dplyr)

params_Sim = read.csv("C:/Users/euble/OneDrive/Master/ACP5/ACP5_PondModel_Main/ACP5_PondModel_Main/cali/params_Simstrat_LHC_n500.csv")
em_Sim = read.csv("C:/Users/euble/OneDrive/Master/ACP5/ACP5_PondModel_Main/ACP5_PondModel_Main/cali/Simstrat_LHC_n500.csv")

params_GOTM = read.csv("C:/Users/euble/OneDrive/Master/ACP5/ACP5_PondModel_Main/ACP5_PondModel_Main/cali/params_GOTM_LHC_n500.csv")
em_GOTM = read.csv("C:/Users/euble/OneDrive/Master/ACP5/ACP5_PondModel_Main/ACP5_PondModel_Main/cali/GOTM_LHC_n500.csv")

params_FLake = read.csv("C:/Users/euble/OneDrive/Master/ACP5/ACP5_PondModel_Main/ACP5_PondModel_Main/cali/params_FLake_LHC_n100.csv")
em_FLake = read.csv("C:/Users/euble/OneDrive/Master/ACP5/ACP5_PondModel_Main/ACP5_PondModel_Main/cali/FLake_LHC_n100.csv")

join_Sim = full_join(params_Sim, em_Sim, by = "par_id")
join_GOTM = full_join(params_GOTM, em_GOTM, by = "par_id")
join_FLake = full_join(params_FLake, em_FLake, by = "par_id")


## Merge data sets for calibration of shared parameters
Sim_red = join_Sim %>% select(-a_seiche, -par_id)
GOTM_red = join_GOTM %>% select(-turb_param.k_min, -par_id)
FLake_red = join_FLake %>% select(-c_relax_C, -par_id)

shared = rbind(Sim_red, GOTM_red, FLake_red)



#### Calculate composite error metric score ####

  ## Weights (scaled with standard deviation to approximate equal influence of individual metrics)
  w1 = 1/sd(join_Sim$rmse) # RMSE
  w2 = 1/sd((1-join_Sim$nse)) # NSE
  w3 = 1/sd(abs(join_Sim$bias)) # Bias
  w4 = 1/sd(1- join_Sim$r) # R^2
  
  ## Compisite error metric
  join_Sim$score = w1*join_Sim$rmse + w2*(1-join_Sim$nse) + w3*abs(join_Sim$bias) + w4*(1- join_Sim$r)

  
  ## Weights (scaled with standard deviation to approximate equal influence of individual metrics)
  w5 = 1/sd(join_GOTM$rmse) # RMSE
  w6 = 1/sd((1-join_GOTM$nse)) # NSE
  w7 = 1/sd(abs(join_GOTM$bias)) # Bias
  w8 = 1/sd(1- join_GOTM$r) # R^2
  
  ## Compisite error metric
  join_GOTM$score = w5*join_GOTM$rmse + w6*(1-join_GOTM$nse) + w7*abs(join_GOTM$bias) + w8*(1- join_GOTM$r)
  

  ## Weights (scaled with standard deviation to approximate equal influence of individual metrics)
  w9 = 1/sd(join_FLake$rmse) # RMSE
  w10 = 1/sd((1-join_FLake$nse)) # NSE
  w11 = 1/sd(abs(join_FLake$bias)) # Bias
  w12 = 1/sd(1- join_FLake$r) # R^2
  
  ## Compisite error metric
  join_FLake$score = w9*join_FLake$rmse + w10*(1-join_FLake$nse) + w11*abs(join_FLake$bias) + w12*(1- join_FLake$r)
  
  
  ## Weights (scaled with standard deviation to approximate equal influence of individual metrics)
  w13 = 1/sd(shared$rmse) # RMSE
  w14 = 1/sd((1-shared$nse)) # NSE
  w15 = 1/sd(abs(shared$bias)) # Bias
  w16 = 1/sd(1- shared$r) # R^2
  
  ## Compisite error metric
  shared$score = w13*shared$rmse + w14*(1-shared$nse) + w15*abs(shared$bias) + w16*(1- shared$r)
  
  
  
 
#### Simstrat ####
##### Visual assessment #####

  ## Wind speed

  # Extract the computed density data from ggplot
  dens_data <- ggplot_build(
    ggplot(join_Sim, aes(x = wind_speed)) +
      geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4)
  )$data[[1]]
  
  # Find wind speed where density is maximized
  optimal_wind_speed <- dens_data$x[which.max(dens_data$y)]
  
  
  
  
  a_ws = ggplot(join_Sim, aes(x = wind_speed, y = score)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", color = "blue") +
    geom_vline(xintercept = optimal_wind_speed, color = "red", linetype = "dashed") +
    labs(title = "Composite Score vs Wind Speed",
         x = "Wind Speed",
         y = "Composite Score")
  
  b_ws = ggplot(join_Sim, aes(x = wind_speed)) +
    geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4) +
    geom_vline(xintercept = optimal_wind_speed, color = "red", linetype = "dashed") +
    labs(title = "Density of Wind Speed for Low Composite Score",
         x = "Wind Speed",
         y = "Weighted Density")
  
  a_ws+b_ws
  
  ## Optimal wind speed
  print(optimal_wind_speed)
  
  
  ## SWR
  # Extract the computed density data from ggplot
  dens_data <- ggplot_build(
    ggplot(join_Sim, aes(x = swr)) +
      geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4)
  )$data[[1]]
  
  # Find wind speed where density is maximized
  optimal_swr <- dens_data$x[which.max(dens_data$y)]
  
  
  a_swr = ggplot(join_Sim, aes(x = swr, y = score)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", color = "blue") +
    geom_vline(xintercept = optimal_swr, color = "red", linetype = "dashed") +
    labs(title = "Composite Score vs swr",
         x = "swr",
         y = "Composite Score")
  
  b_swr = ggplot(join_Sim, aes(x = swr)) +
    geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4) +
    geom_vline(xintercept = optimal_swr, color = "red", linetype = "dashed") +
    labs(title = "Density of swr for Low Composite Score",
         x = "swr",
         y = "Weighted Density")
  
  a_swr+b_swr
  
  print(optimal_swr)
  
  ## Kw
  # Extract the computed density data from ggplot
  dens_data <- ggplot_build(
    ggplot(join_Sim, aes(x = Kw)) +
      geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4)
  )$data[[1]]
  
  # Find wind speed where density is maximized
  optimal_Kw <- dens_data$x[which.max(dens_data$y)]
  
  
  a_Kw = ggplot(join_Sim, aes(x = Kw, y = score)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", color = "blue") +
    geom_vline(xintercept = optimal_Kw, color = "red", linetype = "dashed") +
    labs(title = "Composite Score vs Kw",
         x = "Kw",
         y = "Composite Score")
  
  b_Kw = ggplot(join_Sim, aes(x = Kw)) +
    geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4) +
    geom_vline(xintercept = optimal_Kw, color = "red", linetype = "dashed") +
    labs(title = "Density of Kw for Low Composite Score",
         x = "Kw",
         y = "Weighted Density")
  
  a_Kw+b_Kw
  
  print(optimal_Kw)
  
  ## a_seiche
  # Extract the computed density data from ggplot
  dens_data <- ggplot_build(
    ggplot(join_Sim, aes(x = a_seiche)) +
      geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4)
  )$data[[1]]
  
  # Find wind speed where density is maximized
  optimal_a_seiche <- dens_data$x[which.max(dens_data$y)]
  
  
  a_a_seiche = ggplot(join_Sim, aes(x = a_seiche, y = score)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", color = "blue") +
    geom_vline(xintercept = optimal_a_seiche, color = "red", linetype = "dashed") +
    labs(title = "Composite Score vs a_seiche",
         x = "a_seiche",
         y = "Composite Score")
  
  b_a_seiche = ggplot(join_Sim, aes(x = a_seiche)) +
    geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4) +
    geom_vline(xintercept = optimal_a_seiche, color = "red", linetype = "dashed") +
    labs(title = "Density of a_seiche for Low Composite Score",
         x = "a_seiche",
         y = "Weighted Density")
  
  a_a_seiche+b_a_seiche
  
  print(optimal_a_seiche)
  
  
  ##### Random Forest predictions #####
  ## Use a random forest model to predict the best set of parameters (using brute force grid search) 
  rf_model <- randomForest(score ~ wind_speed + swr + Kw + a_seiche, 
                           data = join_Sim, 
                           ntree = 500, importance = TRUE)
  
  importance(rf_model)  
  varImpPlot(rf_model)  # Plot variable importance
  

  # Define ranges of parameters
  wind_range <- seq(min(join_Sim$wind_speed), max(join_Sim$wind_speed), length.out = 50)
  swr_range <- seq(min(join_Sim$swr), max(join_Sim$swr), length.out = 50)
  Kw_range <- seq(min(join_Sim$Kw), max(join_Sim$Kw), length.out = 50)
  a_seiche_range <- seq(min(join_Sim$a_seiche), max(join_Sim$a_seiche), length.out = 50)
  
  # Create all combinations
  param_grid <- expand.grid(wind_speed = wind_range, swr = swr_range, Kw = Kw_range, a_seiche = a_seiche_range)
  
  # Predict scores
  param_grid$predicted_score <- predict(rf_model, newdata = param_grid)
  
  # Determine the 10% quantile threshold
  quantile_10 <- quantile(param_grid$predicted_score, probs = 0.0001)
  
  # Filter the dataset for the lowest 10% of predicted_score values
  filtered_data <- param_grid[param_grid$predicted_score <= quantile_10, ]
  
  # Calculate the mean for each parameter in the filtered data
  mean_values <- colMeans(filtered_data[, c("wind_speed", "swr", "Kw", "a_seiche")])
  
  # Print the results
  print(mean_values)
  
  

  

  ##### Rank #####
  ## Ranking parameters sets on composite error score and describing parameter ranges of best sets
  
  ## Finding optimal parameters by ranking parameter sets by composite error score
  join_Sim$rank <- rank(join_Sim$score)
  
  # Select the best 20% of models
  best_params_rank <- join_Sim[join_Sim$rank <= quantile(join_Sim$rank, 0.1), ]
  
  # Summarize the best parameter ranges
  #summary(best_params_rank[, c("wind_speed", "swr", "Kw", "a_seiche")])
  rank = sapply((best_params_rank[, c("wind_speed", "swr", "Kw", "a_seiche")]), mean)
  
  
  #### Comparing results
  density_Simstrat = c(optimal_wind_speed, optimal_swr, optimal_Kw, optimal_a_seiche, "Optimal density")
  rf_Simstrat = c(mean_values, "Random Forest")
  rank_Simstrat = c(rank, "Rank")  
  
  comparison_Simstrat = rbind(density_Simstrat, rf_Simstrat, rank_Simstrat)
  comparison_Simstrat = as.data.frame(comparison_Simstrat)
  names(comparison_Simstrat) = c("wind_speed", "swr", "Kw", "a_seiche", "Calibration method")
  
  
  
  
  
#### GOTM ####
  ##### Visual assessment #####
  
  ## Wind speed
  
  # Extract the computed density data from ggplot
  dens_data <- ggplot_build(
    ggplot(join_GOTM, aes(x = wind_speed)) +
      geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4)
  )$data[[1]]
  
  # Find wind speed where density is maximized
  optimal_wind_speed <- dens_data$x[which.max(dens_data$y)]
  
  
  
  
  a_ws = ggplot(join_GOTM, aes(x = wind_speed, y = score)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", color = "blue") +
    geom_vline(xintercept = optimal_wind_speed, color = "red", linetype = "dashed") +
    labs(title = "Composite Score vs Wind Speed",
         x = "Wind Speed",
         y = "Composite Score")
  
  b_ws = ggplot(join_GOTM, aes(x = wind_speed)) +
    geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4) +
    geom_vline(xintercept = optimal_wind_speed, color = "red", linetype = "dashed") +
    labs(title = "Density of Wind Speed for Low Composite Score",
         x = "Wind Speed",
         y = "Weighted Density")
  
  a_ws+b_ws
  
  ## Optimal wind speed
  print(optimal_wind_speed)
  
  
  ## SWR
  # Extract the computed density data from ggplot
  dens_data <- ggplot_build(
    ggplot(join_GOTM, aes(x = swr)) +
      geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4)
  )$data[[1]]
  
  # Find wind speed where density is maximized
  optimal_swr <- dens_data$x[which.max(dens_data$y)]
  
  
  a_swr = ggplot(join_GOTM, aes(x = swr, y = score)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", color = "blue") +
    geom_vline(xintercept = optimal_swr, color = "red", linetype = "dashed") +
    labs(title = "Composite Score vs swr",
         x = "swr",
         y = "Composite Score")
  
  b_swr = ggplot(join_GOTM, aes(x = swr)) +
    geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4) +
    geom_vline(xintercept = optimal_swr, color = "red", linetype = "dashed") +
    labs(title = "Density of swr for Low Composite Score",
         x = "swr",
         y = "Weighted Density")
  
  a_swr+b_swr
  
  print(optimal_swr)
  
  ## Kw
  # Extract the computed density data from ggplot
  dens_data <- ggplot_build(
    ggplot(join_GOTM, aes(x = Kw)) +
      geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4)
  )$data[[1]]
  
  # Find wind speed where density is maximized
  optimal_Kw <- dens_data$x[which.max(dens_data$y)]
  
  
  a_Kw = ggplot(join_GOTM, aes(x = Kw, y = score)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", color = "blue") +
    geom_vline(xintercept = optimal_Kw, color = "red", linetype = "dashed") +
    labs(title = "Composite Score vs Kw",
         x = "Kw",
         y = "Composite Score")
  
  b_Kw = ggplot(join_GOTM, aes(x = Kw)) +
    geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4) +
    geom_vline(xintercept = optimal_Kw, color = "red", linetype = "dashed") +
    labs(title = "Density of Kw for Low Composite Score",
         x = "Kw",
         y = "Weighted Density")
  
  a_Kw+b_Kw
  
  print(optimal_Kw)
  
  ## turb_param.k_min
  # Extract the computed density data from ggplot
  dens_data <- ggplot_build(
    ggplot(join_GOTM, aes(x = turb_param.k_min)) +
      geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4)
  )$data[[1]]
  
  # Find wind speed where density is maximized
  optimal_turb_param.k_min <- dens_data$x[which.max(dens_data$y)]
  
  
  a_turb_param.k_min = ggplot(join_GOTM, aes(x = turb_param.k_min, y = score)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", color = "blue") +
    geom_vline(xintercept = optimal_turb_param.k_min, color = "red", linetype = "dashed") +
    labs(title = "Composite Score vs turb_param.k_min",
         x = "turb_param.k_min",
         y = "Composite Score")
  
  b_turb_param.k_min = ggplot(join_GOTM, aes(x = turb_param.k_min)) +
    geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4) +
    geom_vline(xintercept = optimal_turb_param.k_min, color = "red", linetype = "dashed") +
    labs(title = "Density of turb_param.k_min for Low Composite Score",
         x = "turb_param.k_min",
         y = "Weighted Density")
  
  a_turb_param.k_min+b_turb_param.k_min
  
  print(optimal_turb_param.k_min)
  
  
  ##### Random Forest predictions #####
  ## Use a random forest model to predict the best set of parameters (using brute force grid search) 
  rf_model <- randomForest(score ~ wind_speed + swr + Kw + turb_param.k_min, 
                           data = join_GOTM, 
                           ntree = 500, importance = TRUE)
  
  importance(rf_model)  
  varImpPlot(rf_model)  # Plot variable importance
  
  
  # Define ranges of parameters
  wind_range <- seq(min(join_GOTM$wind_speed), max(join_GOTM$wind_speed), length.out = 50)
  swr_range <- seq(min(join_GOTM$swr), max(join_GOTM$swr), length.out = 50)
  Kw_range <- seq(min(join_GOTM$Kw), max(join_GOTM$Kw), length.out = 50)
  turb_param.k_min_range <- seq(min(join_GOTM$turb_param.k_min), max(join_GOTM$turb_param.k_min), length.out = 50)
  
  # Create all combinations
  param_grid <- expand.grid(wind_speed = wind_range, swr = swr_range, Kw = Kw_range, turb_param.k_min = turb_param.k_min_range)
  
  # Predict scores
  param_grid$predicted_score <- predict(rf_model, newdata = param_grid)
  
  # Determine the 0.01% quantile threshold
  quantile_10 <- quantile(param_grid$predicted_score, probs = 0.0001)
  
  # Filter the dataset for the lowest 10% of predicted_score values
  filtered_data <- param_grid[param_grid$predicted_score <= quantile_10, ]
  
  # Calculate the mean for each parameter in the filtered data
  mean_values <- colMeans(filtered_data[, c("wind_speed", "swr", "Kw", "turb_param.k_min")])
  
  # Print the results
  print(mean_values)
  
  
  
  
  
  ##### Rank #####
  ## Ranking parameters sets on composite error score and describing parameter ranges of best sets
  
  ## Finding optimal parameters by ranking parameter sets by composite error score
  join_GOTM$rank <- rank(join_GOTM$score)
  
  # Select the best 20% of models
  best_params_rank <- join_GOTM[join_GOTM$rank <= quantile(join_GOTM$rank, 0.1), ]
  
  # Summarize the best parameter ranges
  #summary(best_params_rank[, c("wind_speed", "swr", "Kw", "turb_param.k_min")])
  rank = sapply((best_params_rank[, c("wind_speed", "swr", "Kw", "turb_param.k_min")]), mean)
  
  
  #### Comparing results
  density_GOTM = c(optimal_wind_speed, optimal_swr, optimal_Kw, optimal_turb_param.k_min, "Optimal density")
  rf_GOTM = c(mean_values, "Random Forest")
  rank_GOTM = c(rank, "Rank")  
  
  comparison_GOTM = rbind(density_GOTM, rf_GOTM, rank_GOTM)
  comparison_GOTM = as.data.frame(comparison_GOTM)
  names(comparison_GOTM) = c("wind_speed", "swr", "Kw", "turb_param.k_min", "Calibration method")
  
  
  
  
  
  
  
  
#### FLake ####
  ##### Visual assessment #####
  
  ## Wind speed
  
  # Extract the computed density data from ggplot
  dens_data <- ggplot_build(
    ggplot(join_FLake, aes(x = wind_speed)) +
      geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4)
  )$data[[1]]
  
  # Find wind speed where density is maximized
  optimal_wind_speed <- dens_data$x[which.max(dens_data$y)]
  
  
  
  
  a_ws = ggplot(join_FLake, aes(x = wind_speed, y = score)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", color = "blue") +
    geom_vline(xintercept = optimal_wind_speed, color = "red", linetype = "dashed") +
    labs(title = "Composite Score vs Wind Speed",
         x = "Wind Speed",
         y = "Composite Score")
  
  b_ws = ggplot(join_FLake, aes(x = wind_speed)) +
    geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4) +
    geom_vline(xintercept = optimal_wind_speed, color = "red", linetype = "dashed") +
    labs(title = "Density of Wind Speed for Low Composite Score",
         x = "Wind Speed",
         y = "Weighted Density")
  
  a_ws+b_ws
  
  ## Optimal wind speed
  print(optimal_wind_speed)
  
  
  ## SWR
  # Extract the computed density data from ggplot
  dens_data <- ggplot_build(
    ggplot(join_FLake, aes(x = swr)) +
      geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4)
  )$data[[1]]
  
  # Find wind speed where density is maximized
  optimal_swr <- dens_data$x[which.max(dens_data$y)]
  
  
  a_swr = ggplot(join_FLake, aes(x = swr, y = score)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", color = "blue") +
    geom_vline(xintercept = optimal_swr, color = "red", linetype = "dashed") +
    labs(title = "Composite Score vs swr",
         x = "swr",
         y = "Composite Score")
  
  b_swr = ggplot(join_FLake, aes(x = swr)) +
    geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4) +
    geom_vline(xintercept = optimal_swr, color = "red", linetype = "dashed") +
    labs(title = "Density of swr for Low Composite Score",
         x = "swr",
         y = "Weighted Density")
  
  a_swr+b_swr
  
  print(optimal_swr)
  
  ## Kw
  # Extract the computed density data from ggplot
  dens_data <- ggplot_build(
    ggplot(join_FLake, aes(x = Kw)) +
      geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4)
  )$data[[1]]
  
  # Find wind speed where density is maximized
  optimal_Kw <- dens_data$x[which.max(dens_data$y)]
  
  
  a_Kw = ggplot(join_FLake, aes(x = Kw, y = score)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", color = "blue") +
    geom_vline(xintercept = optimal_Kw, color = "red", linetype = "dashed") +
    labs(title = "Composite Score vs Kw",
         x = "Kw",
         y = "Composite Score")
  
  b_Kw = ggplot(join_FLake, aes(x = Kw)) +
    geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4) +
    geom_vline(xintercept = optimal_Kw, color = "red", linetype = "dashed") +
    labs(title = "Density of Kw for Low Composite Score",
         x = "Kw",
         y = "Weighted Density")
  
  a_Kw+b_Kw
  
  print(optimal_Kw)
  
  ## c_relax_C
  # Extract the computed density data from ggplot
  dens_data <- ggplot_build(
    ggplot(join_FLake, aes(x = c_relax_C)) +
      geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4)
  )$data[[1]]
  
  # Find wind speed where density is maximized
  optimal_c_relax_C <- dens_data$x[which.max(dens_data$y)]
  
  
  a_c_relax_C = ggplot(join_FLake, aes(x = c_relax_C, y = score)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", color = "blue") +
    geom_vline(xintercept = optimal_c_relax_C, color = "red", linetype = "dashed") +
    labs(title = "Composite Score vs c_relax_C",
         x = "c_relax_C",
         y = "Composite Score")
  
  b_c_relax_C = ggplot(join_FLake, aes(x = c_relax_C)) +
    geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4) +
    geom_vline(xintercept = optimal_c_relax_C, color = "red", linetype = "dashed") +
    labs(title = "Density of c_relax_C for Low Composite Score",
         x = "c_relax_C",
         y = "Weighted Density")
  
  a_c_relax_C+b_c_relax_C
  
  print(optimal_c_relax_C)
  
  
  ##### Random Forest predictions #####
  ## Use a random forest model to predict the best set of parameters (using brute force grid search) 
  rf_model <- randomForest(score ~ wind_speed + swr + Kw + c_relax_C, 
                           data = join_FLake, 
                           ntree = 500, importance = TRUE)
  
  importance(rf_model)  
  varImpPlot(rf_model)  # Plot variable importance
  
  
  # Define ranges of parameters
  wind_range <- seq(min(join_FLake$wind_speed), max(join_FLake$wind_speed), length.out = 50)
  swr_range <- seq(min(join_FLake$swr), max(join_FLake$swr), length.out = 50)
  Kw_range <- seq(min(join_FLake$Kw), max(join_FLake$Kw), length.out = 50)
  c_relax_C_range <- seq(min(join_FLake$c_relax_C), max(join_FLake$c_relax_C), length.out = 50)
  
  # Create all combinations
  param_grid <- expand.grid(wind_speed = wind_range, swr = swr_range, Kw = Kw_range, c_relax_C = c_relax_C_range)
  
  # Predict scores
  param_grid$predicted_score <- predict(rf_model, newdata = param_grid)
  
  # Determine the 0.01% quantile threshold
  quantile_10 <- quantile(param_grid$predicted_score, probs = 0.0001)
  
  # Filter the dataset for the lowest 10% of predicted_score values
  filtered_data <- param_grid[param_grid$predicted_score <= quantile_10, ]
  
  # Calculate the mean for each parameter in the filtered data
  mean_values <- colMeans(filtered_data[, c("wind_speed", "swr", "Kw", "c_relax_C")])
  
  # Print the results
  print(mean_values)
  
  
  
  
  
  ##### Rank #####
  ## Ranking parameters sets on composite error score and describing parameter ranges of best sets
  
  ## Finding optimal parameters by ranking parameter sets by composite error score
  join_FLake$rank <- rank(join_FLake$score)
  
  # Select the best 20% of models
  best_params_rank <- join_FLake[join_FLake$rank <= quantile(join_FLake$rank, 0.1), ]
  
  # Summarize the best parameter ranges
  #summary(best_params_rank[, c("wind_speed", "swr", "Kw", "c_relax_C")])
  rank = sapply((best_params_rank[, c("wind_speed", "swr", "Kw", "c_relax_C")]), mean)
  
  
  #### Comparing results
  density_FLake = c(optimal_wind_speed, optimal_swr, optimal_Kw, optimal_c_relax_C, "Optimal density")
  rf_FLake = c(mean_values, "Random Forest")
  rank_FLake = c(rank, "Rank")  
  
  comparison_FLake = rbind(density_FLake, rf_FLake, rank_FLake)
  comparison_FLake = as.data.frame(comparison_FLake)
  names(comparison_FLake) = c("wind_speed", "swr", "Kw", "c_relax_C", "Calibration method")
  
  
  
#### Shared parameters ####
  ##### Visual assessment #####
  
  ## Wind speed
  
  # Extract the computed density data from ggplot
  dens_data <- ggplot_build(
    ggplot(shared, aes(x = wind_speed)) +
      geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4)
  )$data[[1]]
  
  # Find wind speed where density is maximized
  optimal_wind_speed <- dens_data$x[which.max(dens_data$y)]
  
  
  
  
  a_ws = ggplot(shared, aes(x = wind_speed, y = score)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", color = "blue") +
    geom_vline(xintercept = optimal_wind_speed, color = "red", linetype = "dashed") +
    labs(title = "Composite Score vs Wind Speed",
         x = "Wind Speed",
         y = "Composite Score")
  
  b_ws = ggplot(shared, aes(x = wind_speed)) +
    geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4) +
    geom_vline(xintercept = optimal_wind_speed, color = "red", linetype = "dashed") +
    labs(title = "Density of Wind Speed for Low Composite Score",
         x = "Wind Speed",
         y = "Weighted Density")
  
  a_ws+b_ws
  
  ## Optimal wind speed
  print(optimal_wind_speed)
  
  
  ## SWR
  # Extract the computed density data from ggplot
  dens_data <- ggplot_build(
    ggplot(shared, aes(x = swr)) +
      geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4)
  )$data[[1]]
  
  # Find wind speed where density is maximized
  optimal_swr <- dens_data$x[which.max(dens_data$y)]
  
  
  a_swr = ggplot(shared, aes(x = swr, y = score)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", color = "blue") +
    geom_vline(xintercept = optimal_swr, color = "red", linetype = "dashed") +
    labs(title = "Composite Score vs swr",
         x = "swr",
         y = "Composite Score")
  
  b_swr = ggplot(shared, aes(x = swr)) +
    geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4) +
    geom_vline(xintercept = optimal_swr, color = "red", linetype = "dashed") +
    labs(title = "Density of swr for Low Composite Score",
         x = "swr",
         y = "Weighted Density")
  
  a_swr+b_swr
  
  print(optimal_swr)
  
  ## Kw
  # Extract the computed density data from ggplot
  dens_data <- ggplot_build(
    ggplot(shared, aes(x = Kw)) +
      geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4)
  )$data[[1]]
  
  # Find wind speed where density is maximized
  optimal_Kw <- dens_data$x[which.max(dens_data$y)]
  
  
  a_Kw = ggplot(shared, aes(x = Kw, y = score)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", color = "blue") +
    geom_vline(xintercept = optimal_Kw, color = "red", linetype = "dashed") +
    labs(title = "Composite Score vs Kw",
         x = "Kw",
         y = "Composite Score")
  
  b_Kw = ggplot(shared, aes(x = Kw)) +
    geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4) +
    geom_vline(xintercept = optimal_Kw, color = "red", linetype = "dashed") +
    labs(title = "Density of Kw for Low Composite Score",
         x = "Kw",
         y = "Weighted Density")
  
  a_Kw+b_Kw
  
  print(optimal_Kw)
  
 
  
  
  ##### Random Forest predictions #####
  ## Use a random forest model to predict the best set of parameters (using brute force grid search) 
  rf_model <- randomForest(score ~ wind_speed + swr + Kw, 
                           data = shared, 
                           ntree = 500, importance = TRUE)
  
  importance(rf_model)  
  varImpPlot(rf_model)  # Plot variable importance
  
  
  # Define ranges of parameters
  wind_range <- seq(min(shared$wind_speed), max(shared$wind_speed), length.out = 50)
  swr_range <- seq(min(shared$swr), max(shared$swr), length.out = 50)
  Kw_range <- seq(min(shared$Kw), max(shared$Kw), length.out = 50)

  # Create all combinations
  param_grid <- expand.grid(wind_speed = wind_range, swr = swr_range, Kw = Kw_range)
  
  # Predict scores
  param_grid$predicted_score <- predict(rf_model, newdata = param_grid)
  
  # Determine the 0.1% quantile threshold
  quantile_10 <- quantile(param_grid$predicted_score, probs = 0.001)
  
  # Filter the dataset for the lowest 10% of predicted_score values
  filtered_data <- param_grid[param_grid$predicted_score <= quantile_10, ]
  
  # Calculate the mean for each parameter in the filtered data
  mean_values <- colMeans(filtered_data[, c("wind_speed", "swr", "Kw")])
  
  # Print the results
  print(mean_values)
  
  
  
  
  
  ##### Rank #####
  ## Ranking parameters sets on composite error score and describing parameter ranges of best sets
  
  ## Finding optimal parameters by ranking parameter sets by composite error score
  shared$rank <- rank(shared$score)
  
  # Select the best 20% of models
  best_params_rank <- shared[shared$rank <= quantile(shared$rank, 0.1), ]
  
  # Summarize the best parameter ranges
  #summary(best_params_rank[, c("wind_speed", "swr", "Kw", "c_relax_C")])
  rank = sapply((best_params_rank[, c("wind_speed", "swr", "Kw")]), mean)
  
  
  #### Comparing results
  density_shared = c(optimal_wind_speed, optimal_swr, optimal_Kw, "Optimal density")
  rf_shared = c(mean_values, "Random Forest")
  rank_shared = c(rank, "Rank")  
  
  comparison_shared = rbind(density_shared, rf_shared, rank_shared)
  comparison_shared = as.data.frame(comparison_shared)
  names(comparison_shared) = c("wind_speed", "swr", "Kw", "Calibration method")
  
  
  #### Visualisation ####
  rf_Simstrat = as.numeric(rf_Simstrat[-c(5)])
  rf_GOTM = as.numeric(rf_GOTM[-c(5)])
  rf_FLake = as.numeric(rf_Simstrat[-c(5)])
  rf_shared = as.numeric(rf_shared[-c(4)])
  
  rank_Simstrat = as.numeric(rank_Simstrat[-c(5)])
  rank_GOTM = as.numeric(rank_GOTM[-c(5)])
  rank_FLake = as.numeric(rank_Simstrat[-c(5)])
  rank_shared = as.numeric(rank_shared[-c(4)])
  
  overview <- data.frame(
    
    value =   as.numeric(c(density_Simstrat[-c(5)], density_GOTM[-c(5)], density_FLake[-c(5)], density_shared[-c(4)],
                rf_Simstrat, rf_GOTM, rf_FLake, rf_shared,
                rank_Simstrat, rank_GOTM, rank_FLake, rank_shared)),
  
    parameter = c("wind_speed", "swr", "Kw", "a_seiche",
                  "wind_speed", "swr", "Kw", "turb_param.k_min",
                  "wind_speed", "swr", "Kw", "c_relax_C",
                  "wind_speed", "swr", "Kw",
                  "wind_speed", "swr", "Kw", "a_seiche",
                  "wind_speed", "swr", "Kw", "turb_param.k_min",
                  "wind_speed", "swr", "Kw", "c_relax_C",
                  "wind_speed", "swr", "Kw",
                  "wind_speed", "swr", "Kw", "a_seiche",
                  "wind_speed", "swr", "Kw", "turb_param.k_min",
                  "wind_speed", "swr", "Kw", "c_relax_C",
                  "wind_speed", "swr", "Kw"),
  
    model = c(rep("Simstrat", 4), rep("GOTM", 4), rep("FLake", 4), rep("Shared", 3),
              rep("Simstrat", 4), rep("GOTM", 4), rep("FLake", 4), rep("Shared", 3),
              rep("Simstrat", 4), rep("GOTM", 4), rep("FLake", 4), rep("Shared", 3)),
    calibration_method = c(rep("Optimal density", 15), rep("Random Forest", 15), rep("Rank", 15)),
    stringsAsFactors = FALSE
  )
  
  
 
  ## parameter values per model and calibration method
  
    ## Shared parameters
    ggplot(overview[overview$parameter %in% c("Kw", "swr", "wind_speed"),], aes(x = model, y = value, fill = calibration_method)) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_wrap(~ parameter, scales = "free_y") +  # One plot per parameter
      labs(title = "Parameter Values by Model and Calibration Method",
           x = "Model",
           y = "Value") +
      theme_minimal()
  
    ## Unique parameters
    ggplot(overview[overview$parameter %in% c("a_seiche", "turb_param.k_min", "c_relax_C"),], aes(x = model, y = value, fill = calibration_method)) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_wrap(~ parameter, scales = "free_y") +  # One plot per parameter
      labs(title = "Parameter Values by Model and Calibration Method",
           x = "Model",
           y = "Value") +
      theme_minimal()
    
    
  ## Distribution of parameter values
    
    ## Shared
    ggplot(overview[overview$parameter %in% c("Kw", "swr", "wind_speed"),], aes(x = parameter, y = value)) +
      geom_boxplot() + geom_jitter() +
      labs(title = "Distribution of Parameter Values Across Calibration Methods",
           x = "Parameter",
           y = "Value") +
      theme_minimal()
  
    ## Unique
    ggplot(overview[overview$parameter %in% c("a_seiche", "turb_param.k_min", "c_relax_C"),], aes(x = parameter, y = value)) +
      geom_boxplot() + geom_jitter() +
    labs(title = "Distribution of Parameter Values Across Calibration Methods",
         x = "Parameter",
         y = "Value") +
      theme_minimal()
    
    