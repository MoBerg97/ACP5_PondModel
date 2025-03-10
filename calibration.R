## Calibriation

library(ggplot2)
library(patchwork)
library(randomForest)

params_Sim = read.csv("C:/Users/euble/OneDrive/Master/ACP5/ACP5_PondModel_Main/ACP5_PondModel_Main/cali/params_Simstrat_LHC_n500.csv")
em_Sim = read.csv("C:/Users/euble/OneDrive/Master/ACP5/ACP5_PondModel_Main/ACP5_PondModel_Main/cali/Simstrat_LHC_n500.csv")

params_GOTM = read.csv("C:/Users/euble/OneDrive/Master/ACP5/ACP5_PondModel_Main/ACP5_PondModel_Main/cali/params_GOTM_LHC_n500.csv")
em_GOTM = read.csv("C:/Users/euble/OneDrive/Master/ACP5/ACP5_PondModel_Main/ACP5_PondModel_Main/cali/GOTM_LHC_n500.csv")

params_FLake = read.csv("C:/Users/euble/OneDrive/Master/ACP5/ACP5_PondModel_Main/ACP5_PondModel_Main/cali/params_FLake_LHC_n100.csv")
em_FLake = read.csv("C:/Users/euble/OneDrive/Master/ACP5/ACP5_PondModel_Main/ACP5_PondModel_Main/cali/FLake_LHC_n100.csv")

join_Sim = full_join(params_Sim, em_Sim, by = "par_id")
join_GOTM = full_join(params_GOTM, em_GOTM, by = "par_id")
join_FLake = full_join(params_FLake, em_FLake, by = "par_id")



## Calculate composite error metric score

  ## Weights (scaled with standard deviation to approximate equal influence of individual metrics)
  w1 = 1/sd(join_Sim$rmse) # RMSE
  w2 = 1/sd((1-join_Sim$nse)) # NSE
  w3 = 1/sd(abs(join_Sim$bias)) # Bias
  w4 = 1/sd(1- join_Sim$r) # R^2
  
  ## Compisite error metric
  join_Sim$score = w1*join_Sim$rmse + 10*w2*(1-join_Sim$nse) + w3*abs(join_Sim$bias) + w4*(1- join_Sim$r)



#### Visual assessment ####

  ## Wind speed

  # Extract the computed density data from ggplot
  dens_data <- ggplot_build(
    ggplot(join_sim, aes(x = wind_speed)) +
      geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4)
  )$data[[1]]
  
  # Find wind speed where density is maximized
  optimal_wind_speed <- dens_data$x[which.max(dens_data$y)]
  
  
  
  
  a_ws = ggplot(join_sim, aes(x = wind_speed, y = score)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", color = "blue") +
    geom_vline(xintercept = optimal_wind_speed, color = "red", linetype = "dashed") +
    labs(title = "Composite Score vs Wind Speed",
         x = "Wind Speed",
         y = "Composite Score")
  
  b_ws = ggplot(join_sim, aes(x = wind_speed)) +
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
    ggplot(join_sim, aes(x = swr)) +
      geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4)
  )$data[[1]]
  
  # Find wind speed where density is maximized
  optimal_swr <- dens_data$x[which.max(dens_data$y)]
  
  
  a_swr = ggplot(join_sim, aes(x = swr, y = score)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", color = "blue") +
    geom_vline(xintercept = optimal_swr, color = "red", linetype = "dashed") +
    labs(title = "Composite Score vs swr",
         x = "swr",
         y = "Composite Score")
  
  b_swr = ggplot(join_sim, aes(x = swr)) +
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
    ggplot(join_sim, aes(x = Kw)) +
      geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4)
  )$data[[1]]
  
  # Find wind speed where density is maximized
  optimal_Kw <- dens_data$x[which.max(dens_data$y)]
  
  
  a_Kw = ggplot(join_sim, aes(x = Kw, y = score)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", color = "blue") +
    geom_vline(xintercept = optimal_Kw, color = "red", linetype = "dashed") +
    labs(title = "Composite Score vs Kw",
         x = "Kw",
         y = "Composite Score")
  
  b_Kw = ggplot(join_sim, aes(x = Kw)) +
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
    ggplot(join_sim, aes(x = a_seiche)) +
      geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4)
  )$data[[1]]
  
  # Find wind speed where density is maximized
  optimal_a_seiche <- dens_data$x[which.max(dens_data$y)]
  
  
  a_a_seiche = ggplot(join_sim, aes(x = a_seiche, y = score)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", color = "blue") +
    geom_vline(xintercept = optimal_a_seiche, color = "red", linetype = "dashed") +
    labs(title = "Composite Score vs a_seiche",
         x = "a_seiche",
         y = "Composite Score")
  
  b_a_seiche = ggplot(join_sim, aes(x = a_seiche)) +
    geom_density(aes(weight = 1 / score), fill = "blue", alpha = 0.4) +
    geom_vline(xintercept = optimal_a_seiche, color = "red", linetype = "dashed") +
    labs(title = "Density of a_seiche for Low Composite Score",
         x = "a_seiche",
         y = "Weighted Density")
  
  a_a_seiche+b_a_seiche
  
  print(optimal_a_seiche)
  
  
  #### Random Forest predictions ####
  ## Use a random forest model to predict the best set of parameters (using brute force grid search) 
  rf_model <- randomForest(score ~ wind_speed + swr + Kw + a_seiche, 
                           data = join_sim, 
                           ntree = 500, importance = TRUE)
  
  importance(rf_model)  
  varImpPlot(rf_model)  # Plot variable importance
  

  # Define ranges of parameters
  wind_range <- seq(min(join_sim$wind_speed), max(join_sim$wind_speed), length.out = 50)
  swr_range <- seq(min(join_sim$swr), max(join_sim$swr), length.out = 50)
  Kw_range <- seq(min(join_sim$Kw), max(join_sim$Kw), length.out = 50)
  a_seiche_range <- seq(min(join_sim$a_seiche), max(join_sim$a_seiche), length.out = 50)
  
  # Create all combinations
  param_grid <- expand.grid(wind_speed = wind_range, swr = swr_range, Kw = Kw_range, a_seiche = a_seiche_range)
  
  # Predict scores
  param_grid$predicted_score <- predict(rf_model, newdata = param_grid)
  
  # Find parameter set with the lowest predicted score
  best_params_rf <- param_grid[which.min(param_grid$predicted_score), ]
  print(best_params_rf)
  

  #### Rank ####
  ## Ranking parameters sets on composite error score and describing parameter ranges of best sets
  
  ## Finding optimal parameters by ranking parameter sets by composite error score
  join_sim$rank <- rank(join_sim$score)
  
  # Select the best 20% of models
  best_params_rank <- join_sim[join_sim$rank <= quantile(join_sim$rank, 0.1), ]
  
  # Summarize the best parameter ranges
  #summary(best_params_rank[, c("wind_speed", "swr", "Kw", "a_seiche")])
  sapply((best_params_rank[, c("wind_speed", "swr", "Kw", "a_seiche")]), mean)
  
  
  #### Comparing results
  print(c(optimal_wind_speed, optimal_swr, optimal_Kw, optimal_a_seiche))
  print(best_params_rf)
  sapply((best_params_rank[, c("wind_speed", "swr", "Kw", "a_seiche")]), mean)
  