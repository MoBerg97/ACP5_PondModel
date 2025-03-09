## Results

#### Packages ####

# Load libraries for post-processing
library(ggplot2)
library(readr)
library(lubridate)
library(tidyr)
library(dplyr)
library(readxl)
# library(zoo)
# library(data.table)

library(rLakeAnalyzer)
library(GLM3r)
library(glmtools)
library(FLakeR)
library(GOTMr)
library(gotmtools)
library(SimstratR)
library(MyLakeR)
library(LakeEnsemblR)

## Set up
ncdf <- "output_FLake_GOTM_Simstrat/ensemble_output.nc"
model = c("FLake", "GOTM", "Simstrat")

#### Results ####
# Take a look at the model fits to the observed data
error_metrics = calc_fit(ncdf = ncdf,
         model = model,
         var = "temp")



##### Model fits #####

##### Heatmap #####
p1 <- plot_heatmap(ncdf)

# Change the theme and increase text size for saving
p1 <- p1 +
  theme_classic(base_size = 24) + 
  scale_colour_gradientn(limits = c(0, 21),
                         colours = rev(RColorBrewer::brewer.pal(11, "Spectral")))

p1

##### Built-in: Models vs. Oberservations, Residuals, Temperature range boxplots (per depth) #####

# Plot ensemble mean at 0m
plot_ensemble(ncdf = ncdf, model = model, var = "temp", depth = 0, residuals = TRUE, boxwhisker = TRUE)

# Plot ensemble mean at 0.1m
plot_ensemble(ncdf = ncdf, model = model, var = "temp", depth = 0.1, residuals = TRUE, boxwhisker = TRUE)

# Plot ensemble mean at 0.15m
plot_ensemble(ncdf = ncdf, model = model, var = "temp", depth = 0.15, residuals = TRUE, boxwhisker = TRUE)

# Plot ensemble mean at 0.2m
plot_ensemble(ncdf = ncdf, model = model, var = "temp", depth = 0.2, residuals = TRUE, boxwhisker = TRUE)

# Plot ensemble mean at 0.25m
plot_ensemble(ncdf = ncdf, model = model, var = "temp", depth = 0.25, residuals = TRUE, boxwhisker = TRUE)

# Plot ensemble mean at 0.3m
plot_ensemble(ncdf = ncdf, model = model, var = "temp", depth = 0.3, residuals = TRUE, boxwhisker = TRUE)



##### ggplot2: Water temperature predictions of models vs. observed temperatures (at depths in m: 0, 0.05, 0.1, 0.15, 0.2, 0.25) #####

# Load post-processed output data into your workspace
analyse_df <- analyse_ncdf(ncdf = ncdf, spin_up = NULL, drho = 0.1, model = model)


# Example plot the summer stratification period
strat_df <- analyse_df$strat
out_df <- analyse_df$out_df
obs <- analyse_df$obs_temp

obs$model = "Observed"
obs$depths = obs$depths - 0.05
out_df = rbind(out_df, obs)


# Convert date column to POSIXct for proper time handling
out_df <- out_df %>%
  mutate(date = as.POSIXct(date, format = "%Y-%m-%d %H:%M:%S"),
         depth_label = paste0(depths, " m"))  # Create depth label

# Create the plot with labeled facets
ggplot(out_df, aes(x = date, y = temp, color = model)) +
  geom_line() +
  facet_wrap(~depth_label, scales = "free_y") +  # Label depth with "m"
  labs(title = "Water Temperature Predictions by Model",
       x = "Date", y = "Temperature (°C)", color = "Model") +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 12))  # Make facet labels clearer



