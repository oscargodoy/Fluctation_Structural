# plot distance to the edge against low density growth rate 
rm(list=ls())

library(dplyr)
library(jpeg)

source("code/R_toolbox/toolbox_coexistence.R")
source("code/R_toolbox/lemkelcp.R")
source("code/R_toolbox/ISbuild.R")
source("code/R_toolbox/ISgraph.R")
source("code/R_toolbox/figs/circleCones.R")
source("code/R_toolbox/figs/sphereCones.R")
source("code/R_toolbox/figs/sphereSubCones.R")
source("code/R_toolbox/figs/InformationField.R")

source("code/R_toolbox/LDGR.R")
source("code/R_toolbox/calculate_distance_to_border_2sp.R")
source("code/R_toolbox/relative_distances_to_centroid_2sp.R")

# load libraries
#devtools::install_github("MITEcology/feasibility_analysis")
library(feasibilityR)
library(anisoFun)

# load data
A <-matrix(data =NA, nrow = 2, ncol = 2)
A[1,1] <- -0.75
A[2,2] <- -0.75
A[1,2] <- -0.25
A[2,1] <- -0.25

alpha11 <- A[1,1]
alpha22 <- A[2,2]
alpha12 <- A[1,2]
alpha21 <- A[2,1]



# Additional function to study impact of varying vectors of intrinsic growth rate on LDGR and distance to the edge
study_parameter_space <- function(ra_range, rb_range, alpha11, alpha22, alpha12, alpha21) {
  n <- length(ra_range)
  m <- length(rb_range)
  
  results <- matrix(NA, nrow = n, ncol = m)
  lambda1_matrix <- matrix(NA, nrow = n, ncol = m)
  lambda2_matrix <- matrix(NA, nrow = n, ncol = m)
  
  # Define outcome codes
  outcomes <- c("Priority effects" = 0, 
                "Species 1 wins" = 1, 
                "Species 2 wins" = 2, 
                "Coexistence" = 3)
  
  
  # Calculate outcomes for all parameter combinations
  for (i in 1:n) {
    for (j in 1:m) {
      ra <- ra_range[i]
      rb <- rb_range[j]
      
      res <- calculate_growth_rate_when_rare(ra, rb, alpha11, alpha22, alpha12, alpha21)
      lambda1 <- res$lambda1
      lambda2 <- res$lambda2
      
      lambda1_matrix[i,j] <- lambda1
      lambda2_matrix[i,j] <- lambda2
      
      # Determine outcome based on invasion criteria
      if (lambda1 > 1 && lambda2 > 1) {
        results[i,j] <- outcomes["Coexistence"]
      } else if (lambda1 > 1 && lambda2 < 1) {
        results[i,j] <- outcomes["Species 1 wins"]
      } else if (lambda1 < 1 && lambda2 > 1) {
        results[i,j] <- outcomes["Species 2 wins"]
      } else {
        results[i,j] <- outcomes["Priority effects"]
      }
    }
  }
  
  return(list(
    results = results,
    lambda1_matrix = lambda1_matrix,
    lambda2_matrix = lambda2_matrix,
    ra_range = ra_range,
    rb_range = rb_range,
    outcomes = outcomes
  ))
}

# Example parameter space exploration
ra_range <- seq(3.8, 1.0, length.out = 50)
rb_range <- seq(1.0,1.0, length.out = 50)

#calculating the low density growth rates
param_space <- study_parameter_space(ra_range, rb_range, alpha11, alpha22, alpha12, alpha21)

# Create a matrix of outcomes
outcome_matrix <- matrix(NA, nrow = length(param_space$ra_range), ncol = length(param_space$rb_range))

# Create a matrix of outcomes
outcome_matrix_rel <- matrix(NA, nrow = length(param_space$ra_range), ncol = length(param_space$rb_range))

#computing distances to the edge for a given matrix but variability in the intrinsic growth rates
for (i in 1:length(param_space$ra_range)) {
  for (j in 1:length(param_space$rb_range)) {
    X <- calculate_distance_to_border_2sp(A, c(param_space$ra_range[i], param_space$rb_range[j]))
    Y <- relative_distances_to_centroid_2sp(A, c(param_space$ra_range[i], param_space$rb_range[j]), norm = "yes")
    outcome_matrix[i,j] <- X[[1]]
    outcome_matrix_rel[i,j] <- Y[[1]]
  }
}

#Plot the outcome of LDGR versus distance to the edge 
# Set theme for consistent, clean appearance
library(ggplot2)

# Create a data frame for plotting
df <- data.frame(
  distance = outcome_matrix[,1],
  lambda1 = log(param_space$lambda1_matrix[,1]),
  lambda2 = log(param_space$lambda2_matrix[,1])
)

# Convert to long format for easier plotting
df_long <- reshape2::melt(df, id.vars = "distance", 
                          variable.name = "parameter", 
                          value.name = "log_value")

# Create enhanced plot
ggplot(df_long, aes(x = distance, y = log_value, color = parameter, shape = parameter)) +
  # Add points with better size and transparency
  geom_point(size = 3, alpha = 0.7) +
  # Add reference lines with better styling
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray", size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray", size = 0.7) +
  # Add annotations for regions
  annotate("text", x = 0.025, y = 0.55, label = "Mutual invasibility", 
           fontface = "bold",angle = 90, size = 4.5) +
  annotate("text", x = -0.025, y = 0.55, label = "Exclusion", 
           fontface = "bold", angle = 90, size = 4.5) +
  # Add regression lines to help visualize trends
  geom_smooth(method = "loess", se = FALSE, linetype = "solid", alpha = 0.7, size = 1) +
  # Set colors with a better palette
  scale_color_manual(values = c("lambda1" = "#0072B2", "lambda2" = "#D55E00"),
                     labels = c("Superior competitor", "Inferior competitor")) +
  scale_shape_manual(values = c(16, 16),
                     labels = c("Superior competitor", "Inferior competitor")) +
  # Improve labels and title
  labs(
    x = "Distance to Edge of Feasibility Domain",
    y = "Low Density Growth Rate (log scale)",
    title = "",
    color = " ",
    shape = " "
  ) +
  # Set plot theme for better appearance
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "darkgray"),
    legend.position = c(0.9,0.75),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA),
    axis.text = element_text(size = 11)
  ) +
  # Set consistent plot dimensions
  coord_cartesian(xlim = c(-0.05, 0.6), ylim = c(-0.3, 1.5))


# Doing the same for relative distances to the centroid.----

# Set theme for consistent, clean appearance
library(ggplot2)

# Create a data frame for plotting
df <- data.frame(
  distance = outcome_matrix_rel[,1],
  lambda1 = log(param_space$lambda1_matrix[,1]),
  lambda2 = log(param_space$lambda2_matrix[,1])
)

# Convert to long format for easier plotting
df_long <- reshape2::melt(df, id.vars = "distance", 
                          variable.name = "parameter", 
                          value.name = "log_value")

# Create enhanced plot
ggplot(df_long, aes(x = distance, y = log_value, color = parameter, shape = parameter)) +
  # Add points with better size and transparency
  geom_point(size = 3, alpha = 0.7) +
  # Add reference lines with better styling
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkgray", size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray", size = 0.7) +
  # Add annotations for regions
  annotate("text", x = 0.7, y = 0.5, label = "Mutual invasibility", 
           fontface = "bold", size = 4.5) +
  annotate("text", x = 1.05, y = 0.7, label = "Exclusion", 
           fontface = "bold", angle = 90, size = 4.5) +
  # Add regression lines to help visualize trends
  geom_smooth(method = "loess", se = FALSE, linetype = "solid", alpha = 0.7, size = 1) +
  # Set colors with a better palette
  scale_color_manual(values = c("lambda1" = "#0072B2", "lambda2" = "#D55E00"),
                     labels = c("Superior competitor", "Inferior competitor")) +
  scale_shape_manual(values = c(16, 16),
                     labels = c("Superior competitor", "Inferior competitor")) +
  # Improve labels and title
  labs(
    x = "Relative distance to the centroid",
    y = "Low Density Growth Rate (log scale)",
    title = "Structural Stability and Low Density Growth Rate",
    color = "Parameter",
    shape = "Parameter"
  ) +
  # Set plot theme for better appearance
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "darkgray"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 11)
  ) +
  # Set consistent plot dimensions
  coord_cartesian(xlim = c(0.1, 2.5), ylim = c(-0.3, 1.5))






