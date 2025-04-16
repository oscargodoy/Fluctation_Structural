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

# 2. Plot the cone without labels 
drawCircleCones(A, allCones = TRUE, drawLabels = FALSE)
r1 <- c(2.5,1)
ra <- r1[1]
rb <- r1[2]

#normalize the r vectors and add the arrow
point1 <-r1/sqrt(sum(r1^2))
points(point1[1], point1[2], pch=19, col="red", cex=0.75)
arrows(0, 0, x1 = point1[1], y1 = point1[2], length = 0.05, angle = 30,
       code = 2, col = par("fg"), lty = 2,
       lwd = par("lwd"))
text(x=point1[1]+0.05, y=point1[2]+0.04, "A", cex=0.7)

# Additional function to study impact of varying vectors of instrinsic growth rate on LDGR and distance to the edge
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

#computing distances to the edge for a given matrix but variability in the intrinsic growth rates
for (i in 1:length(param_space$ra_range)) {
  for (j in 1:length(param_space$rb_range)) {
    X <- calculate_distance_to_border_2sp(A, c(param_space$ra_range[i], param_space$rb_range[j]))
    outcome_matrix[i,j] <- X[[1]]
  }
}

plot(outcome_matrix[,1],
     log(param_space$lambda1_matrix[,1]), 
     xlab = "Distance to the edge of the feasbility domain", 
     ylab = "Low Density Growth Rate (log scale)", 
     xlim = c(-0.05, 0.5),
     ylim = c(-0.3, 1.5),
     main = " Relating structural stability to low density growth rate",
     pch = 19, col = "blue")

#add a line to the plot with other lambda2 values
points(outcome_matrix[,1], log(param_space$lambda2_matrix[,1]), pch = 19, col = "red")

# add a vertical line at x = 0
abline(v = 0, col = "black", lty = 2)
# add a horizontal line at y = 0
abline(h = 0, col = "black", lty = 2)

text(0.15,0.5,"Mutual invasibility")
text(-0.05, 0.5, "Exclusion", v)





