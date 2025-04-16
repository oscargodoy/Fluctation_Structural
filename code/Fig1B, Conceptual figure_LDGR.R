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

#1st step define some points and draw them in the unit circle----
dev.off()
#Define a matrix with some degree of facilitation (positive values)
A <-matrix(data =NA, nrow = 2, ncol = 2)
A[1,1] <- -0.5
A[2,2] <- -0.5
A[1,2] <- 0.2 #facilitation
A[2,1] <- -0.25

#put several vectors on top of the figure. We define four scenarios
r1 <- c(-0.2,1)
r2 <- c(0.7,1)
r3 <- c(2,1)
r4 <- c(2.2,0.3)

# 2. Plot the cone without labels 
drawCircleCones(A, allCones = FALSE, drawLabels = FALSE)

#normalize the r vectors and add the arrow
point1 <-r1/sqrt(sum(r1^2))
points(point1[1], point1[2], pch=19, col="red", cex=0.75)
arrows(0, 0, x1 = point1[1], y1 = point1[2], length = 0.05, angle = 30,
       code = 2, col = par("fg"), lty = 2,
       lwd = par("lwd"))
text(x=point1[1]+0.05, y=point1[2]+0.04, "A", cex=0.7)


point2 <-r2/sqrt(sum(r2^2))
points(point2[1], point2[2], pch=19, col="red", cex=0.75)
arrows(0, 0, x1 = point2[1], y1 = point2[2], length = 0.05, angle = 30,
       code = 2, col = par("fg"), lty = 2,
       lwd = par("lwd"))
text(x=point2[1]+0.05, y=point2[2]+0.04, "B", cex=0.7)

point3 <-r3/sqrt(sum(r3^2))
points(point3[1], point3[2], pch=19, col="red", cex=0.75)
arrows(0, 0, x1 = point3[1], y1 = point3[2], length = 0.05, angle = 30,
       code = 2, col = par("fg"), lty = 2,
       lwd = par("lwd"))
text(x=point3[1]+0.05, y=point3[2]+0.04, "C", cex=0.7)

point4 <-r4/sqrt(sum(r4^2))
points(point4[1], point4[2], pch=19, col="red", cex=0.75)
arrows(0, 0, x1 = point4[1], y1 = point4[2], length = 0.05, angle = 30,
       code = 2, col = par("fg"), lty = 2,
       lwd = par("lwd"))
text(x=point4[1]+0.05, y=point4[2]+0.04, "D", cex=0.7)

text(x=1, y=-0.05, "r1", cex=0.7)
text(x=0.05, y=1.03, "r2", cex=0.7)

dev.off()

#2nd step calculate LDGR for each scenario----
source("code/R_toolbox/LDGR.R")

# Define parameters for these four different scenarios (A,B,C,D)
# to store the results of log(LDGR)
ldgr<-matrix(NA, nrow = 4, ncol = 2)

#put several vectors on top of the figure. We define four scenarios
r1 <- c(-0.2,1)
r2 <- c(0.7,1)
r3 <- c(2,1)
r4 <- c(2.2,0.3)

#scenario A
ra <-r1[1]   # Intrinsic growth rate of species 1
rb <- r1[2]   # Intrinsic growth rate of species 2
alpha11 <- A[1,1]*-1 # Effect of species 1 on itself
alpha22 <- A[2,2]*-1  # Effect of species 2 on itself
alpha12 <- A[1,2]*-1 # Effect of species 2 on species 1
alpha21 <- A[2,1]*-1  # Effect of species 1 on species 2

# Calculate analytical growth rates when rare
analytical_results <- calculate_growth_rate_when_rare(ra, rb, alpha11, alpha22, alpha12, alpha21)
print("Analytical results:")

print(paste("Growth rate of species 1 when rare (lambda1):", analytical_results$lambda1))
print(paste("Growth rate of species 2 when rare (lambda2):", analytical_results$lambda2))

# Verify with simulations
simulation_results <- simulate_invasion(ra, rb, alpha11, alpha22, alpha12, alpha21)
print("Simulation results:")
print(paste("Empirical growth rate of species 1 when rare:", simulation_results$emp_lambda1))
print(paste("Empirical growth rate of species 2 when rare:", simulation_results$emp_lambda2))

ldgr[1,2] <- log(analytical_results$lambda2)

# Visualize the dynamics
visualize_dynamics(simulation_results)

# Analyze coexistence for our example
coexistence_outcome <- analyze_coexistence(analytical_results$lambda1, analytical_results$lambda2)
print(paste("Coexistence analysis:", coexistence_outcome))

dev.off()

#scenario B

ra <-r2[1]   # Intrinsic growth rate of species 1
rb <- r2[2]   # Intrinsic growth rate of species 2
alpha11 <- A[1,1]*-1 # Effect of species 1 on itself
alpha22 <- A[2,2]*-1  # Effect of species 2 on itself
alpha12 <- A[1,2]*-1 # Effect of species 2 on species 1
alpha21 <- A[2,1]*-1  # Effect of species 1 on species 2

# Calculate analytical growth rates when rare
analytical_results <- calculate_growth_rate_when_rare(ra, rb, alpha11, alpha22, alpha12, alpha21)
print("Analytical results:")

print(paste("Growth rate of species 1 when rare (lambda1):", analytical_results$lambda1))
print(paste("Growth rate of species 2 when rare (lambda2):", analytical_results$lambda2))

# Verify with simulations
simulation_results <- simulate_invasion(ra, rb, alpha11, alpha22, alpha12, alpha21)
print("Simulation results:")
print(paste("Empirical growth rate of species 1 when rare:", simulation_results$emp_lambda1))
print(paste("Empirical growth rate of species 2 when rare:", simulation_results$emp_lambda2))

ldgr[2,1] <- log(analytical_results$lambda1)
ldgr[2,2] <- log(analytical_results$lambda2)

# Visualize the dynamics
visualize_dynamics(simulation_results)

# Analyze coexistence for our example
coexistence_outcome <- analyze_coexistence(analytical_results$lambda1, analytical_results$lambda2)
print(paste("Coexistence analysis:", coexistence_outcome))

dev.off()

#scenario C

ra <-r3[1]   # Intrinsic growth rate of species 1
rb <- r3[2]   # Intrinsic growth rate of species 2
alpha11 <- A[1,1]*-1 # Effect of species 1 on itself
alpha22 <- A[2,2]*-1  # Effect of species 2 on itself
alpha12 <- A[1,2]*-1 # Effect of species 2 on species 1
alpha21 <- A[2,1]*-1  # Effect of species 1 on species 2

# Calculate analytical growth rates when rare
analytical_results <- calculate_growth_rate_when_rare(ra, rb, alpha11, alpha22, alpha12, alpha21)
print("Analytical results:")

print(paste("Growth rate of species 1 when rare (lambda1):", analytical_results$lambda1))
print(paste("Growth rate of species 2 when rare (lambda2):", analytical_results$lambda2))

# Verify with simulations
simulation_results <- simulate_invasion(ra, rb, alpha11, alpha22, alpha12, alpha21)
print("Simulation results:")
print(paste("Empirical growth rate of species 1 when rare:", simulation_results$emp_lambda1))
print(paste("Empirical growth rate of species 2 when rare:", simulation_results$emp_lambda2))

ldgr[3,1] <- log(analytical_results$lambda1)
ldgr[3,2] <- log(analytical_results$lambda2)

# Visualize the dynamics
visualize_dynamics(simulation_results)

# Analyze coexistence for our example
coexistence_outcome <- analyze_coexistence(analytical_results$lambda1, analytical_results$lambda2)
print(paste("Coexistence analysis:", coexistence_outcome))

dev.off()

#scenario D

ra <-r4[1]   # Intrinsic growth rate of species 1
rb <- r4[2] # Intrinsic growth rate of species 2
alpha11 <- A[1,1]*-1 # Effect of species 1 on itself
alpha22 <- A[2,2]*-1  # Effect of species 2 on itself
alpha12 <- A[1,2]*-1 # Effect of species 2 on species 1
alpha21 <- A[2,1]*-1  # Effect of species 1 on species 2

# Calculate analytical growth rates when rare
analytical_results <- calculate_growth_rate_when_rare(ra, rb, alpha11, alpha22, alpha12, alpha21)
print("Analytical results:")

print(paste("Growth rate of species 1 when rare (lambda1):", analytical_results$lambda1))
print(paste("Growth rate of species 2 when rare (lambda2):", analytical_results$lambda2))

# Verify with simulations
simulation_results <- simulate_invasion(ra, rb, alpha11, alpha22, alpha12, alpha21)
print("Simulation results:")
print(paste("Empirical growth rate of species 1 when rare:", simulation_results$emp_lambda1))
print(paste("Empirical growth rate of species 2 when rare:", simulation_results$emp_lambda2))

ldgr[4,1] <- log(analytical_results$lambda1)
ldgr[4,2] <- log(analytical_results$lambda2)

# Visualize the dynamics
visualize_dynamics(simulation_results)

# Analyze coexistence for our example
coexistence_outcome <- analyze_coexistence(analytical_results$lambda1, analytical_results$lambda2)
print(paste("Coexistence analysis:", coexistence_outcome))

dev.off()

# resume all in a single figure. 

#Add Log(LDGR)

#Scenario A
# Define the values for the bar plot
species_values <- c(sp1 = ldgr[1,1], sp2 = ldgr[1,2])

# Create the bar plot
jpeg(file="figures/barplot_scenario_A.jpeg")
barplot(species_values,
        ylab = "LDGR",
        col = c("lightgrey", "lightgrey"),
        ylim = c(0, 1.5),
        border = "black",
        names.arg = c(" ", " "))

# Add a grid for better readability
#grid(nx = NA, ny = NULL, lty = 2, col = "gray")

# Add value labels on top of each bar
text(x = c(0.7, 1.9), 
     y = species_values + 0.10, 
     labels = round(species_values, 2),
     cex = 4)


dev.off()

#Scenario B
# Define the values for the bar plot
species_values <- c(sp1 = ldgr[2,1], sp2 = ldgr[2,2])

# Create the bar plot
jpeg(file="figures/barplot_scenario_B.jpeg")
barplot(species_values,
        ylab = "LDGR",
        col = c("lightgrey", "lightgrey"),
        ylim = c(0, 1.5),
        border = "black",
        names.arg = c(" ", " "))

# Add a grid for better readability
#grid(nx = NA, ny = NULL, lty = 2, col = "gray")

# Add value labels on top of each bar
text(x = c(0.7, 1.9), 
     y = species_values + 0.10, 
     labels = round(species_values, 2),
     cex = 4)

dev.off()

#Scenario C
# Define the values for the bar plot
species_values <- c(sp1 = ldgr[3,1], sp2 = ldgr[3,2])
# Create the bar plot
jpeg(file="figures/barplot_scenario_C.jpeg")
barplot(species_values,
        ylab = "LDGR",
        col = c("lightgrey", "lightgrey"),
        ylim = c(0, 1.5),
        border = "black",
        names.arg = c(" ", " "))
# Add a grid for better readability
#grid(nx = NA, ny = NULL, lty = 2, col = "gray")

# Add value labels on top of each bar
text(x = c(0.7, 1.9), 
     y = species_values + 0.10, 
     labels = round(species_values, 2),
     cex = 4)
dev.off()

#Scenario D
# Define the values for the bar plot
species_values <- c(sp1 = ldgr[4,1], sp2 = ldgr[4,2])
# Create the bar plot
jpeg(file="figures/barplot_scenario_D.jpeg")
barplot(species_values,
        ylab = "LDGR",
        col = c("lightgrey", "black"),
        ylim = c(0, 1.5),
        border = "black",
        names.arg = c(" ", " "))
# Add a grid for better readability
#grid(nx = NA, ny = NULL, lty = 2, col = "gray")

# Add value labels on top of each bar
text(x = c(0.7, 1.9), 
     y = species_values +0.13, 
     labels = round(species_values, 2),
     cex = 4)
dev.off()

A <-matrix(data =NA, nrow = 2, ncol = 2)
A[1,1] <- -0.5
A[2,2] <- -0.5
A[1,2] <- 0.2 #facilitation
A[2,1] <- -0.25

#put several vectors on top of the figure. We define four scenarios
r1 <- c(-0.2,1)
r2 <- c(0.7,1)
r3 <- c(2,1)
r4 <- c(2.2,0.3)

# 2. Plot the cone without labels 
drawCircleCones(A, allCones = FALSE, drawLabels = FALSE)

#normalize the r vectors and add the arrow
point1 <-r1/sqrt(sum(r1^2))
points(point1[1], point1[2], pch=19, col="red", cex=1)
arrows(0, 0, x1 = point1[1], y1 = point1[2], length = 0.05, angle = 30,
       code = 2, col = par("fg"), lty = 2,
       lwd = par("lwd"))
text(x=point1[1]+0.05, y=point1[2]+0.04, "A", cex=1)


point2 <-r2/sqrt(sum(r2^2))
points(point2[1], point2[2], pch=19, col="red", cex=1)
arrows(0, 0, x1 = point2[1], y1 = point2[2], length = 0.05, angle = 30,
       code = 2, col = par("fg"), lty = 2,
       lwd = par("lwd"))
text(x=point2[1]+0.05, y=point2[2]+0.04, "B", cex=1)

point3 <-r3/sqrt(sum(r3^2))
points(point3[1], point3[2], pch=19, col="red", cex=1)
arrows(0, 0, x1 = point3[1], y1 = point3[2], length = 0.05, angle = 30,
       code = 2, col = par("fg"), lty = 2,
       lwd = par("lwd"))
text(x=point3[1]+0.05, y=point3[2]+0.04, "C", cex=1)

point4 <-r4/sqrt(sum(r4^2))
points(point4[1], point4[2], pch=19, col="red", cex=1)
arrows(0, 0, x1 = point4[1], y1 = point4[2], length = 0.05, angle = 30,
       code = 2, col = par("fg"), lty = 2,
       lwd = par("lwd"))
text(x=point4[1]+0.05, y=point4[2]+0.04, "D", cex=1)

text(x=1, y=-0.05, "r1", cex=1)
text(x=0.05, y=1.03, "r2", cex=1)

#add the previously saved jpge files 
img1 <- readJPEG("figures/barplot_scenario_A.jpeg")
img2 <- readJPEG("figures/barplot_scenario_B.jpeg")
img3 <- readJPEG("figures/barplot_scenario_C.jpeg")
img4 <- readJPEG("figures/barplot_scenario_D.jpeg")

# Add the images to the plot
rasterImage(img1, -0.67, 0.9, -0.37, 0.6)
rasterImage(img2, 1, 1.1, 0.7, 0.8)
rasterImage(img3, 0.9, 0.53, 1.2, 0.83)
rasterImage(img4, 1, 0.2, 1.3, 0.5)

# Add a title to the plot
title(main = "Low density growth rate under different scenarios", cex.main = 1)

# Save the plot as a JPEG file
jpeg(file="figures/LDGR_cone.jpeg", width = 800, height = 800)
# save the plot as a PDF file
pdf(file="figures/LDGR_cone.pdf", width = 12, height = 12)
dev.off()

  
