rm(list=ls())

library(dplyr)

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
source("code/LDGR_including facilitation.R")

# Define parameters for these four different scenarios (A,B,C,D)

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

# Run full simulation with equal starting populations
sim_results <- simulate_competition(
  ra, rb, alpha11, alpha22, alpha12, alpha21,
  N1_init = 1.0, 
  N2_init = 1.0,
  timesteps = 200
)

cat("\nFinal populations after simulation:\n")
cat(paste("Species 1:", round(sim_results$final_N1, 3), "\n"))
cat(paste("Species 2:", round(sim_results$final_N2, 3), "\n"))

# Determine the outcome
if(sim_results$final_N1 > 0.01 && sim_results$final_N2 < 0.01) {
  outcome <- "Species 1 wins (Species 2 goes extinct)"
} else if(sim_results$final_N1 < 0.01 && sim_results$final_N2 > 0.01) {
  outcome <- "Species 2 wins (Species 1 goes extinct)"
} else if(sim_results$final_N1 > 0.01 && sim_results$final_N2 > 0.01) {
  outcome <- "Coexistence"
} else {
  outcome <- "Both species go extinct"
}

cat(paste("Simulation outcome:", outcome, "\n"))

# Plot the dynamics
par(mfrow=c(1,1))
visualize_dynamics(sim_results, "Population Dynamics with Equal Starting Populations")

#scenario B

ra <-r2[1]   # Intrinsic growth rate of species 1
rb <- r2[2]   # Intrinsic growth rate of species 2
alpha11 <- A[1,1]*-1 # Effect of species 1 on itself
alpha22 <- A[2,2]*-1  # Effect of species 2 on itself
alpha12 <- A[1,2]*-1 # Effect of species 2 on species 1
alpha21 <- A[2,1]*-1  # Effect of species 1 on species 2

# Run full simulation with equal starting populations
sim_results <- simulate_competition(
  ra, rb, alpha11, alpha22, alpha12, alpha21,
  N1_init = 1.0, 
  N2_init = 1.0,
  timesteps = 200
)

cat("\nFinal populations after simulation:\n")
cat(paste("Species 1:", round(sim_results$final_N1, 3), "\n"))
cat(paste("Species 2:", round(sim_results$final_N2, 3), "\n"))

# Determine the outcome
if(sim_results$final_N1 > 0.01 && sim_results$final_N2 < 0.01) {
  outcome <- "Species 1 wins (Species 2 goes extinct)"
} else if(sim_results$final_N1 < 0.01 && sim_results$final_N2 > 0.01) {
  outcome <- "Species 2 wins (Species 1 goes extinct)"
} else if(sim_results$final_N1 > 0.01 && sim_results$final_N2 > 0.01) {
  outcome <- "Coexistence"
} else {
  outcome <- "Both species go extinct"
}

cat(paste("Simulation outcome:", outcome, "\n"))

# Plot the dynamics
par(mfrow=c(1,1))
visualize_dynamics(sim_results, "Population Dynamics with Equal Starting Populations")

#scenario C

ra <-r3[1]   # Intrinsic growth rate of species 1
rb <- r3[2]   # Intrinsic growth rate of species 2
alpha11 <- A[1,1]*-1 # Effect of species 1 on itself
alpha22 <- A[2,2]*-1  # Effect of species 2 on itself
alpha12 <- A[1,2]*-1 # Effect of species 2 on species 1
alpha21 <- A[2,1]*-1  # Effect of species 1 on species 2

# Run full simulation with equal starting populations
sim_results <- simulate_competition(
  ra, rb, alpha11, alpha22, alpha12, alpha21,
  N1_init = 1.0, 
  N2_init = 1.0,
  timesteps = 200
)

cat("\nFinal populations after simulation:\n")
cat(paste("Species 1:", round(sim_results$final_N1, 3), "\n"))
cat(paste("Species 2:", round(sim_results$final_N2, 3), "\n"))

# Determine the outcome
if(sim_results$final_N1 > 0.01 && sim_results$final_N2 < 0.01) {
  outcome <- "Species 1 wins (Species 2 goes extinct)"
} else if(sim_results$final_N1 < 0.01 && sim_results$final_N2 > 0.01) {
  outcome <- "Species 2 wins (Species 1 goes extinct)"
} else if(sim_results$final_N1 > 0.01 && sim_results$final_N2 > 0.01) {
  outcome <- "Coexistence"
} else {
  outcome <- "Both species go extinct"
}

cat(paste("Simulation outcome:", outcome, "\n"))

# Plot the dynamics
par(mfrow=c(1,1))
visualize_dynamics(sim_results, "Population Dynamics with Equal Starting Populations")

#scenario D

ra <-r4[1]   # Intrinsic growth rate of species 1
rb <- r4[2] # Intrinsic growth rate of species 2
alpha11 <- A[1,1]*-1 # Effect of species 1 on itself
alpha22 <- A[2,2]*-1  # Effect of species 2 on itself
alpha12 <- A[1,2]*-1 # Effect of species 2 on species 1
alpha21 <- A[2,1]*-1  # Effect of species 1 on species 2

cat("\nFinal populations after simulation:\n")
cat(paste("Species 1:", round(sim_results$final_N1, 3), "\n"))
cat(paste("Species 2:", round(sim_results$final_N2, 3), "\n"))

# Determine the outcome
if(sim_results$final_N1 > 0.01 && sim_results$final_N2 < 0.01) {
  outcome <- "Species 1 wins (Species 2 goes extinct)"
} else if(sim_results$final_N1 < 0.01 && sim_results$final_N2 > 0.01) {
  outcome <- "Species 2 wins (Species 1 goes extinct)"
} else if(sim_results$final_N1 > 0.01 && sim_results$final_N2 > 0.01) {
  outcome <- "Coexistence"
} else {
  outcome <- "Both species go extinct"
}

cat(paste("Simulation outcome:", outcome, "\n"))

# Plot the dynamics
par(mfrow=c(1,1))
visualize_dynamics(sim_results, "Population Dynamics with Equal Starting Populations")

  
