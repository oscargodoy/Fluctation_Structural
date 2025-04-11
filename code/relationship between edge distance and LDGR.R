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

# load libraries
#devtools::install_github("MITEcology/feasibility_analysis")
library(feasibilityR)

#change the code to allow both positive and negative distances to the edge.
distances <- c()
fuera <- c()

# load data
A <-matrix(data =NA, nrow = 2, ncol = 2)
A[1,1] <- -0.5
A[2,2] <- -0.5
A[1,2] <- -0.25
A[2,1] <- -0.25

r0 <- c(1,1.2)

# 2. Plot the cone without labels 
drawCircleCones(A, allCones = TRUE, drawLabels = TRUE)



aux <- solve(-A_int,r0)

if(all(aux>0)){     
  distances <- feasibility_resistance_full(A_int, r0, nsample = 500)
  fuera <- 0
}else {        
  distances <- feasibility_resistance_full(A_int, r0, nsample = 500)
  fuera <- 1
  distances <- -distances
}
