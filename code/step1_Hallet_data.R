rm(list=ls())

library(dplyr)
require(devtools)
install_github("RadicalCommEcol/anisoFun")
library(anisoFUN) #this package includes most of the functions for structural stability.

source("code/R_toolbox/toolbox_coexistence.R")
source("code/R_toolbox/lemkelcp.R")
source("code/R_toolbox/ISbuild.R")
source("code/R_toolbox/ISgraph.R")
source("code/R_toolbox/figs/circleCones.R")
source("code/R_toolbox/figs/sphereCones.R")
source("code/R_toolbox/figs/sphereSubCones.R")
source("code/R_toolbox/figs/InformationField.R")

#read data 
#Hallet et al. 2018 ELE (https://onlinelibrary.wiley.com/doi/10.1111/ele.13341)
d1 <- read.table(file = "data/avena_erodium_data.csv", header=T, sep=",")

# calculate the feasiblity domain for the average environment 
# 1st step calculate average of each parameter across environments----
lambda_avena <- d1 %>% filter(species == 'Avena') %>% summarise(avg = mean(lambda_mean))
lambda_erodium <- d1 %>% filter(species == 'Erodium') %>% summarise(avg = mean(lambda_mean))
alphaii_avena <- d1 %>% filter(species == 'Avena') %>% summarise(avg = mean(alpha_ii_mean))
alphaii_erodium <- d1 %>% filter(species == 'Erodium') %>% summarise(avg = mean(alpha_ii_mean))
alphaij_avena <- d1 %>% filter(species == 'Avena') %>% summarise(avg = mean(alpha_ij_mean))
alphaij_erodium <- d1 %>% filter(species == 'Erodium') %>% summarise(avg = mean(alpha_ij_mean))

r <-as.vector(unlist(rbind(lambda_avena, lambda_erodium))) # vector of intrinsic growth rates
A <-matrix(data =NA, nrow = 2, ncol = 2)
A[1,1] <-unlist(alphaii_avena)
A[2,2] <-unlist(alphaii_erodium)
A[1,2] <-unlist(alphaij_avena)
A[2,1] <-unlist(alphaij_erodium)

#we multiply by -1 as competition should be negative 
A <- A*-1

#Compute the feasibility domain and whether species can coexist for this two species case. 

niche <- 10^Omega(A) # this is niche differences
fitness <- theta(A, r) # this is fitness differences
feasibility <- test_feasibility(A,r) #does coexistence occur?
feasibility


#normalize the vector r
point <-r/sqrt(sum(r^2))
points(point[1], point[2], pch=19, col="red", cex=0.5)

#Same as before but for the upper bound

lambda_avena_upper <- d1 %>% filter(species == 'Avena') %>% summarise(avg = mean(lambda_upper))
lambda_erodium_upper <- d1 %>% filter(species == 'Erodium') %>% summarise(avg = mean(lambda_upper))
alphaii_avena_upper <- d1 %>% filter(species == 'Avena') %>% summarise(avg = mean(alpha_ii_upper))
alphaii_erodium_upper <- d1 %>% filter(species == 'Erodium') %>% summarise(avg = mean(alpha_ii_upper))
alphaij_avena_upper <- d1 %>% filter(species == 'Avena') %>% summarise(avg = mean(alpha_ij_upper))
alphaij_erodium_upper <- d1 %>% filter(species == 'Erodium') %>% summarise(avg = mean(alpha_ij_upper))

r <-as.vector(unlist(rbind(lambda_avena_upper, lambda_erodium_upper))) # vector of intrinsic growth rates
A <-matrix(data =NA, nrow = 2, ncol = 2)
A[1,1] <-unlist(alphaii_avena_upper)
A[2,2] <-unlist(alphaii_erodium_upper)
A[1,2] <-unlist(alphaij_avena_upper)
A[2,1] <-unlist(alphaij_erodium_upper)

#Compute the feasibility domain and whether species can coexist for this two species case. 

niche <- 10^Omega(A) # this is niche differences
fitness <- theta(A, r) # this is fitness differences
feasibility <- test_feasibility(A,r) #does coexistence occur?
feasibility

#we multiply by -1 as competition should be negative 
A <- A*-1

#plot the outcome
drawCircleCones(A)
# 1. Draw all cones without labels:  
drawCircleCones(A, allCones = TRUE, drawLabels = FALSE)
# 2. Remove cones out of 11:  
drawCircleCones(A, allCones = FALSE, drawLabels = FALSE)
# 3. Add labels just for the 11 cone: 
drawCircleCones(A, allCones = FALSE, drawLabels = TRUE)

#normalize the vector r
point <-r/sqrt(sum(r^2))
points(point[1], point[2], pch=19, col="red", cex=1)

#Same as before but for the lower bound

lambda_avena_lower <- d1 %>% filter(species == 'Avena') %>% summarise(avg = mean(lambda_lower))
lambda_erodium_lower <- d1 %>% filter(species == 'Erodium') %>% summarise(avg = mean(lambda_lower))
alphaii_avena_lower <- d1 %>% filter(species == 'Avena') %>% summarise(avg = mean(alpha_ii_lower))
alphaii_erodium_lower <- d1 %>% filter(species == 'Erodium') %>% summarise(avg = mean(alpha_ii_lower))
alphaij_avena_lower <- d1 %>% filter(species == 'Avena') %>% summarise(avg = mean(alpha_ij_lower))
alphaij_erodium_lower <- d1 %>% filter(species == 'Erodium') %>% summarise(avg = mean(alpha_ij_lower))

r <-as.vector(unlist(rbind(lambda_avena_lower, lambda_erodium_lower))) # vector of intrinsic growth rates
A <-matrix(data =NA, nrow = 2, ncol = 2)
A[1,1] <-unlist(alphaii_avena_lower)
A[2,2] <-unlist(alphaii_erodium_lower)
A[1,2] <-unlist(alphaij_avena_lower)
A[2,1] <-unlist(alphaij_erodium_lower)

#we multiply by -1 as competition should be negative 
A <- A*-1

#Compute the feasibility domain and whether species can coexist for this two species case. 

niche <- 10^Omega(A) # this is niche differences
fitness <- theta(A, r) # this is fitness differences
feasibility <- test_feasibility(A,r) #does coexistence occur?
feasibility

#plot the outcome
drawCircleCones(A)
# 1. Draw all cones without labels:  
drawCircleCones(A, allCones = TRUE, drawLabels = FALSE)
# 2. Remove cones out of 11:  
drawCircleCones(A, allCones = FALSE, drawLabels = FALSE)
# 3. Add labels just for the 11 cone: 
drawCircleCones(A, allCones = FALSE, drawLabels = TRUE)

#normalize the vector r
point <-r/sqrt(sum(r^2))
points(point[1], point[2], pch=19, col="red", cex=1)

# 2st step calculate the cone for each environment----

d1 <- read.table(file = "data/avena_erodium_data.csv", header=T, sep=",")

c_dry <-subset(d1, treatment=="Consistent dry")
esf_dry <-subset(d1, treatment=="Early season fall dry")
lsp_dry <-subset(d1, treatment=="Late season spring dry")
c_wet <-subset(d1, treatment=="Consistent wet")


#consistent dry treatment
r_c_dry <- as.numeric(c_dry$lambda_mean)
A_c_dry <- as.matrix(c_dry[,6:7])*-1

A_c_dry[1,2] <- 0 #for visualization

niche <- 10^Omega(A_c_dry) # this is niche differences
fitness <- theta(A_c_dry, r_c_dry) # this is fitness differences
feasibility <- test_feasibility(A_c_dry,r_c_dry) #does coexistence occur?
feasibility

plot_c_dry <- drawCircleCones(A_c_dry, allCones = FALSE, drawLabels = TRUE)

#normalize the vector r
point <-r_c_dry/sqrt(sum(r_c_dry^2))
points(point[1], point[2], pch=19, col="red", cex=0.5)

#Early season fall dry

r_esf_dry <- as.numeric(esf_dry$lambda_mean)
A_esf_dry <- as.matrix(esf_dry[,6:7])*-1

niche <- 10^Omega(A_esf_dry) # this is niche differences
fitness <- theta(A_esf_dry, r_esf_dry) # this is fitness differences
feasibility <- test_feasibility(A_esf_dry,r_esf_dry) #does coexistence occur?
feasibility

plot_esf_dry <- drawCircleCones(A_esf_dry, allCones = FALSE, drawLabels = TRUE)

#normalize the vector r
point <-r_esf_dry/sqrt(sum(r_esf_dry^2))
points(point[1], point[2], pch=19, col="red", cex=0.5)

# Late season spring dry treatment

r_lsp_dry <- as.numeric(lsp_dry$lambda_mean)
A_lsp_dry <- as.matrix(lsp_dry[,6:7])*-1

niche <- 10^Omega(A_lsp_dry) # this is niche differences
fitness <- theta(A_lsp_dry, r_lsp_dry) # this is fitness differences
feasibility <- test_feasibility(A_lsp_dry,r_lsp_dry) #does coexistence occur?
feasibility

plot_lsp_dry <- drawCircleCones(A_lsp_dry, allCones = FALSE, drawLabels = TRUE)

#normalize the vector r
point <-r_lsp_dry/sqrt(sum(r_lsp_dry^2))
points(point[1], point[2], pch=19, col="red", cex=0.5)

# Consistent wet treatment

r_c_wet <- as.numeric(c_wet$lambda_mean)
A_c_wet <- as.matrix(c_wet[,6:7])*-1

niche <- 10^Omega(A_c_wet) # this is niche differences
fitness <- theta(A_c_wet, r_c_wet) # this is fitness differences
feasibility <- test_feasibility(A_c_wet,r_c_wet) #does coexistence occur?
feasibility

plot_c_wet <- drawCircleCones(A_c_wet, allCones = FALSE, drawLabels = TRUE)

#normalize the vector r
point <-r_c_wet/sqrt(sum(r_c_wet^2))
points(point[1], point[2], pch=19, col="red", cex=0.5)

#plot the four treatments
par(mfrow = c(2, 2))
plot_c_dry
plot_esf_dry
plot_lsp_dry
plot_c_wet




