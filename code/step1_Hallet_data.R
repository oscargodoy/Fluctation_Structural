rm(list=ls())

library(dplyr)
library(anisoFUN) #this package includes most of the functions for structural stability.

source('code/toolbox_coexistence.R')

#read data 
#Hallet et al. 2018 ELE (https://onlinelibrary.wiley.com/doi/10.1111/ele.13341)
d1 <- read.table(file = "data/avena_erodium_data.csv", header=T, sep=",")

# calculate the feasiblity domain for the average environment 
# 1st step calculate average of each parameter across environments. 
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

#Compute the feasibility domain and whether species can coexist for this two species case. 

niche <- 10^Omega(A) # this is niche differences
fitness <- theta(A, r) # this is fitness differences
feasibility <- test_feasibility(A,r) #does coexistence occur?
feasibility

#plot the outcome
NEED TO DO THIS STILL

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


#plot the outcome
NEED TO DO THIS STILL

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

#Compute the feasibility domain and whether species can coexist for this two species case. 

niche <- 10^Omega(A) # this is niche differences
fitness <- theta(A, r) # this is fitness differences
feasibility <- test_feasibility(A,r) #does coexistence occur?
feasibility




