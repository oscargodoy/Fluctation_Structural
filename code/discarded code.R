





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

