#%%
rm(list = ls())
library(tidyr)
library(purrr)
library(dplyr)
library(ggplot2)
source("code/R_toolbox/multispecies distance to the edge.r")
generate_inte_rand <- function(S, sigma, conne = 1, dist = "norm", mu = 0) {
  if (!is.numeric(S)) {stop("S must be numerical")}
  if (!is.vector(S)) {stop("S must be a single number")}
  if (length(S) > 1) {stop("S must be a single number")}
  if (!is.numeric(sigma)) {stop("sigma must be numerical")}
  if (!is.vector(sigma)) {stop("sigma must be a single number")}
  if (length(sigma) > 1) {stop("sigma must be a single number")}
  if (!is.numeric(conne)) {stop("conne must be numerical")}
  if (!is.vector(conne)) {stop("conne must be a single number")}
  if (length(conne) > 1) {stop("conne must be a single number")}
  if (!dist %in% c("norm", "lnorm", "half-norm")) {stop("dist must be one either norm, lnorm or half-norm")}
  if (!is.numeric(mu)) {stop("mu must be numerical")}
  if (!is.vector(mu)) {stop("mu must be a single number")}
  if (length(mu) > 1) {stop("mu must be a single number")}
  if (dist == "norm") {
    matA <- rnorm(S * S, mean = mu, sd = sigma)
  } else if (dist == "lnorm") {
    matA <- -rlnorm(S * S, mean = mu, sd = sigma)
  } else if (dist == "half-norm") {
    matA <- abs(rnorm(S * S, mean = mu, sd = sigma))
  }
  zeroes <- sample(
    c(rep.int(1, floor(S * S * conne)), rep.int(0, (S * S - floor(S * S * conne))))
  )
  matA[which(zeroes == 0)] <- 0
  matA <- matrix(matA, ncol = S, nrow = S)
  diag(matA) <- -1
  return(matA)
}
#%% load data
A_int <- matrix(nrow = 5,  data = c(
  0.119368637, -0.05263433, -0.1749874, -0.09584908, -0.28214625,
  0.502825014, -0.20448711,  0.3722590, -0.06135142,  0.12308237,
 -0.003493152, -0.10534205, -0.3028440, -0.06804226,  0.60075539,
 -0.112259616, -0.06606926, -0.0526510,  0.07152292,  0.08208773,
  0.002410886, -0.06914852,  0.5126345, -0.11674101, -0.41231919))

r0 <- c( 0.89213202, -0.63922380,  1.05595753, -0.02177712,  0.79042702)

df <- Measure_Distances(A_int, r0, norm = "l2", nsample = 1000, just_min = FALSE)
print(df)
df <- Measure_Distances(A_int, r0, norm = "l2", nsample = 1000, just_min = TRUE)
print(df) #NB the minimum distance does not necessarily correspond to the minimum of the first output, since the distances are recalculated
#%%
#%% test 0 inside the first hyperquadrant
n_sp <- 10
n_samples <- 70
measures_matrix <- matrix(nrow = n_samples, ncol = n_sp)
A_int <- -diag(n_sp) #the first hyperquadrant
#seed:
set.seed(1234)
for (i in 1:n_samples) {
  r0 <- c(runif(n_sp, 0, 1))
  r0 <- - A_int %*% c(runif(n_sp, 0, 1)) ## Generate random location inside the feasibility domain
   
  df <- Measure_Distances(A_int, r0, norm = "l2", nsample = 1000, just_min = FALSE)
  measures_matrix[i,] <-df[,1]

  print(i)
}

hist( matrix(measures_matrix,ncol = 1 ), breaks = 50, main = "Histogram of distances to the border inside the first hyperquadrant")

#%% test 1 random A but in the first hyperquadrant
n_sp <- 10
n_samples <- 70
measures_matrix <- matrix(nrow = n_samples, ncol = n_sp)
#seed:
set.seed(1234)
for (i in 1:n_samples) {
  A_int <- - abs(generate_inte_rand(n_sp, 1, 1, "norm")) #the first hyperquadrant
  r0 <- c(runif(n_sp, 0, 1))
  r0 <- - A_int %*% c(runif(n_sp, 0, 1)) ## Generate random location inside the feasibility domain
   
  df <- Measure_Distances(A_int, r0, norm = "l2", nsample = 1000, just_min = FALSE)
  measures_matrix[i,] <-df[,1]

  print(i)
}

hist( matrix(measures_matrix,ncol = 1 ), breaks = 50, main = "Histogram of distances to the border inside the FD, which is in first hyperquadrant")

#%% test 2 random A and r random, anywhere
n_sp <- 10
n_samples <- 70
measures_matrix <- matrix(nrow = n_samples, ncol = n_sp)
#seed:
set.seed(1234)
for (i in 1:n_samples) {
  A_int <- generate_inte_rand(n_sp, 1, 1, "norm")
  r0 <- c(runif(n_sp, 0, 1))
  random_sign <- sample(c(-1, 1), n_sp, replace = TRUE)
  r0 <- r0 * random_sign
   
  df <- Measure_Distances(A_int, r0, norm = "l2", nsample = 1000, just_min = FALSE)
  measures_matrix[i,] <-df[,1]

  print(i)
}

hist( matrix(measures_matrix,ncol = 1 ), breaks = 50, main = "Histogram of distances to the border inside, random r and A")

#%% test 3 random A and r random, anywhere, by species number
n_sp_vector <- c(3, 4, 5, 10, 15)
n_samples <- 70
measures <- list()
extincts <- list()
min_distances <- list()
set.seed(1234)
for(s in 1:length(n_sp_vector)) {
n_sp <- n_sp_vector[s]
measures_matrix <- matrix(nrow = n_samples, ncol = n_sp)
extinct_matrix <- matrix(nrow = n_samples, ncol = 1)
min_distances_matrix <- matrix(nrow = n_samples, ncol = 1)

    for (i in 1:n_samples) {
    A_int <- generate_inte_rand(n_sp, 1, 1, "norm")
    r0 <- c(runif(n_sp, 0, 1))
    random_sign <- sample(c(-1, 1), n_sp, replace = TRUE)
    r0 <- r0 * random_sign
    
    df <- Measure_Distances(A_int, r0, norm = "l2", nsample = 1000, just_min = FALSE)
    measures_matrix[i,] <-df[,1]
    extinct_matrix[i,] <- sum(df[,2])/n_sp
    min_distances_matrix[i,] <- min(df[,1])

    print(i)
    }
    measures[[s]] <- measures_matrix
    extincts[[s]] <- extinct_matrix
    min_distances[[s]] <- min_distances_matrix
}
#%% A Plot histograms for each species together in the same plot with ggplot2
measures_df <- do.call(rbind, lapply(1:length(measures), function(i) {
  data.frame(distance = as.vector(measures[[i]]), species = rep(n_sp_vector[i], n_samples))
}))
ggplot(measures_df, aes(x = distance, color = factor(species))) +
  geom_freqpoly(bins = 50, linewidth = 1) +
  labs(title = "Histogram of distances to the border",
       x = "Distance to the border",
       y = "Frequency",
       color = "Species") +
  theme_minimal()

#%%  B
extinct_df <- do.call(rbind, lapply(1:length(extincts), function(i) {
  data.frame(extinct = as.vector(extincts[[i]]), species = rep(n_sp_vector[i], n_samples))
}))
ggplot(extinct_df, aes(x = extinct, color = factor(species))) +
facet_wrap(~ species, scales = "free") +
  geom_freqpoly(bins = 10, linewidth = 1) +
  labs(title = "Histogram of extinction rates (1 = all species extinct)",
       x = "Extinction rate",
       y = "Frequency",
       color = "Species") +
  theme_minimal()

#%% C Plot minimum distances for each r vector
min_distances_df <- do.call(rbind, lapply(1:length(min_distances), function(i) {
  data.frame(min_distance = as.vector(min_distances[[i]]), species = rep(n_sp_vector[i], n_samples))
}))
ggplot(min_distances_df, aes(x = min_distance, color = factor(species))) +
  geom_freqpoly(bins = 15, linewidth = 1) +
  labs(title = "Histogram of minimum distances to the border",
       x = "Minimum distance to the border",
       y = "Frequency",
       color = "Species") +
  theme_minimal()
  