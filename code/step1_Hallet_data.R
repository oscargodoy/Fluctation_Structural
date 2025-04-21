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
source("code/R_toolbox/Functions to plot Hallet data.R")

par(mar = c(3, 3, 1.5, 1.5))  # Set margins
# Hallet Average conditions----

#read data 
#Hallet et al. 2018 ELE (https://onlinelibrary.wiley.com/doi/10.1111/ele.13341)
d1 <- read.table(file = "data/avena_erodium_data.csv", header=T, sep=";")

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

#draw the cone
draw_biodiversity_cone(A)
#draw the vector
point <-r/sqrt(sum(r^2))
arrows(0, 0, x1 = point[1], y1 = point[2], length = 0.15, angle = 30,
       code = 2, col = "brown", lty =1, lwd=2)
text(0.8, 1, "A. Average conditions", cex = 1.2, pos = 1)


# Hallet relative non-linearity in the lambdas----

d1 <- read.table(file = "data/avena_erodium_data.csv", header=T, sep=",")

c_dry <-subset(d1, treatment=="Consistent dry")
esf_dry <-subset(d1, treatment=="Early season fall dry")
lsp_dry <-subset(d1, treatment=="Late season spring dry")
c_wet <-subset(d1, treatment=="Consistent wet")

#obtain r and A and their se for each environment
#Consistent dry
r_c_dry <- as.numeric(c_dry$lambda_mean)
r_c_dry_error <- as.numeric(c_dry$lambda_se)
A_c_dry <- as.matrix(c_dry[,6:7])*-1
A_c_dry_error <- as.matrix(c_dry[,8:9])*-1
#A_c_dry[1,2] <- 0 #for visualization
#Early season fall dry
r_esf_dry <- as.numeric(esf_dry$lambda_mean)
r_esf_dry_error <- as.numeric(esf_dry$lambda_se)
A_esf_dry <- as.matrix(esf_dry[,6:7])*-1
A_esf_dry_error <- as.matrix(esf_dry[,8:9])*-1
#Late season spring dry
r_lsp_dry <- as.numeric(lsp_dry$lambda_mean)
r_lsp_dry_error <- as.numeric(lsp_dry$lambda_se)
A_lsp_dry <- as.matrix(lsp_dry[,6:7])*-1
A_lsp_dry_error <- as.matrix(lsp_dry[,8:9])*-1
#Consistent wet
r_c_wet <- as.numeric(c_wet$lambda_mean)
r_c_wet_error <- as.numeric(c_wet$lambda_se)
A_c_wet <- as.matrix(c_wet[,6:7])*-1
A_c_wet_error <- as.matrix(c_wet[,8:9])*-1

#plot it 

#draw the cone
draw_biodiversity_cone(A)
text(0.8, 1, "B. Relative non-linearity in the lambdas", cex = 1.2, pos = 1)
#draw the vector Consistently dry
point <-r_c_dry/sqrt(sum(r_c_dry^2))
arrows(0, 0, x1 = point[1], y1 = point[2], length = 0.15, angle = 30,
       code = 2, col = "brown", lty =1, lwd=2)
text(point[1] + 0.05, point[2] + 0.06, "Consistent dry", cex = 1.2, pos = 1)

#draw the vector Early season fall dry
point <-r_esf_dry/sqrt(sum(r_esf_dry^2))
arrows(0, 0, x1 = point[1], y1 = point[2], length = 0.15, angle = 30,
       code = 2, col = "brown", lty =1, lwd=2)
text(point[1] + 0.08, point[2] + 0.05, "Fall dry", cex = 1.2, pos = 1)

#draw the vector Late season spring dry
point <-r_lsp_dry/sqrt(sum(r_lsp_dry^2))
arrows(0, 0, x1 = point[1], y1 = point[2], length = 0.15, angle = 30,
       code = 2, col = "brown", lty =1, lwd=2)
text(point[1] + 0.05, point[2] + 0.07, "Spring dry", cex = 1.2, pos = 1)

#draw the vector Consistently wet
point <-r_c_wet/sqrt(sum(r_c_wet^2))
arrows(0, 0, x1 = point[1], y1 = point[2], length = 0.15, angle = 30,
       code = 2, col = "brown", lty =1, lwd=2)
text(point[1] + 0.05, point[2] + 0.07, "Consistent wet", cex = 1.2, pos = 1)


# Hallet relative non-linearity in the alphas----
#draw the cone

draw_multiple_cones(
list(A_c_dry, A_esf_dry, A_lsp_dry, A_c_wet),
labels = c("Consistently dry", "Early season fall dry","Late season spring dry", "Consistently wet"))
text(0.8, 1, "C. Relative non-linearity in the alphas", cex = 1.2, pos = 1)
#average r

r <-as.vector(unlist(rbind(lambda_avena, lambda_erodium))) # vector of intrinsic growth rates
point <-r/sqrt(sum(r^2))
arrows(0, 0, x1 = point[1], y1 = point[2], length = 0.15, angle = 30,
       code = 2, col = "brown", lty =1, lwd=2)

text(point[1]-0.07, point[2] + 0.14, "Fall dry", cex = 1.2, pos = 1)
lines (c(0.125, 0.01), c(1.08, 1), col='black')
text(point[1]-0.01, point[2] + 0.1, "Consistent dry", cex = 1.2, pos = 1)
lines (c(0.125, 0.05), c(1.035, 1), col='black')
text(point[1]+0.03, point[2] -0.3, "Spring dry", cex = 1.2, pos = 1, srt=65)
text(point[1]+0.17, point[2] -0.7, "Consistent Wet", cex = 1.2, pos = 1, srt=27)


# The storage effect ----

draw_multiple_cones(
  list(A_c_dry, A_esf_dry, A_lsp_dry, A_c_wet),
  labels = c("Consistently dry", "Early season fall dry","Late season spring dry", "Consistently wet"))
text(0.8, 1, "D. The storage effect", cex = 1.2, pos = 1)
#average r

text(point[1]-0.07, point[2] + 0.14, "Fall dry", cex = 1.2, pos = 1)
lines (c(0.125, 0.01), c(1.08, 1), col='black')
text(point[1]-0.01, point[2] + 0.1, "Consistent dry", cex = 1.2, pos = 1)
lines (c(0.125, 0.05), c(1.035, 1), col='black')
text(point[1]+0.03, point[2] -0.3, "Spring dry", cex = 1.2, pos = 1, srt=65)
text(point[1]+0.17, point[2] -0.7, "Consistent Wet", cex = 1.2, pos = 1, srt=27)

point <-r_c_dry/sqrt(sum(r_c_dry^2))
arrows(0, 0, x1 = point[1], y1 = point[2], length = 0.15, angle = 30,
       code = 2, col = "brown", lty =1, lwd=2)
text(point[1] + 0.05, point[2] + 0.06, "Consistent dry", cex = 1.2, pos = 1)

#draw the vector Early season fall dry
point <-r_esf_dry/sqrt(sum(r_esf_dry^2))
arrows(0, 0, x1 = point[1], y1 = point[2], length = 0.15, angle = 30,
       code = 2, col = "brown", lty =1, lwd=2)
text(point[1] + 0.10, point[2], "Fall dry", cex = 1.2, pos = 1)

#draw the vector Late season spring dry
point <-r_lsp_dry/sqrt(sum(r_lsp_dry^2))
arrows(0, 0, x1 = point[1], y1 = point[2], length = 0.15, angle = 30,
       code = 2, col = "brown", lty =1, lwd=2)
text(point[1] + 0.05, point[2] + 0.07, "Spring dry", cex = 1.2, pos = 1)

#draw the vector Consistently wet
point <-r_c_wet/sqrt(sum(r_c_wet^2))
arrows(0, 0, x1 = point[1], y1 = point[2], length = 0.15, angle = 30,
       code = 2, col = "brown", lty =1, lwd=2)
text(point[1] + 0.05, point[2] + 0.07, "Consistent wet", cex = 1.2, pos = 1)

#compute for each environment the distance to the edge
source("code/R_toolbox/calculate_distance_to_border_2sp.R")

#consistent dry
d_c_dry<- calculate_distance_to_border_2sp(A_c_dry, r_c_dry)
#early season fall dry
d_esf_dry<- calculate_distance_to_border_2sp(A_esf_dry, r_esf_dry)
#late season spring dry
d_lsp_dry<- calculate_distance_to_border_2sp(A_lsp_dry, r_lsp_dry)
#consistent wet
d_c_wet<- calculate_distance_to_border_2sp(A_c_wet, r_c_wet)

distances<- c(d_c_dry[[1]], d_esf_dry[[1]], d_lsp_dry[[1]], d_c_wet[[1]])
distances #in all cases species can not coexist under each particular environment

####
#### now calculate distances for one hundred matrices for each environment
####

#consistent dry----

# Set a seed for reproducibility
set.seed(123)

# Define the mean matrix and standard error matrix
#consistent dry
mean_matrix <- A_c_dry*-1 #back to positive values for calculating the matrices
se_matrix <- A_c_dry_error*-1

# Create a list to store the 100 generated matrices
matrices_list <- vector("list", 100)

# Generate 100 random matrices
for (i in 1:100) {
  # Create a matrix of the same dimensions as mean_matrix
  random_matrix <- matrix(0, nrow = nrow(mean_matrix), ncol = ncol(mean_matrix))
  
  # Fill each position with a random value from a normal distribution
  # with the corresponding mean and standard error
  for (row in 1:nrow(mean_matrix)) {
    for (col in 1:ncol(mean_matrix)) {
      random_matrix[row, col] <- rnorm(
        n = 1,
        mean = mean_matrix[row, col],
        sd = se_matrix[row, col]
      )
    }
  }
  
  # Store the matrix in the list
  matrices_list[[i]] <- random_matrix*-1 # back to being negative competition
}

# Print the first 3 matrices as an example
for (i in 1:3) {
  cat("Matrix", i, ":\n")
  print(matrices_list[[i]])
  cat("\n")
}

# Define the mean vector and standard error vector
mean_vector <- r_c_dry
se_vector <- r_c_dry_error

# Check that the vectors have the same length
if (length(mean_vector) != length(se_vector)) {
  stop("Mean vector and standard error vector must have the same length")
}

# Set a seed for reproducibility
set.seed(123)

# Create a list to store the 100 generated vectors
vectors_list <- vector("list", 100)

# Generate 100 random vectors
for (i in 1:100) {
  # Create a vector of the same length as mean_vector
  random_vector <- numeric(length(mean_vector))
  
  # Fill each position with a random value from a normal distribution
  # with the corresponding mean and standard error
  for (j in 1:length(mean_vector)) {
    random_vector[j] <- rnorm(
      n = 1,
      mean = mean_vector[j],
      sd = se_vector[j]
    )
  }
  
  # Store the vector in the list
  vectors_list[[i]] <- random_vector
}

# Print the first 3 vectors as an example
for (i in 1:3) {
  cat("Vector", i, ":\n")
  print(vectors_list[[i]])
  cat("\n")
}

#compute distance to the edge for these 100 vectors and matrices
distances_list_c_dry <- vector("list", 100)
for (i in 1:100) {
  # Calculate the distance to the edge for each matrix and vector
  distances_list_c_dry[[i]] <- calculate_distance_to_border_2sp(matrices_list[[i]], vectors_list[[i]])
}

#Early season fall dry----

# Set a seed for reproducibility
set.seed(123)

# Define the mean matrix and standard error matrix
mean_matrix <- A_esf_dry*-1 #back to positive values for calculating the matrices
se_matrix <- A_esf_dry_error*-1

# Create a list to store the 100 generated matrices
matrices_list <- vector("list", 100)

# Generate 100 random matrices
for (i in 1:100) {
  # Create a matrix of the same dimensions as mean_matrix
  random_matrix <- matrix(0, nrow = nrow(mean_matrix), ncol = ncol(mean_matrix))
  
  # Fill each position with a random value from a normal distribution
  # with the corresponding mean and standard error
  for (row in 1:nrow(mean_matrix)) {
    for (col in 1:ncol(mean_matrix)) {
      random_matrix[row, col] <- rnorm(
        n = 1,
        mean = mean_matrix[row, col],
        sd = se_matrix[row, col]
      )
    }
  }
  
  # Store the matrix in the list
  matrices_list[[i]] <- random_matrix*-1 # back to being negative competition
}

# Print the first 3 matrices as an example
for (i in 1:3) {
  cat("Matrix", i, ":\n")
  print(matrices_list[[i]])
  cat("\n")
}

# Define the mean vector and standard error vector
mean_vector <- r_esf_dry
se_vector <- r_esf_dry_error

# Check that the vectors have the same length
if (length(mean_vector) != length(se_vector)) {
  stop("Mean vector and standard error vector must have the same length")
}

# Set a seed for reproducibility
set.seed(123)

# Create a list to store the 100 generated vectors
vectors_list <- vector("list", 100)

# Generate 100 random vectors
for (i in 1:100) {
  # Create a vector of the same length as mean_vector
  random_vector <- numeric(length(mean_vector))
  
  # Fill each position with a random value from a normal distribution
  # with the corresponding mean and standard error
  for (j in 1:length(mean_vector)) {
    random_vector[j] <- rnorm(
      n = 1,
      mean = mean_vector[j],
      sd = se_vector[j]
    )
  }
  
  # Store the vector in the list
  vectors_list[[i]] <- random_vector
}

# Print the first 3 vectors as an example
for (i in 1:3) {
  cat("Vector", i, ":\n")
  print(vectors_list[[i]])
  cat("\n")
}

#compute distance to the edge for these 100 vectors and matrices
distances_list_esf_dry <- vector("list", 100)
for (i in 1:100) {
  # Calculate the distance to the edge for each matrix and vector
  distances_list_esf_dry[[i]] <- calculate_distance_to_border_2sp(matrices_list[[i]], vectors_list[[i]])
}

#late season spring dry----
# Set a seed for reproducibility
set.seed(123)
# Define the mean matrix and standard error matrix
mean_matrix <- A_lsp_dry*-1 #back to positive values for calculating the matrices
se_matrix <- A_lsp_dry_error*-1
# Create a list to store the 100 generated matrices
matrices_list <- vector("list", 100)
# Generate 100 random matrices
for (i in 1:100) {
  # Create a matrix of the same dimensions as mean_matrix
  random_matrix <- matrix(0, nrow = nrow(mean_matrix), ncol = ncol(mean_matrix))
  
  # Fill each position with a random value from a normal distribution
  # with the corresponding mean and standard error
  for (row in 1:nrow(mean_matrix)) {
    for (col in 1:ncol(mean_matrix)) {
      random_matrix[row, col] <- rnorm(
        n = 1,
        mean = mean_matrix[row, col],
        sd = se_matrix[row, col]
      )
    }
  }
  
  # Store the matrix in the list
  matrices_list[[i]] <- random_matrix*-1 # back to being negative competition
}

# Print the first 3 matrices as an example
for (i in 1:3) {
  cat("Matrix", i, ":\n")
  print(matrices_list[[i]])
  cat("\n")
}

# Define the mean vector and standard error vector
mean_vector <- r_lsp_dry
se_vector <- r_lsp_dry_error
# Check that the vectors have the same length
if (length(mean_vector) != length(se_vector)) {
  stop("Mean vector and standard error vector must have the same length")
}

# Print the first 3 vectors as an example
for (i in 1:3) {
  cat("Vector", i, ":\n")
  print(vectors_list[[i]])
  cat("\n")
}

#compute distance to the edge for these 100 vectors and matrices
distances_list_lsp_dry <- vector("list", 100)
for (i in 1:100) {
  # Calculate the distance to the edge for each matrix and vector
  distances_list_lsp_dry[[i]] <- calculate_distance_to_border_2sp(matrices_list[[i]], vectors_list[[i]])
}

#consistent wet----
# Set a seed for reproducibility
set.seed(123)
# Define the mean matrix and standard error matrix
mean_matrix <- A_c_wet*-1 #back to positive values for calculating the matrices
se_matrix <- A_c_wet_error*-1
# Create a list to store the 100 generated matrices
matrices_list <- vector("list", 100)
# Generate 100 random matrices
for (i in 1:100) {
  # Create a matrix of the same dimensions as mean_matrix
  random_matrix <- matrix(0, nrow = nrow(mean_matrix), ncol = ncol(mean_matrix))
  
  # Fill each position with a random value from a normal distribution
  # with the corresponding mean and standard error
  for (row in 1:nrow(mean_matrix)) {
    for (col in 1:ncol(mean_matrix)) {
      random_matrix[row, col] <- rnorm(
        n = 1,
        mean = mean_matrix[row, col],
        sd = se_matrix[row, col]
      )
    }
  }
  
  # Store the matrix in the list
  matrices_list[[i]] <- random_matrix*-1 # back to being negative competition
}

# Print the first 3 matrices as an example
for (i in 1:3) {
  cat("Matrix", i, ":\n")
  print(matrices_list[[i]])
  cat("\n")
}

# Define the mean vector and standard error vector
mean_vector <- r_c_wet
se_vector <- r_c_wet_error
# Check that the vectors have the same length
if (length(mean_vector) != length(se_vector)) {
  stop("Mean vector and standard error vector must have the same length")
}

# Print the first 3 vectors as an example
for (i in 1:3) {
  cat("Vector", i, ":\n")
  print(vectors_list[[i]])
  cat("\n")
}

#compute distance to the edge for these 100 vectors and matrices
distances_list_c_wet <- vector("list", 100)
for (i in 1:100) {
  # Calculate the distance to the edge for each matrix and vector
  distances_list_c_wet[[i]] <- calculate_distance_to_border_2sp(matrices_list[[i]], vectors_list[[i]])
}


# Helper function for tidy data
make_df <- function(dist_list, label) {
  data.frame(
    distance = sapply(dist_list, function(x) x[[1]]),
    group = label
  )
}

# Build tidy dataframe
df_all <- bind_rows(
  make_df(distances_list_c_dry, "Consistent Dry"),
  make_df(distances_list_esf_dry, "Early season fall dry"),
  make_df(distances_list_lsp_dry, "Late season spring dry"),
  make_df(distances_list_c_wet, "Consistent wet")
)

ggplot(df_all, aes(x = group, y = distance, fill = group)) +
  geom_boxplot(outlier.size = 1, outlier.alpha = 0.7) +
  labs(
    title = "Structural stability of each environment",
    x = " ",
    y = "Distance to the Edge"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_brewer(palette = "Set2") +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )





#code below in case needed. 

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


# Weighted average----
#the paper says
# Across the time series 
#50% of years were consistent dry, 
#12% of years were early-season dry, 
#12% of years were late-season dry, 
#26% of years were consistent wet. 

weights <- c(0.49586778, 0.1157025, 0.1239669, 0.2644628)  # The weights for each value

# average versus relative non-linearity in the lambdas----
#draw the cone
draw_biodiversity_cone(A)
#draw the vector
point <-r/sqrt(sum(r^2))
arrows(0, 0, x1 = point[1], y1 = point[2], length = 0.15, angle = 30,
       code = 2, col = "brown", lty =1, lwd=2)
text(0.35, 1.05, labels = expression(Delta[i]^lambda), cex = 2, pos = 1)

#compute the weighted average
# Create a matrix from the vectors (each column is one of our vectors)
data_matrix <- cbind(r_c_dry, r_esf_dry, r_lsp_dry, r_c_wet)

# Compute weighted average for each row
# This gives us a vector with the same length as each input vector
weighted_avg_vector <- (data_matrix %*% weights) / sum(weights)

# plot the weighted average vector in the cone
point <-weighted_avg_vector/sqrt(sum(weighted_avg_vector^2))
arrows(0, 0, x1 = point[1], y1 = point[2], length = 0.15, angle = 30,
       code = 2, col = "brown", lty =1, lwd=2)
# arrow for point out the change in relative non linearity of the lambda
arrows(x0 = 0.35, y0 = 0.8, x1 = 0.23, y1 = 0.85,
       length = 0.15, angle = 30,
       code = 1, col = "black", lty =1, lwd=2)

# average versus relative non-linearity in the alphas----

#Function to calculate weighted average of matrices
weighted_average_matrices <- function(matrices_list, weights) {
  # Check if dimensions are consistent
  dims <- dim(matrices_list[[1]])
  for (i in 2:length(matrices_list)) {
    if (!all(dim(matrices_list[[i]]) == dims)) {
      stop("All matrices must have the same dimensions")
    }
  }
  
  # Initialize result matrix with zeros
  result <- matrix(0, nrow = dims[1], ncol = dims[2])
  
  # Calculate weighted sum
  for (i in 1:length(matrices_list)) {
    result <- result + weights[i] * matrices_list[[i]]
  }
  
  # Normalize by sum of weights
  result <- result / sum(weights)
  
  return(result)
}

# Put matrices in a list
matrices_list <- list(A_c_dry, A_esf_dry, A_lsp_dry, A_c_wet)

# Calculate weighted average
weighted_matrix <- weighted_average_matrices(matrices_list, weights)
#weighted_matrix[1,2]<-0

# Display the result
print(result_matrix)
#plot both the weighted average and the original
draw_multiple_cones(list(weighted_matrix, A))
text(0.20, 1.09, labels = expression(Delta[i]^alpha), cex = 2, pos = 1)
# arrow for point out the change in relative non linearity of the alphas
arrows(x0 = 0.28, y0 = 0.88, x1 = 0.07, y1 = 0.92,
       length = 0.15, angle = 30,
       code = 2, col = "black", lty =1, lwd=2)

point <-r/sqrt(sum(r^2))
arrows(0, 0, x1 = point[1], y1 = point[2], length = 0.15, angle = 30,
       code = 2, col = "brown", lty =1, lwd=2)

#Storage effect----

#plot both the weighted average and the original
draw_multiple_cones(list(weighted_matrix, A))
text(0.33, 1.08, labels = expression(Delta[i]^{lambda * alpha}), cex = 2, pos = 1)
# arrow for point out the change in relative non linearity of the alphas
arrows(x0 = 0.28, y0 = 0.88, x1 = 0.07, y1 = 0.92,
       length = 0.15, angle = 30,
       code = 2, col = "black", lty =1, lwd=2)

# arrow for point out the change in relative non linearity of the lambdas
point <-r/sqrt(sum(r^2))
arrows(0, 0, x1 = point[1], y1 = point[2], length = 0.15, angle = 30,
       code = 2, col = "brown", lty =1, lwd=2)

point <-weighted_avg_vector/sqrt(sum(weighted_avg_vector^2))
arrows(0, 0, x1 = point[1], y1 = point[2], length = 0.15, angle = 30,
       code = 2, col = "brown", lty =1, lwd=2)
# arrow for point out the change in relative non linearity of the lambda
arrows(x0 = 0.35, y0 = 0.82, x1 = 0.23, y1 = 0.85,
       length = 0.15, angle = 30,
       code = 1, col = "black", lty =1, lwd=2)











