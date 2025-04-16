# Function to calculate low-density growth rate when rare for two competing species
# using a discrete Lotka-Volterra competition model with intraspecific competition

calculate_growth_rate_when_rare <- function(ra, rb, alpha11, alpha22, alpha12, alpha21) {
  # Parameters:
  # ra, rb: intrinsic growth rates of species 1 and 2
  # alpha11: intraspecific competition coefficient for species 1
  # alpha22: intraspecific competition coefficient for species 2
  # alpha12: competitive effect of species 2 on species 1 (interspecific)
  # alpha21: competitive effect of species 1 on species 2 (interspecific)
  
  # Calculate equilibrium density of species 1 in absence of species 2
  N1_eq <- ra/alpha11
  
  # Calculate equilibrium density of species 2 in absence of species 1
  N2_eq <- rb/alpha22
  
  # Calculate invasion fitness (growth rate when rare) for species 1
  # when species 2 is at equilibrium
  lambda1 <- 1 + ra * (1 - alpha12 * N2_eq)
  
  # Calculate invasion fitness (growth rate when rare) for species 2
  # when species 1 is at equilibrium
  lambda2 <- 1 + rb * (1 - alpha21 * N1_eq)
  
  # Return the results
  return(list(
    lambda1 = lambda1,
    lambda2 = lambda2,
    ra = ra,
    rb = rb,
    alpha11 = alpha11,
    alpha22 = alpha22,
    alpha12 = alpha12,
    alpha21 = alpha21,
    N1_eq = N1_eq,
    N2_eq = N2_eq
  ))
}

# Function to run simulations to verify growth rates
simulate_invasion <- function(ra, rb, alpha11, alpha22, alpha12, alpha21, 
                              initial_rare = 0.01, timesteps = 1000) {
  # Parameters for simulation
  # initial_rare: initial population size of the invader
  # timesteps: number of time steps to simulate
  
  # Calculate equilibrium densities
  N1_eq <- ra/alpha11
  N2_eq <- rb/alpha22
  
  # Run two simulations: 
  # 1. Species 1 invading Species 2 at equilibrium
  # 2. Species 2 invading Species 1 at equilibrium
  
  # Simulation 1: Species 1 invading
  N1 <- numeric(timesteps)
  N2 <- numeric(timesteps)
  N1[1] <- initial_rare
  N2[1] <- N2_eq
  
  for (t in 1:(timesteps-1)) {
    N1[t+1] <- N1[t] * (1 + ra * (1 - alpha11 * N1[t] - alpha12 * N2[t]))
    N2[t+1] <- N2[t] * (1 + rb * (1 - alpha22 * N2[t] - alpha21 * N1[t]))
  }
  
  # Calculate empirical growth rate for species 1 when rare
  emp_lambda1 <- N1[2] / N1[1]
  
  # Simulation 2: Species 2 invading
  M1 <- numeric(timesteps)
  M2 <- numeric(timesteps)
  M1[1] <- N1_eq
  M2[1] <- initial_rare
  
  for (t in 1:(timesteps-1)) {
    M1[t+1] <- M1[t] * (1 + ra * (1 - alpha11 * M1[t] - alpha12 * M2[t]))
    M2[t+1] <- M2[t] * (1 + rb * (1 - alpha22 * M2[t] - alpha21 * M1[t]))
  }
  
  # Calculate empirical growth rate for species 2 when rare
  emp_lambda2 <- M2[2] / M2[1]
  
  # Return the results
  return(list(
    emp_lambda1 = emp_lambda1,
    emp_lambda2 = emp_lambda2,
    N1_trajectory = N1,
    N2_trajectory = N2,
    M1_trajectory = M1,
    M2_trajectory = M2
  ))
}

# Function to analyze coexistence criteria based on invasion growth rates
analyze_coexistence <- function(lambda1, lambda2) {
  if (lambda1 > 1 && lambda2 > 1) {
    return("Stable coexistence: Both species can invade when rare")
  } else if (lambda1 > 1 && lambda2 < 1) {
    return("Competitive exclusion: Species 1 excludes Species 2")
  } else if (lambda1 < 1 && lambda2 > 1) {
    return("Competitive exclusion: Species 2 excludes Species 1")
  } else {
    return("Priority effects: Outcome depends on initial conditions")
  }
}

# Function to visualize dynamics over time for both invasion scenarios
visualize_dynamics <- function(sim_results) {
  par(mfrow = c(2, 1))
  
  # Species 1 invading
  plot(1:length(sim_results$N1_trajectory), 
       sim_results$N1_trajectory, 
       type = "l", col = "blue", 
       xlab = "Time", ylab = "Population Size", 
       main = "Species 1 (blue) invading Species 2 (red)")
  lines(1:length(sim_results$N2_trajectory), 
        sim_results$N2_trajectory, 
        col = "red")
  legend("topright", legend = c("Species 1 (invader)", "Species 2 (resident)"), 
         col = c("blue", "red"), lty = 1)
  
  # Species 2 invading
  plot(1:length(sim_results$M2_trajectory), 
       sim_results$M2_trajectory, 
       type = "l", col = "red", 
       xlab = "Time", ylab = "Population Size", 
       main = "Species 2 (red) invading Species 1 (blue)")
  lines(1:length(sim_results$M1_trajectory), 
        sim_results$M1_trajectory, 
        col = "blue")
  legend("topright", legend = c("Species 1 (resident)", "Species 2 (invader)"), 
         col = c("blue", "red"), lty = 1)
}

# Example usage with intraspecific competition
# Define parameters
ra <- 3.5        # Intrinsic growth rate of species 1
rb <- 1        # Intrinsic growth rate of species 2
alpha11 <- 0.75   # Intraspecific competition for species 1
alpha22 <- 0.75   # Intraspecific competition for species 2
alpha12 <- 0.25  # Competitive effect of species 2 on species 1
alpha21 <- 0.25  # Competitive effect of species 1 on species 2

# Calculate analytical growth rates when rare
analytical_results <- calculate_growth_rate_when_rare(ra, rb, alpha11, alpha22, alpha12, alpha21)
print("Analytical results:")
print(paste("Equilibrium density of species 1 alone:", analytical_results$N1_eq))
print(paste("Equilibrium density of species 2 alone:", analytical_results$N2_eq))
print(paste("Growth rate of species 1 when rare (lambda1):", analytical_results$lambda1))
print(paste("Growth rate of species 2 when rare (lambda2):", analytical_results$lambda2))

# Verify with simulations
simulation_results <- simulate_invasion(ra, rb, alpha11, alpha22, alpha12, alpha21)
print("Simulation results:")
print(paste("Empirical growth rate of species 1 when rare:", simulation_results$emp_lambda1))
print(paste("Empirical growth rate of species 2 when rare:", simulation_results$emp_lambda2))

# Visualize the dynamics
visualize_dynamics(simulation_results)

# Analyze coexistence for our example
coexistence_outcome <- analyze_coexistence(analytical_results$lambda1, analytical_results$lambda2)
print(paste("Coexistence analysis:", coexistence_outcome))

dev.off()

# Additional function to study impact of varying competition strengths on coexistence
study_parameter_space <- function(ra, rb, alpha11, alpha22, alpha12_range, alpha21_range) {
  n <- length(alpha12_range)
  m <- length(alpha21_range)
  
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
      alpha12 <- alpha12_range[i]
      alpha21 <- alpha21_range[j]
      
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
    alpha12_range = alpha12_range,
    alpha21_range = alpha21_range,
    outcomes = outcomes
  ))
}

# Example parameter space exploration
alpha12_range <- seq(-0.3, 0.8, length.out = 20)
alpha21_range <- seq(-0.3, 0.8, length.out = 20)

param_space <- study_parameter_space(ra, rb, alpha11, alpha22, alpha12_range, alpha21_range)

# Visualize parameter space
image(param_space$alpha12_range, 
      param_space$alpha21_range, 
      param_space$results, 
      col = c("gray", "blue", "red", "green"),
      xlab = "alpha12 (effect of species 2 on 1)",
      ylab = "alpha21 (effect of species 1 on 2)",
      main = "Competition outcomes")

# Add a legend
outcome_names <- names(param_space$outcomes)
outcome_colors <- c("gray", "blue", "red", "green")
legend("topright", 
       legend = outcome_names, 
       fill = outcome_colors, 
       title = "Outcomes")

