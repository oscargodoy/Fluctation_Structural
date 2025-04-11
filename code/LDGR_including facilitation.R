# Function to calculate low-density growth rate when rare for two competing species
# using a discrete Lotka-Volterra competition model with intraspecific competition
# and the possibility of facilitation (negative alpha values)

calculate_growth_rate_when_rare <- function(r1, r2, alpha11, alpha22, alpha12, alpha21) {
  # Parameters:
  # r1, r2: intrinsic growth rates of species 1 and 2
  # alpha11: intraspecific competition coefficient for species 1 (must be positive)
  # alpha22: intraspecific competition coefficient for species 2 (must be positive)
  # alpha12: effect of species 2 on species 1 (negative = facilitation, positive = competition)
  # alpha21: effect of species 1 on species 2 (negative = facilitation, positive = competition)
  
  # Calculate equilibrium density of species 1 in absence of species 2
  N1_eq <- 1/alpha11
  
  # Calculate equilibrium density of species 2 in absence of species 1
  N2_eq <- 1/alpha22
  
  # Calculate equilibrium when both species are present
  # This is complex with facilitation, but we'll compute it for reference
  denom <- alpha11*alpha22 - alpha12*alpha21
  
  if(denom != 0) {
    N1_coex <- (1 - alpha12/alpha22)/denom
    N2_coex <- (1 - alpha21/alpha11)/denom
  } else {
    N1_coex <- NA
    N2_coex <- NA
  }
  
  # Calculate invasion fitness (growth rate when rare) for species 1
  # when species 2 is at equilibrium
  lambda1 <- 1 + r1 * (1 - alpha11 * 0 - alpha12 * N2_eq)
  
  # Calculate invasion fitness (growth rate when rare) for species 2
  # when species 1 is at equilibrium
  lambda2 <- 1 + r2 * (1 - alpha22 * 0 - alpha21 * N1_eq)
  
  # Return the results
  return(list(
    lambda1 = lambda1,
    lambda2 = lambda2,
    r1 = r1,
    r2 = r2,
    alpha11 = alpha11,
    alpha22 = alpha22,
    alpha12 = alpha12,
    alpha21 = alpha21,
    N1_eq = N1_eq,
    N2_eq = N2_eq,
    N1_coex = N1_coex,
    N2_coex = N2_coex
  ))
}

# Function to run simulations with actual population dynamics
simulate_competition <- function(r1, r2, alpha11, alpha22, alpha12, alpha21, 
                                 N1_init, N2_init, timesteps = 200) {
  
  # Initialize population vectors
  N1 <- numeric(timesteps)
  N2 <- numeric(timesteps)
  
  # Set initial conditions
  N1[1] <- N1_init
  N2[1] <- N2_init
  
  # Run the simulation
  for (t in 1:(timesteps-1)) {
    N1[t+1] <- N1[t] * (1 + r1 * (1 - alpha11 * N1[t] - alpha12 * N2[t]))
    N2[t+1] <- N2[t] * (1 + r2 * (1 - alpha22 * N2[t] - alpha21 * N1[t]))
    
    # Prevent negative populations
    if(N1[t+1] < 0) N1[t+1] <- 0
    if(N2[t+1] < 0) N2[t+1] <- 0
    
    # Check for extinction (very small populations)
    if(N1[t+1] < 1e-6) N1[t+1] <- 0
    if(N2[t+1] < 1e-6) N2[t+1] <- 0
  }
  
  return(list(
    N1_trajectory = N1,
    N2_trajectory = N2,
    final_N1 = N1[timesteps],
    final_N2 = N2[timesteps]
  ))
}

# Function to run invasion analysis
simulate_invasion <- function(r1, r2, alpha11, alpha22, alpha12, alpha21, 
                              initial_rare = 0.001, timesteps = 100) {
  # Calculate equilibrium densities
  N1_eq <- 1/alpha11
  N2_eq <- 1/alpha22
  
  # Simulation 1: Species 1 invading
  N1 <- numeric(timesteps)
  N2 <- numeric(timesteps)
  N1[1] <- initial_rare
  N2[1] <- N2_eq
  
  for (t in 1:(timesteps-1)) {
    N1[t+1] <- N1[t] * (1 + r1 * (1 - alpha11 * N1[t] - alpha12 * N2[t]))
    N2[t+1] <- N2[t] * (1 + r2 * (1 - alpha22 * N2[t] - alpha21 * N1[t]))
    
    # Prevent negative populations
    if(N1[t+1] < 0) N1[t+1] <- 0
    if(N2[t+1] < 0) N2[t+1] <- 0
  }
  
  # Calculate empirical growth rate for species 1 when rare
  emp_lambda1 <- N1[2] / N1[1]
  
  # Simulation 2: Species 2 invading
  M1 <- numeric(timesteps)
  M2 <- numeric(timesteps)
  M1[1] <- N1_eq
  M2[1] <- initial_rare
  
  for (t in 1:(timesteps-1)) {
    M1[t+1] <- M1[t] * (1 + r1 * (1 - alpha11 * M1[t] - alpha12 * M2[t]))
    M2[t+1] <- M2[t] * (1 + r2 * (1 - alpha22 * M2[t] - alpha21 * M1[t]))
    
    # Prevent negative populations
    if(M1[t+1] < 0) M1[t+1] <- 0
    if(M2[t+1] < 0) M2[t+1] <- 0
  }
  
  # Calculate empirical growth rate for species 2 when rare
  emp_lambda2 <- M2[2] / M2[1]
  
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
analyze_coexistence <- function(lambda1, lambda2, N1_coex, N2_coex) {
  # Basic invasion analysis
  if (lambda1 > 1 && lambda2 > 1) {
    coex_type <- "Stable coexistence: Both species can invade when rare"
  } else if (lambda1 > 1 && lambda2 < 1) {
    coex_type <- "Competitive exclusion: Species 1 excludes Species 2"
  } else if (lambda1 < 1 && lambda2 > 1) {
    coex_type <- "Competitive exclusion: Species 2 excludes Species 1"
  } else {
    coex_type <- "Priority effects: Outcome depends on initial conditions"
  }
  
  # Check feasibility of coexistence equilibrium
  if (!is.na(N1_coex) && !is.na(N2_coex)) {
    if (N1_coex <= 0 || N2_coex <= 0) {
      coex_note <- "Note: Theoretical coexistence equilibrium is not biologically feasible."
      if (lambda1 > 1 && lambda2 > 1) {
        coex_type <- "Despite mutual invasibility, coexistence is not stable due to negative equilibrium densities."
      }
    } else {
      coex_note <- paste("Theoretical coexistence equilibrium: Species 1 =", round(N1_coex, 3), 
                         ", Species 2 =", round(N2_coex, 3))
    }
  } else {
    coex_note <- "Cannot determine coexistence equilibrium (singular matrix)."
  }
  
  return(list(
    outcome = coex_type,
    note = coex_note
  ))
}

# Function to visualize dynamics over time
visualize_dynamics <- function(sim_results, title="Population Dynamics") {
  plot(1:length(sim_results$N1_trajectory), 
       sim_results$N1_trajectory, 
       type = "l", col = "blue", lwd=2,
       xlab = "Time", ylab = "Population Size", 
       main = title,
       ylim = c(0, max(c(sim_results$N1_trajectory, sim_results$N2_trajectory))*1.1))
  
  lines(1:length(sim_results$N2_trajectory), 
        sim_results$N2_trajectory, 
        col = "red", lwd=2)
  
  legend("topright", 
         legend = c("Species 1", "Species 2"), 
         col = c("blue", "red"), 
         lty = 1, lwd=2)
}

# Run analysis with the given parameters
r1 <- 2.2        # Intrinsic growth rate of species 1
r2 <- 0.3        # Intrinsic growth rate of species 2
alpha11 <- 0.5   # Intraspecific competition for species 1
alpha22 <- 0.5   # Intraspecific competition for species 2
alpha12 <- 0.2  # Effect of species 2 on species 1 (negative = facilitation)
alpha21 <- 0.25  # Effect of species 1 on species 2

# Print parameters
cat("Model Parameters:\n")
cat(paste("r1 =", r1, ", r2 =", r2, "\n"))
cat(paste("alpha11 =", alpha11, ", alpha22 =", alpha22, "\n"))
cat(paste("alpha12 =", alpha12, " (", ifelse(alpha12 < 0, "facilitation", "competition"), ")\n", sep=""))
cat(paste("alpha21 =", alpha21, " (", ifelse(alpha21 < 0, "facilitation", "competition"), ")\n", sep=""))

# Calculate analytical growth rates when rare
analytical_results <- calculate_growth_rate_when_rare(r1, r2, alpha11, alpha22, alpha12, alpha21)

cat("\nEquilibrium densities in monoculture:\n")
cat(paste("Species 1 alone:", round(analytical_results$N1_eq, 3), "\n"))
cat(paste("Species 2 alone:", round(analytical_results$N2_eq, 3), "\n"))

cat("\nInvasion analysis:\n")
cat(paste("Growth rate of species 1 when rare (lambda1):", round(analytical_results$lambda1, 3), 
          ifelse(analytical_results$lambda1 > 1, " (can invade)", " (cannot invade)"), "\n"))
cat(paste("Growth rate of species 2 when rare (lambda2):", round(analytical_results$lambda2, 3),
          ifelse(analytical_results$lambda2 > 1, " (can invade)", " (cannot invade)"), "\n"))

# Analyze coexistence
coexistence_analysis <- analyze_coexistence(
  analytical_results$lambda1, 
  analytical_results$lambda2,
  analytical_results$N1_coex,
  analytical_results$N2_coex
)

cat("\nCoexistence analysis:\n")
cat(paste(coexistence_analysis$outcome, "\n"))
cat(paste(coexistence_analysis$note, "\n"))

# Run full simulation with equal starting populations
sim_results <- simulate_competition(
  r1, r2, alpha11, alpha22, alpha12, alpha21,
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

# Additional simulation with species 2 having an advantage
sim_results2 <- simulate_competition(
  r1, r2, alpha11, alpha22, alpha12, alpha21,
  N1_init = 0.1, 
  N2_init = 2.0,
  timesteps = 200
)

# Determine the outcome of second simulation
if(sim_results2$final_N1 > 0.01 && sim_results2$final_N2 < 0.01) {
  outcome2 <- "Species 1 wins (Species 2 goes extinct)"
} else if(sim_results2$final_N1 < 0.01 && sim_results2$final_N2 > 0.01) {
  outcome2 <- "Species 2 wins (Species 1 goes extinct)"
} else if(sim_results2$final_N1 > 0.01 && sim_results2$final_N2 > 0.01) {
  outcome2 <- "Coexistence"
} else {
  outcome2 <- "Both species go extinct"
}

cat("\nSimulation with species 2 starting advantage:\n")
cat(paste("Final species 1:", round(sim_results2$final_N1, 3), "\n"))
cat(paste("Final species 2:", round(sim_results2$final_N2, 3), "\n"))
cat(paste("Outcome:", outcome2, "\n"))

# Plot second scenario
visualize_dynamics(sim_results2, "Population Dynamics (Species 2 Starting Advantage)")

