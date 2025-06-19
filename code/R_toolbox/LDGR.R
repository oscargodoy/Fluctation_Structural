calculate_growth_rate_when_rare_corrected <- function(ra, rb, alpha11, alpha22, alpha12, alpha21) {
  # Parameters:
  # ra, rb: intrinsic growth rates of species 1 and 2
  # alpha11: intraspecific competition coefficient for species 1
  # alpha22: intraspecific competition coefficient for species 2
  # alpha12: competitive effect of species 2 on species 1 (interspecific)
  # alpha21: competitive effect of species 1 on species 2 (interspecific)
  
  # Calculate equilibrium density of species 1 in absence of species 2
  N1_eq <- - ra/alpha11
  
  # Calculate equilibrium density of species 2 in absence of species 1
  N2_eq <- - rb/alpha22
  
  # Calculate invasion fitness (growth rate when rare) for species 1
  # when species 2 is at equilibrium
  lambda1 <- ra + alpha12 * N2_eq
  
  # Calculate invasion fitness (growth rate when rare) for species 2
  # when species 1 is at equilibrium
  lambda2 <- rb + alpha21 * N1_eq
  
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

# Additional function to study impact of varying vectors of intrinsic growth rate on LDGR and distance to the edge
study_parameter_space <- function(ra_range, rb_range, alpha11, alpha22, alpha12, alpha21) {
  n <- length(ra_range)
  m <- length(rb_range)
  
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
      ra <- ra_range[i]
      rb <- rb_range[j]
      
      res <- calculate_growth_rate_when_rare_corrected(ra, rb, alpha11, alpha22, alpha12, alpha21)
      lambda1 <- res$lambda1
      lambda2 <- res$lambda2
      
      lambda1_matrix[i,j] <- lambda1
      lambda2_matrix[i,j] <- lambda2
      
      # Determine outcome based on invasion criteria
      if (lambda1 > 0 && lambda2 > 0) {
        results[i,j] <- outcomes["Coexistence"]
      } else if (lambda1 > 0 && lambda2 < 0) {
        results[i,j] <- outcomes["Species 1 wins"]
      } else if (lambda1 < 0 && lambda2 > 0) {
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
    ra_range = ra_range,
    rb_range = rb_range,
    outcomes = outcomes
  ))
}

norm_vec <- function(x){
  return(x/sqrt(sum(x^2)))
}