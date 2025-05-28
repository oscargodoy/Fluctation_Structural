# Load necessary libraries
library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)

# Function implementations
# First method - calculate_distance_to_border_2sp
# This is Oscar's initial version - focused on the angle
# prueba commit 

calculate_distance_to_border_2sp <- function(A_int, r, norm = "no"){
  norm_vec <- function(x){
    return(x/sqrt(sum(x^2)))
  }
  
  arc_length <- function(a, b) {
    acos(sum(a * b))
  }
  
  if(length(r) != 2 || length(A_int) != 4){
    cat("This function is only for communities of 2 species. \n")
  } else {
    r_norm <- norm_vec(r)
    A1 <- - norm_vec(A_int[,1])
    A2 <- - norm_vec(A_int[,2])
    aux <- solve(-A_int,r)
    
    if(all(aux>0)){
      distance <- min(arc_length(r_norm,A1),arc_length(r_norm,A2))
      outside <- 0
    }else {
      distance <- - min(arc_length(r_norm,A1),arc_length(r_norm,A2))
      outside <- 1
    }
    
    if(norm == "yes"){
      distance_norm <- distance / pi
      return(list(distance_norm,outside))
    } else if(norm == "no"){
      return(list(distance,outside))
    } else {
      cat("Norm only takes values yes and no. \n")
    }
  }
}

# Second method - calculate_distance_to_edge
# This is Lauren and Claude's version
# Focused on Eucidean distances
# Need to confirm that distance is being calculated as the shortest path to the feasibility domain

calculate_distance_to_edge <- function(r1, r2, alpha11, alpha12, alpha21, alpha22) {
  # Calculate the slopes of the boundaries
  left_slope <- alpha22 / alpha12
  right_slope <- alpha21 / alpha11
  
  # Calculate the equilibrium values (N1*, N2*)
  det <- alpha11 * alpha22 - alpha12 * alpha21
  N2_star <- (r1 * alpha21 - r2 * alpha11) / det
  N1_star <- (-r1 - N2_star * alpha12) /  alpha11
  
  # Check if point is in feasibility domain (both N1* and N2* > 0)
  if (N1_star > 0 && N2_star > 0) {
    # Inside the feasibility domain
    # Calculate distances to both boundaries and return the minimum
    # We'll use the perpendicular distance from point (r1, r2) to each boundary line
    
    # For the left boundary (with slope left_slope)
    # Line equation: r2 = left_slope * r1
    # Distance formula: |r2 - left_slope * r1| / sqrt(1 + left_slope^2)
    dist_left <- abs(r2 - left_slope * r1) / sqrt(1 + left_slope^2)
    
    # For the right boundary (with slope right_slope)
    # Line equation: r2 = right_slope * r1
    # Distance formula: |r2 - right_slope * r1| / sqrt(1 + right_slope^2)
    dist_right <- abs(r2 - right_slope * r1) / sqrt(1 + right_slope^2)
    
    # Return the minimum distance to any boundary
    return(min(dist_left, dist_right))
  } else {
    # Outside the feasibility domain
    # Calculate the distance as negative
    # First, determine which boundary is closest
    
    # Calculate closest point on left boundary
    if (left_slope == 0) {
      # If left boundary is horizontal
      closest_left <- c(r1, 0)
    } else if (is.infinite(left_slope)) {
      # If left boundary is vertical
      closest_left <- c(0, r2)
    } else {
      # General case
      # Closest point on line r2 = left_slope * r1 is given by:
      x <- (r1 + left_slope * r2) / (1 + left_slope^2)
      y <- left_slope * x
      closest_left <- c(x, y)
    }
    
    # Calculate closest point on right boundary
    if (right_slope == 0) {
      # If right boundary is horizontal
      closest_right <- c(r1, 0)
    } else if (is.infinite(right_slope)) {
      # If right boundary is vertical
      closest_right <- c(0, r2)
    } else {
      # General case
      # Closest point on line r2 = right_slope * r1 is given by:
      x <- (r1 + right_slope * r2) / (1 + right_slope^2)
      y <- right_slope * x
      closest_right <- c(x, y)
    }
    
    # Calculate distances to both closest points
    dist_left <- sqrt((r1 - closest_left[1])^2 + (r2 - closest_left[2])^2)
    dist_right <- sqrt((r1 - closest_right[1])^2 + (r2 - closest_right[2])^2)
    
    # Return the negative of the minimum distance
    return(-min(dist_left, dist_right))
  }
}

# Function to wrap the first method so it returns a single numeric value
calculate_distance_method1 <- function(r1, r2, alpha11, alpha12, alpha21, alpha22) {
  A_int <- matrix(c(alpha11, alpha21, alpha12, alpha22), nrow=2)
  r <- c(r1, r2)
  result <- calculate_distance_to_border_2sp(A_int, r)
  return(result[[1]])
}

# LDGR calculation function
calculate_ldgr <- function(r1, r2, alpha11, alpha12, alpha21, alpha22) {
  # Calculate equilibrium densities when each species is alone
  N1_alone <- r1 / alpha11  # Equilibrium of species 1 when alone
  N2_alone <- r2 / alpha22  # Equilibrium of species 2 when alone
  
  # Calculate LDGR for species 1 (when introduced to equilibrium of species 2)
  LDGR1 <- r1 - alpha12 * N2_alone
  
  # Calculate LDGR for species 2 (when introduced to equilibrium of species 1)
  LDGR2 <- r2 - alpha21 * N1_alone
  
  return(list(LDGR1 = LDGR1, LDGR2 = LDGR2))
}


# Function to analyze coexistence metrics for given parameters
analyze_coexistence <- function(r1, r2, alpha11, alpha12, alpha21, alpha22) {
  # Calculate distances using both methods
  distance1 <- calculate_distance_method1(r1, r2, alpha11, alpha12, alpha21, alpha22)
  distance2 <- calculate_distance_to_edge(r1, r2, alpha11, alpha12, alpha21, alpha22)
  
  # Calculate LDGR
  ldgr_results <- calculate_ldgr(r1, r2, alpha11, alpha12, alpha21, alpha22)
  LDGR1 <- ldgr_results$LDGR1
  LDGR2 <- ldgr_results$LDGR2
  min_LDGR <- min(LDGR1, LDGR2)
  
  # Calculate the slopes of the feasibility domain boundaries
  left_slope <- alpha22 / alpha12
  right_slope <- alpha21 / alpha11
  
  # Return results
  return(data.frame(
    r1 = r1,
    r2 = r2,
    alpha11 = alpha11,
    alpha22 = alpha22,
    alpha12 = alpha12,
    alpha21 = alpha21,
    distance1 = distance1,
    distance2 = distance2,
    LDGR1 = LDGR1,
    LDGR2 = LDGR2,
    min_LDGR = min_LDGR,
    left_slope = left_slope,
    right_slope = right_slope,
    coexistence = (LDGR1 > 0 & LDGR2 > 0)
  ))
}

#################################################
# Scenario 1: Fixed parameters except varying r1
#################################################
scenario1_fixed_r2 <- function() {
  # Fixed parameters
  alpha11 <- 0.075
  alpha22 <- 0.075
  alpha12 <- 0.05
  alpha21 <- 0.05
  r2 <- 1.5
  
  # Varying r1
  r1_values <- seq(-1, 3.8, by = 0.05)
  
  results <- data.frame()
  for (r1 in r1_values) {
    result <- analyze_coexistence(r1, r2, alpha11, alpha12, alpha21, alpha22)
    results <- rbind(results, result)
  }
  
  results$scenario <- "varying_r1_fixed_r2"
  return(results)
}

#################################################
# Scenario 2: Fixed parameters except varying r2
#################################################
scenario2_fixed_r1 <- function() {
  # Fixed parameters
  alpha11 <- 0.075
  alpha22 <- 0.075
  alpha12 <- 0.05
  alpha21 <- 0.05
  r1 <- 1.5
  
  # Varying r2
  r2_values <- seq(-1, 3.8, by = 0.05)
  
  results <- data.frame()
  for (r2 in r2_values) {
    result <- analyze_coexistence(r1, r2, alpha11, alpha12, alpha21, alpha22)
    results <- rbind(results, result)
  }
  
  results$scenario <- "varying_r2_fixed_r1"
  return(results)
}

#################################################
# Scenario 3: Fixed parameters except varying alpha12
#################################################
scenario3_varying_alpha12 <- function() {
  # Fixed parameters
  alpha11 <- 0.075
  alpha22 <- 0.075
  alpha21 <- 0.05
  r1 <- 1.5
  r2 <- 1.5
  
  # Varying alpha12
  alpha12_values <- seq(-0.15, 0.15, by = 0.01)
  
  results <- data.frame()
  for (alpha12 in alpha12_values) {
    result <- analyze_coexistence(r1, r2, alpha11, alpha12, alpha21, alpha22)
    results <- rbind(results, result)
  }
  
  results$scenario <- "varying_alpha12"
  return(results)
}

#################################################
# Scenario 4: Fixed parameters except varying alpha21
#################################################
scenario4_varying_alpha21 <- function() {
  # Fixed parameters
  alpha11 <- 0.075
  alpha22 <- 0.075
  alpha12 <- 0.05
  r1 <- 1.5
  r2 <- 1.5
  
  # Varying alpha21
  alpha21_values <- seq(-0.15, 0.15, by = 0.01)
  
  results <- data.frame()
  for (alpha21 in alpha21_values) {
    result <- analyze_coexistence(r1, r2, alpha11, alpha12, alpha21, alpha22)
    results <- rbind(results, result)
  }
  
  results$scenario <- "varying_alpha21"
  return(results)
}

#################################################
# Scenario 5: Fixed parameters except varying alpha11
#################################################
scenario5_varying_alpha11 <- function() {
  # Fixed parameters
  alpha22 <- 0.075
  alpha12 <- 0.05
  alpha21 <- 0.05
  r1 <- 1.5
  r2 <- 1.5
  
  # Varying alpha11
  alpha11_values <- seq(0.01, 0.15, by = 0.01)  # Start at 0.01 to avoid division by zero
  
  results <- data.frame()
  for (alpha11 in alpha11_values) {
    result <- analyze_coexistence(r1, r2, alpha11, alpha12, alpha21, alpha22)
    results <- rbind(results, result)
  }
  
  results$scenario <- "varying_alpha11"
  return(results)
}

#################################################
# Scenario 6: Fixed parameters except varying alpha22
#################################################
scenario6_varying_alpha22 <- function() {
  # Fixed parameters
  alpha11 <- 0.075
  alpha12 <- 0.05
  alpha21 <- 0.05
  r1 <- 1.5
  r2 <- 1.5
  
  # Varying alpha22
  alpha22_values <- seq(0.01, 0.15, by = 0.01)  # Start at 0.01 to avoid division by zero
  
  results <- data.frame()
  for (alpha22 in alpha22_values) {
    result <- analyze_coexistence(r1, r2, alpha11, alpha12, alpha21, alpha22)
    results <- rbind(results, result)
  }
  
  results$scenario <- "varying_alpha22"
  return(results)
}

#################################################
# Run all scenarios and combine results
#################################################
run_parameter_sweep <- function() {
  cat("Running scenario 1: Varying r1 with fixed r2...\n")
  results_scenario1 <- scenario1_fixed_r2()
  
  cat("Running scenario 2: Varying r2 with fixed r1...\n")
  results_scenario2 <- scenario2_fixed_r1()
  
  cat("Running scenario 3: Varying interspecific competition alpha12...\n")
  results_scenario3 <- scenario3_varying_alpha12()
  
  cat("Running scenario 4: Varying interspecific competition alpha21...\n")
  results_scenario4 <- scenario4_varying_alpha21()
  
  cat("Running scenario 5: Varying intraspecific competition alpha11...\n")
  results_scenario5 <- scenario5_varying_alpha11()
  
  cat("Running scenario 6: Varying intraspecific competition alpha22...\n")
  results_scenario6 <- scenario6_varying_alpha22()
  
  # Combine all results
  all_results <- rbind(
    results_scenario1,
    results_scenario2,
    results_scenario3,
    results_scenario4,
    results_scenario5,
    results_scenario6
  )
  
  return(all_results)
}

#################################################
# Run parameter sweep and analyze results
#################################################
set.seed(42)
all_results <- run_parameter_sweep()

#################################################
# Create plots for each scenario
#################################################

# Function to create plot pairs for each scenario
create_scenario_plots <- function(data, scenario_name, varying_param) {
  scenario_data <- data %>% filter(scenario == scenario_name)
  
  # Plot for Method 1 vs. LDGR
  p1 <- ggplot(scenario_data, aes(x=min_LDGR, y=distance1, color=!!sym(varying_param))) +
    geom_point(size=3) +
    geom_path(aes(group=1), color="grey50", alpha=0.5) +
    geom_hline(yintercept=0, linetype="dashed", color="red") +
    geom_vline(xintercept=0, linetype="dashed", color="red") +
    labs(title=paste("Arc Distance vs. min LDGR -", gsub("_", " ", scenario_name)),
         x="Minimum LDGR",
         y="Arc Distance (Method 1)",
         color=gsub("_", " ", varying_param)) +
    theme_minimal()  + geom_vline(xintercept = 0)
  
  # Plot for Method 2 vs. LDGR
  p2 <- ggplot(scenario_data, aes(x=min_LDGR, y=distance2, color=!!sym(varying_param))) +
    geom_point(size=3) +
    geom_path(aes(group=1), color="grey50", alpha=0.5) +
    geom_hline(yintercept=0, linetype="dashed", color="red") +
    geom_vline(xintercept=0, linetype="dashed", color="red") +
    labs(title=paste("Euclidean Distance vs. min LDGR -", gsub("_", " ", scenario_name)),
         x="Minimum LDGR",
         y="Euclidean Distance (Method 2)",
         color=gsub("_", " ", varying_param)) +
    theme_minimal()  + geom_vline(xintercept = 0)
  
  # Create parameter vs. metrics plot
  p3 <- scenario_data %>%
    select(!!sym(varying_param), min_LDGR, distance1, distance2) %>%
    reshape2::melt(id.vars = varying_param) %>%
    ggplot(aes(x=!!sym(varying_param), y=value, color=variable)) +
    geom_line(size=1) +
    geom_hline(yintercept=0, linetype="dashed", color="red") +
    labs(title=paste("Metrics vs.", gsub("_", " ", varying_param)),
         x=gsub("_", " ", varying_param),
         y="Metric Value",
         color="Metric") +
    scale_color_manual(values=c("min_LDGR"="blue", "distance1"="red", "distance2"="green"),
                       labels=c("min_LDGR"="Min LDGR", 
                                "distance1"="Arc Distance", 
                                "distance2"="Euclidean Distance")) +
    theme_minimal() + geom_vline(xintercept = 0)
  
  return(list(p1, p2, p3))
}

# Create plots for each scenario
plots_scenario1 <- create_scenario_plots(all_results, "varying_r1_fixed_r2", "r1")
plots_scenario2 <- create_scenario_plots(all_results, "varying_r2_fixed_r1", "r2")
plots_scenario3 <- create_scenario_plots(all_results, "varying_alpha12", "alpha12")
plots_scenario4 <- create_scenario_plots(all_results, "varying_alpha21", "alpha21")
plots_scenario5 <- create_scenario_plots(all_results, "varying_alpha11", "alpha11")
plots_scenario6 <- create_scenario_plots(all_results, "varying_alpha22", "alpha22")

plots_scenario1 + geom_vline(xintercept = 0)

#################################################
# Add correlation and regression analysis
#################################################

# Function to analyze relationships for each scenario
analyze_scenario_relationships <- function(data) {
  scenarios <- unique(data$scenario)
  
  # Initialize result dataframe
  relationship_stats <- data.frame(
    scenario = character(),
    varying_param = character(),
    cor_minLDGR_dist1 = numeric(),
    cor_minLDGR_dist2 = numeric(),
    r2_minLDGR_dist1 = numeric(),
    r2_minLDGR_dist2 = numeric(),
    linearity_dist1 = character(),
    linearity_dist2 = character(),
    stringsAsFactors = FALSE
  )
  
  for (s in scenarios) {
    scenario_data <- data %>% filter(scenario == s)
    
    # Determine which parameter is varying
    varying_param <- switch(s,
                            "varying_r1_fixed_r2" = "r1",
                            "varying_r2_fixed_r1" = "r2",
                            "varying_alpha12" = "alpha12",
                            "varying_alpha21" = "alpha21",
                            "varying_alpha11" = "alpha11",
                            "varying_alpha22" = "alpha22")
    
    # Calculate correlations
    cor_minLDGR_dist1 <- cor(scenario_data$min_LDGR, scenario_data$distance1)
    cor_minLDGR_dist2 <- cor(scenario_data$min_LDGR, scenario_data$distance2)
    
    # Fit linear models
    lm1 <- lm(distance1 ~ min_LDGR, data=scenario_data)
    lm2 <- lm(distance2 ~ min_LDGR, data=scenario_data)
    
    r2_1 <- summary(lm1)$r.squared
    r2_2 <- summary(lm2)$r.squared
    
    # Check for linearity
    # We'll call R² > 0.95 as "strongly linear"
    linearity_dist1 <- ifelse(r2_1 > 0.95, "Strongly Linear", 
                              ifelse(r2_1 > 0.8, "Moderately Linear", "Non-linear"))
    linearity_dist2 <- ifelse(r2_2 > 0.95, "Strongly Linear", 
                              ifelse(r2_2 > 0.8, "Moderately Linear", "Non-linear"))
    
    # Add to results
    relationship_stats <- rbind(relationship_stats, data.frame(
      scenario = s,
      varying_param = varying_param,
      cor_minLDGR_dist1 = cor_minLDGR_dist1,
      cor_minLDGR_dist2 = cor_minLDGR_dist2,
      r2_minLDGR_dist1 = r2_1,
      r2_minLDGR_dist2 = r2_2,
      linearity_dist1 = linearity_dist1,
      linearity_dist2 = linearity_dist2,
      stringsAsFactors = FALSE
    ))
  }
  
  return(relationship_stats)
}

# Analyze relationships for all scenarios
relationship_stats <- analyze_scenario_relationships(all_results)

#################################################
# Generate comprehensive visualizations
#################################################

# Function to create feasibility domain plots with LDGR and distance contours
create_feasibility_plot <- function(alpha11, alpha22, alpha12, alpha21) {
  # Create grid of r1, r2 values
  r1_seq <- seq(-1, 3.8, by = 0.1)
  r2_seq <- seq(-1, 3.8, by = 0.1)
  grid <- expand.grid(r1=r1_seq, r2=r2_seq)
  
  # Calculate metrics for each point
  results <- data.frame(
    r1 = grid$r1,
    r2 = grid$r2,
    distance1 = NA,
    distance2 = NA,
    min_LDGR = NA,
    coexistence = NA
  )
  
  for (i in 1:nrow(grid)) {
    r1 <- grid$r1[i]
    r2 <- grid$r2[i]
    
    # Calculate distances
    results$distance1[i] <- calculate_distance_method1(r1, r2, alpha11, alpha12, alpha21, alpha22)
    results$distance2[i] <- calculate_distance_to_edge(r1, r2, alpha11, alpha12, alpha21, alpha22)
    
    # Calculate LDGR
    ldgr_res <- calculate_ldgr(r1, r2, alpha11, alpha12, alpha21, alpha22)
    results$min_LDGR[i] <- min(ldgr_res$LDGR1, ldgr_res$LDGR2)
    
    # Determine coexistence
    results$coexistence[i] <- ldgr_res$LDGR1 > 0 && ldgr_res$LDGR2 > 0
  }
  
  # Calculate the slopes of the boundaries
  left_slope <- alpha22 / alpha12
  right_slope <- alpha21 / alpha11
  
  # Create plots
  # Plot 1: Feasibility domain with LDGR contours
  p1 <- ggplot(results, aes(x=r1, y=r2)) +
    geom_tile(aes(fill=min_LDGR)) +
    geom_contour(aes(z=min_LDGR), breaks=0, color="black", linewidth=1) +
    geom_abline(slope=left_slope, intercept=0, linetype="dashed", color="white") +
    geom_abline(slope=right_slope, intercept=0, linetype="dashed", color="white") +
    scale_fill_gradient2(low="red", mid="white", high="blue", midpoint=0) +
    labs(title="Minimum LDGR in r1-r2 Space",
         subtitle=paste("α11=", alpha11, ", α12=", alpha12, ", α21=", alpha21, ", α22=", alpha22),
         x="r1", y="r2", fill="Min LDGR") +
    theme_minimal() +
    coord_cartesian(xlim=c(-1, 3.8), ylim=c(-1, 3.8))
  
  # Plot 2: Feasibility domain with Method 1 (Arc Distance) contours
  p2 <- ggplot(results, aes(x=r1, y=r2)) +
    geom_tile(aes(fill=distance1)) +
    geom_contour(aes(z=distance1), breaks=0, color="black", linewidth=1) +
    geom_abline(slope=left_slope, intercept=0, linetype="dashed", color="white") +
    geom_abline(slope=right_slope, intercept=0, linetype="dashed", color="white") +
    scale_fill_gradient2(low="red", mid="white", high="blue", midpoint=0) +
    labs(title="Arc Distance (Method 1) in r1-r2 Space",
         subtitle=paste("α11=", alpha11, ", α12=", alpha12, ", α21=", alpha21, ", α22=", alpha22),
         x="r1", y="r2", fill="Arc Distance") +
    theme_minimal() +
    coord_cartesian(xlim=c(-1, 3.8), ylim=c(-1, 3.8))
  
  # Plot 3: Feasibility domain with Method 2 (Euclidean Distance) contours
  p3 <- ggplot(results, aes(x=r1, y=r2)) +
    geom_tile(aes(fill=distance2)) +
    geom_contour(aes(z=distance2), breaks=0, color="black", linewidth=1) +
    geom_abline(slope=left_slope, intercept=0, linetype="dashed", color="white") +
    geom_abline(slope=right_slope, intercept=0, linetype="dashed", color="white") +
    scale_fill_gradient2(low="red", mid="white", high="blue", midpoint=0) +
    labs(title="Euclidean Distance (Method 2) in r1-r2 Space",
         subtitle=paste("α11=", alpha11, ", α12=", alpha12, ", α21=", alpha21, ", α22=", alpha22),
         x="r1", y="r2", fill="Euclidean Distance") +
    theme_minimal() +
    coord_cartesian(xlim=c(-1, 3.8), ylim=c(-1, 3.8))
  
  # Plot 4: Difference between contours
  results$diff_dist <- results$distance1 - results$distance2
  p4 <- ggplot(results, aes(x=r1, y=r2)) +
    geom_tile(aes(fill=diff_dist)) +
    geom_contour(aes(z=min_LDGR), breaks=0, color="black", linetype="solid", linewidth=1) +
    geom_contour(aes(z=distance1), breaks=0, color="red", linetype="dashed", linewidth=1) +
    geom_contour(aes(z=distance2), breaks=0, color="blue", linetype="dotted", linewidth=1) +
    scale_fill_gradient2(low="purple", mid="white", high="orange", midpoint=0) +
    labs(title="Difference Between Distance Metrics",
         subtitle="Black=LDGR=0, Red=Arc Dist=0, Blue=Euclidean Dist=0",
         x="r1", y="r2", fill="Dist1 - Dist2") +
    theme_minimal() +
    coord_cartesian(xlim=c(-1, 3.8), ylim=c(-1, 3.8))
  
  return(list(p1, p2, p3, p4))
}

# Create feasibility domain plots for a few parameter sets
default_plots <- create_feasibility_plot(0.075, 0.075, 0.05, 0.05)
competitive_plots <- create_feasibility_plot(0.075, 0.075, 0.10, 0.10)
facilitative_plots <- create_feasibility_plot(0.075, 0.075, -0.05, -0.05)
asymmetric_plots <- create_feasibility_plot(0.075, 0.075, 0.10, -0.05)

#################################################
# Print summary results and insights
#################################################

# Print relationship statistics
print(relationship_stats)

# Calculate overall correlations
overall_cor_minLDGR_dist1 <- cor(all_results$min_LDGR, all_results$distance1)
overall_cor_minLDGR_dist2 <- cor(all_results$min_LDGR, all_results$distance2)

cat("\n\nOverall Correlation between min_LDGR and Arc Distance (Method 1):", overall_cor_minLDGR_dist1)
cat("\nOverall Correlation between min_LDGR and Euclidean Distance (Method 2):", overall_cor_minLDGR_dist2)

# Summary insights
cat("\n\nKey Insights from Parameter Sweep Analysis:\n")
cat("1. Effect of varying growth rates (r1, r2):\n")
cat("   - When varying r1 with fixed r2, Method 2 shows a", 
    relationship_stats$linearity_dist2[relationship_stats$scenario == "varying_r1_fixed_r2"], 
    "relationship with min_LDGR\n")
cat("   - When varying r2 with fixed r1, Method 2 shows a", 
    relationship_stats$linearity_dist2[relationship_stats$scenario == "varying_r2_fixed_r1"], 
    "relationship with min_LDGR\n")

cat("\n2. Effect of varying interspecific competition (alpha12, alpha21):\n")
cat("   - When varying alpha12, Method 2 shows a", 
    relationship_stats$linearity_dist2[relationship_stats$scenario == "varying_alpha12"], 
    "relationship with min_LDGR\n")
cat("   - When varying alpha21, Method 2 shows a", 
    relationship_stats$linearity_dist2[relationship_stats$scenario == "varying_alpha21"], 
    "relationship with min_LDGR\n")

cat("\n3. Effect of varying intraspecific competition (alpha11, alpha22):\n")

                                       
                                       