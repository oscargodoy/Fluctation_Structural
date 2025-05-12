# Structural Stability Feasibility Domain Visualization (Revised)
# Load required libraries
library(ggplot2)
library(shape)

# Define alpha values (competition coefficients)
# Using values that create a wider feasibility domain
alpha11 <- 0.7  # intraspecific competition of species 1
alpha12 <- 0.3  # effect of species 2 on species 1  
alpha21 <- 0.4  # effect of species 1 on species 2
alpha22 <- 0.8  # intraspecific competition of species 2

# Calculate the slopes of the boundaries
slope_right <- alpha21/alpha11  # Right boundary slope
slope_left <- alpha22/alpha12   # Left boundary slope

# Create the plot
par(mfrow=c(1,1), mar=c(5,5,4,2))

# Set up the plot area (focusing on positive quadrant)
plot(0, 0, xlim=c(0, 1), ylim=c(0, 1), type="n",
     xlab="r1", ylab="r2",
     cex.lab=1.2, cex.main=1.4, asp=1) +
scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1))

# Add the unit circle
circle <- seq(0, 2*pi, length.out=100)
x_circle <- cos(circle)
y_circle <- sin(circle)
lines(x_circle, y_circle, lwd=2)

# Add axes
abline(h=0, v=0, col="gray", lty=2)

# Draw the feasibility domain boundaries (only in positive quadrant)
# Right boundary (slope = alpha21/alpha11)
x_vals <- seq(0, 1, length.out=100)
y_vals_right <- slope_right * x_vals

# Left boundary (slope = alpha22/alpha12)  
y_vals_left <- slope_left * x_vals

# Clip to the circle in the positive quadrant
# For right boundary
valid_right <- x_vals^2 + y_vals_right^2 <= 1 & x_vals >= 0 & y_vals_right >= 0
x_right <- x_vals[valid_right]
y_right <- y_vals_right[valid_right]

# For left boundary
valid_left <- x_vals^2 + y_vals_left^2 <= 1 & x_vals >= 0 & y_vals_left >= 0
x_left <- x_vals[valid_left]
y_left <- y_vals_left[valid_left]

# # Draw the boundaries with dashed lines
# lines(x_right, y_right, col="black", lty=2, lwd=1.5)
# lines(x_left, y_left, col="black", lty=2, lwd=1.5)

# Create the polygon for the feasibility domain
# Find intersection points with circle
# Right boundary: y = slope_right * x, x^2 + y^2 = 1
x_intersect_right <- 1/sqrt(1 + slope_right^2)
y_intersect_right <- slope_right * x_intersect_right

# Left boundary: y = slope_left * x, x^2 + y^2 = 1
x_intersect_left <- 1/sqrt(1 + slope_left^2)
y_intersect_left <- slope_left * x_intersect_left

# Create arc between the two intersection points
arc_angle_right <- atan2(y_intersect_right, x_intersect_right)
arc_angle_left <- atan2(y_intersect_left, x_intersect_left)
arc_angles <- seq(arc_angle_right, arc_angle_left, length.out=50)
arc_x <- cos(arc_angles)
arc_y <- sin(arc_angles)

# Combine to create polygon
polygon_x <- c(0, x_intersect_right, arc_x, x_intersect_left)
polygon_y <- c(0, y_intersect_right, arc_y, y_intersect_left)

# Fill the feasibility domain
polygon(polygon_x, polygon_y, col=rgb(0.5, 1, 0.5, alpha=0.5), border=NA)

# Calculate centroid of the feasibility domain
centroid_x <- mean(polygon_x)
centroid_y <- mean(polygon_y)

# Define the vectors
# 1. Long vector at the centroid
vector1_x <- centroid_x
vector1_y <- centroid_y
length_centroid <- sqrt(vector1_x^2 + vector1_y^2)

# 3. Vector near the edge (same length as centroid vector)
# Choose a point on the right boundary
t <- 0.75  # Parameter along the boundary  
edge_point_x <- t * x_intersect_right
edge_point_y <- t * y_intersect_right

# Calculate the perpendicular distance from this point to the boundary
dist_to_right_edge <- abs(-alpha21 * edge_point_x + alpha11 * edge_point_y) / sqrt(alpha21^2 + alpha11^2)

# Normal vector to right boundary (perpendicular, pointing inward)
normal_x <- -slope_right / sqrt(1 + slope_right^2)
normal_y <- 1 / sqrt(1 + slope_right^2)

# Move the edge point inward by a small distance
inward_distance <- 0.05
temp_x <- edge_point_x + inward_distance * normal_x
temp_y <- edge_point_y + inward_distance * normal_y

# Now scale this vector to have the same length as the centroid vector
current_length <- sqrt(temp_x^2 + temp_y^2)
scale_to_match <- length_centroid / current_length
vector3_x <- temp_x * scale_to_match
vector3_y <- temp_y * scale_to_match

# Calculate the actual distance from vector3 to the edge
dist_edge_vector <- abs(-alpha21 * vector3_x + alpha11 * vector3_y) / sqrt(alpha21^2 + alpha11^2)

# 2. Short vector with the same distance to edge as the edge vector
# We need to find a point that has the same distance to the edge as vector3
# Start with a point on the same ray as the centroid but closer to origin
trial_scale <- 0.3
trial_x <- centroid_x * trial_scale
trial_y <- centroid_y * trial_scale

# Calculate its distance to the edges
dist_trial_right <- abs(-alpha21 * trial_x + alpha11 * trial_y) / sqrt(alpha21^2 + alpha11^2)
dist_trial_left <- abs(-alpha22 * trial_x + alpha12 * trial_y) / sqrt(alpha22^2 + alpha12^2)
dist_trial <- min(dist_trial_right, dist_trial_left)

# Adjust the scale to match the distance
# Use iterative approach to find the right scale
target_distance <- dist_edge_vector
best_scale <- trial_scale
for (i in 1:20) {
  test_x <- centroid_x * best_scale
  test_y <- centroid_y * best_scale
  dist_right <- abs(-alpha21 * test_x + alpha11 * test_y) / sqrt(alpha21^2 + alpha11^2)
  dist_left <- abs(-alpha22 * test_x + alpha12 * test_y) / sqrt(alpha22^2 + alpha12^2)
  current_dist <- min(dist_right, dist_left)
  
  if (abs(current_dist - target_distance) < 0.001) break
  
  if (current_dist > target_distance) {
    best_scale <- best_scale * 0.95
  } else {
    best_scale <- best_scale * 1.05
  }
}

vector2_x <- centroid_x * best_scale
vector2_y <- centroid_y * best_scale

# Recalculate to ensure accuracy
dist_B_right <- abs(-alpha21 * vector2_x + alpha11 * vector2_y) / sqrt(alpha21^2 + alpha11^2)
dist_B_left <- abs(-alpha22 * vector2_x + alpha12 * vector2_y) / sqrt(alpha22^2 + alpha12^2)
dist_B <- min(dist_B_right, dist_B_left)

# Draw the vectors
arrows(0, 0, vector1_x, vector1_y, col="black", lwd=2, length=0.1)
arrows(0, 0, vector2_x, vector2_y, col="black", lwd=2, length=0.1)
arrows(0, 0, vector3_x, vector3_y, col="black", lwd=2, length=0.1)

# # Add points at vector tips
# points(c(vector1_x, vector2_x, vector3_x), 
#        c(vector1_y, vector2_y, vector3_y), 
#        pch=21, cex=2, col="black", bg="red")

# Add labels
text(vector1_x - 0.05, vector1_y, "A", col="black", cex=1.2, font=2)
text(vector2_x - 0.05, vector2_y, "B", col="black", cex=1.2, font=2)
text(vector3_x - 0.05, vector3_y , "C", col="black", cex=1.2, font=2)

# # Add a fourth point D on the circle boundary
# vector4_x <- 0.9
# vector4_y <- 0.2
# # points(vector4_x, vector4_y, pch=21, cex=2, col="black", bg="red")
# text(vector4_x + 0.05, vector4_y, "D", col="black", cex=1.2, font=2)

# Add dashed lines from origin to all points
segments(0, 0, vector1_x, vector1_y, lty=2, col="black", lwd=1)
segments(0, 0, vector2_x, vector2_y, lty=2, col="black", lwd=1)
segments(0, 0, vector3_x, vector3_y, lty=2, col="black", lwd=1)
#segments(0, 0, vector4_x, vector4_y, lty=2, col="black", lwd=1)

# Print the values used
cat("Alpha values used:\n")
cat(sprintf("α₁₁ = %.2f, α₁₂ = %.2f\n", alpha11, alpha12))
cat(sprintf("α₂₁ = %.2f, α₂₂ = %.2f\n", alpha21, alpha22))
cat("\nBoundary slopes:\n")
cat(sprintf("Right boundary (α₂₁/α₁₁) = %.3f\n", slope_right))
cat(sprintf("Left boundary (α₂₂/α₁₂) = %.3f\n", slope_left))
cat("\nGrowth rate vectors:\n")
cat(sprintf("A - Long vector (centroid): r₁ = %.3f, r₂ = %.3f\n", vector1_x, vector1_y))
cat(sprintf("B - Short vector: r₁ = %.3f, r₂ = %.3f\n", vector2_x, vector2_y))
cat(sprintf("C - Edge vector: r₁ = %.3f, r₂ = %.3f\n", vector3_x, vector3_y))
#cat(sprintf("D - Outside vector: r₁ = %.3f, r₂ = %.3f\n", vector4_x, vector4_y))

# Calculate vector lengths
length_A <- sqrt(vector1_x^2 + vector1_y^2)
length_B <- sqrt(vector2_x^2 + vector2_y^2)
length_C <- sqrt(vector3_x^2 + vector3_y^2)

cat("\nVector lengths:\n")
cat(sprintf("A - Long vector: %.3f\n", length_A))
cat(sprintf("B - Short vector: %.3f\n", length_B))
cat(sprintf("C - Edge vector: %.3f\n", length_C))

# Calculate and report distances to edge
calc_distance_to_edge <- function(x, y, alpha11, alpha12, alpha21, alpha22) {
  # Distance to right boundary
  dist_right <- abs(-alpha21 * x + alpha11 * y) / sqrt(alpha21^2 + alpha11^2)
  # Distance to left boundary
  dist_left <- abs(-alpha22 * x + alpha12 * y) / sqrt(alpha22^2 + alpha12^2)
  
  # Check if inside feasibility domain
  if (y >= alpha21/alpha11 * x && y <= alpha22/alpha12 * x && x >= 0 && y >= 0) {
    return(min(dist_right, dist_left))
  } else {
    return(-min(dist_right, dist_left))
  }
}

cat("\nDistances to feasibility domain edge:\n")
dist_A <- calc_distance_to_edge(vector1_x, vector1_y, alpha11, alpha12, alpha21, alpha22)
dist_B <- calc_distance_to_edge(vector2_x, vector2_y, alpha11, alpha12, alpha21, alpha22)
dist_C <- calc_distance_to_edge(vector3_x, vector3_y, alpha11, alpha12, alpha21, alpha22)
#dist_D <- calc_distance_to_edge(vector4_x, vector4_y, alpha11, alpha12, alpha21, alpha22)

cat(sprintf("A - Long vector: %.3f\n", dist_A))
cat(sprintf("B - Short vector: %.3f\n", dist_B))
cat(sprintf("C - Edge vector: %.3f\n", dist_C))
#cat(sprintf("D - Outside vector: %.3f\n", dist_D))

# Verify conditions
cat("\nVerification:\n")
cat(sprintf("Length A = %.3f, Length C = %.3f (should be equal)\n", length_A, length_C))
cat(sprintf("Distance B = %.3f, Distance C = %.3f (should be equal)\n", dist_B, dist_C))

