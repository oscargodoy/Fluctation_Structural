# Load required libraries
library(ggplot2)
library(dplyr)

# Define alpha values (competition coefficients)
alpha11 <- 0.7  # intraspecific competition of species 1
alpha12 <- 0.3  # effect of species 2 on species 1  
alpha21 <- 0.4  # effect of species 1 on species 2
alpha22 <- 0.8  # intraspecific competition of species 2

# Calculate slopes for feasibility domain boundaries
slope_left <- alpha22/alpha12  # 0.8/0.3 = 2.667
slope_right <- alpha21/alpha11  # 0.4/0.7 = 0.571

# Define colors for edges
color_left <- "blue"
color_right <- "orange"

# Create circular boundary
radius <- 1.0
n_points <- 100
angles <- seq(0, 2*pi, length.out = n_points)
circle_df <- data.frame(
  x = radius * cos(angles),
  y = radius * sin(angles)
)

# Create feasibility domain extending to the circle
theta_right <- atan(slope_right)
theta_left <- atan(slope_left)

# Create arc for the curved part of the feasibility domain
arc_angles <- seq(theta_right, theta_left, length.out = 50)
arc_x <- radius * cos(arc_angles)
arc_y <- radius * sin(arc_angles)

# Create feasibility domain polygon including the arc
fd_vertices <- data.frame(
  r1 = c(0, arc_x, 0),
  r2 = c(0, arc_y, 0)
)

# Separate the edges for different coloring
left_edge <- data.frame(
  x = c(0, radius * cos(theta_left)),
  y = c(0, radius * sin(theta_left))
)

right_edge <- data.frame(
  x = c(0, radius * cos(theta_right)),
  y = c(0, radius * sin(theta_right))
)

# Calculate centroid of feasibility domain  
centroid_r1 <- mean(c(0, radius * cos(theta_left), radius * cos(theta_right)))
centroid_r2 <- mean(c(0, radius * sin(theta_left), radius * sin(theta_right)))

# Define vectors
# Vector 1: Inside feasibility domain (closer to centroid)
v1_length <- 0.7  # Increased length to prevent overlap
v1_angle <- (theta_left + theta_right) / 2  # Middle of feasibility domain
v1_start <- c(0, 0)
v1_end <- c(v1_length * cos(v1_angle), v1_length * sin(v1_angle))

# Label for inside vector
v1_label <- "A"

# Vector 2: Outside feasibility domain (to the right of right boundary)
v2_length <- 0.7
v2_angle <- 0.2  # Small angle, outside the feasibility domain
v2_start <- c(0, 0)
v2_end <- c(v2_length * cos(v2_angle), v2_length * sin(v2_angle))

# Label for outside vector
v2_label <- "B"

# Function to find perpendicular point on a line from origin
find_perpendicular_point <- function(point, slope) {
  # For a line y = slope * x through origin
  # Perpendicular line has slope -1/slope
  # Find intersection point
  x_perp <- (point[1] + slope * point[2]) / (1 + slope^2)
  y_perp <- x_perp * slope
  return(c(x_perp, y_perp))
}

# Find perpendicular points on both boundaries for each vector
# For vector 1 (inside)
v1_perp_left <- find_perpendicular_point(v1_end, slope_left)
v1_perp_right <- find_perpendicular_point(v1_end, slope_right)

# For vector 2 (outside)
v2_perp_left <- find_perpendicular_point(v2_end, slope_left)
v2_perp_right <- find_perpendicular_point(v2_end, slope_right)

# Constrain perpendicular points to be within the circle
constrain_to_circle <- function(point) {
  dist <- sqrt(point[1]^2 + point[2]^2)
  if (dist > radius) {
    return(point * radius / dist)
  }
  return(point)
}

v1_perp_left <- constrain_to_circle(v1_perp_left)
v1_perp_right <- constrain_to_circle(v1_perp_right)
v2_perp_left <- constrain_to_circle(v2_perp_left)
v2_perp_right <- constrain_to_circle(v2_perp_right)

# Create data frames for plotting
vectors_df <- data.frame(
  x_start = c(v1_start[1], v2_start[1]),
  y_start = c(v1_start[2], v2_start[2]),
  x_end = c(v1_end[1], v2_end[1]),
  y_end = c(v1_end[2], v2_end[2]),
  label = c(v1_label, v2_label)
)

# Create distance lines for both edges with colors
distance_lines_df <- data.frame(
  x_start = c(v1_end[1], v1_end[1], v2_end[1], v2_end[1]),
  y_start = c(v1_end[2], v1_end[2], v2_end[2], v2_end[2]),
  x_end = c(v1_perp_left[1], v1_perp_right[1], v2_perp_left[1], v2_perp_right[1]),
  y_end = c(v1_perp_left[2], v1_perp_right[2], v2_perp_left[2], v2_perp_right[2]),
  vector_type = c("v1", "v1", "v2", "v2"),
  edge_type = c("left", "right", "left", "right"),
  color = c(color_left, color_right, color_left, color_right),
  pn = c("apositive", "apositive", "apositive", "negative")
)



# Create the plot
p <- ggplot() +
  # Background circle
  geom_path(data = circle_df, aes(x = x, y = y), 
            color = "black", size = 1.2) +
  
  # Feasibility domain (filled)
  geom_polygon(data = fd_vertices, aes(x = r1, y = r2), 
               fill = "lightgreen", alpha = 0.6) +
  
  # Colored edges of feasibility domain
  geom_segment(data = left_edge, 
               aes(x = x[1], y = y[1], xend = x[2], yend = y[2]),
               color = color_left, size = 1.5) +
  geom_segment(data = right_edge, 
               aes(x = x[1], y = y[1], xend = x[2], yend = y[2]),
               color = color_right, size = 1.5) +
  
  # Arc portion of the feasibility domain edge
  geom_path(data = data.frame(x = arc_x, y = arc_y), 
            aes(x = x, y = y), 
            color = "darkgreen", size = 1.5) +
  
  # Vectors
  geom_segment(data = vectors_df, 
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(0.3, "cm")), 
               size = 1.2, color = "black") +
  
  # Distance lines (dashed) - colored to match edges
  geom_segment(data = distance_lines_df,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = color, linetype = pn),
              size = 1) +
  scale_color_identity() +

  # Labels for vectors
  geom_point(data = vectors_df, 
             aes(x = x_end, y = y_end),
             color = "black", size = 3) +
  geom_text(data = vectors_df, 
            aes(x = x_end, y = y_end, label = label),
            vjust = -1.5, hjust = 0.5, size = 5, fontface = "bold") +
  
  # Axes
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed", alpha = 0.5) +
  
  # Labels
  labs(x = expression(r[1]),
       y = expression(r[2])) +
  
  # Theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(color = "black"),
    aspect.ratio = 1,
    legend.position = "none"
  ) +
  
  # Set axis limits
  coord_fixed(xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1))

# Calculate distances to both edges
dist_v1_left <- sqrt((v1_end[1] - v1_perp_left[1])^2 + (v1_end[2] - v1_perp_left[2])^2)
dist_v1_right <- sqrt((v1_end[1] - v1_perp_right[1])^2 + (v1_end[2] - v1_perp_right[2])^2)
dist_v2_left <- sqrt((v2_end[1] - v2_perp_left[1])^2 + (v2_end[2] - v2_perp_left[2])^2)
dist_v2_right <- sqrt((v2_end[1] - v2_perp_right[1])^2 + (v2_end[2] - v2_perp_right[2])^2)

# Determine which distances to show (closest for inside, both for outside)
# For vector 1 (inside), show the minimum distance
min_dist_v1 <- min(dist_v1_left, dist_v1_right)

# # Add distance annotations
# p <- p +
#   # Distance for vector A (inside) - just show the minimum
#   annotate("text", 
#            x = v1_end[1] + 0.1, 
#            y = v1_end[2] + 0.05,
#            label = paste0("d = ", round(min_dist_v1, 2)), 
#            size = 3.5, color = "darkgreen") +
#   # Distance for vector D (outside) - show as negative
#   annotate("text", 
#            x = v2_end[1] - 0.05, 
#            y = v2_end[2] - 0.1,
#            label = paste0("d = -", round(dist_v2_right, 2)), 
#            size = 3.5, color = "darkred")

# Add small annotations for the boundary distances with matching colors
p <- p +
  # Annotate distances on the perpendicular lines for vector A
  annotate("text", 
           x = mean(c(v1_end[1], v1_perp_left[1])) - 0.03, 
           y = mean(c(v1_end[2], v1_perp_left[2])) + 0.05,
           label = round(dist_v1_left, 2), 
           size = 4, color = color_left) +
  annotate("text", 
           x = mean(c(v1_end[1], v1_perp_right[1])) + 0.05, 
           y = mean(c(v1_end[2], v1_perp_right[2])) - 0.01,
           label = round(dist_v1_right, 2), 
           size = 4, color = color_right) +
  # Annotate distances on the perpendicular lines for vector D
  annotate("text", 
           x = mean(c(v2_end[1], v2_perp_left[1])) - 0.06, 
           y = mean(c(v2_end[2], v2_perp_left[2])) + 0.06,
           label = round(dist_v2_left, 2), 
           size = 4, color = color_left) +
  annotate("text", 
           x = mean(c(v2_end[1], v2_perp_right[1])) + 0.03, 
           y = mean(c(v2_end[2], v2_perp_right[2])) +.05,
           label = -round(dist_v2_right, 2), 
           size = 4, color = color_right)

# Display the plot
print(p)

# Save the plot
ggsave("feasibility_domain_colored_edges.png", p, width = 8, height = 8, dpi = 300)