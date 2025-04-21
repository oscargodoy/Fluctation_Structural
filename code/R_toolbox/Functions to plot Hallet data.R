
# simplified version for a single cone----

library(circlize)

# Calculate the cone angles from interaction matrix
get_cone_angles <- function(gamma) {
  # Normalize vectors and get their angles
  normalize <- function(v) v/sqrt(sum(v^2))
  
  # Get angle in degrees for a point
  point_to_angle <- function(p) {
    angle <- ifelse(p[2] > 0,
                    acos(p[1]),
                    2*pi - acos(p[1]))
    return(360 * angle/(2*pi))
  }
  
  # Get the two vectors that define the cone
  m <- -gamma
  v1 <- normalize(m[,1])
  v2 <- normalize(m[,2])
  
  return(c(point_to_angle(v1), point_to_angle(v2)))
}

# Calculate cone size as fraction of circle
get_cone_size <- function(cone) {
  ((cone[2] - cone[1]) %% 360) / 360
}

# Draw only the main biodiversity cone
draw_biodiversity_cone <- function(gamma) {
  # Get cone boundaries
  cone <- get_cone_angles(gamma)
  
  # Setup plotting area
  plot(c(0, 1.1), c(0, 1.1), type = "n", 
       axes = FALSE, ann = FALSE, asp = 1)
  
  
  # Determine start and end angles of the main cone
  start_angle <- cone[1]
  end_angle <- cone[2]
  
  # Handle cases where the cone crosses the positive x-axis
  if (start_angle > end_angle) {
    end_angle <- end_angle + 360
  }
  
  # Draw the main biodiversity cone
  draw.sector(start_angle, end_angle, clock.wise = FALSE, 
              col = adjustcolor("green3", alpha = 0.3), border = "orange3", lwd=2)
  
  # Draw axes
  lines(c(-0.2, 1), c(0, 0), lwd = 1)
  lines(c(0, 0), c(-0.2, 1), lwd = 1)
  text(1.2, 0, "Avena", cex = 1.2, pos = 2)
  text(0, 1.1, "Erodium", cex = 1.2, pos = 1)
  axis(1,at=c(0, 0.2, 0.4, 0.6, 0.8 ,1), labels=c("0","0.2", "0.4", "0.6", "0.8","1.0"), pos=c(0,0) ,cex.axis=1.2)
  axis(2,at=c(0, 0.2, 0.4, 0.6, 0.8 ,1), labels=c("0","0.2", "0.4", "0.6", "0.8","1.0"), pos=c(0,0) ,cex.axis=1.2)
  
  # Draw unit circle outline for reference
  symbols(0, 0, circles = 1, inches = FALSE, add = TRUE, fg = "gray", lwd = 1)
}

# Example usage:
#gamma <- -1*matrix(c(1, 0.5, 0.2, 1), nrow = 2)
#draw_biodiversity_cone(gamma)


#mulitple matrices----
library(circlize)

# Calculate the cone angles from interaction matrix
get_cone_angles <- function(gamma) {
  # Normalize vectors and get their angles
  normalize <- function(v) v/sqrt(sum(v^2))
  
  # Get angle in degrees for a point
  point_to_angle <- function(p) {
    angle <- ifelse(p[2] > 0,
                    acos(p[1]),
                    2*pi - acos(p[1]))
    return(360 * angle/(2*pi))
  }
  
  # Get the two vectors that define the cone
  m <- -gamma
  v1 <- normalize(m[,1])
  v2 <- normalize(m[,2])
  
  return(c(point_to_angle(v1), point_to_angle(v2)))
}

# Calculate cone size as fraction of circle
get_cone_size <- function(cone) {
  size <- ((cone[2] - cone[1]) %% 360) / 360
  return(size)
}

# Draw multiple biodiversity cones from a list of gamma matrices
draw_multiple_cones <- function(gamma_list, colors = NULL, labels = NULL, show_sizes = FALSE) {
  # Setup plotting area
  plot(c(0, 1.1), c(0, 1.1), type = "n", 
       axes = FALSE, ann = FALSE, asp = 1, xlab = "", ylab = "")
  
  # Draw unit circle outline for reference
  symbols(0, 0, circles = 1, inches = FALSE, add = TRUE, fg = "gray", lwd = 1)
  
  # Draw axes
  lines(c(-0.2, 1), c(0, 0), lwd = 1)
  lines(c(0, 0), c(-0.2, 1), lwd = 1)
  text(1.2, 0, "Avena", cex = 1.2, pos = 2)
  text(0, 1.05, "Erodium", cex = 1.2, pos = 3, srt=0)
  axis(1,at=c(0, 0.2, 0.4, 0.6, 0.8 ,1), labels=c("0","0.2", "0.4", "0.6", "0.8","1.0"), pos=c(0,0) ,cex.axis=1.2)
  axis(2,at=c(0, 0.2, 0.4, 0.6, 0.8 ,1), labels=c("0","0.2", "0.4", "0.6", "0.8","1.0"), pos=c(0,0) ,cex.axis=1.2)
  
  
  # Set default colors if not provided
  if(is.null(colors)) {
    colors <- rep("green3", times=length(gamma_list))
  }
  
  # Set default labels if not provided
  if(is.null(labels)) {
    labels <- paste("Gamma", 1:length(gamma_list))
  }
  
  # Initialize legend data
  legend_colors <- c()
  legend_labels <- c()
  
  # Draw each cone
  for(i in 1:length(gamma_list)) {
    gamma <- gamma_list[[i]]
    cone <- get_cone_angles(gamma)
    
    # Handle cases where the cone crosses the positive x-axis
    start_angle <- cone[1]
    end_angle <- cone[2]
    if(start_angle > end_angle) {
      end_angle <- end_angle + 360
    }
    
    # Draw this cone
    draw.sector(start_angle, end_angle, clock.wise = FALSE, 
                col = adjustcolor(colors[i], alpha = 0.3), border = randomColor(), lwd=2)
  
    # Add cone size label if requested
    if(show_sizes) {
      cone_size <- get_cone_size(cone)
      mid_angle <- 2*pi*(start_angle + (end_angle - start_angle)/2)/360
      label_radius <- 1.15 + (0.1 * ((i-1) %% 3))  # Stagger labels at different radii
      text(label_radius*cos(mid_angle), label_radius*sin(mid_angle), 
           round(cone_size, digits = 2), cex = 0.9, col = colors[i])
    }
    
    # Add to legend data
    #legend_colors <- c(legend_colors, colors[i])
    #legend_labels <- c(legend_labels, labels[i])
    
  }
  
  # Add legend
  #legend("topright", legend = legend_labels, fill = legend_colors, 
         #cex = 0.8, bty = "n", border = NA)
}
    

# Example usage:
# Create a list of gamma matrices to compare
#gamma1 <- matrix(c(1, 0.5, 0.5, 1)*-1, nrow = 2)
#gamma2 <- matrix(c(1, 0.8, 0.4, 1)*-1, nrow = 2)
#gamma3 <- matrix(c(1, 0.3, 0.7, 1)*-1, nrow = 2)
# 
# # Draw multiple cones
#draw_multiple_cones(
#  list(gamma1, gamma2, gamma3),
# labels = c("Low competition", "Medium competition", "High competition"))


