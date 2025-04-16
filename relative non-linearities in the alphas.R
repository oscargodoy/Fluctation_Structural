theta <- seq(0, 2*pi, length.out = 1000)

# Calculate x and y coordinates on the unit circle
x <- cos(theta)
y <- sin(theta)

# Create a plot
plot(x, y, type = "l", col = "blue", lwd = 2,
     asp = 1, # Set aspect ratio to 1 to make circle appear circular
     xlab = "x", ylab = "y",
     main = "Unit Circle",
     xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5))

# Add axes
abline(h = 0, v = 0, col = "gray", lty = 2)

# Add a point at (1,0) to mark the starting point
points(1, 0, col = "red", pch = 16)

# Add text labels for key points
text(1.1, 0, "(1,0)", pos = 4, col = "red")
text(0, 1.1, "(0,1)", pos = 3, col = "darkgreen")
text(-1.1, 0, "(-1,0)", pos = 2, col = "darkgreen")
text(0, -1.1, "(0,-1)", pos = 1, col = "darkgreen")

# Draw radius to illustrate r=1
segments(0, 0, 1, 0, col = "red", lty = 2)



# R code to draw cones from origin to unit circle

# Set up the plot
par(mar = c(3, 3, 2, 2))  # Set margins
plot(NULL, xlim = c(-1, 1), ylim = c(-1, 1), asp = 1,
     xlab = "x", ylab = "y", main = "Cones from Origin to Unit Circle")

# Draw the unit circle
theta <- seq(0, 2*pi, length.out = 1000)
x_circle <- cos(theta)
y_circle <- sin(theta)
lines(x_circle, y_circle, col = "black", lwd = 2)

# Add coordinate axes
abline(h = 0, v = 0, col = "gray", lty = 2)

# Number of cones to draw
num_cones <- 4

# Draw cones with different angles
cone_angles <- seq(0, 2*pi, length.out = num_cones + 1)[-(num_cones + 1)]
cone_colors <- rainbow(num_cones)

cone_angles <- c(0, -0.3, 0.6, -0.9, 1.2)  # Example angles for cones
cone_colors <- c("green","green","green","green","green")

# Draw each cone
for (i in 1:num_cones) {
  # Cone starting angle and cone width (in radians)
  start_angle <- cone_angles[i]
  cone_width <- jitter(pi/4, 0.2)  # width of each cone
  
  # Calculate arc on the unit circle for this cone
  arc_angles <- seq(start_angle - cone_width/2, start_angle + cone_width/2, length.out = 100)
  x_arc <- cos(arc_angles)
  y_arc <- sin(arc_angles)
  
  # Draw filled cone
  polygon(c(0, x_arc, 0), c(0, y_arc, 0), col = adjustcolor(cone_colors[i], alpha = 0.3), border = NA)
  
  # Draw cone borders
  segments(0, 0, cos(start_angle - cone_width/2), sin(start_angle - cone_width/2), col = cone_colors[i], lwd = 2)
  segments(0, 0, cos(start_angle + cone_width/2), sin(start_angle + cone_width/2), col = cone_colors[i], lwd = 2)
  
  # Draw arc on unit circle
  lines(x_arc, y_arc, col = cone_colors[i], lwd = 2)
  
  # Add label in the middle of each cone
  label_angle <- start_angle
  label_r <- 0.7  # position of label (distance from origin)
  text(label_r * cos(label_angle), label_r * sin(label_angle), 
       paste("Cone", i), col = cone_colors[i], cex = 0.8)
}

# Add origin label
text(0.1, 0.1, "(0,0)", cex = 0.8)

# Add title at the bottom for explanation
mtext("Each cone projects from the origin (0,0) to an arc on the unit circle", side = 1, line = 2.5)



# R code to draw relative non-linearities in the alpha

# Set up the plot
par(mar = c(4, 4, 3, 3))  # Set margins
plot(NULL, xlim = c(-0.99, 0.99), ylim = c(-0.99, 0.99), asp = 1,
     xlab = " ", ylab = " ", main = "Relative non-linearity in the alphas")

# Draw the unit circle
theta <- seq(0, 2*pi, length.out = 1000)
x_circle <- cos(theta)
y_circle <- sin(theta)
lines(x_circle, y_circle, col = "black", lwd = 2)

# Add coordinate axes
abline(h = 0, v = 0, col = "gray", lty = 2)

# Number of cones to draw
num_cones <- 4

# Starting positions for cones (in radians)
start_positions <- c(-0.1, -0.3, 1, 0.8)


# Different widths for each cone (in radians)
cone_widths <- c(pi/6, pi/12, pi/3, pi/4)

# Colors for the cones
cone_colors <- rep("lightgreen", times = num_cones)

env <- c("dry", "very dry","very wet", "wet")

# Draw each cone
for (i in 1:num_cones) {
  # Cone center angle and width
  center_angle <- start_positions[i]
  width <- cone_widths[i]
  
  # Calculate arc on the unit circle for this cone
  arc_angles <- seq(center_angle - width/2, center_angle + width/2, length.out = 100)
  x_arc <- cos(arc_angles)
  y_arc <- sin(arc_angles)
  
  # Draw filled cone
  polygon(c(0, x_arc, 0), c(0, y_arc, 0), col = adjustcolor(cone_colors[i], alpha = 0.3), border = NA)
  
  # Draw cone borders
  segments(0, 0, cos(center_angle - width/2), sin(center_angle - width/2), col = cone_colors[i], lwd = 2)
  segments(0, 0, cos(center_angle + width/2), sin(center_angle + width/2), col = cone_colors[i], lwd = 2)
  
  # Draw arc on unit circle
  lines(x_arc, y_arc, col = cone_colors[i], lwd = 2)
  
  # Add label in the middle of each cone
  label_angle <- center_angle
  label_r <- 0.65  # position of label (distance from origin)
  text(label_r * cos(label_angle), label_r * sin(label_angle), 
       paste0("Env. ", env[i]), col = "black", cex = 0.8)
}

lines(x_circle, y_circle, col = "black", lwd = 2)

#draw a single vector that starts from origin zero, and its normalized
r <- c(1,0.7)
point1 <-r/sqrt(sum(r^2))
arrows(0, 0, x1 = point1[1], y1 = point1[2], length = 0.15, angle = 30,
       code = 2, col = par("fg"), lty = 3,
       lwd = par("lwd"))
text(x=point1[1]+0.1, y=point1[2]+0.04, "Average r", cex=0.8)


# Add title at the bottom for explanation
mtext("Changes in environmental conditons lead to changes in alphas but not in lambdas", side = 1, line = 2.5)




