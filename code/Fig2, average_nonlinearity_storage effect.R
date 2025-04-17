library(jpeg)

par(mfrow = c(2, 2))

# Average conditions ----

# Set up the plot
#jpeg("figures/Average conditions.jpg")


par(mar = c(4, 4, 3, 3))  # Set margins
plot(NULL, xlim = c(-0.99, 0.99), ylim = c(-0.99, 0.99), asp = 1,
     xlab = " ", ylab = " ", main = " A. Average conditions")

# Draw the unit circle
theta <- seq(0, 2*pi, length.out = 1000)
x_circle <- cos(theta)
y_circle <- sin(theta)
lines(x_circle, y_circle, col = "black", lwd = 2)

# Add coordinate axes
abline(h = 0, v = 0, col = "gray", lty = 2)

# Number of cones to draw
num_cones <- 1

# Starting positions for cones (in radians)
start_positions <- 0.4


# Different widths for each cone (in radians)
cone_widths <- c(pi/3)

# Colors for the cones
cone_colors <- "lightgreen"

env <- c("Average environ.")

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
r <- c(1,1.05)
point1 <-r/sqrt(sum(r^2))
arrows(0, 0, x1 = point1[1], y1 = point1[2], length = 0.15, angle = 30,
       code = 2, col = "brown", lty =1, lwd=2)
text(x=point1[1]+0.12, y=point1[2]+0.05, "Average r", cex=0.8)


# Add title at the bottom for explanation
mtext("Average environmental conditons implies neither variation in the lambdas nor in the alphas", side = 1, line = 2.5)
#dev.off()

# Relative non-linearity in the lambdas ----

jpeg("figures/Relative non linearity lambdas.jpg")

# Set up the plot
par(mar = c(4, 4, 3, 3))  # Set margins
plot(NULL, xlim = c(-0.99, 0.99), ylim = c(-0.99, 0.99), asp = 1,
     xlab = " ", ylab = " ", main = " B. Relative non-linearity in the lambdas")

# Draw the unit circle
theta <- seq(0, 2*pi, length.out = 1000)
x_circle <- cos(theta)
y_circle <- sin(theta)
lines(x_circle, y_circle, col = "black", lwd = 2)

# Add coordinate axes
abline(h = 0, v = 0, col = "gray", lty = 2)

# Number of cones to draw
num_cones <- 1

# Starting positions for cones (in radians)
start_positions <- 0.4


# Different widths for each cone (in radians)
cone_widths <- c(pi/3)

# Colors for the cones
cone_colors <- "lightgreen"

env <- c("Average environ.")

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
r <- c(1,1)
point1 <-r/sqrt(sum(r^2))
arrows(0, 0, x1 = point1[1], y1 = point1[2], length = 0.15, angle = 30,
       code = 2, col = "brown", lty =1, lwd=2)
text(x=point1[1]+0.05, y=point1[2]+0.04, "dry", cex=0.8)

r <- c(1,1.8)
point1 <-r/sqrt(sum(r^2))
arrows(0, 0, x1 = point1[1], y1 = point1[2], length = 0.15, angle = 30,
       code = 2, col = "brown", lty =1, lwd=2)
text(x=point1[1]+0.05, y=point1[2]+0.04, "very dry", cex=0.8)

r <- c(2.5,0.6)
point1 <-r/sqrt(sum(r^2))
arrows(0, 0, x1 = point1[1], y1 = point1[2], length = 0.15, angle = 30,
       code = 2, col = "brown", lty =1, lwd=2)
text(x=point1[1]+0.05, y=point1[2]+0.04, "wet", cex=0.8)

r <- c(2.5,-0.7)
point1 <-r/sqrt(sum(r^2))
arrows(0, 0, x1 = point1[1], y1 = point1[2], length = 0.15, angle = 30,
       code = 2, col = "brown", lty =1, lwd=2)
text(x=point1[1]-0.17, y=point1[2]-0.02, "very wet", cex=0.8)

# Add title at the bottom for explanation
mtext("Changes in environmental conditons lead to changes in lambdas but not in alphas", side = 1, line = 2.5)
dev.off()

# Relative non-linearity in the alphas ----

# Set up the plot
jpeg("figures/Relative non linearity alphas.jpg")
par(mar = c(4, 4, 3, 3))  # Set margins
plot(NULL, xlim = c(-0.99, 0.99), ylim = c(-0.99, 0.99), asp = 1,
     xlab = " ", ylab = " ", main = " C. Relative non-linearity in the alphas")

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

env <- c("wet", "very wet","very dry", "dry")

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
       code = 2, col = "brown", lty =1, lwd=2)
text(x=point1[1]+0.12, y=point1[2]+0.04, "Average λ", cex=0.8)


# Add title at the bottom for explanation
mtext("Changes in environmental conditons lead to changes in alphas but not in lambdas", side = 1, line = 2.5)

dev.off()
# Storage effect ----

# Set up the plot
jpeg("figures/The storage effect.jpg")
par(mar = c(4, 4, 3, 3))  # Set margins
plot(NULL, xlim = c(-0.99, 0.99), ylim = c(-0.99, 0.99), asp = 1,
     xlab = " ", ylab = " ", main = " D. The storage effect")

# Draw the unit circle
theta <- seq(0, 2*pi, length.out = 1000)
x_circle <- cos(theta)
y_circle <- sin(theta)
lines(x_circle, y_circle, col = "black", lwd = 2)

# Add coordinate axes
abline(h = 0, v = 0, col = "gray", lty = 2)

# Number of cones to draw
num_cones <- 1

# Starting positions for cones (in radians)
start_positions <- 0.5

# Different widths for each cone (in radians)
cone_widths <- c(pi/2)

# Colors for the cones
cone_colors <- "lightgreen"


# Draw each cone
for (i in 1:num_cones) {
  # Cone center angle and width
  center_angle <- start_positions[i]
  width <- cone_widths[i]
  
  # Calculate arc on the unit circle for this cone
  arc_angles <- seq(center_angle - width/2, center_angle + width/2, length.out = 100)
  x_arc <- cos(arc_angles)
  y_arc <- sin(arc_angles)
  
  # Create curved edges by using Bezier curves
  # For left edge
  left_angle <- center_angle - width/2
  left_control_x <- 0.5 * cos(left_angle + 0.2)  # Control point for curve
  left_control_y <- 0.5 * sin(left_angle + 0.2)
  left_edge_x <- seq(0, cos(left_angle), length.out = 50)
  # Create curved path using quadratic Bezier function
  t <- seq(0, 1, length.out = 50)
  left_edge_y <- (1-t)^2 * 0 + 2*(1-t)*t * left_control_y + t^2 * sin(left_angle)
  left_edge_x <- (1-t)^2 * 0 + 2*(1-t)*t * left_control_x + t^2 * cos(left_angle)
  
  # For right edge
  right_angle <- center_angle + width/2
  right_control_x <- 0.5 * cos(right_angle - 0.2)  # Control point for curve
  right_control_y <- 0.5 * sin(right_angle - 0.2)
  right_edge_x <- seq(0, cos(right_angle), length.out = 50)
  # Create curved path using quadratic Bezier function
  right_edge_y <- (1-t)^2 * 0 + 2*(1-t)*t * right_control_y + t^2 * sin(right_angle)
  right_edge_x <- (1-t)^2 * 0 + 2*(1-t)*t * right_control_x + t^2 * cos(right_angle)
  
  # Draw filled cone with curved edges
  polygon(c(0, left_edge_x, x_arc, right_edge_x[50:1]), 
          c(0, left_edge_y, y_arc, right_edge_y[50:1]), 
          col = adjustcolor(cone_colors[i], alpha = 0.3), border = NA)
  
  # Draw curved cone borders
  lines(left_edge_x, left_edge_y, col = cone_colors[i], lwd = 2)
  lines(right_edge_x, right_edge_y, col = cone_colors[i], lwd = 2)
  
  # Draw arc on unit circle
  lines(x_arc, y_arc, col = cone_colors[i], lwd = 2)
  
  # Add label in the middle of each cone
  label_angle <- center_angle
  label_r <- 0.65  # position of label (distance from origin)
  text(label_r * cos(label_angle), label_r * sin(label_angle), 
       paste0(" "), col = "black", cex = 0.8)
}
lines(x_circle, y_circle, col = "black", lwd = 2)

#draw a single vector that starts from origin zero, and its normalized
r <- c(1,0.5)
point1 <-r/sqrt(sum(r^2))
arrows(0, 0, x1 = point1[1], y1 = point1[2], length = 0.15, angle = 30,
       code = 2, col = "brown", lty =1, lwd=2)
text(x=point1[1]-0.15, y=point1[2]+0.04, "λ incenter", cex=0.8)

# Add title at the bottom for explanation
mtext("The storage effect implies a positive covariance between lambdas and alphas", side = 1, line = 2.5)

dev.off()

# Put all together ----
#add the previously saved jpge files 
img1 <- readJPEG("figures/Average conditions.jpg")
img2 <- readJPEG("figures/Relative non linearity lambdas.jpg")
img3 <- readJPEG("figures/Relative non linearity alphas.jpg")
img4 <- readJPEG("figures/The storage effect.jpg")

# Plot the images
par(mfrow = c(2, 2))
plot(img1)
plot(img2)
plot(img3)
plot(img4)


