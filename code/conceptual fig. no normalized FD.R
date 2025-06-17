# Species interaction matrix
A <- matrix(data = NA, nrow = 2, ncol = 2)
A[1,1] <- -0.75
A[2,2] <- -0.75
A[1,2] <- -0.35
A[2,1] <- -0.35

# Calculate slopes for cone borders
slope1 <- A[1,2]/A[1,1]  # -0.35/-0.75 = 0.4667
slope2 <- A[2,2]/A[2,1]  # -0.75/-0.35 = 2.1429

# Create a simple plot with gray dashed lines at x=1 and y=1
plot(0:1.1, 0:1.1, type = "n", 
     xlim = c(0, 1.1), ylim = c(0, 1.1),
     xlab = "r1", ylab = "r2", main = "")

# Add vertical line at x=1 (gray dashed)
abline(v = 1, col = "gray", lty = 2)

# Add horizontal line at y=1 (gray dashed)
abline(h = 1, col = "gray", lty = 2)

# Draw the cone borders (lines from origin)
x_vals <- c(0, 1.1)
y1_vals <- c(0, 1.1 * slope1)

# For slope2, stop at y=1.1
x_end_slope2 <- 1.1 / slope2  # x value where y=1.1
y2_vals <- c(0, 1.1)
x2_vals <- c(0, x_end_slope2)

lines(x_vals, y1_vals, col = "orange", lwd = 3)
lines(x2_vals, y2_vals, col = "blue", lwd = 3)

# Fill the cone area
x_fill <- c(0, 1.1, 1.1/slope2, 0)
y_fill <- c(0, 1.1 * slope1, 1.1, 0)
polygon(x_fill, y_fill, col = "lightgreen", border = NA)

# Add vector arrow from (0,0) to (0.6, 0.6)
arrows(0, 0, 0.6, 0.6, col = "black", lwd = 2, length = 0.1)
# Add letter "A" at the tip of the arrow
text(0.6, 0.6, "A", pos = 4, col = "black", cex = 1.2, font = 2)

# Add vector arrow B from (0,0) to (0.7, 0.1)
arrows(0, 0, 0.7, 0.1, col = "black", lwd = 2, length = 0.1)

# Add letter "B" at the tip of arrow B
text(0.7, 0.1, "B", pos = 4, col = "black", cex = 1.2, font = 2)

# Calculate shortest distances from A(0.6, 0.6) to each slope line
# Point A coordinates
Ax <- 0.6
Ay <- 0.6

# Distance to slope1 line
dist1_A <- abs(slope1 * Ax - Ay) / sqrt(slope1^2 + 1)
# Distance to slope2 line  
dist2_A <- abs(slope2 * Ax - Ay) / sqrt(slope2^2 + 1)

# Find perpendicular points on each line for point A
x1_perp_A <- (Ax + slope1 * Ay) / (slope1^2 + 1)
y1_perp_A <- slope1 * x1_perp_A

x2_perp_A <- (Ax + slope2 * Ay) / (slope2^2 + 1)
y2_perp_A <- slope2 * x2_perp_A

# Draw perpendicular lines from A
lines(c(Ax, x1_perp_A), c(Ay, y1_perp_A), col = "orange", lwd = 2, lty = )
lines(c(Ax, x2_perp_A), c(Ay, y2_perp_A), col = "blue", lwd = 2, lty = 1)

# Add distance labels on the lines for A
mid1_x_A <- (Ax + x1_perp_A) / 2
mid1_y_A <- (Ay + y1_perp_A) / 2
mid2_x_A <- (Ax + x2_perp_A) / 2
mid2_y_A <- (Ay + y2_perp_A) / 2

text(mid1_x_A, mid1_y_A, sprintf("%.2f", dist1_A), pos = 2, col = "orange", cex = 0.9, font = 2)
text(mid2_x_A, mid2_y_A, sprintf("%.2f", dist2_A), pos = 3, col = "blue", cex = 0.9, font = 2)

# Add vector arrow B from (0,0) to (0.7, 0.1)
arrows(0, 0, 0.7, 0.1, col = "black", lwd = 2, length = 0.1)

# Add letter "B" at the tip of arrow B
text(0.7, 0.1, "B", pos = 4, col = "black", cex = 1.2, font = 2)

# Calculate shortest distances from B(0.7, 0.1) to each slope line
# Point B coordinates
Bx <- 0.7
By <- 0.1

# Distance from point to line ax + by + c = 0 is |ax + by + c|/sqrt(a^2 + b^2)
# Line 1: y = slope1*x => slope1*x - y = 0 => a=slope1, b=-1, c=0
# Line 2: y = slope2*x => slope2*x - y = 0 => a=slope2, b=-1, c=0

# Distance to slope1 line
dist1 <- abs(slope1 * Bx - By) / sqrt(slope1^2 + 1)
# Distance to slope2 line  
dist2 <- abs(slope2 * Bx - By) / sqrt(slope2^2 + 1)

# Find perpendicular points on each line
# For line y = slope1*x, perpendicular from (Bx, By)
# Perpendicular slope = -1/slope1
# Point on line: solve system of equations
x1_perp <- (Bx + slope1 * By) / (slope1^2 + 1)
y1_perp <- slope1 * x1_perp

# For line y = slope2*x, perpendicular from (Bx, By)
x2_perp <- (Bx + slope2 * By) / (slope2^2 + 1)
y2_perp <- slope2 * x2_perp

# Draw perpendicular lines
lines(c(Bx, x1_perp), c(By, y1_perp), col = "orange", lwd = 2, lty = 3)
lines(c(Bx, x2_perp), c(By, y2_perp), col = "blue", lwd = 2, lty = 1)

# Add distance labels on the lines
mid1_x <- (Bx + x1_perp) / 2
mid1_y <- (By + y1_perp) / 2
mid2_x <- (Bx + x2_perp) / 2
mid2_y <- (By + y2_perp) / 2

text(mid1_x, mid1_y, sprintf("%.2f", dist1), pos = 2, col = "orange", cex = 0.9, font = 2)
text(mid2_x, mid2_y, sprintf("%.2f", dist2), pos = 3, col = "blue", cex = 0.9, font = 2)

# Add the text
text(0.2, 0.5, labels = "species 1 excluded", srt = 55
     , adj = 0, col = "blue")

# Add the text
text(0.73, 0.30, labels = "species 2 excluded", srt = 16
     , adj = 0, col = "orange")


#### ahora con dos vectores con misma distancia

# Species interaction matrix
A <- matrix(data = NA, nrow = 2, ncol = 2)
A[1,1] <- -0.75
A[2,2] <- -0.75
A[1,2] <- -0.35
A[2,1] <- -0.35

# Calculate slopes for cone borders
slope1 <- A[1,2]/A[1,1]  # -0.35/-0.75 = 0.4667
slope2 <- A[2,2]/A[2,1]  # -0.75/-0.35 = 2.1429

# Create a simple plot with gray dashed lines at x=1 and y=1
plot(0:1.1, 0:1.1, type = "n", 
     xlim = c(0, 1.1), ylim = c(0, 1.1),
     xlab = "r1", ylab = "r2", main = "")

# Add vertical line at x=1 (gray dashed)
abline(v = 1, col = "gray", lty = 2)

# Add horizontal line at y=1 (gray dashed)
abline(h = 1, col = "gray", lty = 2)

# Draw the cone borders (lines from origin)
x_vals <- c(0, 1.1)
y1_vals <- c(0, 1.1 * slope1)

# For slope2, stop at y=1.1
x_end_slope2 <- 1.1 / slope2  # x value where y=1.1
y2_vals <- c(0, 1.1)
x2_vals <- c(0, x_end_slope2)

lines(x_vals, y1_vals, col = "orange", lwd = 3)
lines(x2_vals, y2_vals, col = "blue", lwd = 3)

# Fill the cone area
x_fill <- c(0, 1.1, 1.1/slope2, 0)
y_fill <- c(0, 1.1 * slope1, 1.1, 0)
polygon(x_fill, y_fill, col = "lightgreen", border = NA)

# Add vector arrow from (0,0) to (0.6, 0.6)
arrows(0, 0, 0.6, 0.6, col = "black", lwd = 2, length = 0.1)
# Add letter "A" at the tip of the arrow
text(0.6, 0.6, "A", pos = 1, col = "black", cex = 1.2, font = 2)

# Add vector arrow B from (0,0) to (0.3, 0.3)
arrows(0, 0, 0.3, 0.3, col = "black", lwd = 2, length = 0.1)
# Add letter "B" at the tip of arrow B
text(0.3, 0.3, "B", pos = 3, col = "black", cex = 1.2, font = 2)

# Calculate shortest distances from A(0.6, 0.6) to each slope line
# Point A coordinates
Ax <- 0.6
Ay <- 0.6

# Distance to slope1 line
dist1_A <- abs(slope1 * Ax - Ay) / sqrt(slope1^2 + 1)
# Distance to slope2 line  
dist2_A <- abs(slope2 * Ax - Ay) / sqrt(slope2^2 + 1)

# Find perpendicular points on each line for point A
x1_perp_A <- (Ax + slope1 * Ay) / (slope1^2 + 1)
y1_perp_A <- slope1 * x1_perp_A

x2_perp_A <- (Ax + slope2 * Ay) / (slope2^2 + 1)
y2_perp_A <- slope2 * x2_perp_A

# Draw perpendicular lines from A
lines(c(Ax, x1_perp_A), c(Ay, y1_perp_A), col = "orange", lwd = 2, lty = 1)
lines(c(Ax, x2_perp_A), c(Ay, y2_perp_A), col = "blue", lwd = 2, lty = 1)

# Add distance labels on the lines for A
mid1_x_A <- (Ax + x1_perp_A) / 2
mid1_y_A <- (Ay + y1_perp_A) / 2
mid2_x_A <- (Ax + x2_perp_A) / 2
mid2_y_A <- (Ay + y2_perp_A) / 2

text(mid1_x_A, mid1_y_A, sprintf("%.2f", dist1_A), pos = 2, col = "orange", cex = 0.9, font = 2)
text(mid2_x_A, mid2_y_A, sprintf("%.2f", dist2_A), pos = 1, col = "blue", cex = 0.9, font = 2)


# Calculate shortest distances from B(0.7, 0.1) to each slope line
# Point B coordinates
Bx <- 0.3
By <- 0.3

# Distance from point to line ax + by + c = 0 is |ax + by + c|/sqrt(a^2 + b^2)
# Line 1: y = slope1*x => slope1*x - y = 0 => a=slope1, b=-1, c=0
# Line 2: y = slope2*x => slope2*x - y = 0 => a=slope2, b=-1, c=0

# Distance to slope1 line
dist1 <- abs(slope1 * Bx - By) / sqrt(slope1^2 + 1)
# Distance to slope2 line  
dist2 <- abs(slope2 * Bx - By) / sqrt(slope2^2 + 1)

# Find perpendicular points on each line
# For line y = slope1*x, perpendicular from (Bx, By)
# Perpendicular slope = -1/slope1
# Point on line: solve system of equations
x1_perp <- (Bx + slope1 * By) / (slope1^2 + 1)
y1_perp <- slope1 * x1_perp

# For line y = slope2*x, perpendicular from (Bx, By)
x2_perp <- (Bx + slope2 * By) / (slope2^2 + 1)
y2_perp <- slope2 * x2_perp

# Draw perpendicular lines
lines(c(Bx, x1_perp), c(By, y1_perp), col = "orange", lwd = 2, lty = 1)
lines(c(Bx, x2_perp), c(By, y2_perp), col = "blue", lwd = 2, lty = 1)

# Add distance labels on the lines
mid1_x <- (Bx + x1_perp) / 2
mid1_y <- (By + y1_perp) / 2
mid2_x <- (Bx + x2_perp) / 2
mid2_y <- (By + y2_perp) / 2

text(mid1_x, mid1_y, sprintf("%.2f", dist1), pos = 2, col = "orange", cex = 0.9, font = 2)
text(mid2_x, mid2_y, sprintf("%.2f", dist2), pos = 3, col = "blue", cex = 0.9, font = 2)

# Add vector arrow C from (0,0) to (0.8, 0.55)
arrows(0, 0, 0.8, 0.53, col = "black", lwd = 2, length = 0.1)
# Add letter "B" at the tip of arrow B
text(0.8, 0.53, "C", pos = 4, col = "black", cex = 1.2, font = 2)


# Calculate shortest distances from B(0.8, 0.55) to each slope line
# Point B coordinates
Cx <- 0.8
Cy <- 0.53

# Distance from point to line ax + by + c = 0 is |ax + by + c|/sqrt(a^2 + b^2)
# Line 1: y = slope1*x => slope1*x - y = 0 => a=slope1, b=-1, c=0
# Line 2: y = slope2*x => slope2*x - y = 0 => a=slope2, b=-1, c=0

# Distance to slope1 line
dist1 <- abs(slope1 * Cx - Cy) / sqrt(slope1^2 + 1)
# Distance to slope2 line  
dist2 <- abs(slope2 * Cx - Cy) / sqrt(slope2^2 + 1)

# Find perpendicular points on each line
# For line y = slope1*x, perpendicular from (Bx, By)
# Perpendicular slope = -1/slope1
# Point on line: solve system of equations
x1_perp <- (Cx + slope1 * Cy) / (slope1^2 + 1)
y1_perp <- slope1 * x1_perp

# For line y = slope2*x, perpendicular from (Bx, By)
x2_perp <- (Cx + slope2 * Cy) / (slope2^2 + 1)
y2_perp <- slope2 * x2_perp

# Draw perpendicular lines
lines(c(Cx, x1_perp), c(Cy, y1_perp), col = "orange", lwd = 2, lty = 1)
lines(c(Cx, x2_perp), c(Cy, y2_perp), col = "blue", lwd = 2, lty = 1)

# Add distance labels on the lines
mid1_x <- (Cx + x1_perp) / 2
mid1_y <- (Cy + y1_perp) / 2
mid2_x <- (Cx + x2_perp) / 2
mid2_y <- (Cy + y2_perp) / 2

text(mid1_x, mid1_y, sprintf("%.2f", dist1), pos = 2, col = "orange", cex = 0.9, font = 2)
text(mid2_x, mid2_y, sprintf("%.2f", dist2), pos = 3, col = "blue", cex = 0.9, font = 2)







