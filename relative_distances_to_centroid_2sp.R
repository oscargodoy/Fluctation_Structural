#' Calculate the (arc) relative distance to the centroid of the feasibility domain
#' 
#'
#' @description  the relative distance to the centroid is calculated as the fraction between two distance
#' One distance is between the vector of intrinsic growth rates and the centroid 
#' Another distance is calculated between the edge and the centroid. 
#' 
#'
#' @param A_int the interaction matrix
#' @param r the intrinsic growth rate
#' @param norm whether the arclength is normalised by the maximun possible length (pi).
#' 
#' @return a list with: A numeric value of division between arc lengths. 
#' 
#' 

relative_distances_to_centroid_2sp <- function(A_int, r, norm = "no"){
  norm_vec <- function(x){
    return(x/sqrt(sum(x^2)))
  }
  arc_length <- function(a, b) {
    acos(sum(a * b))
  }
  if(length(r) != 2 || ncol(A_int) != 2){
    cat("This function is only for communities of 2 species. \n")
    
  } else {
    # Calculate centroid of A_int
    centroid <- rowMeans(A_int)
    centroid_norm <- - norm_vec(centroid)
    r_norm <- norm_vec(r)
    A1 <- - norm_vec(A_int[,1])
    A2 <- - norm_vec(A_int[,2])
    aux <- solve(-A_int,r)
    
  
    if(all(aux>0)){     
      distance_r_rc <- arc_length(r_norm, centroid_norm)
      distance_edge_rc <- min(arc_length(centroid_norm, A1), arc_length(centroid_norm, A2))
      outside <- 0
    } else {        
      distance_r_rc <- arc_length(r_norm, centroid_norm)
      distance_edge_rc <- min(arc_length(centroid_norm, A1), arc_length(centroid_norm, A2))
      outside <- 1
    }
    
    if(norm == "yes"){
      distance_r_rc_norm <- distance_r_rc / pi
      distance_edge_rc_norm <- distance_edge_rc / pi
      relative_distance <- distance_r_rc_norm/distance_edge_rc_norm
      
      return(list(relative_distance=relative_distance, distance_r_rc=distance_r_rc_norm, 
                  distance_edge_rc=distance_edge_rc_norm, outside=outside, centroid=centroid))
    } else if(norm == "no"){
      relative_distance <- distance_r_rc/distance_edge_rc
      return(list(relative_distance=relative_distance, distance_r_rc=distance_r_rc, 
                  distance_edge_rc=distance_edge_rc, outside=outside, centroid=centroid))
    } else {
      cat("Norm only takes values yes and no. \n")
    }
  }
}

# Function to plot the results on a unit circle
plot_unit_circle_results <- function(A_int, r, result) {
  # Extract data
  centroid <- result$centroid
  outside <- result$outside
  
  # Normalize vectors to plot on unit circle
  norm_vec <- function(x) x/sqrt(sum(x^2))
  
  r_norm <- norm_vec(r)
  A1_norm <- norm_vec(A_int[,1])
  A2_norm <- norm_vec(A_int[,2])
  centroid_norm <- norm_vec(centroid)
  
  # Set up the plot
  plot(0, 0, type="n", xlim=c(-1, 1), ylim=c(-1, 1), 
       xlab="x", ylab="y", main="Unit Circle with Vectors", asp=1)
  
  # Draw the unit circle
  theta <- seq(0, 2*pi, length.out=100)
  lines(cos(theta), sin(theta))
  
  # Draw vectors
  arrows(0, 0, r_norm[1], r_norm[2], col="blue", lwd=2)
  arrows(0, 0, A1_norm[1], A1_norm[2], col="red", lwd=2)
  arrows(0, 0, A2_norm[1], A2_norm[2], col="red", lwd=2)
  arrows(0, 0, centroid_norm[1], centroid_norm[2], col="green", lwd=2)
  
  # Add points to highlight the tips of the vectors
  points(r_norm[1], r_norm[2], col="blue", pch=16)
  points(A1_norm[1], A1_norm[2], col="red", pch=16)
  points(A2_norm[1], A2_norm[2], col="red", pch=16)
  points(centroid_norm[1], centroid_norm[2], col="green", pch=16)
  
  # Add a legend
  legend("topright", legend=c("r vector", "A1 vector", "A2 vector", "Centroid"), 
         col=c("blue", "red", "red", "green"), lwd=2, pch=16)
  
  # Display status and distances
  status <- ifelse(outside == 0, "Inside", "Outside")
  title(sub=paste("Status:", status, 
                  "| Relative distance:", round(result$relative_distance, 3)))
}

# Example inputs
A_int <- matrix(c(1, 0.5, 
                  0.5, 1), nrow=2, byrow=TRUE)
A_int_minus <-A_int*-1
r <- c(0.2, 0.1)

# Run the function
result <- relative_distances_to_centroid_2sp(A_int_minus, r, norm="no")
result
# Plot the results
plot_unit_circle_results(A_int, r, result)

# Try another example with a different r vector that might be outside
r2 <- c(1, 2)
result2 <- relative_distances_to_centroid_2sp(A_int_minus, r2, norm="no")
result2
plot_unit_circle_results(A_int, r2, result2)

