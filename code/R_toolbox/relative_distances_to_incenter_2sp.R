#' Calculate the (arc) relative distance to the incenter of the feasibility domain
#' 
#'
#' @description  the relative distance to the incenter is calculated as the fraction between two distance
#' One distance is between the vector of intrinsic growth rates and the incenter 
#' Another distance is calculated between the edge and the incenter. 
#' 
#'
#' @param A_int the interaction matrix
#' @param r the intrinsic growth rate
#' @param norm whether the arclength is normalised by the maximun possible length (pi).
#' 
#' @return a list with: A numeric value of division between arc lengths. 
#' 
#' 

relative_distances_to_incenter_2sp <- function(A_int, r, norm = "no"){
  norm_vec <- function(x){
    return(x/sqrt(sum(x^2)))
  }
  arc_length <- function(a, b) {
    acos(sum(a * b))
  }
  if(length(r) != 2 || ncol(A_int) != 2){
    cat("This function is only for communities of 2 species. \n")
    
  } else {
    # Calculate incenter of A_int
    incenter <- rowMeans(A_int)
    incenter_norm <- - norm_vec(incenter)
    r_norm <- norm_vec(r)
    A1 <- - norm_vec(A_int[,1])
    A2 <- - norm_vec(A_int[,2])
    aux <- solve(-A_int,r)
    
  
    if(all(aux>0)){     
      distance_r_rc <- arc_length(r_norm, incenter_norm)
      distance_edge_rc <- min(arc_length(incenter_norm, A1), arc_length(incenter_norm, A2))
      outside <- 0
    } else {        
      distance_r_rc <- arc_length(r_norm, incenter_norm)
      distance_edge_rc <- min(arc_length(incenter_norm, A1), arc_length(incenter_norm, A2))
      outside <- 1
    }
    
    if(norm == "yes"){
      distance_r_rc_norm <- distance_r_rc / pi
      distance_edge_rc_norm <- distance_edge_rc / pi
      relative_distance <- distance_r_rc_norm/distance_edge_rc_norm
      
      return(list(relative_distance=relative_distance, distance_r_rc=distance_r_rc_norm, 
                  distance_edge_rc=distance_edge_rc_norm, outside=outside, incenter=incenter))
    } else if(norm == "no"){
      relative_distance <- distance_r_rc/distance_edge_rc
      return(list(relative_distance=relative_distance, distance_r_rc=distance_r_rc, 
                  distance_edge_rc=distance_edge_rc, outside=outside, incenter=incenter))
    } else {
      cat("Norm only takes values yes and no. \n")
    }
  }
}
