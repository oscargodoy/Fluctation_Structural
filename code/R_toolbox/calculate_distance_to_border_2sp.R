#' Calculate the (arc) distance to the border of the feasibility domain
#'
#' @description  the distance between an intrinsic growth rate vector inside OR outside the feasibility domain and a border of the domain. It can be calculated by the nearest arc length from \eqn{r} to a border as: \eqn{d_b = \arccos <r(N^*), r(\mathrm{border})>}. 
#'
#' @param A_int the interaction matrix
#' @param r the intrinsic growth rate
#' @param norm whether the arclength is normalised by the maximun possible length (pi).
#' 
#' @return a list with: A numeric value of the arc length. A numeric value indicating wheter r is inside (0 and positive distance) or outside (1 and negative distance) the feasibility domain
#'
#' @examples
#' A_int <- matrix(rnorm(4),nrow=2) ## Generate random interaction matrix with n=2 species
#' r <- A_int %*% c(runif(2, 0, 1)) ## Generate random location
#' myresult <- calculate_distance_to_border_2sp(A_int, r) 
#' cat("distance = ",myresult[[1]],"\n is outside?",myresult[[2]])
#' @export

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
        inferior_comp <- paste("sp", which(aux == min(aux)))
    }else {        
        distance <- - min(arc_length(r_norm,A1),arc_length(r_norm,A2))
        outside <- 1
        inferior_comp <- paste("sp", which(aux == min(aux)))
    }

    if(norm == "yes"){
        distance_norm <- distance / pi
        return(list(distance_norm,outside, inferior_comp))
    } else if(norm == "no"){
        return(list(distance,outside, inferior_comp))
    } else {
        cat("Norm only takes values yes and no. \n")
    }

}
}
