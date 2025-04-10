# Building Information Structures (IS)

# The code in this file builds Information Structures as described in:
# Pablo Almaraz, Piotr Kalita, José A. Langa, Fernando Soler-Toscano,
# "Structural stability of invasion graphs for generalized 
# Lotka--Volterra systems". arXiv:2209.09802
# https://arxiv.org/abs/2209.09802

# Author of the code: Fernando Soler-Toscano - fsoler@us.es

source("code/R_toolbox/lemkelcp.R")
library(quadprog)

# function to get the GASS with the package LemkeLCP
getGASS_LCP_Lemke <- function(alphas, gammas, maxIter = 100){
  sol <- unlist(lemkelcp(matrix(-gammas,nrow = length(alphas)), t(-alphas), maxIter)[1])
  if(is.na(sol[1])){
    stop("LCP has no solution")
  }
  round(sol, digits=12)
}

# function to check if a given vector (gass) satisfied the LCP condition
# to be the GASS for alphas and gammas
checkGASS_LCP <- function(alphas, gammas, gass){
  abs(t(gass) %*% (gammas %*% gass + t(alphas))) < 1e-10
}

# function to get the GASS of given alphas and gammas by using
# quadratic programming
getGASS_QP <- function(alphas, gammas){
  G   <- t(rbind(-gammas,diag(length(alphas))))
  h <- matrix(c(alphas,rep(0,length(alphas))))
  sol <- solve.QP(-t(gammas+t(gammas)), alphas, G, bvec = h)
  round(sol$solution, digits=12)
}

# numerical encoding (1 to 2^n) of an equilibrium point in a system of
# n species. The encoding considers the non-zero positions
getSolPosition <- function(sol){
  sum(2**((1:length(sol))[sol>0]-1))+1
}

# numerical encoding (1 to 2^n) of a (sub)community in a system of n species.
# The encoding considers the present species.
getSubPosition <- function(sub){
  sum(2**(sub-1))+1
}

# Given disjunct (sub)communities 'list1' and 'list2', this function returns
# a vector with the numerical encodings (see 'getSubPosition') of all subcommunities
# having all specied of 'list1' and any combination of species of 'list2'.
getAllSuperPos <- function(list1,list2){
  l2 <- length(list2)
  v2 <- getSubPosition(list1)-1
  sapply(1:(2**l2-1),function(x){
    prov <- as.integer(intToBits(x)[1:l2])*list2
    getSubPosition(prov[prov>0])+v2
  })
}

# Building the Information Structure
# Input:
# - alphas: data frame with n alpha values
# - gammas: data matrix of n*n values
# Output: List with the following elements
# - $points: matrix with all stationary points, one point per row
# - $subsetGASS: data frame indicating the index of the GASS (2n column)
# of each subsystem (1st column)
# - $connectivity: connectivity matrix. Rows and columns are the indexes of
# the points in $points
# - $gassInd: index of the point which is the GASS of the whole system
ISbuild <- function(alphas0, gammas){
  alphas <- as.data.frame(matrix(alphas0,nrow = 1))
  # STEP 1: Finding the GASSes for all subsystems

  # Size of the system (number of species)
  size <- length(alphas)

  # list to store all stationary points
  allStatPoints <- list()

  # Number of stationary points found (trival solution)
  foundSols <- 1

  # Trivial solution (0 for all species) is added to the list
  allStatPoints[[foundSols]] <- as.data.frame(matrix(rep(0.0, size),nrow=1))

  # list indicating which (sub)communities are equilibrium points. This will be
  # used in step 2 to set the connectivity. When a new stationary point is found
  # with species 'subc', the value of listSubsGASS[n] (being n the number
  # encoding subc) is set to the index of allStatPoints containing that point.
  listSubsGASS <- rep(0,2**size)
  listSubsGASS[1] <- 1

  # list indicating index of the stationary point which if the gass of each
  # subsystem. After the first step of the algorithm is completed, each
  # subsystem will have a gass.
  # Value of listGASSindex[scom] is the number of the solution (in the list
  # allStatPoints) which is the gass of the (sub)community scom (en number
  # encoding it).
  listGASSindex <- rep(0,2**size)
  listGASSindex[1] <- 1

  # Data frame to store the GASS in a readable way
  subsetGASS <- data.frame(subset=character(), ind=integer(), stringsAsFactors=FALSE)

  # loop to get the GASSes for all subsystems
  for(s in size:1){
    allSubS  <-combn(1:size,s)    # All subsystems of size s (from size to 1)
    for(i in ncol(allSubS):1){
      ss     <- allSubS[,i]           # 'ss': current subsystem to look for the GASS
      subsCode <- getSubPosition(ss)  # Code of the community (listGASSindex)

      # In case that the GASS of this subcommunity was not found before
      if(listGASSindex[subsCode] == 0){
        # Obtaining the GASS
        gass   <- tryCatch(getGASS_QP(alphas[ss], gammas[ss,ss]), error = function(cond){getGASS_LCP_Lemke(alphas[ss], gammas[ss,ss])})
        v <- rep(0.0, size)
        v[ss] <- gass                   # v contains the GASS with 0s
        gassCode <- getSolPosition(v)   # Code of the obtained gass

        # Check if the stable point was found before
        lookAt <- listSubsGASS[gassCode]  # 0 if the stationary point is new
        if(lookAt == 0){             # New point found
          foundSols <- foundSols+1   # Increase the number of points found
          lookAt <- foundSols
          # Add the new stationary point at the new position:
          allStatPoints[[lookAt]] <- as.data.frame(matrix(v,nrow=1))
          # Indicate that the community represented by 'gassCode' is stationary
          listSubsGASS[gassCode] <- lookAt

          listGASSindex[subsCode] <- lookAt
        }
        # Indicate the GASS for the current community
        listGASSindex[subsCode] <- lookAt
        subsetGASS[nrow(subsetGASS)+1,]<-list(subset=toString(ss), ind=lookAt)

        # Now look for subcommunities that may have the same GASS
        # recall that 'ss' is the current community for which subcommunities
        # with the same GASS are searched
        eqZeroGASS <- ss[gass == 0]     # positions of 'ss' with 0
        if(length(eqZeroGASS) > 0){     # in case there is some 0
          grZeroGASS <- ss[gass > 0]    # positions with a value >0

          # Set the GASS for the subcommunity with all > 0
          if(listGASSindex[getSubPosition(grZeroGASS)] == 0){
            if(checkGASS_LCP(alphas[grZeroGASS], gammas[grZeroGASS,grZeroGASS], v[grZeroGASS])){
              listGASSindex[getSubPosition(grZeroGASS)] <- lookAt
              subsetGASS[nrow(subsetGASS)+1,]<-list(subset=toString(grZeroGASS), ind=lookAt)
            }
          }
          # Set the GASS for intermediate communities in which was not
          # previously set and chechGASS_LCP succeeds
          l2 <- length(eqZeroGASS)
          # Look at all possible subcommunities (subC) between 'grZeroGASS' and
          # 'grZeroGASS + eqZeroGASS'
          if(l2>1){
            for(n in 1:(2**l2-2)){
              subC <- sort(unique(c(grZeroGASS,
                                    as.integer(intToBits(n)[1:l2])*eqZeroGASS)))[-1]
              subCNum <- getSubPosition(subC)  # position of the subcomm to look at
              # In case it was not previously set
              if(listGASSindex[subCNum] == 0){
                # and the found gass is also the stationary point of subC
                if(checkGASS_LCP(alphas[subC], gammas[subC,subC], v[subC])){
                  # register the gass for that subcommunity
                  listGASSindex[subCNum] <- lookAt
                  subsetGASS[nrow(subsetGASS)+1,]<-list(subset=toString(subC), ind=lookAt)
                }
              }
            }
          }
        }
      }
    }
  }
  subsetGASS[nrow(subsetGASS)+1,]<-list(subset="0", ind=1)
  # number of stationary points
  nPoints <- foundSols
  # convert 'allStatPoints' into a matrix
  allStatPoints <- as.matrix(do.call(rbind,allStatPoints))
  rownames(allStatPoints) <- 1:nPoints  # enumerate rows
  colnames(allStatPoints) <- c()        # columns

  # STEP 2: Connecting stationary points
  # Connectivity matrix
  connectivity <- matrix(0, nrow = nPoints, ncol = nPoints)
  # loop to connect all points
  for(p in 1:nPoints){  # p is the index of the point to look for its connections
    point <- allStatPoints[p,]    # stationary point
    I  <- (1:size)[point > 0]     # positions with a value  > 0
    NI <- (1:size)[point == 0]    # positions with a value == 0
    if(length(NI) != 0){          # if there are 0 values
      if(length(I) == 0){   # when all values are 0, we only look at alphas
        R <- alphas[NI]
      } else {              # otherwise, compute r_i values
        R  <- alphas[NI] + gammas[NI,I] %*% as.matrix(point[I])
      }
      J <- NI[(1:length(R))[R>0]] # J is the i's with r_i > 0
      if(length(J) > 0){          # if there are values in J look for connections
        connectivity[p,unique(listGASSindex[getAllSuperPos(I,J)])] <- 1
      }
    }
  }
  # Index of the GASS
  gassInd <- listGASSindex[getSubPosition(1:size)]

  list(
    # The matrix with all stationary points: one point per row
    points = allStatPoints,
    # The dataframe indicating the index of the GASS (2n column)
    # of each subsystem (1st column)
    subsetGASS   = subsetGASS,
    # The connectivity matrix. Rows and columns are the indexes of
    # the points in $points
    connectivity = connectivity,
    # Index of the point which is the GASS of the whole system
    gassInd      = gassInd
  )
}
