rm(list=ls())

library(dplyr)

source("code/R_toolbox/toolbox_coexistence.R")
source("code/R_toolbox/lemkelcp.R")
source("code/R_toolbox/ISbuild.R")
source("code/R_toolbox/ISgraph.R")
source("code/R_toolbox/figs/circleCones.R")
source("code/R_toolbox/figs/sphereCones.R")
source("code/R_toolbox/figs/sphereSubCones.R")
source("code/R_toolbox/figs/InformationField.R")

#Define a matrix with some degree of facilitation (positive values)
A <-matrix(data =NA, nrow = 2, ncol = 2)
A[1,1] <- -1
A[2,2] <- -1
A[1,2] <- 0.3 #facilitation
A[2,1] <- -0.5

#put several vector on top of the figure 
r1 <- c(-0.2,1)
r2 <- c(0.7,1)
r3 <- c(2,1)
r4 <- c(2.2,0.6)

# 2. Plot the cone without labels 
drawCircleCones(A, allCones = FALSE, drawLabels = FALSE)

#normalize the vector r
point1 <-r1/sqrt(sum(r1^2))
points(point1[1], point1[2], pch=19, col="red", cex=0.75)
arrows(0, 0, x1 = point1[1], y1 = point1[2], length = 0.05, angle = 30,
       code = 2, col = par("fg"), lty = 2,
       lwd = par("lwd"))

point2 <-r2/sqrt(sum(r2^2))
points(point2[1], point2[2], pch=19, col="red", cex=0.75)
arrows(0, 0, x1 = point2[1], y1 = point2[2], length = 0.05, angle = 30,
       code = 2, col = par("fg"), lty = 2,
       lwd = par("lwd"))

point3 <-r3/sqrt(sum(r3^2))
points(point3[1], point3[2], pch=19, col="red", cex=0.75)
arrows(0, 0, x1 = point3[1], y1 = point3[2], length = 0.05, angle = 30,
       code = 2, col = par("fg"), lty = 2,
       lwd = par("lwd"))

point4 <-r4/sqrt(sum(r4^2))
points(point4[1], point4[2], pch=19, col="red", cex=0.75)
arrows(0, 0, x1 = point4[1], y1 = point4[2], length = 0.05, angle = 30,
       code = 2, col = par("fg"), lty = 2,
       lwd = par("lwd"))


  
