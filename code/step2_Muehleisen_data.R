rm(list=ls())

library(dplyr)
require(devtools)
#install_github("RadicalCommEcol/anisoFun")
library(anisoFun) #this package includes most of the functions for structural stability.

source("code/R_toolbox/toolbox_coexistence.R")
source("code/R_toolbox/lemkelcp.R")
source("code/R_toolbox/ISbuild.R")
source("code/R_toolbox/ISgraph.R")
source("code/R_toolbox/figs/circleCones.R")
source("code/R_toolbox/figs/sphereCones.R")
source("code/R_toolbox/figs/sphereSubCones.R")
source("code/R_toolbox/figs/InformationField.R")
source("code/R_toolbox/Functions to plot Muehleisen data.R")
source("code/R_toolbox/calculate_distance_to_border_2sp.R")
source("code/R_toolbox/multispecies distance to the edge.r")

par(mar = c(3, 3, 1.5, 1.5))  # Set margins
# Muehleisen Average conditions----


d2vr <- read.table(file = "data/california_vital_rates.csv", header=T, sep=",")
d2int <- read.table(file = "data/california_interactions.csv", header = T, sep = ",") %>%
  select(species, treatment, AVFA, BRHO, ESCA, LACA, VUMY)

## Multispecies coexistence under average conditions
vrmean <- d2vr %>%
  select(-treatment) %>%
  group_by(species) %>%
  summarize_all(list(mean = mean))

r <- vrmean$nu_mean

intmean <- d2int %>%
  select(-treatment) %>%
  group_by(species) %>%
  summarize_all(list(mean = mean)) %>%
  select(species, AVFA_mean, BRHO_mean, ESCA_mean, LACA_mean, VUMY_mean)

A <- as.matrix(intmean[,2:6])

#we multiply by -1 as competition should be negative
A <- A*-1

distout <- Measure_Distances(A, r)
species <- unique(d2vr$species)
avgdistout <- data.frame(distout, splist)


## Calculate for every pairwise combination of Muehleisen

pairwise_average <- data.frame(species1=character(), species2=character(), de=numeric(), extinct=numeric(), inferior = character())
# create vector of species
splist <- unique(d2vr$species)

#create matrix of every unique species combo
spcombos <- combn(splist, 2)

for (i in 1:ncol(spcombos)) {
  #subset to focal species pairs
  spcomboi <- spcombos[,i]
  vrsub <- subset(d2vr, species%in%c(spcomboi))
  intsub <- subset(d2int, species%in%c(spcomboi))
  
  intsub2 <- intsub %>%
    pivot_longer(names_to = "competitor", values_to ="value", AVFA:VUMY) %>%
    filter(competitor%in%c(spcombos[,i])) %>%
    pivot_wider(names_from = competitor, values_from = value)
  
# take the mean vital rates across conditions (should this be weighted by treatment frequency??)
vrmean <- vrsub %>%
  select(-treatment) %>%
  group_by(species) %>%
  summarize_all(list(mean = mean))

# take the mean interaction strengths across conditions (should this be weighted by treatment frequency??)
intmean <- intsub2 %>%
  select(-treatment) %>%
  group_by(species) %>%
  summarize_all(list(mean = mean)) 

# create growth rate vector and interaction matrix
r <- vrmean$nu_mean
A <- as.matrix(intmean[,2:3])

#we multiply by -1 as competition should be negative
A <- A*-1

de <- calculate_distance_to_border_2sp(A, r)

tempout <- data.frame(species1= spcomboi[1], species2 = spcomboi[2], de = de[[1]][1], extinct = de[[2]][1], inferior = de[[3]][1])
pairwise_average <- rbind( pairwise_average, tempout)


#draw the cone
draw_biodiversity_cone(A, spcomboi[1], spcomboi[2])
#draw the vector
point <-r/sqrt(sum(r^2))
arrows(0, 0, x1 = point[1], y1 = point[2], length = 0.15, angle = 30,
       code = 2, col = "brown", lty =1, lwd=2)
text(0.8, 1, "A. Average conditions", cex = 1.2, pos = 1)


}

pairwise_average
