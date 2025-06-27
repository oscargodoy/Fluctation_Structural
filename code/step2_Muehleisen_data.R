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


## READ IN THE MUEHLEISEN DATA
d2vr <- read.table(file = "data/california_vital_rates.csv", header=T, sep=",") %>%
  mutate(prop = ifelse(treatment == "dry", .62, .38))
d2int <- read.table(file = "data/california_interactions.csv", header = T, sep = ",") %>%
  select(species, treatment, AVFA, BRHO, ESCA, LACA, VUMY) %>%
  mutate(prop = ifelse(treatment == "dry", .62, .38))

# NOTE: dry occurs 62% of the time and wet 38% over the last 50 years (61.2 versus 38.8 in the 100+ timeseries)


######## MULTISPECIES ANALYSIS ######
## Multispecies - average conditions
vrmean <- d2vr %>%
  mutate(nu_weighted = nu*prop) %>%
  group_by(species) %>%
  summarize(nu_weightedmean = sum(nu_weighted))

r_mean <- vrmean$nu_weightedmean

intmean <- d2int %>%
  mutate(AVFA = AVFA*prop, BRHO = BRHO*prop, ESCA = ESCA*prop, LACA = LACA*prop, VUMY = VUMY*prop) %>% 
  select(-treatment, -prop) %>%
  group_by(species) %>%
  summarize_all(list(sum = sum)) 

A_mean <- as.matrix(intmean[,2:6])

#we multiply by -1 as competition should be negative
A_mean <- A_mean*-1

distout_m <- Measure_Distances(A_mean, r_mean)
species <- unique(d2vr$species)
distout_mean <- data.frame(distout_m, species)
rm(distout_m)

## Multispecies - dry conditions
vrdry <- d2vr %>%
  filter(treatment == "dry")

r_dry <- vrdry$nu

intdry <- d2int %>%
  filter(treatment == "dry") %>%
  select(-treatment, -prop)

A_dry <- as.matrix(intdry[,2:6])

#we multiply by -1 as competition should be negative
A_dry <- A_dry*-1

distout_d <- Measure_Distances(A_dry, r_dry)
distout_dry <- data.frame(distout_d, species)
rm(distout_d)

## Multispecies - wet conditions
vrwet <- d2vr %>%
  filter(treatment == "wet")

r_wet <- vrwet$nu

intwet <- d2int %>%
  filter(treatment == "wet") %>%
  select(-treatment, -prop)

A_wet <- as.matrix(intwet[,2:6])

#we multiply by -1 as competition should be negative
A_wet <- A_wet*-1

distout_w <- Measure_Distances(A_wet, r_wet)
distout_wet <- data.frame(distout_w, species)
rm(distout_w)






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
