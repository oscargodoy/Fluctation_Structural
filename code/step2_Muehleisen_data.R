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


# NOTE: dry occurs 62% of the time and wet 38% over the last 50 years (61.2 versus 38.8 in the 100+ timeseries)
pw <- .38
pd <- .62

## READ IN THE MUEHLEISEN DATA
d2vr <- read.table(file = "data/california_vital_rates.csv", header=T, sep=",") %>%
  mutate(prop = ifelse(treatment == "dry", pd, pw))
d2int <- read.table(file = "data/california_interactions.csv", header = T, sep = ",") %>%
  select(species, treatment, AVFA, BRHO, ESCA, LACA, VUMY) %>%
  mutate(prop = ifelse(treatment == "dry", pd, pw))

######## MULTISPECIES ANALYSIS ##########
### Multispecies - average conditions ###
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

species <- unique(d2vr$species)
distout_mean <- data.frame(Measure_Distances(A_mean, r_mean), species) %>%
  mutate(distance = ifelse(extinct == 1, -distance, distance))

### Multispecies - dry conditions ###
vrdry <- d2vr %>%
  filter(treatment == "dry")

r_dry <- vrdry$nu

intdry <- d2int %>%
  filter(treatment == "dry") %>%
  select(-treatment, -prop)

A_dry <- as.matrix(intdry[,2:6])

A_dry <- A_dry*-1

distout_dry <- data.frame(Measure_Distances(A_dry, r_dry), species) %>%
  mutate(distance = ifelse(extinct == 1, -distance, distance)) %>%
  mutate(treatment = "dry", prop = pd)

distout_dry

## Multispecies - wet conditions
vrwet <- d2vr %>%
  filter(treatment == "wet")

r_wet <- vrwet$nu

intwet <- d2int %>%
  filter(treatment == "wet") %>%
  select(-treatment, -prop)

A_wet <- as.matrix(intwet[,2:6])

A_wet <- A_wet*-1

distout_wet <- data.frame(Measure_Distances(A_wet, r_wet), species) %>%
  mutate(distance = ifelse(extinct == 1, -distance, distance)) %>%
  mutate(treatment = "wet", prop = pw)

distout_wet

### COEXISTENCE WITH VARIABILITY ###
distvar <- rbind(distout_wet, distout_dry) %>%
  mutate(distance_weighted = distance * prop) %>%
  group_by(species) %>%
  summarize(dist_withvariability = sum(distance_weighted))


### RELATIVE NONLINEARITY IN LAMBDA ###
# calculate distance of the wet and dry r from the average feasibility domain
distout_wr_ai <- data.frame(Measure_Distances(A_mean, r_wet), species) %>%
  mutate(distance = ifelse(extinct == 1, -distance, distance), 
         propdist = pw*distance)
names(distout_wr_ai) <- c("distance-wr-ai", "extinct-wr-ai","species", "propdist_wr_ai")

distout_dr_ai <- data.frame(Measure_Distances(A_mean, r_dry), species) %>%
  mutate(distance = ifelse(extinct == 1, -distance, distance),
         propdist = pd*distance)
names(distout_dr_ai) <- c("distance-dr-ai", "extinct-dr-ai","species", "propdist_dr_ai")

wr_dr <- left_join(distout_wr_ai, distout_dr_ai) %>%
  mutate(dist_onlylambda = propdist_wr_ai + propdist_dr_ai)


distout_lambda <- left_join(distout_mean, wr_dr) %>%
  select(species, distance, dist_onlylambda) %>%
  mutate(relnonlambda_abs = abs(distance - dist_onlylambda),
         relnonlambda_directional = ifelse(dist_onlylambda > distance, relnonlambda_abs, -relnonlambda_abs))


### RELATIVE NONLINEARITY IN ALPHA ###
# calculate distance of the wet and dry r from the average feasibility domain
distout_mr_wi <- data.frame(Measure_Distances(A_wet, r_mean), species) %>%
  mutate(distance = ifelse(extinct == 1, -distance, distance), 
         propdist = pw*distance)
names(distout_mr_wi) <- c("distance-mr-wi", "extinct-mr-wi","species", "propdist_mr_wi")

distout_mr_di <- data.frame(Measure_Distances(A_dry, r_mean), species) %>%
  mutate(distance = ifelse(extinct == 1, -distance, distance),
         propdist = pd*distance)
names(distout_mr_di) <- c("distance-mr-di", "extinct-mr-di","species", "propdist_mr_di")

wi_di <- left_join(distout_mr_wi, distout_mr_di) %>%
  mutate(dist_onlyalpha = propdist_mr_wi + propdist_mr_di)


distout_alpha <- left_join(distout_mean, wi_di) %>%
  select(species, distance, dist_onlyalpha) %>%
  mutate(relnonalpha_abs = abs(distance - dist_onlyalpha),
         relnonalpha_directional = ifelse(dist_onlyalpha > distance, relnonalpha_abs, -relnonalpha_abs))

distout_relnon <- left_join(distout_lambda, distout_alpha) %>%
  select(-relnonlambda_abs, -relnonalpha_abs)

distout_all <- left_join(distout_relnon, distvar) %>%
  mutate(storage_abs = abs(distance - dist_withvariability),
         storage_directional = ifelse(dist_withvariability > distance, storage_abs, -storage_abs)) %>%
  select(-storage_abs) %>%
  mutate(check = distance + relnonlambda_directional + relnonalpha_directional + storage_directional,
         check2 = relnonlambda_directional + relnonalpha_directional + storage_directional)
  select(-relnonlambda_directional, -relnonalpha_directional)


  
# ## Calculate for every pairwise combination of Muehleisen
  
  pairwise_outcomes <- data.frame(species1=character(), species2=character(), de=numeric(), extinct=numeric(), inferior = character(), scenario = character(), prop = numeric())
  
  pairwise_calculations <- data.frame(species1=character(), species2=character(), de=numeric(), de_rellambda = numeric(), de_relalpha = numeric(),
                                      de_storage = numeric(), rellambda = numeric(), relalpha = numeric(), storage = numeric(), superior = character(), swap_rellambda = numeric(),
                                      superior_relnonlambda = numeric(), swap_relalpha = numeric(), superior_relnonalpha = numeric(), swap_storage = numeric())
  # create vector of species
  splist <- unique(d2vr$species)
  
  #create matrix of every unique species combo
  spcombos <- combn(splist, 2)
  
  for (i in 1:ncol(spcombos)) {
    #subset to focal species pairs
    spcomboi <- spcombos[,i]
    
    ### mean conditions - create vector and matrix ###
    vrsubmean <- subset(vrmean, species%in%c(spcomboi))
    intsubmean <- subset(intmean, species%in%c(spcomboi)) %>%
      pivot_longer(names_to = "competitor", values_to ="value", AVFA_sum:VUMY_sum) %>%
      separate(competitor, into = c("competitor", "todelete")) %>%
      select(-todelete) %>%
      filter(competitor%in%c(spcombos[,i])) %>%
      pivot_wider(names_from = competitor, values_from = value)
    
    # create growth rate vector and interaction matrix
    r_mean <- vrsubmean$nu_weightedmean
    A_mean <- as.matrix(intsubmean[,2:3])
    
    #we multiply by -1 as competition should be negative
    A_mean <- A_mean*-1
    
    ### dry conditions - create vector and matrix ###
    vrsubdry <- subset(vrdry, species%in%c(spcomboi))
    intsubdry <- subset(intdry, species%in%c(spcomboi)) %>%
      pivot_longer(names_to = "competitor", values_to ="value", AVFA:VUMY) %>%
      filter(competitor%in%c(spcombos[,i])) %>%
      pivot_wider(names_from = competitor, values_from = value)
    
    r_dry <- vrsubdry$nu
    A_dry <- as.matrix(intsubdry[,2:3])*-1
      

    ### wet conditions - create vector and matrix ###
    vrsubwet <- subset(vrwet, species%in%c(spcomboi))
    intsubwet <- subset(intwet, species%in%c(spcomboi)) %>%
      pivot_longer(names_to = "competitor", values_to ="value", AVFA:VUMY) %>%
      filter(competitor%in%c(spcombos[,i])) %>%
      pivot_wider(names_from = competitor, values_from = value)
    
    r_wet <- vrsubwet$nu
    A_wet <- as.matrix(intsubwet[,2:3])*-1


    ## coexistence in average and static conditions
    de_mean <- calculate_distance_to_border_2sp(A_mean, r_mean)
    de_dry <- calculate_distance_to_border_2sp(A_dry, r_dry)
    de_wet <- calculate_distance_to_border_2sp(A_wet, r_wet)
    
    ## relative nonlinearity in lambdas
    de_wr_ai <- calculate_distance_to_border_2sp(A_mean, r_wet)
    de_dr_ai <- calculate_distance_to_border_2sp(A_mean, r_dry)
    
    ## relative nonlinearity in alphas
    de_ar_wi <- calculate_distance_to_border_2sp(A_wet, r_mean)
    de_ar_di <- calculate_distance_to_border_2sp(A_dry, r_mean)
    
    
    

    de_mean_out <- data.frame(species1= spcomboi[1], species2 = spcomboi[2], de = de_mean[[1]][1], extinct = de_mean[[2]][1], superior = de_mean[[3]][1], scenario = "average", prop = 1)
   
    de_dry_out <- data.frame(species1= spcomboi[1], species2 = spcomboi[2], de = de_dry[[1]][1], extinct = de_dry[[2]][1], superior = de_dry[[3]][1], scenario = "dry", prop = pd)
    de_wet_out <- data.frame(species1= spcomboi[1], species2 = spcomboi[2], de = de_wet[[1]][1], extinct = de_wet[[2]][1], superior = de_wet[[3]][1], scenario = "wet", prop = pw)
   
    sharedwinner_storage <- ifelse(identical(de_dry_out$superior, de_wet_out$superior) == T, 1, 0)
     storage_out <- rbind(de_dry_out, de_wet_out) %>%
       mutate(de_weighted = de*prop) %>%
       mutate(de_storage = ifelse(identical(de_dry_out$superior, de_wet_out$superior) & sum(de_dry_out$extinct, de_wet_out$extinct) == 2, sum(de_weighted), sum(abs(de_weighted))),
              swap_storage = ifelse(!identical(de_dry_out$superior, de_wet_out$superior) & (sum(de_dry_out$extinct, de_wet_out$extinct) == 2), 1, 0)) %>%
       mutate(superior_storage = superior) %>%
       select(-(de:de_weighted))
       
    de_wr_ai_out <- data.frame(species1= spcomboi[1], species2 = spcomboi[2], de = de_wr_ai[[1]][1], extinct = de_wr_ai[[2]][1], superior = de_wr_ai[[3]][1], scenario = "wr_ai", prop = pw)
    de_dr_ai_out <- data.frame(species1= spcomboi[1], species2 = spcomboi[2], de = de_dr_ai[[1]][1], extinct = de_dr_ai[[2]][1], superior = de_dr_ai[[3]][1], scenario = "dr_ai", prop = pd)
    relnonlambda_out <- rbind(de_wr_ai_out, de_dr_ai_out) %>%
      mutate(de_weighted = de*prop) %>%
      mutate(de_rellambda = ifelse(identical(de_wr_ai_out$superior, de_dr_ai_out$superior) & sum(de_wr_ai_out$extinct, de_dr_ai_out$extinct) == 2, sum(de_weighted), sum(abs(de_weighted))),
             swap_rellambda = ifelse(!identical(de_wr_ai_out$superior, de_dr_ai_out$superior) & (sum(de_wr_ai_out$extinct, de_dr_ai_out$extinct) == 2), 1, 0)) %>%
      mutate(superior_relnonlambda = superior) %>%
      select(-(de:de_weighted))
    
    de_ar_wi_out <- data.frame(species1= spcomboi[1], species2 = spcomboi[2], de = de_ar_wi[[1]][1], extinct = de_ar_wi[[2]][1], superior = de_ar_wi[[3]][1], scenario = "ar_wi", prop = pw)
    de_ar_di_out <- data.frame(species1= spcomboi[1], species2 = spcomboi[2], de = de_ar_di[[1]][1], extinct = de_ar_di[[2]][1], superior = de_ar_di[[3]][1], scenario = "ar_di", prop = pd)
    relnonalpha_out <- rbind(de_ar_wi_out, de_ar_di_out) %>%
      mutate(de_weighted = de*prop) %>%
      mutate(de_relalpha = ifelse(identical(de_ar_wi_out$superior, de_ar_di_out$superior) & sum(de_ar_wi_out$extinct, de_ar_di_out$extinct) == 2, sum(de_weighted), sum(abs(de_weighted))),
             swap_relalpha = ifelse(!identical(de_ar_wi_out$superior, de_ar_di_out$superior) & (sum(de_ar_wi_out$extinct, de_ar_di_out$extinct) == 2), 1, 0)) %>%
      mutate(superior_relnonalpha = superior) %>%
      select(-(de:de_weighted))
    
    
    
    calculations_temp <- full_join(full_join(de_mean_out %>% select(-extinct, -scenario, -prop), full_join(storage_out[1,], relnonlambda_out[1,])), relnonalpha_out[1,]) %>%
      mutate(rellambda = de_rellambda - de, relalpha = de_relalpha - de, storage = de_storage - de) %>%
      select(species1, species2, de, de_rellambda, de_relalpha, de_storage, rellambda, relalpha, storage,  swap_rellambda, superior_relnonlambda, swap_relalpha, superior, superior_relnonalpha, swap_storage, superior_storage) 
    
    pairwise_calculations <- rbind(pairwise_calculations, calculations_temp)
    pairwise_outcomes <- rbind( pairwise_outcomes, de_mean_out, de_dry_out, de_wet_out, de_wr_ai_out, de_dr_ai_out, de_ar_wi_out, de_ar_di_out)
    
  }
    
    
    #draw the cone - mean
    draw_biodiversity_cone(A_mean, spcomboi[1], spcomboi[2])
    #draw the vector
    point <-r_mean/sqrt(sum(r_mean^2))
    arrows(0, 0, x1 = point[1], y1 = point[2], length = 0.15, angle = 30,
           code = 2, col = "brown", lty =1, lwd=2)
    text(0.8, 1, "A. Average conditions", cex = 1.2, pos = 1)
    
    
    
    
    
    
    
    
    
  }
  
  pairwise_average
  
  
  
  
  
# ## Calculate for every pairwise combination of Muehleisen
# 
# pairwise_average <- data.frame(species1=character(), species2=character(), de=numeric(), extinct=numeric(), inferior = character())
# # create vector of species
# splist <- unique(d2vr$species)
# 
# #create matrix of every unique species combo
# spcombos <- combn(splist, 2)
# 
# for (i in 1:ncol(spcombos)) {
#   #subset to focal species pairs
#   spcomboi <- spcombos[,i]
#   vrsub <- subset(d2vr, species%in%c(spcomboi))
#   intsub <- subset(d2int, species%in%c(spcomboi))
#   
#   intsub2 <- intsub %>%
#     pivot_longer(names_to = "competitor", values_to ="value", AVFA:VUMY) %>%
#     filter(competitor%in%c(spcombos[,i])) %>%
#     pivot_wider(names_from = competitor, values_from = value)
#   
# # take the mean vital rates across conditions (should this be weighted by treatment frequency??)
# vrmeanx <- vrsub %>%
#   select(-treatment) %>%
#   group_by(species) %>%
#   summarize_all(list(mean = mean))
# 
# # take the mean interaction strengths across conditions (should this be weighted by treatment frequency??)
# intmeanx <- intsub2 %>%
#   select(-treatment) %>%
#   group_by(species) %>%
#   summarize_all(list(mean = mean)) 
# 
# # create growth rate vector and interaction matrix
# r <- vrmeanx$nu_mean
# A <- as.matrix(intmeanx[,2:3])
# 
# #we multiply by -1 as competition should be negative
# A <- A*-1
# 
# de <- calculate_distance_to_border_2sp(A, r)
# 
# tempout <- data.frame(species1= spcomboi[1], species2 = spcomboi[2], de = de[[1]][1], extinct = de[[2]][1], inferior = de[[3]][1])
# pairwise_average <- rbind( pairwise_average, tempout)
# 
# 
# #draw the cone
# draw_biodiversity_cone(A, spcomboi[1], spcomboi[2])
# #draw the vector
# point <-r/sqrt(sum(r^2))
# arrows(0, 0, x1 = point[1], y1 = point[2], length = 0.15, angle = 30,
#        code = 2, col = "brown", lty =1, lwd=2)
# text(0.8, 1, "A. Average conditions", cex = 1.2, pos = 1)
# 
# 
# }
# 
# pairwise_average
