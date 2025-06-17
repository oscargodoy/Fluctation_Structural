#code to transform vital rates in order to put Beverton-Holt into LV units. 
#following appendix S1 Godoy and Levine 2014

#create the function, nu is the key parameter which is equivalente to the r in the LV
#nu=(lambda*germination)/beta where beta=(1-(1-g)*s), these are the seeds that do not stay
#in the soil bank due to death or germination. 

nu <- function(lambda, germination, survival) {
  beta <- 1 - (1 - germination) * survival
  nu_value <- (lambda * germination) / beta
  return(nu_value)
}

#load data
lambda <- read.csv("data/california_data.csv", header = TRUE)
lambda <- lambda[lambda$species != "TRHI", ] #remove species TRHI
species <- unique(lambda$species) #get the species names
g_s <- read.csv("data/california_germination_survival.csv", header = TRUE, sep=";")
#remove column names from file g_s
g_s <- g_s[,-1]


#calculate nu for each species and environment wet and dry using the nu function
#select the rows with wet as as treatment remove species TRHI
lambda_wet <- lambda[lambda$treatment == "wet", "lambda"]
lambda_dry <- lambda[lambda$treatment == "dry", "lambda"]
germination_wet <- g_s[g_s$treatment == "wet", "germination"]
germination_dry <- g_s[g_s$treatment == "dry", "germination"]
survival_wet <- g_s[g_s$treatment == "wet", "survival"]
survival_dry <- g_s[g_s$treatment == "dry", "survival"]

nu_wet <- nu(lambda_wet, germination_wet, survival_wet)
nu_dry <- nu(lambda_dry, germination_dry, survival_dry)
#add nu to the lambda data frame with species code

nu_wet <- cbind(species, treatment = "wet", nu = nu_wet)
nu_dry <- cbind(species, treatment = "dry", nu = nu_dry)
#combine the two data frames
nu_combined <- rbind(nu_wet, nu_dry)
#join with the data_frame with lambda
lambda <- merge(lambda, nu_combined, by = c("species", "treatment"), all.x = TRUE)
#add to lambda the germination and survival rates
lambda <- merge(lambda, g_s, by = c("species", "treatment"), all.x = TRUE)
#reorder the columns
vital_rates <- lambda[, c("species", "treatment", "lambda", "germination", "survival", "nu")]
#save the data frame with vital rates
write.csv(vital_rates, "data/california_vital_rates.csv", row.names = FALSE)

#save interactions as a csv file
#remove from lambda the columns = order, lambda, nu, survival, and germination
interactions <- lambda[, c("species", "treatment", species)]
write.csv(interactions, "data/california_interactions.csv", row.names = FALSE)



