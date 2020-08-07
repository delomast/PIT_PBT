# combined models
# let's first make it work for ONE stratum

library(tidyverse)
library(lubridate)
library(rjags)

detect1 <- rbinom(1, 100, .05)
detect2 <- rbinom(1, 100, .1)

data <- list(obs = c(detect1, detect2),
				 t = c(.05, .1))

# model
modFile <- "combined_models/test.txt"
cat("
    model{
   	# likelihood
   	for (g in 1:2){
   		obs[g] ~ dbin(t[g], round(exp(logtotalG)))
   	}

		# priors
		logtotalG ~  dnorm(0,1.0E-3)
		

   	totalG <- round(exp(logtotalG))
   	

    }

    ",
	 file = modFile)

# iterations <- 20000
iterations <- 2000
chains <- 1

inits <- list(logtotalG = log(1000))

model <- jags.model(file = modFile, data=data, inits = inits, n.chains = chains)
samples <- coda.samples(model, c("totalG"), n.iter=iterations)

summary(samples)

modFile <- "combined_models/test.txt"
cat("
    model{
   	# likelihood
   	for (g in 1:2){
   		obs[g] ~ dbin(t[g], round(exp(logtotalG)))
   	}

		# priors
		logtotalG ~  dnorm(0,1.0E-3)
		

   	totalG <- round(exp(logtotalG))
   	

    }

    ",
	 file = modFile)
