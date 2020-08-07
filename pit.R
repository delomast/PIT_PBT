# models with just PIT tag data

library(tidyverse)
library(lubridate)
library(rjags)

if(!dir.exists("PIT_models")) dir.create("PIT_models")

d <- read_csv("data/raw_data.csv", guess_max = 1e6)
# 2016 is first year both PBT and PIT can be used for all release groups
d <- d %>% rename(RAL_RTR = `RAL/RTR`, Expansion_Ref = `Expansion Reference`)

e <- read_csv("data/expansions.csv")
colnames(e) <- gsub(" ", "_", colnames(e))

# filtering and pulling out return year
aPit <- d %>% filter(!is.na(GraLastDate)) %>% 
	select(MY, TagID, RelName, Expansion_Ref, RelSite, GraLastDate) %>% 
	mutate(GraLastDate = mdy(GraLastDate), rYear = year(GraLastDate)) %>% 
	mutate(tagRate = 1 / e$RAL_Expansion[match(Expansion_Ref, e$Expansion_Reference)])
aPit

tagRates <- aPit %>% select(Expansion_Ref, tagRate) %>% distinct
# summary of RAL PIT tag rates
tagRates %>% filter(tagRate < 1) %>% pull(tagRate) %>% summary
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.004225 0.014319 0.022004 0.028286 0.032472 0.140845

# model
modFile <- "PIT_models/simple.txt"
cat("
    model{
   	# likelihood
   	for (g in 1:G){
   		obsPIT[g] ~ dbin(t[g], round(exp(logtotalG[g])))
   	}

		# priors
		for (g in 1:G){
   		logtotalG[g] ~  dnorm(0,1.0E-3)
		}

   	# summary
		for (g in 1:G){
   		totalG[g] <- round(exp(logtotalG[g]))
   	}

    }

    ",
	 file = modFile)

iterations <- 20000
# iterations <- 2000
chains <- 1

# fit for all years, for entire run
estimates <- list()
for(y in 2016:2019){
	# now subset data by return (spawn) year and prep for estimates
	dy <- aPit %>% filter(rYear == y) %>% count(Expansion_Ref) %>% 
		left_join(tagRates, by = "Expansion_Ref")
	inits <- list(	.RNG.name = "base::Mersenne-Twister",
		.RNG.seed = 7,
		logtotalG = log(dy$n / dy$tagRate)
	)

	data <- list(obsPIT = dy$n, t = dy$tagRate, G = nrow(dy))
	
	# fit model
	model <- jags.model(file = modFile, data=data, inits=inits, n.chains = chains)
	samples <- coda.samples(model, c("totalG"), n.iter=iterations)
	
	estimates[[as.character(y)]] <- samples
	
}

save(estimates, file = "pit_only_results.rda")


