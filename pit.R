# models with just PIT tag data

library(tidyverse)
library(lubridate)
library(rjags)

if(!dir.exists("PIT_models")) dir.create("PIT_models")

d <- read_csv("data/raw_data.csv", guess_max = 1e6)
# 2013 is first year both PBT and PIT can be used for all release groups
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

# looking at distribution of returns by group x year
return_by_year <- aPit %>% group_by(rYear) %>% count(Expansion_Ref) %>% 
	left_join(tagRates, by = "Expansion_Ref") %>%
	mutate(t = n / tagRate) %>% pull(t)
hist(log(return_by_year), freq = FALSE)
points(seq(0,12,.1), dnorm(seq(0,12,.1),sd = sqrt(1/.001)))

plot(seq(0,12,.1), dnorm(seq(0,12,.1),sd = sqrt(1/.001)))


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
   		logtotalG[g] ~  dnorm(0,1.0E-3) T(0,)
		}

   	# summary
		for (g in 1:G){
   		totalG[g] <- round(exp(logtotalG[g]))
   	}

    }

    ",
	 file = modFile)

iterations <- 10000
chains <- 1

# fit for all years, for entire run
estimates <- list()
for(y in 2013:2019){
	# now subset data by return (spawn) year and prep for estimates
	dy <- aPit %>% filter(rYear == y) %>% count(Expansion_Ref) %>% 
		left_join(tagRates, by = "Expansion_Ref")
	inits <- list(logtotalG = log(dy$n / dy$tagRate))
	data <- list(obsPIT = dy$n, t = dy$tagRate, G = nrow(dy))
	
	# fit model
	model <- jags.model(file = modFile, data=data, inits=inits, n.chains = chains)
	samples <- coda.samples(model, c("totalG"), n.iter=iterations)
	
	estimates[[as.character(y)]] <- samples
	
}

save(estimates, file = "pit_only_results.rda")

# look at some results

plot(samples[[1]][,2])
m <- samples[[1]] %>% as.data.frame %>% select(contains("totalG")) %>%
	apply(2, mean)
cbind(m, exp(inits[[1]]))

plot(m, exp(inits[[1]]))
abline(0,1)
hist(m - exp(inits[[1]]), breaks = 20)
hist((m - exp(inits[[1]]))/exp(inits[[1]]), breaks = 20)
names(m) <- dy$Expansion_Ref
m
sum(m)




