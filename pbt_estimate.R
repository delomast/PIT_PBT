# pbt models

library(tidyverse)
library(lubridate)
library(rjags)

if(!dir.exists("PBT_models")) dir.create("PBT_models")

# the three "really_final_..." files have the filtered data to use for estimates
trap <- read_tsv("filtered_data/really_final_filtered_trap_data.txt", guess_max = 1e4) # these are genotyped samples
trap_LGT <- read_csv("data/LGTrapping_16_19.csv", guess_max = 1e6) %>% # trap data from LGTrapping - 2016-19 for AD / AI
	mutate(CollectionDate = mdy(CollectionDate))
pools <- read_tsv("filtered_data/really_final_pools.txt")
tag_rates <- read_tsv("filtered_data/really_final_tag_rates.txt") %>% 
	bind_rows(tibble(PBT_RELEASE_GROUP = "Unassigned", PBT_RELEASE_GROUP_TAGRATE = 0)) # adding Unassigned
wc <- read_csv("data/window_counts.csv") %>% mutate(OperationsDate = mdy(OperationsDate)) %>% arrange(OperationsDate)
# add standardized week numbers to trap and window counts
wc$WeekNumber <- week(wc$OperationsDate)
wc$Year <- year(wc$OperationsDate)
trap$weekSampled <- week(trap$DateSampled)
trap_LGT$weekSampled <- week(trap_LGT$CollectionDate)
# sum window counts by week and year
wc <- wc %>% group_by(Year, WeekNumber) %>% summarise(count = sum(Total), count = replace_na(count, 0))
wc %>% group_by(Year) %>% summarise(t = sum(count))




######################################################################
######################################################################

y <- 2016 # will be for y in year
# organize data
# totals trapped and counted at window by week
# filter by year and make NA as "Unassigned"
temp_trap <- trap %>% filter(year(DateSampled) == y) %>% mutate(relGroup = replace_na(relGroup, "Unassigned"))
temp_trap_LGT <- trap_LGT %>% filter(year(CollectionDate) == y)
tot_trap <- temp_trap %>% count(weekSampled) %>% rename(t_trap = n)
w_count <- wc %>% filter(Year == y) %>% rename(weekSampled = WeekNumber) %>% left_join(tot_trap, by = "weekSampled") %>%
	mutate(t_trap = replace_na(t_trap, 0)) %>% mutate(diff = count < t_trap) %>% arrange(weekSampled)
# make sure window count is not less than trap count
if(any(w_count$diff)) stop("error wc and trap counts in ", y)
# remove weeks with zero and outside the sp/su managment period
w_count <- w_count %>% filter(count > 0, weekSampled <= max(trap$weekSampled))
# counts of AD and total trapped
w_count <- temp_trap_LGT %>% group_by(weekSampled) %>% summarise(AD_count = sum(LGDMarkAD == "AD", na.rm = TRUE),
																	  AI_count = sum(LGDMarkAD == "AI", na.rm = TRUE)) %>%
	right_join(w_count, by = c("weekSampled")) %>% select(weekSampled, count, t_trap, AD_count, AI_count)

# number PBT tagged fish in each group by week
# matrix with weeks as columns and groups as rows
# AD-clipped
template <- temp_trap %>% filter(grepl("OtsLGRA", temp_trap$Pedigree_name)) %>%
	pull(relGroup) %>% unique %>% expand_grid(relGroup = ., weekSampled = unique(w_count$weekSampled))
AD_PBT_counts <- temp_trap %>% filter(grepl("OtsLGRA", temp_trap$Pedigree_name)) %>% 
	group_by(weekSampled) %>% count(relGroup) %>% right_join(template, by = c("relGroup", "weekSampled")) %>%
	arrange(relGroup, weekSampled) %>% mutate(n = replace_na(n,0)) %>% spread(weekSampled, n)
if(any(as.numeric(colnames(AD_PBT_counts)[2:ncol(AD_PBT_counts)]) != w_count$weekSampled)) stop("error AD colnames in ", y)
AD_tag <- AD_PBT_counts %>% select(relGroup) %>% rename(PBT_RELEASE_GROUP = relGroup) %>% 
	left_join(tag_rates, by = "PBT_RELEASE_GROUP")
if(any(AD_tag$PBT_RELEASE_GROUP != AD_PBT_counts$relGroup)) stop("error AD rownames in ", y)

# AD-intact
template <- temp_trap %>% filter(grepl("OtsLGRU", temp_trap$Pedigree_name)) %>%
	pull(relGroup) %>% unique %>% expand_grid(relGroup = ., weekSampled = unique(w_count$weekSampled))
AI_PBT_counts <- temp_trap %>% filter(grepl("OtsLGRU", temp_trap$Pedigree_name)) %>% 
	group_by(weekSampled) %>% count(relGroup) %>% right_join(template, by = c("relGroup", "weekSampled")) %>%
	arrange(relGroup, weekSampled) %>% mutate(n = replace_na(n,0)) %>% spread(weekSampled, n)
if(any(as.numeric(colnames(AI_PBT_counts)[2:ncol(AI_PBT_counts)]) != w_count$weekSampled)) stop("error AI colnames in ", y)
AI_tag <- AI_PBT_counts %>% select(relGroup) %>% rename(PBT_RELEASE_GROUP = relGroup) %>% 
	left_join(tag_rates, by = "PBT_RELEASE_GROUP")
if(any(AI_tag$PBT_RELEASE_GROUP != AI_PBT_counts$relGroup)) stop("error AI rownames in ", y)
# turn into list to pass to JAGS
all_data <- list(w_count = w_count$count,
					  AD_PBT = as.matrix(AD_PBT_counts[2:ncol(AD_PBT_counts)]),
					  AD_tag = AD_tag$PBT_RELEASE_GROUP_TAGRATE,
					  AI_PBT = as.matrix(AI_PBT_counts[2:ncol(AI_PBT_counts)]),
					  AI_tag = AI_tag$PBT_RELEASE_GROUP_TAGRATE,
					  n_AD = w_count$AD_count,
					  n_AI = w_count$AI_count
					  )


# trying to incorporate window counts
modFile <- "PBT_models/model_pbt_wc.txt"
cat("
	data{
		nWeeks <- length(w_count)
		n_AD_groups <- length(AD_tag) - 1 # -1 to ignore Unassigned group
		n_AI_groups <- length(AI_tag) - 1
		for (w in 1:length(w_count)){
			n_AD_genotype[w] <- sum(AD_PBT[,w])
			n_AI_genotype[w] <- sum(AI_PBT[,w])
			n_clip_trap[w] <- n_AD[w] + n_AI[w]
		}
	}

    model{
   	# likelihood
   	
   	for (w in 1:nWeeks){
   		# proportion clipped
   		n_AD[w] ~ dbin(p_clip[w], n_clip_trap[w])
   		# proportion of each AD group
   		true_AD_PBT[1:(n_AD_groups + 1),w] ~ dmulti(pi_AD[,w], n_AD_genotype[w])
   		for (g in 1:n_AD_groups){
   			AD_PBT[g,w] ~ dbin(AD_tag[g], true_AD_PBT[g,w])
   		}
   		# proportion of each AI group
   		true_AI_PBT[1:(n_AI_groups + 1),w] ~ dmulti(pi_AI[,w], n_AI_genotype[w])
   		for (g in 1:n_AI_groups){
   			AI_PBT[g,w] ~ dbin(AI_tag[g], true_AI_PBT[g,w])
   		}
   		
   	}
   	
   	# walk across weeks
   	for (w in 2:nWeeks){
   		for(g in 1:(n_AD_groups + 1)){
   			log_abund_AD[g,w] ~ dnorm(log_abund_AD[g,w-1], AD_tau[g])
   		}
   		for(g in 1:(n_AI_groups + 1)){
   			log_abund_AI[g,w] ~ dnorm(log_abund_AI[g,w-1], AI_tau[g])
   		}
   	}
   	
   	# conversion between abundance and proportions
   	for (w in 1:nWeeks){
   		for (g in 1:(n_AD_groups + 1)){
   			tempEXP_AD[g,w] <- exp(log_abund_AD[g,w])
   		}
   		for (g in 1:(n_AI_groups + 1)){
   			tempEXP_AI[g,w] <- exp(log_abund_AI[g,w])
   		}
   		tot_AD[w] <- sum(tempEXP_AD[,w])
   		tot_AI[w] <- sum(tempEXP_AI[,w])
   		w_count[w] ~ sum(tot_AD[w], tot_AI[w])
   		
   		p_clip[w] <- tot_AD[w] / w_count[w]
   	}

   	for (w in 1:nWeeks){
   		for (g in 1:(n_AD_groups + 1)){
   			pi_AD[g,w] <- exp(log_abund_AD[g,w]) / tot_AD[w]
   		}
   		for (g in 1:(n_AI_groups + 1)){
   			pi_AI[g,w] <- exp(log_abund_AI[g,w]) / tot_AI[w]
   		}
   	}
   	

		# priors on first week
		for(g in 1:(n_AD_groups + 1)){
			log_abund_AD[g,1] ~ dnorm(0,1.0E-3)
		}
		for(g in 1:(n_AI_groups + 1)){
			log_abund_AI[g,1] ~ dnorm(0,1.0E-3)
		}

   	# priors on walks
	   for(g in 1:(n_AD_groups + 1)){
	   	AD_tau[g] ~ dgamma(1.0E-3, 1.0E-3)
	   	# sig[g] ~ dunif(0,10)
	   	# AD_tau[g] <- 1 / pow(sig[g], 2)
	   	# AD_tau[g] ~ dunif(999,1000)
	   }
		
		for(g in 1:(n_AI_groups + 1)){
			AI_tau[g] ~ dgamma(1.0E-3, 1.0E-3)
		}
    }
    ",
	 file = modFile)
iterations <- 500
chains <- 1
log_abund_AD <- matrix(NA, nrow = length(all_data$AD_tag), ncol = length(all_data$w_count))
log_abund_AI <- matrix(NA, nrow = length(all_data$AI_tag), ncol = length(all_data$w_count))
pClip <- rep(NA, length(all_data$w_count))
for(w in 1:length(all_data$w_count)){
	pClip[w] <- (all_data$n_AD[w] + .01) / (all_data$n_AD[w] + all_data$n_AI[w] + .02)
	log_abund_AD[,w] <- log((all_data$w_count[w] * pClip[w]) * ((all_data$AD_PBT[,w] + .001) / sum(all_data$AD_PBT[,w] + .001)))
	log_abund_AI[,w] <- log((all_data$w_count[w] * (1 - pClip[w])) * ((all_data$AI_PBT[,w] + .001) / sum(all_data$AI_PBT[,w] + .001)))
}
pClip <- log(pClip / (1 - pClip))

inits <- list(logit_p_clip = pClip,
				  log_abund_AD = log_abund_AD,
				  log_abund_AI = log_abund_AI
				  )
inits <- list(log_abund_AD = log_abund_AD,
				  log_abund_AI = log_abund_AI
				  )
log(
(all_data$w_count[18] * (exp(inits$logit_p_clip[18]) / (1 + exp(inits$logit_p_clip[18])))) - sum(exp(inits$log_abund_AD[,18]), na.rm = TRUE)
)

model <- jags.model(file = modFile, data=all_data, n.chains = chains, inits = inits)
samples <- coda.samples(model, c("p_clip"), n.iter=iterations)
round(cbind((all_data$n_AD + all_data$n_AI),
all_data$n_AD / (all_data$n_AD + all_data$n_AI),
apply(samples[[1]],2, mean)
),4)

samples <- coda.samples(model, c("pi_AD"), n.iter=iterations)
temp <- round(apply(samples[[1]],2, mean),4)
cbind(temp[grepl("\\[37,", names(temp))],
		all_data$AD_PBT[37,] / colSums(all_data$AD_PBT),
		colSums(all_data$AD_PBT), all_data$w_count, colSums(all_data$AI_PBT)
)


summary(samples[[1]])

plot(samples[[1]][,1])
plot(samples[[1]][,5])
plot(samples[[1]][,725])
plot(samples[[1]][,"pi_AD[22,18]"])

samples <- coda.samples(model, c("AD_tau"), n.iter=iterations)
plot(samples)
hist(apply(samples[[1]],2,mean), breaks = 15)
# need to mess with priors and possibly grouping early and late weeks where there are no trapped fish

samples <- coda.samples(model, c("pi_AD"), n.iter=10000)

# this performed very poorly. mixing was terrible, only very very small steps
# stayed at initial values. Seemed to be b/c of w_count[w] ~ sum(tot_AD[w], tot_AI[w])
# making it difficult to sample
# one option is to write own sampler...
# other option is to just stratify and forget random walk




