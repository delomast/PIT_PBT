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

# stratifying
modFile <- "PBT_models/model_pbt_strat.txt"
cat("
	# note: everything is set up as weeks, but really stratum = week
	data{
		nWeeks <- length(w_count)
		n_AD_groups <- length(AD_tag) - 1 # -1 to ignore Unassigned group
		n_AI_groups <- length(AI_tag) - 1
		for (w in 1:length(w_count)){
			n_AD_genotype[w] <- sum(AD_PBT[,w])
			n_AI_genotype[w] <- sum(AI_PBT[,w])
			n_clip_trap[w] <- n_AD[w] + n_AI[w]
			alpha_AD[1:length(AD_tag),w] <- ifelse(AD_PBT[,w] > 0, .1, 1e-20)
   		alpha_AI[1:length(AI_tag),w] <- ifelse(AI_PBT[,w] > 0, .1, 1e-20)
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
   	
   	# priors
   	for (w in 1:nWeeks){
   		p_clip[w] ~ dbeta(.01,.01)
   	
   		pi_AD[1:(n_AD_groups + 1),w] ~ ddirich(alpha_AD[,w])
   		
			pi_AI[1:(n_AI_groups + 1),w] ~ ddirich(alpha_AI[,w])
   	}

   	# summary - mult by window counts and add up
   	for (w in 1:nWeeks){
   		tot_AD[w] <- p_clip[w] * w_count[w]
   		tot_AI[w] <- w_count[w] - tot_AD[w]
   		tot_AD_groups[1:(n_AD_groups + 1),w] <- pi_AD[,w] * tot_AD[w]
   		tot_AI_groups[1:(n_AI_groups + 1),w] <- pi_AI[,w] * tot_AI[w]
   		
   	}
   	for(g in 1:(n_AD_groups + 1)){
   		grand_tot_AD_groups[g] <- sum(tot_AD_groups[g,])
   	}
   	for(g in 1:(n_AI_groups + 1)){
   		grand_tot_AI_groups[g] <- sum(tot_AI_groups[g,])
   	}
   	grand_tot_AD <- sum(tot_AD)
   	grand_tot_AI <- sum(tot_AI)
   	
   	
    }

    ",
	 file = modFile)


PBT_estimates <- list()
# y <- 2019
for(y in 2016:2019){
	# organize data
	# totals trapped and counted at window by week
	# filter by year and make NA as "Unassigned"
	temp_trap <- trap %>% filter(year(DateSampled) == y) %>% mutate(relGroup = replace_na(relGroup, "Unassigned"))
	temp_trap_LGT <- trap_LGT %>% filter(year(CollectionDate) == y)
	tot_trap <- temp_trap %>% count(weekSampled) %>% rename(t_trap = n)
	w_count <- wc %>% filter(Year == y) %>% rename(weekSampled = WeekNumber) %>% left_join(tot_trap, by = "weekSampled") %>%
		mutate(t_trap = replace_na(t_trap, 0)) %>% mutate(diff = count < t_trap) %>% arrange(weekSampled)
	# make sure window count is not less than trap count
	if(any(w_count$diff & w_count$t_trap > 0)) stop("error wc and trap counts in ", y)
	# remove weeks with zero and outside the sp/su managment period
	w_count <- w_count %>% filter(count != 0, weekSampled <= max(trap$weekSampled))
	# counts of AD and total trapped
	w_count <- temp_trap_LGT %>% group_by(weekSampled) %>% summarise(AD_count = sum(LGDMarkAD == "AD", na.rm = TRUE),
																		  AI_count = sum(LGDMarkAD == "AI", na.rm = TRUE)) %>%
		right_join(w_count, by = c("weekSampled")) %>% mutate(AD_count = replace_na(AD_count, 0), AI_count = replace_na(AI_count, 0)) %>%
		select(weekSampled, count, t_trap, AD_count, AI_count) %>% arrange(weekSampled)
	# choosing strata
	# as.data.frame(w_count)
	# grouping early and late weeks into strata as needed
	# manually determined targeting minimum number genotyped ~50
	if( y == 2016){
		w_count$weekSampled[w_count$weekSampled < 18] <- 17
		temp_trap$weekSampled[temp_trap$weekSampled < 18] <- 17
	} else if (y == 2017){
		w_count$weekSampled[w_count$weekSampled < 21] <- 20
		temp_trap$weekSampled[temp_trap$weekSampled < 21] <- 20
		w_count$weekSampled[w_count$weekSampled > 30] <- 31
		temp_trap$weekSampled[temp_trap$weekSampled > 30] <- 31
	} else if (y == 2018){
		w_count$weekSampled[w_count$weekSampled < 20] <- 19
		temp_trap$weekSampled[temp_trap$weekSampled < 20] <- 19
		w_count$weekSampled[w_count$weekSampled == 31] <- 30
		temp_trap$weekSampled[temp_trap$weekSampled == 31] <- 30
		w_count$weekSampled[w_count$weekSampled > 30] <- 31
		temp_trap$weekSampled[temp_trap$weekSampled > 30] <- 31
	} else if (y == 2019){
		w_count$weekSampled[w_count$weekSampled < 20] <- 19
		temp_trap$weekSampled[temp_trap$weekSampled < 20] <- 19
		w_count$weekSampled[w_count$weekSampled > 31] <- 32
		temp_trap$weekSampled[temp_trap$weekSampled > 31] <- 32
	}
	
	w_count <- w_count %>% group_by(weekSampled) %>% 
		summarize(count = sum(count), t_trap = sum(t_trap), AD_count = sum(AD_count), AI_count = sum(AI_count))
	
	
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

	inits <- list(
		.RNG.name = "base::Mersenne-Twister",
		.RNG.seed = 7,
		pi_AD = t(t(all_data$AD_PBT) / colSums(all_data$AD_PBT)),
		pi_AI = t(t(all_data$AI_PBT) / colSums(all_data$AI_PBT))
	)

	iterations <- 20000
	chains <- 1
	
	model <- jags.model(file = modFile, data=all_data, n.chains = chains, inits = inits)
	samples <- coda.samples(model, c("grand_tot_AD_groups", "grand_tot_AI_groups", "grand_tot_AD", "grand_tot_AI"), n.iter=iterations)
	
	PBT_estimates[[as.character(y)]] <- samples
}

save(PBT_estimates, file = "pbt_only_results.rda")
