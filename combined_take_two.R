# combined models
# let's first make it work for ONE stratum

library(tidyverse)
library(lubridate)
library(rjags)

if(!dir.exists("combined_models")) dir.create("combined_models")

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

# pools
pit_pools <- tibble()
pbt_pools <- tibble()
for(i in 1:nrow(pools)){
	pit_pools <- pit_pools %>% bind_rows(tibble(group = strsplit(pools$PIT_groups[i], ",")[[1]],
					  pool = pools$Pool[i]))
	pbt_pools <- pbt_pools %>% bind_rows(tibble(group = strsplit(pools$PBT_groups[i], ",")[[1]],
					  pool = pools$Pool[i]))
}

# PIT data
d <- read_csv("data/raw_data.csv", guess_max = 1e6) %>% rename(RAL_RTR = `RAL/RTR`, Expansion_Ref = `Expansion Reference`)
e <- read_csv("data/expansions.csv")
colnames(e) <- gsub(" ", "_", colnames(e))

# filtering and pulling out return year
aPit <- d %>% filter(!is.na(GraLastDate)) %>% 
	select(MY, TagID, RelName, Expansion_Ref, RelSite, GraLastDate) %>% 
	mutate(GraLastDate = mdy(GraLastDate), rYear = year(GraLastDate), weekSampled = week(GraLastDate)) %>% 
	mutate(tagRate = 1 / e$RAL_Expansion[match(Expansion_Ref, e$Expansion_Reference)])

PIT_tagRates <- aPit %>% select(Expansion_Ref, tagRate) %>% distinct



######################################################################
######################################################################


modFile <- "combined_models/model_no_adjust_wc.txt"
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
		n_pools <- length(indexes)
   	for(p in 1:length(indexes)){
			alpha_p[p,1:3] <- ifelse(obsPIT_both[p,] > 0, .1, 1e-20)
		}

	}

    model{
    
   	# likelihood
   	
   	# PBT
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
   	
   	# PIT
   	for (p in 1:n_pools){ # for each pool
   		for(s in 1:indexes[p]){ # for each PIT group in the pool
   			# obsPIT_both[p,s] ~ dbin(pit_tag_both[p,s] * (1 - tag_loss[p,s]), round(pool_abundance[p] * pit_pi[p,s]))
   			obsPIT_both[p,s] ~ dbin(pit_tag_both[p,s] * .9999, round(pool_abundance[p] * pit_pi[p,s]))
   		}
   	}
   	
   	# mult by window counts and add up
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
   	
   	# abundance of pools
   	for(p in 1:n_pools){
   		pool_AD_total[p] <- grand_tot_AD_groups %*% AD_col_vectors[,p]
   		pool_AI_total[p] <- grand_tot_AI_groups %*% AI_col_vectors[,p]
   		pool_abundance[p] ~ sum(pool_AD_total[p], pool_AI_total[p])
   	}
   	
   	# priors
   	# PBT
   	for (w in 1:nWeeks){
   		p_clip[w] ~ dbeta(.01,.01)
   		pi_AD[1:(n_AD_groups + 1),w] ~ ddirich(alpha_AD[,w])
   		pi_AI[1:(n_AI_groups + 1),w] ~ ddirich(alpha_AI[,w])
   	}
   	# PIT
   	for (p in 1:n_pools){
   		pit_pi[p,1:3] ~ ddirich(alpha_p[p,])
   	}

    }

    ",
	 file = modFile)


PBT_estimates <- list()
# y <- 2019
for(y in 2016:2019){
	
	# just make it work for one year and one stratum
	y = 2016
	s_test = 20
	
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

	w_count <- w_count %>% filter(weekSampled == s_test)
	temp_trap <- temp_trap %>% filter(weekSampled == s_test)
	temp_pit <- aPit %>% filter(rYear == y, weekSampled == s_test)
	
	######
	# now, need, for each pool, the counts of the subgroups and the tagging rates of the subgroups
	######
	AD_pbt_counts <- temp_trap %>% filter(grepl("OtsLGRA", Individual_name)) %>% count(relGroup)
	AI_pbt_counts <- temp_trap %>% filter(grepl("OtsLGRU", Individual_name)) %>% count(relGroup)
	
	temp_pbt_pools <- pbt_pools %>% filter(group %in% c(temp_trap$relGroup)) %>% 
		pull(pool) %>% unique
	temp_pit_pools <- pit_pools %>% filter(group %in% temp_pit$Expansion_Ref) %>% pull(pool) %>%
		unique
	all_pools <- union(temp_pbt_pools, temp_pit_pools)
	AD_pool_count <- array(0, dim = c(length(all_pools) + 1, nrow(AD_pbt_counts)))
	AI_pool_count <- array(0, dim = c(length(all_pools) + 1, nrow(AI_pbt_counts)))
	for(p in 1:length(all_pools)){
		t_pools <- pbt_pools %>% filter(pool == all_pools[p]) %>% pull(group)
		temp <- AD_pbt_counts$n
		temp[!AD_pbt_counts$relGroup %in% t_pools] <- 0
		AD_pool_count[p,] <- temp
		temp <- AI_pbt_counts$n
		temp[!AI_pbt_counts$relGroup %in% t_pools] <- 0
		AI_pool_count[p,] <- temp
	}
	# last row and last column is "Unassigned" not exactly sure how to use this right now
	AD_pool_count[length(all_pools) + 1,nrow(AD_pbt_counts)] <- AD_pbt_counts$n[AD_pbt_counts$relGroup == "Unassigned"]
	AI_pool_count[length(all_pools) + 1,nrow(AI_pbt_counts)] <- AI_pbt_counts$n[AI_pbt_counts$relGroup == "Unassigned"]

	
		
	inits <- list(
		.RNG.name = "base::Mersenne-Twister",
		.RNG.seed = 7,
		pi_AD = t(t(all_data$AD_PBT) / colSums(all_data$AD_PBT)),
		pi_AI = t(t(all_data$AI_PBT) / colSums(all_data$AI_PBT)),
		pit_pi = (obsPIT_both / pit_tag_both) / rowSums(obsPIT_both / pit_tag_both, na.rm = TRUE)
	)
	inits$pit_pi[is.na(inits$pit_pi)] <- 0
	iterations <- 20000
	chains <- 1	
	
	model <- jags.model(file = modFile, data=all_data, n.chains = chains, inits = inits)
	samples <- coda.samples(model, c("grand_tot_AD_groups", "grand_tot_AI_groups", "grand_tot_AD", "grand_tot_AI"), n.iter=iterations)
	
	comb_estimates[[as.character(y)]] <- samples
}

save(comb_estimates, file = "combined_results.rda")

obsPIT_both / pit_tag_both
