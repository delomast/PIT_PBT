# combined models
###########
## this file is current active attempt: stratifying PIT estimates and performing all inference on proportions
###########

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


# it will drastically simplify programming if PBT release groups are either all AI or all AD
# so just splitting each group
add_suffix <- rep("_ad", nrow(trap))
add_suffix[grepl("LGRU", trap$Pedigree_name)] <- "_ai"
trap$relGroup[!is.na(trap$relGroup)] <- paste0(trap$relGroup[!is.na(trap$relGroup)], add_suffix[!is.na(trap$relGroup)])

tag_rates <- bind_rows(tag_rates %>% mutate(PBT_RELEASE_GROUP = paste0(PBT_RELEASE_GROUP, "_ad")),
						 tag_rates %>% mutate(PBT_RELEASE_GROUP = paste0(PBT_RELEASE_GROUP, "_ai"))
						 )
pbt_pools <- bind_rows(pbt_pools %>% mutate(group = paste0(group, "_ad")),
						 pbt_pools %>% mutate(group = paste0(group, "_ai"))
						 )
# table(sort(unique(trap$relGroup)) %in% tag_rates$PBT_RELEASE_GROUP, useNA = "if")

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

## I think this model could be written in stan
## because parameters are all proportions
## may have to try it
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
		}
		n_PIT_groups <- length(pit_tag) - 1
	}

    model{

   	
   	for (w in 1:nWeeks){
	   	# prior
	   	# this is the prior for the top level pool composition for the week
	   	pi_pools[1:(n_pools + 1),w] <- ddirich(pool_priors[1:(n_pools + 1),w]) # +1 pool is fish not in a relevant pool
	    
	    
	   	for(p in 1:n_pools){
	   		# these are priors splitting each pool into pbt groups
	   		split_pools_pbt[p,w,1:(n_AD_groups + n_AI_groups + 1)] ~ ddirich()
	   		# these are priors splitting each pool into pit groups
	   		split_pools_pit[p,w,1:(n_PIT_groups + 1)] ~ ddirich()
	   	}
   	
   	
   		# likelihood
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
   		
   		# proportion of each PIT group
   		true_PIT[1:(n_PIT_groups + 1),w] ~ dmulti(pi_PIT[,w], n_PIT_scanned[w])
   		for (g in 1:n_PIT_groups){
   			# pit_detect[g,w] ~ dbin(pit_tag[g] - tag_loss[g], true_PIT[g,w])
   			pit_detect[g,w] ~ dbin(pit_tag[g], true_PIT[g,w])
   		}
   		
   		# link of pi_PIT, pi_AD, and pi_AI all back to pi_pools
   		for (g in 1:n_AD_groups){
   			pi_AD[g,w] <- pi_pools[,w] %*% split_pools_pbt[,w,g]
   		}
   		# need to think about if this is really correct
   		pi_AD[(n_AD_groups + 1),w] <- 1 - sum(pi_AD[1:n_AD_groups,w])
   		
   		##
   		# Need to split unassigned into unassigned AD and AI
   		# also, split_pools_pbt is all pbt groups, ad and ai, need to reconcile
   		##
   		
   		for (g in (n_AD_groups + 1):(n_AD_groups + n_AI_groups)){
   			# so split pools is AD followed by AI and then unassigned (for pit only groups)
   			pi_AI[g - n_AD_groups,w] <- pi_pools[,w] %*% split_pools_pbt[,w,g]
   		}
   		# need to think about if this is really correct
   		pi_AI[(n_AI_groups + 1),w] <- 1 - sum(pi_AI[1:n_AI_groups,w])
   		
   		for (g in 1:n_PIT_groups){
   			pi_PIT[g,w] <- pi_pools[,w] %*% split_pools_pit[,w,g]
   		}
   		# need to think about if this is really correct
   		pi_AI[(n_PIT_groups + 1),w] <- 1 - sum(pi_PIT[1:n_PIT_groups,w])
   	}
   	
    }

    ",
	 file = modFile)


PBT_estimates <- list()
# y <- 2016
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
	if(y == 2016){
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

	# PIT data
	template <- aPit %>% filter(rYear == y) %>% pull(Expansion_Ref) %>% unique %>% 
		expand_grid(Expansion_Ref = ., weekSampled = unique(w_count$weekSampled))
	dy <- aPit %>% filter(rYear == y) %>% group_by(weekSampled) %>% count(Expansion_Ref) %>% 
		right_join(template, by = c("Expansion_Ref", "weekSampled")) %>%
		arrange(Expansion_Ref, weekSampled) %>% mutate(n = replace_na(n,0)) %>% spread(weekSampled, n)
	temp <- slice(dy,1) %>% mutate(Expansion_Ref = "Unassigned") %>% as.data.frame
	temp[,2:ncol(temp)] <- w_count$count - colSums(dy[,2:ncol(dy)])
	if(any(temp[,2:ncol(temp)] < 0)) stop("Error in pit tag window count compilation")
	dy <- bind_rows(dy, temp)
	dy_tag <- dy %>% select(Expansion_Ref) %>% left_join(PIT_tagRates, by = "Expansion_Ref")
	if(any(dy_tag$Expansion_Ref != dy$Expansion_Ref)) stop("error pit rownames in ", y)
	
	
	# pools
	temp_pbt_pools <- pbt_pools %>% filter(group %in% c(AD_tag$PBT_RELEASE_GROUP, AI_tag$PBT_RELEASE_GROUP)) %>% 
		pull(pool) %>% unique
	temp_pit_pools <- pit_pools %>% filter(group %in% dy$Expansion_Ref) %>% pull(pool) %>%
		unique
	union_pools <- sort(union(temp_pbt_pools, temp_pit_pools))
	
	# each column is a 1-column matrix with a 1 meaning that (AD/AI/pit) group is in the pool
	# it is assumed that either a pool is only represented by tagged groups, or is only represented by untagged groups
	# within each data set (PBT or pit)
	AD_col_vectors <- matrix(0, nrow = nrow(AD_PBT_counts), ncol = length(union_pools))
	AI_col_vectors <- matrix(0, nrow = nrow(AI_PBT_counts), ncol = length(union_pools))
	PIT_col_vectors <- matrix(0, nrow = nrow(dy), ncol = length(union_pools))
	for(i in 1:length(union_pools)){
		t_PBT <- pbt_pools %>% filter(pool == union_pools[i]) %>% pull(group)
		AD_col_vectors[which(AD_PBT_counts$relGroup %in% t_PBT),i] <- 1
		AI_col_vectors[which(AI_PBT_counts$relGroup %in% t_PBT),i] <- 1
		t_PIT <- pit_pools %>% filter(pool == union_pools[i]) %>% pull(group)
		PIT_col_vectors[which(dy$Expansion_Ref %in% t_PIT),i] <- 1
	}
	colSums(AD_col_vectors)
	colSums(AI_col_vectors)
	colSums(PIT_col_vectors)
	
	# last row should be zero, all others should be 1
	rowSums(AD_col_vectors)
	rowSums(AI_col_vectors)
	rowSums(PIT_col_vectors)
	
	# now need priors for the top level pools
	pool_priors <- matrix(0, nrow = length(union_pools), ncol = nrow(w_count))
	for(i in 1:length(union_pools)){
		pool_priors[i,] <- ifelse(AD_col_vectors[,i] %*% as.matrix(AD_PBT_counts[2:ncol(AD_PBT_counts)]) +
				AI_col_vectors[,i] %*% as.matrix(AI_PBT_counts[2:ncol(AI_PBT_counts)]) +
				PIT_col_vectors[,i] %*% as.matrix(dy[2:ncol(dy)]) > 0, 1, 0)
	}
	pool_priors <- rbind(pool_priors, rep(1, nrow(pool_priors))) # adding fish not in a relevant pool
	pool_priors <- t(t(pool_priors) / colSums(pool_priors)) # normalizing so each prior has 1 psuedocount
	
	# now need priors splitting each pool into PBT groups IN EACH STRATA
	# pool, week, group(AD,AI, unassigned)
	# -1 is b/c boht AD and AI have unassigned rows
	split_pools_pbt <- array(0, dim = c(length(union_pools), nrow(w_count), nrow(AD_PBT_counts) + nrow(AI_PBT_counts) - 1))
	# first put 1's in each cell with observed tagged fish
	# I'm sure there is a quicker vectorized way to do this, but it is hard enough to follow as it is
	for(i in 1:length(union_pools)){
		temp_p <- pbt_pools %>% filter(pool == union_pools[i]) %>% pull(group)
		for(w in 1:nrow(w_count)){
			for(g in 1:(nrow(AD_PBT_counts)-1)){
				if(!AD_PBT_counts$relGroup[g] %in% temp_p) next
				split_pools_pbt[i,w,g] <- as.numeric(AD_PBT_counts[[w + 1]][g] > 0)
			}
			for(g in 1:(nrow(AI_PBT_counts)-1)){
				if(!AI_PBT_counts$relGroup[g] %in% temp_p) next
				split_pools_pbt[i,w,g + nrow(AD_PBT_counts) - 1] <- as.numeric(AI_PBT_counts[[w + 1]][g] > 0)
			}
		}
	}
	# now normalize split_pools_pbt so psuedocounts sum to 1
	# and for those that were not observed, put 1 for the "Unassigned" group
	
	
	
	
	
	
	
	
	# turn into list to pass to JAGS
	all_data <- list(w_count = w_count$count,
						  AD_PBT = as.matrix(AD_PBT_counts[2:ncol(AD_PBT_counts)]),
						  AD_tag = AD_tag$PBT_RELEASE_GROUP_TAGRATE,
						  AI_PBT = as.matrix(AI_PBT_counts[2:ncol(AI_PBT_counts)]),
						  AI_tag = AI_tag$PBT_RELEASE_GROUP_TAGRATE,
						  n_AD = w_count$AD_count,
						  n_AI = w_count$AI_count,
						  AD_col_vectors = AD_col_vectors,
						  AI_col_vectors = AI_col_vectors,
						  PIT_col_vectors = PIT_col_vectors,
						  pit_detect = as.matrix(dy[2:ncol(dy)]),
						  pit_tag = dy_tag$tagRate,
						  n_pools = length(union_pools),
						  pool_priors = pool_priors
						  )
	
	
	
	
	
	
	
	
	
	pool_abund <- rep(0, length(intersect_pools))
	for(p in 1:length(intersect_pools)){
		temp <- dy$n[PIT_col_vectors[,p] > 0]
		temp_tag <- dy$tagRate[PIT_col_vectors[,p] > 0]
		pool_abund[p] <- sum(temp / temp_tag)
	}
	
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
	# all_data$obsPIT_both <- ceiling(all_data$obsPIT_both / 10000)
	model <- jags.model(file = modFile, data=all_data, n.chains = chains, inits = inits)
	samples <- coda.samples(model, c("grand_tot_AD_groups", "grand_tot_AI_groups", "grand_tot_AD", "grand_tot_AI"), n.iter=iterations)
	
	comb_estimates[[as.character(y)]] <- samples
}

save(comb_estimates, file = "combined_results.rda")

obsPIT_both / pit_tag_both
