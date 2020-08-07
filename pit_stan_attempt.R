# PIT with stan

library(tidyverse)

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

dy <- aPit %>% filter(rYear == y) %>% count(Expansion_Ref) %>% 
	left_join(tagRates, by = "Expansion_Ref")

data <- list(obsPIT = dy$n, t = dy$tagRate, G = nrow(dy))


data <- list(N = nrow(dy),
				 obs = dy$n,
				 tag = dy$tagRate)

fit <- stan(file = "pit.stan", data = data)

# so, Stan cannot handle discrete parameters
# and there is really no way to reparameterize this to marginalize out
# the very parameter we are trying to infer!
