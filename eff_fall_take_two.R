# TAKE TWO: detection efficiency and fallback-reascension exploration

library(tidyverse)
library(lubridate)

# all PIT tag fish in analysis
d <- read_csv("data/raw_data.csv", guess_max = 1e6)
# 2013 is first year both PBT and PIT can be used for all release groups
d <- d %>% rename(RAL_RTR = `RAL/RTR`, Expansion_Ref = `Expansion Reference`)

# filtering and pulling out return year
aPit <- d %>% filter(!is.na(GraLastDate)) %>% 
	select(MY, TagID, RelName, Expansion_Ref, RelSite, GraLastDate) %>% 
	mutate(GraLastDate = mdy(GraLastDate), rYear = year(GraLastDate))

# PIT tag data at granite - all detections by time
pd <- read_csv("data/LGR_2013_2019_CHNK.csv")
colnames(pd) <- gsub(" ", "_", colnames(pd))
pd <- pd %>% select(-Obs_Count, -Unique_Tags)

# select only PITs in analysis... may not be necessary...
pd <- pd %>% filter(Tag_Code %in% d$TagID)

# group antennae
antGroups <- tibble(antenna = c("01",	"02",	"03",	"04",	"05",	"06",	"07",	"08",
										  "12",	"14",	"16",	"18",	"22",	"24",	"26",	"28",
										  "A1",	"A2",	"B1",	"B2",	"B3",	"B4"),
						  group = c("postTrap",	"postTrap",	"postTrap",	"postTrap",	
						  			 "postTrap",	"postTrap",	"postTrap",	"postTrap",	
						  			 "Trap",	"Trap",	"Trap",	"Trap",	"Trap",	
						  			 "Trap",	"Trap",	"Trap",	"Exit",	"Exit",	
						  			 "Window",	"Window",	"Window",	"Window"))
pd <- pd %>% mutate(antG = antGroups$group[match(Antenna_ID, antGroups$antenna)]) %>%
	mutate(obsDT = mdy_hms(Obs_Time_Value))

cpd <- pd %>% arrange(year(obsDT), Tag_Code, obsDT)

# rough detection efficiency
# of post-trap
# of fish detected at trap, how many were detected at postTrap on the same day
cpd_2013 <- cpd %>% filter(year(obsDT) == 2013)
days <- sort(unique(date(cpd_2013$obsDT)))
trap_tags <- cpd_2013 %>% filter(antG == "Trap")
postTrap_tags <- cpd_2013 %>% filter(antG == "postTrap")
detect <- 0
notDetect <- 0
for(d in days){
	# fish detected on that day
	tags_t <- trap_tags %>% filter(date(obsDT) == d) %>% pull(Tag_Code) %>% unique
	postTags_t <- postTrap_tags %>% filter(date(obsDT) == d) %>% pull(Tag_Code) %>% unique
	detect <- detect + sum(tags_t %in% postTags_t)
	notDetect <- notDetect + sum(!tags_t %in% postTags_t)
}
detect
notDetect
detect / (detect + notDetect)
# 89 %

# of trap
# of fish detected post-trap, how many were detected at trap on the same day
detect <- 0
notDetect <- 0
for(d in days){
	# fish detected on that day
	tags_t <- trap_tags %>% filter(date(obsDT) == d) %>% pull(Tag_Code) %>% unique
	postTags_t <- postTrap_tags %>% filter(date(obsDT) == d) %>% pull(Tag_Code) %>% unique
	detect <- detect + sum(postTags_t %in% tags_t)
	notDetect <- notDetect + sum(!postTags_t %in% tags_t)
}
detect
notDetect
detect / (detect + notDetect)
# 78.6 %
# combined
1 - (.11 * .214)
# .976

# known reascensions
r <- rep(FALSE, nrow(cpd))
for(i in 2:nrow(cpd)){
	if(cpd$Tag_Code[i] == cpd$Tag_Code[i-1] &&
		(cpd$antG[i] %in% c("Trap", "Window") && cpd$antG[i-1] %in% c("postTrap", "Exit"))) r[i] <- TRUE
}
table(r)
rdt <- rep(as.duration(seconds(1)), sum(r))
for(i in 1:sum(r)){
	rdt[i] <- int_length(interval(cpd$obsDT[which(r)[i] - 1], cpd$obsDT[which(r)[i]]))
}
summary(rdt)
quantile(rdt, seq(0,1,.05))
length(rdt)
hist(log(rdt), breaks = 50)
hist(sort(rdt)[1:300])
# can probably model time to reascension as a lognormal

cpd %>% group_by(Tag_Code) %>% summarise(nDays = n_distinct(date(obsDT))) %>%
	count(nDays)
# nDays     n
# <int> <int>
# 	1     1  5357
# 2     2  1112
# 3     3   114
# 4     4    33
# 5     5    15
# 6     6     5
# 7     7     2
# 8     8     1
# 9     9     2
# 10    12     1
# 11    14     1
# 12    19     1
# 13    26     1

# if detectors have high efficiency, then there should be more than 454
# known reascensions, or it's not uncommon for fish to spend more than a day 
# in the ladder

# Two related variables
# 1. estimate detection efficiency
# 2. Estimate reascension rate using detection efficiency and known reascensions
# 
# data
# 1. known reascensions and time period
# 2. known detections (for detection efficiency) and time difference
# 
# Can we express each in terms of the other?
# and then EM or MCMC to find some values?


######
## ok, taking simpler approach
######

# first, determine minimum fallback and reascension time
#   using all years combined

summary(rdt)
min_time <- dseconds(min(rdt))
min_time <- dseconds(quantile(rdt, .05))
# min_time <- dhours(100)

# now, calculate detection probability using fish that were detected at
# two surrounding detectors in less time than the minimum reasecension
# looks like window and exit were only added in 2016, so not available for
# 2013 - 2015

### first, detection prob for the trap
deData <- cpd %>% filter(year(obsDT) == 2016) %>% select(Tag_Code, antG, obsDT) %>%
	arrange(Tag_Code, obsDT)

# find transitions within time limit, then determine if fish was detected in between
fish_at_window <- deData %>% filter(antG == "Window") %>% pull(Tag_Code) %>% unique
detect <- 0
notDetect <- 0

for(f in fish_at_window){
	# only data for one fish
	td <- deData %>% filter(Tag_Code == f) %>% arrange(obsDT)
	if(nrow(td) < 2) next
	trap_detect <- td %>% filter(antG == "Trap")
	# get only window detections followed by a non-window detection
	r <- td$antG[1:(nrow(td) - 1)] == "Window" & td$antG[2:nrow(td)] != "Window"
	if(!any(r)) next
	# for each one, determine if detected at post trap or exit within minimum time
	for(i in which(r)){
		timespan <- interval(td$obsDT[i], td$obsDT[i] + min_time)
		end <- NA
		for(j in (i+1):nrow(td)){
			if(!td$obsDT[j] %within% timespan) break # if outside of time span, stop
			if(td$antG[j] == "Window") break # avoid potential for double counting
			if(td$antG[j] %in% c("postTrap", "Exit")) {
				end <- j
				break
			}
		}
		if(is.na(end)) next # no valid detection within timespan
		# determine if detected at trap between i and end
		if(any(trap_detect$obsDT %within% interval(td$obsDT[i], td$obsDT[end]))){
			detect <- detect + 1
		} else {
			notDetect <- notDetect + 1
		}
	}
}

detect / (detect + notDetect)
# apart from 2016, trap seems to have high detection probability


# now, for the post-trap
deData <- cpd %>% filter(year(obsDT) == 2019) %>% select(Tag_Code, antG, obsDT) %>%
	arrange(Tag_Code, obsDT)

# find transitions within time limit, then determine if fish was detected in between
fish_at_exit <- deData %>% filter(antG == "Exit") %>% pull(Tag_Code) %>% unique
detect <- 0
notDetect <- 0

for(f in fish_at_exit){
	# only data for one fish
	td <- deData %>% filter(Tag_Code == f) %>% arrange(obsDT)
	if(nrow(td) < 2) next
	posttrap_detect <- td %>% filter(antG == "postTrap")
	# get only window detections followed by a non-window detection
	r <- td$antG[1:(nrow(td) - 1)] %in% c("Window", "Trap") & !td$antG[2:nrow(td)] %in% c("Window", "Trap")
	if(!any(r)) next
	# for each one, determine if detected at post trap or exit within minimum time
	for(i in which(r)){
		timespan <- interval(td$obsDT[i], td$obsDT[i] + min_time)
		end <- NA
		for(j in (i+1):nrow(td)){
			if(!td$obsDT[j] %within% timespan) break # if outside of time span, stop
			if(td$antG[j] %in% c("Window", "Trap")) break # avoid potential for double counting
			if(td$antG[j] == "Exit") {
				end <- j
				break
			}
		}
		if(is.na(end)) next # no valid detection within timespan
		# determine if detected at trap between i and end
		if(any(posttrap_detect$obsDT %within% interval(td$obsDT[i], td$obsDT[end]))){
			detect <- detect + 1
		} else {
			notDetect <- notDetect + 1
		}
	}
}

detect / (detect + notDetect)

temp <- cpd %>% filter(year(obsDT) >= 2016) %>% select(Tag_Code, antG, obsDT) %>%
	arrange(Tag_Code, obsDT)
tBool <- temp$antG[1:(nrow(temp) - 1)] %in% c("Window", "Trap") & temp$antG[2:nrow(temp)] == "Exit"
table(tBool)

temp[which(tBool),]
temp[which(tBool) + 1,]
# different fish
trapFish <- cpd %>% filter(antG == "Trap") %>% pull(Tag_Code) %>% unique
postTrapFish <- cpd %>% filter(antG == "postTrap") %>% pull(Tag_Code) %>% unique
table(trapFish %in% postTrapFish)
32 / 5772
# so yes, postTrap does appear to have very high probability of detecting fish
#   from 2016-2019
# maybe don't need to worry so much about detection efficiency
