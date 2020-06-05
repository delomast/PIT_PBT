# detection efficiency and fallback-reascension exploration
######
## most of this is junk, computation times to long to be useful if
##     implemented properly
######

library(tidyverse)
library(lubridate)
library(rjags)

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

# collapse consecutive detections at the same antennae group
cpd <- pd %>% arrange(year(obsDT), Tag_Code, obsDT)

rows_to_remove <- rep(FALSE, nrow(cpd))
for(i in 2:nrow(cpd)){
	if(cpd$Tag_Code[i] == cpd$Tag_Code[i-1] &&
		cpd$antG[i] == cpd$antG[i-1]) rows_to_remove[i] <- TRUE
}
table(rows_to_remove)
hist(as.duration(cpd$obsDT[rows_to_remove] - cpd$obsDT[which(rows_to_remove) - 1]))
diffSec <- as.duration(cpd$obsDT[rows_to_remove] - cpd$obsDT[which(rows_to_remove) - 1])
cpd %>% filter(Tag_Code == cpd$Tag_Code[rows_to_remove][which(as.numeric(diffSec) == 3017466)]) %>% 
	as.data.frame()
hist(diffSec[as.numeric(diffSec) < 3600])
hist(diffSec[as.numeric(diffSec) < 100])
summary(diffSec)
quantile(diffSec, seq(.01,.99, .04))
sum(diffSec > minutes(5))

# movement back through trap should be impossible?
r <- rep(FALSE, nrow(cpd))
for(i in 2:nrow(cpd)){
	if(cpd$Tag_Code[i] == cpd$Tag_Code[i-1] &&
		(cpd$antG[i] %in% c("Trap", "Window") & cpd$antG[i-1] %in% c("postTrap", "Exit"))) r[i] <- TRUE
}
table(r)
rdt <- as.duration(cpd$obsDT[r] - cpd$obsDT[which(r) - 1])
summary(rdt)
# yes, looks like it hasn't happened
sum(rdt < hours(1))
rdt[rdt < hours(1)]

# this is proving to be a bit of a mess (ok maybe not all that bad) trying to 
#  collapse heuristically
# maybe better approach is to model movements and detection efficiency together?
# don't forget that time spent in the trap will show up in the movements...

# rough detection efficiency
# of post-trap
# of fish detected at trap, how many were detected at postTrap on the same day
cpd_2013 <- cpd %>% filter(year(obsDT) == 2013)
days <- sort(unique(date(cpd_2013$obsDT)))
trap_tags <- cpd_2013 %>% filter(antG == "Trap")
postTrap_tags <- cpd_2013 %>% filter(antG == "postTrap")
detect <- 0
notDetect <- 0
temp <- c()
for(d in days){
	# fish detected on that day
	tags_t <- trap_tags %>% filter(date(obsDT) == d) %>% pull(Tag_Code) %>% unique
	postTags_t <- postTrap_tags %>% filter(date(obsDT) == d) %>% pull(Tag_Code) %>% unique
	detect <- detect + sum(tags_t %in% postTags_t)
	notDetect <- notDetect + sum(!tags_t %in% postTags_t)
	temp <- c(temp, tags_t[!tags_t %in% postTags_t])
}
detect
notDetect
detect / (detect + notDetect)
# 89 %
cpd_2013 %>% filter(Tag_Code %in% temp) %>% select(Tag_Code, antG, obsDT) %>%
	mutate(d = date(obsDT), h = hour(obsDT)) %>% select(-obsDT) %>% distinct %>%
	as.data.frame

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

# model
modFile <- "PIT_models/detections_and_movements.txt"
cat("
	data{
	   # dimension parameters
   	dataDims <- dim(d) 
   	F <- dataDims[1]# number of fish
   	T <- dataDims[2] # number of time steps
	}
    model{
   	# all fish start below dam
   	for(f in 1:F){
   			z[f,1] <- 1
   	}
    
   	# likelihood
   	for (t in 2:T){
   		for(f in 1:F){
   			# state of fish f
   			z[f,t] ~ dcat(p[z[f,t-1],1:9])
   			
   			# detections of fish f at each state
   			for(s in 1:4){
   				d[f,t,s] ~ dbern(dp[s] * ((z[f,t]/2) == s) + ((z[f,t]/2) != s) * .000001)
   			}
   		}
   	}
   	
   	# priors
   	# detection probs
   	for(s in 1:4){
   		dp[s] ~ dbeta(1,1)
   	}
   	
   	# transition probs
   	# rows are current pos, columns are new pos
   	p[1,1:2] ~ ddirich(c(1,1))
   	p[1,3:9] <- rep(.000001,7)
   	p[2,1:3] ~ ddirich(c(1,1,1))
   	p[2,4:9] <- rep(.000001,6)
   	p[3,1] <- .000001
   	p[3,2:4] ~ ddirich(c(1,1,1))
   	p[3,5:9] <- rep(.000001,5)
   	p[4,1:3] <- rep(.000001,3)
   	p[4,4:5] ~ ddirich(c(1,1))
   	p[4,6:9] <- rep(.000001,4)
   	p[5,1:4] <- rep(.000001,4)
   	p[5,5:6] ~ ddirich(c(1,1))
   	p[5,7:9] <- rep(.000001,3)
   	p[6,1:4] <- rep(.000001,4)
   	p[6,5:7] ~ ddirich(c(1,1,1))
   	p[6,8:9] <- rep(.000001,2)
   	p[7,1:5] <- rep(.000001,5)
   	p[7,6:8] ~ ddirich(c(1,1,1))
   	p[7,9] <- .000001
   	p[8,1:6] <- rep(.000001,6)
   	p[8,7:9] ~ ddirich(c(1,1,1))
   	p[9,2:8] <- rep(.000001,7)
   	p[9,1] ~ dbeta(1,1)
   	p[9,9] <- 1.0 - p[9,1]

    }

    ",
	 file = modFile)

iterations <- 10
chains <- 1

# input data is three dimensional array
# dims: fish, time step, detector
data <- cpd %>% filter(year(obsDT) == 2013)
tstepsize <- seconds(hours(12)) # using hours just to see what computation times might be like
start <- min(data$obsDT) - (tstepsize * 10)
end <- max(data$obsDT) + (tstepsize * 10)
steps <- seq(start, end, tstepsize)
# steps_int <- list()
steps_int <- rep(interval(min(data$obsDT), max(data$obsDT)), length(steps) -1)
for(s in 2:length(steps)){
	# steps_int[[s - 1]] <- interval(steps[s-1], steps[s])
	steps_int[s - 1] <- interval(steps[s-1], steps[s] - seconds(0.1))
}
fish <- unique(data$Tag_Code)
# this isn't quite right, but will show comp times
dataJags <- array(0, dim = c(length(fish), length(steps_int), 4))
for(f in 1:length(fish)){
	dt <- data %>% filter(Tag_Code == fish[f]) %>% arrange(obsDT)
	
	for(i in 1:nrow(dt)){
		pos <- which(dt$obsDT[i] %within% steps_int)
		if(dataJags[f, pos, which(c("Window", "Trap", "postTrap", "Exit") == dt$antG[i])] == 1) next
		while(1 %in% dataJags[f, pos, ]) pos <- pos + 1
		dataJags[f, pos, which(c("Window", "Trap", "postTrap", "Exit") == dt$antG[i])] <- 1
	}
}



subset_data <- dataJags[1:3,,]
# fit model
system.time(
model <- jags.model(file = modFile, data=list(d = subset_data), n.chains = chains,
						  n.adapt = 2)
)
system.time(
samples <- coda.samples(model, c("p[9,1]", "dp"), n.iter=100)
)

# too long
