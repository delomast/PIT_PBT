library(tidyverse)
library(lubridate)

# PIT tag data summarizing only unique detections crossing Granite

d <- read_csv("data/raw_data.csv", guess_max = 1e6)
d <- d %>% rename(RAL_RTR = `RAL/RTR`, Expansion_Ref = `Expansion Reference`)

# make sure no replicate tags
length(unique(d$TagID))
nrow(d)

d %>% pull(GraLastDate)
d %>% filter(is.na(GraLastDate)) %>% as.data.frame %>% head
d %>% count(is.na(GraLastDate))

# age breakdown of detections
d %>% filter(!is.na(GraLastDate)) %>% mutate(diff = mdy(GraLastDate) - mdy(RelDate)) %>%
	pull(diff) %>% as.numeric() %>% (function(x){x/365}) %>% hist()

# time of detections
d %>% filter(!is.na(GraLastDate)) %>% mutate(d = mdy(GraLastDate)) %>%
	filter(year(d) == 2013) %>%
	pull(d) %>% summary

d %>% filter(!is.na(GraLastDate)) %>% mutate(d = mdy(GraLastDate)) %>%
	pull(d) %>% year %>% summary


# PIT tag data at granite - all detections by time
pd <- read_csv("data/LGR_2013_2019_CHNK.csv")
colnames(pd) <- gsub(" ", "_", colnames(pd))
for(i in 1:ncol(pd)){
	u <- unique(pull(pd,i))
	if(length(u) < 20){
		print(colnames(pd)[i])
		print(u)
	}
}
pd <- pd %>% select(-Obs_Count, -Unique_Tags)
table(d %>% filter(!is.na(GraLastDate)) %>% pull(TagID) %in% pd$Tag_Code)

d %>% filter(!is.na(GraLastDate)) %>% mutate(b = TagID %in% pd$Tag_Code) %>% group_by(b) %>%
	count(year(mdy(GraLastDate))) %>% as.data.frame
# all the fish for the years we are interested in are present in the LGR file

pd <- pd %>% filter(Tag_Code %in% d$TagID)

# look at detections
ant <- sort(unique(pd$Antenna_ID))
head(pd$Tag_Code)
pd %>% filter(Tag_Code == "384.3B2394899C") %>% as.data.frame()
pd %>% filter(Tag_Code == "384.3B2394899C") %>% as.data.frame()

table(pd$Antenna_ID)

# window counts
wc <- read_csv("data/window_counts.csv") %>% mutate(OperationsDate = mdy(OperationsDate))

