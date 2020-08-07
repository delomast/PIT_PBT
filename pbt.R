# some PBT data sorting and exploration
library(tidyverse)
library(lubridate)

# window counts
wc <- read_csv("data/window_counts.csv") %>% mutate(OperationsDate = mdy(OperationsDate))

# trap data
trap <- read_tsv("data/lgru_lgra_data.txt", guess_max = 1e6) %>% 
	mutate(DateSampled = mdy(DateSampled))
colnames(trap) <- gsub(" ", "_", colnames(trap))

# ma info
maInfo <- read_tsv("data/Ots_ma_data.txt",guess_max = 1e6)
colnames(maInfo) <- gsub(" ", "_", colnames(maInfo))

# tag rates
tagRates <- maInfo %>% select(PBT_RELEASE_GROUP, PBT_RELEASE_GROUP_TAGRATE) %>%
	distinct() %>% filter(!is.na(PBT_RELEASE_GROUP))

# add Ma pedigree and release group
trap <- trap %>% mutate(maPed = maInfo$Pedigree_name[match(GenMa, maInfo$Individual_name)],
					 relGroup = maInfo$PBT_RELEASE_GROUP[match(GenMa, maInfo$Individual_name)],
					 relGroup = na_if(relGroup, "0"),
					 relStrat = maInfo$OffspringReleaseStrategy[match(GenMa, maInfo$Individual_name)])
trap %>% filter(is.na(maPed), !is.na(GenMa)) %>% pull(GenMa) %>% unique
# NA's are fall, lower Columbia, and ungenotyped given "NG": all can be treated as unassigned
trap$GenMa[is.na(trap$maPed)] <- NA
trap$GenPa[is.na(trap$maPed)] <- NA
trap %>% filter(!is.na(GenMa), is.na(relGroup)) %>% pull(maPed) %>% unique %>% sort
# hatchery codes to ignore - lower columbia
toIgnore <- c("^OtsSWWH[0-9]{2}S$", "^OtsKLFH[0-9]{2}S$")
for(s in toIgnore){
	trap$GenMa[grepl(s, trap$maPed)] <- NA
	trap$GenPa[grepl(s, trap$maPed)] <- NA
	trap$maPed[grepl(s, trap$maPed)] <- NA
}
trap %>% filter(!is.na(GenMa), is.na(relGroup)) %>% count(relStrat) %>% as.data.frame()



trap %>% filter(!is.na(GenMa), is.na(relGroup)) %>% pull(maPed) %>% unique %>% sort
# make list of mas with no release groups for Matt
# Matt only has trackign data for BY 2011 onward, so analysis limited to SY2016-19
trap <- trap %>% filter(year(DateSampled) >= 2016)
noReleaseGroup <- trap %>% filter(!is.na(GenMa), is.na(relGroup)) %>% select(maPed, GenMa, relStrat) %>% 
	distinct %>% arrange(maPed, GenMa) %>% rename(OffspringReleaseStrategy = relStrat)
write.table(noReleaseGroup, "Mas_no_release.txt", row.names = FALSE,
				col.names = TRUE, quote = FALSE, sep = "\t")
pas <- trap %>% filter(!is.na(GenMa), is.na(relGroup)) %>% pull(GenPa)
maInfo %>% filter(Individual_name %in% pas) %>% select(contains("PBT"))
# these two just have the pas and mas switched so just need to swap the mas and pas


# tag rates
tagRates <- maInfo %>% select(PBT_RELEASE_GROUP, PBT_RELEASE_GROUP_TAGRATE) %>%
	mutate(PBT_RELEASE_GROUP_TAGRATE = as.numeric(PBT_RELEASE_GROUP_TAGRATE)) %>%
	distinct() %>% filter(!is.na(PBT_RELEASE_GROUP), !is.na(PBT_RELEASE_GROUP_TAGRATE),
								 PBT_RELEASE_GROUP_TAGRATE > 0, 
								 # higher one is old
								 !(PBT_RELEASE_GROUP == "2016-KOOS-KOOS-ClearCr" & PBT_RELEASE_GROUP_TAGRATE == 0.9849))
n_distinct(tagRates$PBT_RELEASE_GROUP) == nrow(tagRates)
dups <- tagRates %>% arrange(PBT_RELEASE_GROUP) %>% count(PBT_RELEASE_GROUP) %>% filter(n > 1) %>%
	pull(PBT_RELEASE_GROUP)
tagRates %>% filter(PBT_RELEASE_GROUP %in% dups) %>% arrange(PBT_RELEASE_GROUP) %>% as.data.frame


# correspondence between PBT and PIT tag release groups


# model
# Maybe a walk could work for the numbers of each group, and the 
# lognormal prior work b/c it puts a good amount of weight close to 0?


