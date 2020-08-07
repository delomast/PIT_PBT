# PBT data sorting and exploration
library(tidyverse)
library(lubridate)
library(EFGLmh)

# window counts
wc <- read_csv("data/window_counts.csv") %>% mutate(OperationsDate = mdy(OperationsDate))

# from progeny with genotypes to exclude failed samples
genoData <- readInData("data/lgra_lgru_16_19.txt", guess_max = 1e6)
gs <- genoSuccess(genoData)
gs %>% filter(success < .90) %>% count(Pop)
# need to take into account 2019 at Ots343
genoData$metadata %>% group_by(Pop) %>% count(grepl("[Ff]ail", GenComments))
# 18 and 19 can use gencomments, 16 and 17 calculate using 95 loci
p95 <- genoData %>% removePops(pops = getPops(.)[getPops(.) != "OtsLGRA16S"]) %>%
	lociSuccess %>% filter(success > 0) %>% pull(locus)
gs <- genoSuccess(genoData, loci = p95)
failInds <- gs %>% filter(success < .90, grepl("1[67]S_", Ind)) %>% pull(Ind)
failInds <- c(failInds, genoData$metadata %>% filter(grepl("[Ff]ail", GenComments)) %>% pull(Ind))
genoData <- genoData %>% removeInds(failInds)
goodInds <- getInds(genoData)

# trap data - only 2016 onward
trap <- read_tsv("data/lgru_lgra_data.txt", guess_max = 1e6) %>% 
	mutate(DateSampled = mdy(DateSampled)) %>% filter(year(DateSampled) >= 2016)
colnames(trap) <- gsub(" ", "_", colnames(trap))
# add biosamples ID
trap <- trap %>% mutate(BiosamplesID = genoData$metadata$BiosamplesID[match(Individual_name, genoData$metadata$Ind)]) %>%
	filter(!is.na(BiosamplesID))

# trap data from LGTrapping - 2016-19
trap2 <- read_csv("data/LGTrapping_16_19.csv", guess_max = 1e6) %>%
	mutate(CollectionDate = mdy(CollectionDate))
trap2 <- trap2 %>% filter(BioSamplesID %in% trap$BiosamplesID)

# ma info
maInfo <- read_tsv("data/Ots_ma_data.txt",guess_max = 1e6)
colnames(maInfo) <- gsub(" ", "_", colnames(maInfo))
maInfo <- maInfo %>% mutate(PBT_RELEASE_GROUP = na_if(PBT_RELEASE_GROUP, "0")) %>%
	filter(!is.na(PBT_RELEASE_GROUP))

# swap ma and pa if needed
swap <- (!trap$GenMa %in% maInfo$Individual_name) & (trap$GenPa %in% maInfo$Individual_name)
sum(swap)
temp <- trap$GenMa[swap]
trap$GenMa[swap] <- trap$GenPa[swap]
trap$GenPa[swap] <- temp
rm(temp)

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

# look for fish with release group in trapping but not in progeny
norelgroup <- trap %>% filter(is.na(relGroup)) %>% pull(BiosamplesID)
trap2 <- trap2 %>% mutate(GenPBT_RGroup = na_if(GenPBT_RGroup, "#N/A"),
								  GenPBT_RGroup = na_if(GenPBT_RGroup, "NG"),
								  GenPBT_RGroup = na_if(GenPBT_RGroup, "Unassigned-H"),
								  GenPBT_RGroup = na_if(GenPBT_RGroup, "Unassigned-W")
								  )
relGrouptrap2 <- trap2 %>% filter(!is.na(GenPBT_RGroup)) %>% pull(BioSamplesID)

table(norelgroup %in% relGrouptrap2, useNA = "if")
trap %>% filter(BiosamplesID %in% norelgroup[norelgroup %in% relGrouptrap2]) %>%
	count(Pedigree_name)

trap2 %>% filter(BioSamplesID %in% norelgroup[norelgroup %in% relGrouptrap2]) %>%
	count(GenPBT_RGroup) %>% as.data.frame
# those can be ignored, they aren't relevant, except the odd format that will be fixed

# correspondence between PBT and PIT tag release groups
# tag rates
tagRates <- maInfo %>% select(PBT_RELEASE_GROUP, PBT_RELEASE_GROUP_TAGRATE) %>%
	distinct() %>% filter(!is.na(PBT_RELEASE_GROUP), PBT_RELEASE_GROUP != "UNK") %>% arrange(PBT_RELEASE_GROUP)
tagRates <- tagRates %>% filter(PBT_RELEASE_GROUP %in% trap$relGroup)
n_distinct(tagRates$PBT_RELEASE_GROUP) == nrow(tagRates) 
tagRates %>% filter(PBT_RELEASE_GROUP %in% (tagRates %>% filter(duplicated(PBT_RELEASE_GROUP)) %>% pull(PBT_RELEASE_GROUP)))
# will have to resolve duplicated groups later

tagRates <- tagRates %>% select(PBT_RELEASE_GROUP) %>% distinct

pit_pbt_corres <- read_csv("data/pit_pbt_codes.csv")
e <- read_csv("data/expansions.csv")
colnames(e) <- gsub(" ", "_", colnames(e))


# multiple PITcodes correspond to the same PBT code
pit_pbt_corres
e <- e %>% mutate(PBTcode = pit_pbt_corres$PBTcode[match(PTAGIS_RelSite, pit_pbt_corres$PITcode)])
unique(pit_pbt_corres$PBTcode)

# these are just the pbt groups that seem to fit a uniform style
tagRates_uniform <- tagRates %>% filter(grepl("^[0-9]{4}-", PBT_RELEASE_GROUP))
# those that don't
tagRates_odd <- tagRates %>% filter(!PBT_RELEASE_GROUP %in% tagRates_uniform$PBT_RELEASE_GROUP)
# now trying to automate matching the uniform ones
tagRates_uniform_split <- tibble()
for(i in 1:nrow(tagRates_uniform)){
	s <- strsplit(tagRates_uniform$PBT_RELEASE_GROUP[i], "-")[[1]]
	if(length(s) < 5) s <- c(s, rep(NA, 5 - length(s)))
	tagRates_uniform_split <- rbind(tagRates_uniform_split, s)
}
names(tagRates_uniform_split) <- c("y", "c1", "c2", "c3", "c4")
tagRates_uniform_split <- as_tibble(tagRates_uniform_split) %>% mutate(y = as.numeric(y))
table(tagRates_uniform_split$c3 %in% pit_pbt_corres$PBTcode, useNA = "if")
tagRates_uniform_split %>% filter(!tagRates_uniform_split$c3 %in% pit_pbt_corres$PBTcode) %>% as.data.frame()
tagRates_uniform <- tagRates_uniform %>% bind_cols(tagRates_uniform_split)

# filtering expansions to be just those observed in 16-19
d <- read_csv("data/raw_data.csv", guess_max = 1e6) %>% 
	rename(RAL_RTR = `RAL/RTR`, Expansion_Ref = `Expansion Reference`)
aPit <- d %>% filter(!is.na(GraLastDate)) %>% 
	select(MY, TagID, RelName, Expansion_Ref, RelSite, GraLastDate) %>% 
	mutate(GraLastDate = mdy(GraLastDate), rYear = year(GraLastDate)) %>% 
	mutate(tagRate = 1 / e$RAL_Expansion[match(Expansion_Ref, e$Expansion_Reference)]) %>%
	filter(rYear >= 2016)
e <- e %>% filter(e$Expansion_Reference %in% aPit$Expansion_Ref)
write.table(tagRates_uniform, "uniform_tag_rates2.txt", sep = "\t", row.names = FALSE, col.names = TRUE,
				quote = FALSE)
write.table(tagRates_odd, "nonuniform_tag_rates2.txt", sep = "\t", row.names = FALSE, col.names = TRUE,
				quote = FALSE)


maInfo %>% filter(PBT_RELEASE_GROUP %in% tagRates_odd$PBT_RELEASE_GROUP,
						Individual_name %in% trap$GenMa) %>%
	group_by(PBT_RELEASE_GROUP) %>% count(Pedigree_name) %>% dumpTable("odd_tags_by_ped.txt")

maInfo %>% filter(PBT_RELEASE_GROUP %in% tagRates_odd$PBT_RELEASE_GROUP,
						Individual_name %in% trap$GenMa) %>% filter(PBT_RELEASE_GROUP == "RapidR-Smolt", grepl("DWOR", Pedigree_name)) %>% as.data.frame

maInfo %>% filter(PBT_RELEASE_GROUP %in% tagRates_odd$PBT_RELEASE_GROUP,
						Individual_name %in% trap$GenMa) %>%
	group_by(PBT_RELEASE_GROUP, PBT_BY_HAT_STOCK) %>% count(Pedigree_name) %>% dumpTable("odd_tags_by_ped_take2.txt")


sink("problems_matching2.txt")
for(i in 1:nrow(e)){
	cat(unlist(e[i,]), sep = "\t")
	# +2 to y b/c y is BY, and MY is in e
	s <- which(tagRates_uniform$c3 == e$PBTcode[i] & (tagRates_uniform$y + 2) == e$Migr._Year[i])
	if(length(s) == 0){
		cat("\tno match\n")
		next
	}
	cat("\t")
	cat(tagRates_uniform$PBT_RELEASE_GROUP[s], sep = "\t")
	cat("\n")
}
sink()

eAll <- read_csv("data/expansions.csv")
colnames(eAll) <- gsub(" ", "_", colnames(eAll))
eAll %>% filter(!eAll$Expansion_Reference %in% aPit$Expansion_Ref) %>%
	write.table("notDetect_pit.txt", sep = "\t", row.names = FALSE, col.names = TRUE,
				quote = FALSE)


# ok, now the "odd" groups, most of them should be PBT_BY_HAT_STOCK + PBT_RELEASE_GROUP
# so editing the Mas and retrying
editMa <- maInfo %>% filter(PBT_RELEASE_GROUP %in% tagRates_odd$PBT_RELEASE_GROUP) %>%
	mutate(PBT_RELEASE_GROUP = paste0(PBT_BY_HAT_STOCK, "-", PBT_RELEASE_GROUP))

trap$relGroup[trap$GenMa %in% editMa$Individual_name] <- 
	editMa$PBT_RELEASE_GROUP[match(trap$GenMa[trap$GenMa %in% editMa$Individual_name], editMa$Individual_name)]

pit_pbt_matchup <- read_tsv("data/pit_pbt_match_ups.txt")
colnames(pit_pbt_matchup) <- gsub(" ", "_", colnames(pit_pbt_matchup))
pit_pbt_matchup %>% select(contains("MJB")) %>% gather(k, v, 1:4) %>% select(v) %>%
	distinct %>% filter(!v %in% trap$relGroup)
# all are in the data!

# now constructing distinguishable groups

# first assign each PBT group a "pool" number
pool_key <- pit_pbt_matchup %>% select(contains("MJB")) %>% gather(k, PBT_relGroup, 1:4) %>% select(PBT_relGroup) %>%
	distinct %>% filter(!is.na(PBT_relGroup)) %>% add_column(pool = 1:nrow(.))

# then combine pools based on co-occurrence of PBT groups
while(TRUE){
	print("pooling")
	nComb <- 0
	for(i in 1:nrow(pool_key)){
		f <- pool_key$PBT_relGroup[i]
		p <- pool_key$pool[i]
		other <- pit_pbt_matchup %>% filter(MJB_PBT_Release_Group_1 %in% f |
												MJB_PBT_Release_Group_2 %in% f |
												MJB_PBT_Release_Group_3 %in% f |
												MJB_PBT_Release_Group_4 %in% f) %>%
			select(contains("MJB")) %>% gather(k, PBT_relGroup, 1:4) %>% 
			select(PBT_relGroup) %>% distinct %>% filter(!is.na(PBT_relGroup), PBT_relGroup != f) %>%
			unlist
		other_p <- pool_key %>% filter(PBT_relGroup %in% other) %>% select(pool) %>%
			distinct %>% unlist
		if(any(other_p != p)){
			nComb <- nComb + 1
			pool_key$pool[pool_key$pool %in% other_p] <- p
		}
	}
	if(nComb == 0) break
}
pool_key %>% as.data.frame
n_distinct(pool_key$pool)

# then assign PIT tag groups to pools
# if pooling done correctly, can just use first PBT group name
pit_pbt_matchup <- pit_pbt_matchup %>% mutate(pool = pool_key$pool[match(MJB_PBT_Release_Group_1, pool_key$PBT_relGroup)])

# now outputting all pools for confirmation
allPITs <- c()
allPBTs <- c()
sink("confirm_pools.txt")
for(i in unique(pool_key$pool)){
	pit <- pit_pbt_matchup %>% filter(pool == i) %>% pull(Expansion_Reference)
	pbt <- pool_key %>% filter(pool == i) %>% pull(PBT_relGroup)
	allPITs <- c(allPITs, pit)# to check if any duplicates
	allPBTs <- c(allPBTs, pbt)
	cat(i, "\t", sep = "")
	cat(pit, sep = ",")
	cat("\t")
	cat(pbt, sep = ",")
	cat("\n")
}
sink()
n_distinct(allPITs)
nrow(pit_pbt_matchup)
n_distinct(allPBTs)
nrow(pool_key)

# MB manually corrected some pools, now loading in final
final_pools <- read_tsv("data/final_pools.txt")
t_pbt <- final_pools %>% select(PBT_groups) %>% unlist %>% strsplit(",") %>% unlist
t_pbt[!t_pbt %in% trap$relGroup]
table(t_pbt %in% trap$relGroup, useNA = "if")
# all are in the data!
# make all the ones not in this list NA (they are irrelevant for comparing to the PIT estimates)
trap$relGroup[!trap$relGroup %in% t_pbt] <- NA
trap %>% count(relGroup) %>% arrange(n) %>% as.data.frame

# now tag rates
tagRates <- editMa %>% bind_rows(maInfo) %>% select(PBT_RELEASE_GROUP, PBT_RELEASE_GROUP_TAGRATE) %>%
	distinct() %>% filter(!is.na(PBT_RELEASE_GROUP), PBT_RELEASE_GROUP %in% trap$relGroup) %>% arrange(PBT_RELEASE_GROUP)
tagRates %>% group_by(PBT_RELEASE_GROUP) %>% summarise(dups = n_distinct(PBT_RELEASE_GROUP_TAGRATE)) %>%
	filter(dups > 1)
dumpTable(tagRates, "tagRates.txt")

dumpTable(trap, "filtered_trap_data.txt")

trap %>% filter(relGroup == "2013-CLWH-SFSW-Powell") %>% count(Pedigree_name)
trap %>% filter(relGroup == "2013-CLWH-SFSW-Powell") %>% count(GenMa) %>% pull(GenMa) %>% writeClipboard()
trap %>% filter(relGroup == "2013-CLWH-SFSW-Powell") %>% pull(Individual_name) %>% writeClipboard()
trap %>% filter(relGroup == "2013-CLWH-SFSW-Powell") %>% select(Pedigree_name, Individual_name, GenMa, GenPa, Marks) %>% dumpTable("2013-CLWH-SFSW-Powell.txt")

trap %>% filter(relGroup == "2013-CLWH-SFSW-Powell") %>% count(Marks)
trap$relGroup[!is.na(trap$relGroup) & !is.na(trap$Marks) & trap$relGroup == "2013-CLWH-SFSW-Powell" & trap$Marks == "AD"] <- "2013-CLWH-SFSW-Powell-AD"
trap$relGroup[!is.na(trap$relGroup) & !is.na(trap$Marks) & trap$relGroup == "2013-CLWH-SFSW-Powell" & trap$Marks == "AI"] <- "2013-CLWH-SFSW-Powell-Adint"

trap %>% filter(relGroup == "2014-CLWH-Powell-SelwayR-AD/Adint") %>% count(Marks)
trap$relGroup[!is.na(trap$relGroup) & !is.na(trap$Marks) & trap$relGroup == "2014-CLWH-Powell-SelwayR-AD/Adint" & trap$Marks == "AD"] <- 
	"2014-CLWH-Powell-SelwayR-AD"
trap$relGroup[!is.na(trap$relGroup) & !is.na(trap$Marks) & trap$relGroup == "2014-CLWH-Powell-SelwayR-AD/Adint" & trap$Marks == "AI"] <- 
	"2014-CLWH-Powell-SelwayR-Adint"

trap %>% filter(relGroup == "2015-CLWH-SFSAL-Powell-AD/Adint") %>% count(Marks)
trap$relGroup[!is.na(trap$relGroup) & !is.na(trap$Marks) & trap$relGroup == "2015-CLWH-SFSAL-Powell-AD/Adint" & trap$Marks == "AD"] <- 
	"2015-CLWH-SFSAL-Powell-AD"
trap$relGroup[!is.na(trap$relGroup) & !is.na(trap$Marks) & trap$relGroup == "2015-CLWH-SFSAL-Powell-AD/Adint" & trap$Marks == "AI"] <- 
	"2015-CLWH-SFSAL-Powell-Adint"

trap %>% filter(relGroup == "2016-CLWH-DWOR-SelwayR-AD/Adint") %>% count(Marks)
trap$relGroup[!is.na(trap$relGroup) & !is.na(trap$Marks) & trap$relGroup == "2016-CLWH-DWOR-SelwayR-AD/Adint" & trap$Marks == "AD"] <- 
	"2016-CLWH-DWOR-SelwayR-AD"
trap$relGroup[!is.na(trap$relGroup) & !is.na(trap$Marks) & trap$relGroup == "2016-CLWH-DWOR-SelwayR-AD/Adint" & trap$Marks == "AI"] <- 
	"2016-CLWH-DWOR-SelwayR-Adint"

trap %>% filter(relGroup == "2016-CLWH-Powell-Powell-AD/Adint") %>% count(Marks)
trap$relGroup[!is.na(trap$relGroup) & !is.na(trap$Marks) & trap$relGroup == "2016-CLWH-Powell-Powell-AD/Adint" & trap$Marks == "AD"] <- 
	"2016-CLWH-Powell-Powell-AD"
trap$relGroup[!is.na(trap$relGroup) & !is.na(trap$Marks) & trap$relGroup == "2016-CLWH-Powell-Powell-AD/Adint" & trap$Marks == "AI"] <- 
	"2016-CLWH-Powell-Powell-Adint"

final_tag_rates <- read_tsv("data/final_tag_rates.txt")

final_tag_rates %>% filter(!PBT_RELEASE_GROUP %in% trap$relGroup)
table(t_pbt %in% final_tag_rates$PBT_RELEASE_GROUP, useNA = "if")
t_pbt[!t_pbt %in% final_tag_rates$PBT_RELEASE_GROUP]

# ok, updating release groups for 2011 and 2012 mas with new groups/data from Forrest
earlyMa <- read_tsv("data/maInfo_11_12.txt", col_types = "ccccc")
trap$relGroup[trap$GenMa %in% earlyMa$`Individual name`] <- earlyMa$`PBT_RELEASE GROUP`[match(trap$GenMa[trap$GenMa %in% earlyMa$`Individual name`], earlyMa$`Individual name`)]

final_tag_rates %>% filter(!PBT_RELEASE_GROUP %in% trap$relGroup)
# all the groups with tag rates are sampled
trap %>% select(relGroup) %>% filter(!is.na(relGroup), !relGroup %in% final_tag_rates$PBT_RELEASE_GROUP) %>%
	count(relGroup) %>% as.data.frame()

trap %>% filter(relGroup == "2012-CLWH-Powell-Selway-AD/Adint") %>% count(Marks)
trap$relGroup[!is.na(trap$relGroup) & !is.na(trap$Marks) & trap$relGroup == "2012-CLWH-Powell-Selway-AD/Adint" & trap$Marks == "AD"] <- 
	"2012-CLWH-Powell-Selway-AD"
trap$relGroup[!is.na(trap$relGroup) & !is.na(trap$Marks) & trap$relGroup == "2012-CLWH-Powell-Selway-AD/Adint" & trap$Marks == "AI"] <- 
	"2012-CLWH-Powell-Selway-Adint"

final_tag_rates %>% filter(!PBT_RELEASE_GROUP %in% trap$relGroup)
# all the groups with tag rates are sampled

# now let's look at the groups that don't have tag rates
trap %>% select(relGroup) %>% filter(!is.na(relGroup), !relGroup %in% final_tag_rates$PBT_RELEASE_GROUP) %>%
	count(relGroup) %>% as.data.frame() %>% dumpTable("no_tag_rate_groups.txt")

rm(genoData)
# asked Matt and Forrest for tag rates for relevant groups
trap %>% filter(!is.na(relGroup), !relGroup %in% final_tag_rates$PBT_RELEASE_GROUP) %>% 
	select(Pedigree_name, Individual_name, GenMa, GenPa, relGroup) %>% dumpTable("individual data for Forrest.txt")

trap %>% filter(relGroup == "SFSalmonR-IntSmolt") %>% pull(GenMa)

# updating relGroups and tag rates with info from Forrest
updateRel <- tibble(old = c("2011-DWOR-NFClearwaterR-Smolt",
	"2012-CLWH-ClearCr-Smolt",
	"2012-CLWH-Powell-Smolt",
	"2012-DWOR-NFClearwaterR-Smolt",
	"2012-KOOS-ClearCr-Smolt",
	"2012-MCCA-SFSalmonR-IntSmolt",
	"2012-RAPH-RapidR-Smolt",
	"SFSalmonR-IntSmolt"), PBT_RELEASE_GROUP = c("2011-DWOR-DWOR-NFClearwaterR",
	"2012-CLWH-POWELL-ClearCr",
	"2012-CLWH-SFSW-Powell",
	"2012-DWOR-DWOR-NFClearwaterR",
	"2012-KOOS-KOOS-ClearCr",
	"2012-MCCA-SFSW-SFSalmonR-Int",
	"2012-RAPH-RAPH-RapidR",
	"2011-MCCA-SF SAL-SFSalmonR-IntSmolt"
	), PBT_RELEASE_GROUP_TAGRATE = c(0.986,
	0.9934,
	0.9917,
	0.9937,
	0.996,
	1,
	0.9536,
	0.9846)
)
for(i in 1:nrow(updateRel)) trap$relGroup[!is.na(trap$relGroup) & trap$relGroup == updateRel$old[i]] <- updateRel$PBT_RELEASE_GROUP[i]
final_tag_rates <- bind_rows(final_tag_rates, updateRel[,2:3]) %>% distinct()
nrow(final_tag_rates) == n_distinct(final_tag_rates$PBT_RELEASE_GROUP)
trap %>% select(relGroup) %>% filter(!is.na(relGroup), !relGroup %in% final_tag_rates$PBT_RELEASE_GROUP) %>%
	count(relGroup) %>% as.data.frame()
# these remaining groups are all irrelevant, so mark as NA
trap$relGroup[!is.na(trap$relGroup) & !trap$relGroup %in% final_tag_rates$PBT_RELEASE_GROUP] <- NA
trap %>% filter(!relGroup %in% final_tag_rates$PBT_RELEASE_GROUP) %>% select(relGroup) %>% distinct
# beautiful
dumpTable(trap, "final_filtered_trap_data.txt")
dumpTable(final_tag_rates, "really_final_tag_rates.txt")


# and now... redo the pooling for SY11 and 12!!!!!

final_pools %>% select(PBT_groups) %>% unlist %>% strsplit(",") %>% unlist %>%
	as_tibble %>% rename(PBT_relGroup = value) %>%
	filter(!PBT_relGroup %in% final_tag_rates$PBT_RELEASE_GROUP) %>% as.data.frame() %>% pull(PBT_relGroup) %>% writeClipboard()
sum(trap$relGroup == "2016-CLWH-DWOR-SelwayR-Adint", na.rm = T)

final_tag_rates %>% filter(PBT_RELEASE_GROUP == "2011-SAWT-UpperSalmonR-SegSmolt")

rf_pools <- read_tsv("really_final_pools.txt")

rf_pools %>% select(PBT_groups) %>% unlist %>% strsplit(",") %>% unlist %>%
	as_tibble %>% rename(PBT_relGroup = value) %>%
	filter(!PBT_relGroup %in% final_tag_rates$PBT_RELEASE_GROUP) %>% as.data.frame() 

table(trap$relGroup[!is.na(trap$relGroup)] %in% (rf_pools %>% select(PBT_groups) %>% unlist %>% strsplit(",") %>% unlist), useNA = "if")

trap %>% filter(!is.na(relGroup), !relGroup %in% (rf_pools %>% select(PBT_groups) %>% unlist %>% strsplit(",") %>% unlist)) %>%
	count(relGroup) %>% pull(relGroup) %>% writeClipboard

# these groups aren't relevant, so zeroing out
trap$relGroup[!is.na(trap$relGroup) & !trap$relGroup %in% (rf_pools %>% select(PBT_groups) %>% unlist %>% strsplit(",") %>% unlist)] <- NA
trap %>% filter(!is.na(relGroup), !relGroup %in% (rf_pools %>% select(PBT_groups) %>% unlist %>% strsplit(",") %>% unlist)) %>%
	count(relGroup)

dumpTable(trap, "really_final_filtered_trap_data.txt")


# the three "really_final_..." files have the filtered data to use for estimates
# Let's perform the modeling in a different file

# model
# Maybe a walk could work for the numbers of each group, and the 
# lognormal prior work b/c it puts a good amount of weight close to 0?
