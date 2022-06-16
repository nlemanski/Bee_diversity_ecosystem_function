### clean up BEF raw data where each row is a single bee specimen, convert to correct format, then save as rdata files ###

## load packages ##
library(tidyverse)

## file paths and file names ##
path <- "C:/Documents/Bee_diversity_ecosystem_function/SQL data/"
figpath <- "C:/Documents/Bee_diversity_ecosystem_function/figures"
rpath <- "C:/Documents/Bee_diversity_ecosystem_function/R data/"
outpath <- "C:/Documents/Bee_diversity_ecosystem_function/BEFresults/"
path2 <- "C:/Documents/Bee_diversity_ecosystem_function/Williams lab data/"

# load nj sql data files
crops = c("blue", "njwat")
crop = "njwat"
for (crop in c("blue","njwat")) {
  fname = paste(path,crop,'_obs.csv',sep='')
  df_visits = read.csv(fname, header=T, stringsAsFactors=F)
  glimpse(df_visits)
  summary(df_visits)
  
  # convert to correct data types and clean up data
  # remove  obs with no year or date listed
  df_visits = filter(df_visits, year!="0")
  
  # convert to correct data type
  df_visits$round[df_visits$round == "NULL"] = "1"  # rename "null" rounds to "1"
  if (crop == "blue") { df_visits$SV_pollen_tubes = as.numeric(df_visits$SV_pollen_tubes) }
  if (crop == "njwat") {df_visits$SV_pollen = as.numeric(df_visits$SV_pollen)}
  df_visits = mutate(df_visits, site=factor(site), year=factor(year), round=factor(round))
  
  # standardize species names
  df_visits = mutate(df_visits, gen_sp = gsub('_seeTN','',gen_sp))
  df_visits = mutate(df_visits, gen_sp = gsub('calcarata_dupla_mikmaqi','calcarata',gen_sp))
  df_visits = mutate(df_visits, gen_sp = gsub('calcarata_or_dupla','calcarata', gen_sp)) # assigned to c. calcarata
  df_visits = mutate(df_visits, gen_sp = gsub('dup_calc','calcarata', gen_sp))
  df_visits = mutate(df_visits, gen_sp = gsub('hitchensi_weemsi','weemsi',gen_sp)) # assigned to l. weemsi b/c morphogroup is small dark
  
  # filter out apis mellifera
  df_visits = filter(df_visits, gen_sp != 'Apis_mellifera')
  
  # fill in NA values with mean function for all species groups
  if (crop == "blue") {
    meantubes  = mean(df_visits$SV_pollen_tubes, na.rm=TRUE)
    df_visits$SV_pollen_tubes[is.na(df_visits$SV_pollen_tubes)] = meantubes }
  if (crop == "njwat") {
    meantubes  = mean(df_visits$SV_pollen, na.rm=TRUE)
    df_visits$SV_pollen[is.na(df_visits$SV_pollen)] = meantubes
  }
  
  if (crop == "blue") {df_obs_blue = df_visits}
  if (crop == "njwat") {df_obs_njwat = df_visits}
  
}
summary(df_obs_blue)
summary(df_obs_njwat)

# save as rdata files
save(df_obs_blue,file=paste(rpath,'obs_blue.Rdata',sep=''))
save(df_obs_njwat,file=paste(rpath,'obs_njwat.Rdata',sep=''))

### clean up CA data and convert to rdata for null model analysis ###
crop = "cawat"
df_obs_cawat = read.csv(paste(path2,"3yrs_wbees_oncrop.csv",sep=''))
summary(df_obs_cawat)

# select watermelon on the crop transect only
df_obs_cawat = filter(df_obs_cawat, flower_sps == 'Citrullus_lanatus', ON_OFF=='ON')

## combine equivalent watermelon sites, rename rounds
# sites  that are actually the second three rounds get relabeled round 4-6
df_obs_cawat = df_obs_cawat %>% 
  mutate(round = ifelse(site %in% c('PAC2','RIVCAC2','TERWIN2','FULOFF2','DURST2','EAT2','FREE2'), (round + 3), round))

# site TERWIN has a third set of 3 rounds, relabel as 7-9
df_obs_cawat = df_obs_cawat %>% 
  mutate(round = ifelse(site %in% c('TERWIN3'), (round + 6), round))

# site PAC has rounds with different labels in 2010, need to relabel site 98 and 99 to 1 and 2
df_obs_cawat = df_obs_cawat %>% 
  mutate(round = ifelse(site == 'PAC', 
                        ifelse(year == '2010', 
                               ifelse(round > 90, 
                                      (round-97), 
                                      (round+2)),
                               round), 
                        round))

# rename sites so that sites sampled multiple seasons have the same name for each round
df_obs_cawat = df_obs_cawat %>% 
  mutate(site = ifelse(str_detect(site,'PAC'), 'PAC',
                ifelse(str_detect(site,'RIVCAC'), 'RIVCAC',
                ifelse(str_detect(site,'TERWIN'), 'TERWIN',
                ifelse(str_detect(site,'FULOFF'), 'FULOFF',
                ifelse(str_detect(site,'DURST'), 'DURST',
                ifelse(str_detect(site,'EAT'), 'EAT',
                ifelse(str_detect(site,'FREE'), 'FREE',
                       site))))))))

# get rid of unnecessary columns
df_obs_cawat = df_obs_cawat %>% 
  select(year, date, round, site, bee_genus, bee_sps, sex)

## import single-visit function data
fname1 = paste(path2,"CAL2001Melon_vis_depFINAL_meandep.csv",sep='')
df_function = read.csv(fname1,header=T, stringsAsFactors=T)
summary(df_function)

# rename columns
df_function = df_function %>% rename(SV_group = ï..Species)
df_function = df_function %>% rename(sex = Sex)
df_function = df_function %>% rename(SV_pollen = SIGmean)

# get rid of unnecessary columns
df_function = df_function %>% select(SV_group, sex, SV_pollen)

## add functional group names for each species/genus
df_obs_cawat = df_obs_cawat %>%
  mutate(SV_group = ifelse(bee_genus == 'Agapostemon', 'Agapostemon', 
                    ifelse(bee_genus == 'Anthophora', 'Anthophora urbana',
                    ifelse(bee_sps == 'Apis_mellifera', 'Apis mellifera',
                    ifelse(bee_sps == 'Bombus_californicus', 'Bombus californicus',
                    ifelse(bee_sps == 'Bombus_vosnesenskii', 'Bombus vosnesenskii',
                    ifelse(bee_genus == 'Bombus', 'Bombus other',
                    ifelse(bee_genus == 'Dialictus', 'Dialictus',
                    ifelse(bee_genus == 'Evylaeus', 'Evylaeus',
                    ifelse(bee_sps == 'Halictus_farinosus', 'Halictus farinosus',
                    ifelse(bee_sps == 'Halictus_ligatus', 'Halictus ligatus',
                    ifelse(bee_sps == 'Halictus_tripartitus', 'Halictus tripartitus',
                    ifelse(bee_genus == 'Hylaeus', 'Hylaeus',
                    ifelse(bee_genus == 'Lasioglossum_Dialictus', 'Dialictus',
                    ifelse(bee_genus ==  'Lasioglossum_Evylaeus', 'Evylaeus',
                    ifelse(bee_genus == 'Lasioglossum_Lasioglossum' , 'Lasioglossum',
                    ifelse(bee_genus == 'Melissodes', 'Melissodes',
                    ifelse(bee_genus == 'Peponapis', 'Peponapis pruinosa', 
                           'other'))))))))))))))))))

# rename columns to match nj data
df_obs_cawat = df_obs_cawat %>% 
  rename(genus = bee_genus, gen_sp = bee_sps)

# fill in "U" for unknown sex
df_obs_cawat = df_obs_cawat %>% replace_na(list(sex='U'))

# add single visit function values to CA watermelon bee observations
df_obs_cawat = left_join(df_obs_cawat, df_function, by=c("SV_group", "sex"))

# convert columns to correct data type
df_obs_cawat = df_obs_cawat %>% 
  mutate(year = factor(year), round = factor(round), site = factor(site), sex = factor(sex), SV_group = factor(SV_group))

# assign the average of all svgroups as pollen function value to "other" group
meanfunall = mean(df_obs_cawat$SV_pollen, na.rm=T)
df_obs_cawat = df_obs_cawat %>% replace_na(list(SV_pollen=meanfunall))

# assign average of bombus svpollen values to "bombus other"
meanb1 = filter(df_obs_cawat, SV_group == 'Bombus californicus')$SV_pollen %>% mean()
meanb2 = filter(df_obs_cawat, SV_group == 'Bombus vosnesenskii')$SV_pollen %>% mean()
meanfunbombus = mean(c(meanb1,meanb2))
df_obs_cawat = df_obs_cawat %>% 
  mutate(SV_pollen = ifelse(SV_group == 'Bombus other', meanfunbombus, SV_pollen))

# drop rows with no species IDs
df_obs_cawat = df_obs_cawat %>% drop_na()

summary(df_obs_cawat)

# save as rdata file
save(df_obs_cawat,file=paste(rpath,'obs_cawat.Rdata',sep=''))
