# import the western watermelon data set and clean data for use in the analysis
# each row is a bee specimen collected off a watermelon flower in the crop transect

## import packages ##
library(tidyverse)

## import CA bee visit data ##
path = "Bee_diversity_ecosystem_function/" # change to filepath where the raw data is located


# import California watermelon bee observations and process for analysis
fname1 = paste(path,'3yrs_wbees_oncrop.csv',sep='')
df_ca_wat = read.csv(fname1,header=T, stringsAsFactors=F)

## select only watermelon data ##
df_ca_wat1 = df_ca %>% filter(flower_sps=='Citrullus_lanatus')
# select only on-crop transect data
df_ca_wat = df_ca_wat1 %>% filter(ON_OFF=='ON')
head(as_tibble(df_ca_wat))

## find unique sites, years, and sampling rounds ##
sites = unique(df_ca_wat$site)
years = unique(df_ca_wat$year)
rounds = unique(df_ca_wat$round)
plants = unique(df_ca_wat$flower_sps)
dates = sort(unique(df_ca_wat$date))
unique(df_ca_wat$sex)

df_ca_wat$sex <- df_ca_wat$sex %>%
  replace_na('unknown')

## combine equivalent watermelon sites, rename rounds ##
siteswat = sort(unique(df_ca_wat$site))
sitedays = df_ca_wat %>% group_by(site,year,round) %>%
  summarize(date = unique(date))

# assign new round numbers 1-6 for the sites sampled multiple seasons per year
# sites that are actually the second three rounds get relabeled round 4-6
df_ca_wat = df_ca_wat %>% 
  mutate(round_original = round, site_original = site) # to record original round and site names

df_ca_wat = df_ca_wat %>% 
  mutate(round = ifelse(site %in% c('PAC2','RIVCAC2','TERWIN2','FULOFF2','DURST2','EAT2','FREE2'), (round + 3), round))

# site TERWIN has a third set of 3 rounds, relabel as 7-9
df_ca_wat = df_ca_wat %>% 
  mutate(round = ifelse(site %in% c('TERWIN3'), (round + 6), round))

# site PAC has rounds with different labels in 2010, need to relabel site 98 and 99 to 1 and 2
df_ca_wat = df_ca_wat %>% 
  mutate(round = ifelse(site == 'PAC', 
                        ifelse(year == '2010', 
                               ifelse(round > 90, 
                                      (round-97), 
                                      (round+2)),
                               round), 
                        round))

# rename sites so that sites sampled multiple seasons have the same name for each round
df_ca_wat = df_ca_wat %>% 
  mutate(site = ifelse(str_detect(site,'PAC'), 'PAC',
                       ifelse(str_detect(site,'RIVCAC'), 'RIVCAC',
                              ifelse(str_detect(site,'TERWIN'), 'TERWIN',
                                     ifelse(str_detect(site,'FULOFF'), 'FULOFF',
                                            ifelse(str_detect(site,'DURST'), 'DURST',
                                                   ifelse(str_detect(site,'EAT'), 'EAT',
                                                          ifelse(str_detect(site,'FREE'), 'FREE',
                                                                 site))))))))


## find number of visits of each bee species for each site-date in watermelon
df_cawatvis = df_ca_wat %>% 
  group_by(year,site,round,bee_genus,bee_sps,sex) %>%
  summarise(visits = length(bee_sps))

## import single-visit function data ##
fname1 = paste(path,"CAL2001Melon_vis_depFINAL_meandep.csv",sep='')
df_function = read.csv(fname1)
df_function = df_function %>% select(1:5)
df_function = df_function %>% rename(fungroup=ï..Species)
fungroups = unique(df_function$fungroup)

# fill in female function values for groups where male value missing
hyl = filter(df_function, fungroup == 'Hylaeus', Sex == 'F')$SIGmean
df_function$SIGmean = ifelse(df_function$fungroup == 'Hylaeus', hyl, df_function$SIGmean)
las = filter(df_function, fungroup == 'Lasioglossum', Sex == 'F')$SIGmean
df_function$SIGmean = ifelse(df_function$fungroup == 'Lasioglossum', las, df_function$SIGmean)

## add functional group names for each species/genus
df_cawatvis = df_cawatvis %>%
  mutate(fungroup = ifelse(bee_genus == 'Agapostemon', 'Agapostemon', 
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

## add pollen function values to the visitation data 
# separate into male and female
dffemale = filter(df_cawatvis,sex=='F')
dfmale = filter(df_cawatvis, sex=='M')
dfnan = filter(df_cawatvis, sex=='unknown')

funfemale = filter(df_function,Sex=='F')
funmale = filter(df_function,Sex=='M')

# add function values
dffemale = left_join(dffemale,funfemale, by="fungroup")
dfmale = left_join(dfmale,funmale, by="fungroup")
dfnan = left_join(dfnan,funfemale, by="fungroup")

# rejoin into one dataframe
df_cawatvis = rbind(dffemale,dfmale,dfnan)
# df_cawatvis = pd.concat([dffemale,dfmale,dfnan])

# assign the average of all svgroups as pollen function value to "other" group
meanfunall = mean(df_function$SIGmean)
df_cawatvis$SIGmean = ifelse(df_cawatvis$fungroup == 'other', meanfunall, df_cawatvis$SIGmean)

meanb1 = mean(filter(df_function,fungroup=='Bombus californicus')$SIGmean)
meanb2 = mean(filter(df_function, fungroup=='Bombus vosnesenskii')$SIGmean)
meanfunbombus = mean(meanb1,meanb2)
df_cawatvis$SIGmean = ifelse(df_cawatvis$fungroup == 'Bombus other', meanfunbombus, df_cawatvis$SIGmean)

# drop rows with no species IDs
df_cawatvis = df_cawatvis %>% filter(!is.na(bee_sps))

# save watermelon visit data to file
fname = paste(path,'ca_wat_visits.csv',sep='')
write.csv(df_cawatvis, file=fname, row.names=F)

ca_wat_visits = read.csv(paste(path,'ca_wat_visits.csv',sep=''))

