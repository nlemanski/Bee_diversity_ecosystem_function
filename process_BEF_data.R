## process crop visitation data for analysis ##

## load packages ##
library(tidyverse)

## file paths and file names ##
path <- "D:/natal/D_Documents/winfree lab/SQL data/"
path2 <- "D:/natal/D_Documents/winfree lab/Williams lab data/"
figpath <- "D:/natal/D_Documents/winfree lab/figures"
rpath <- "D:/natal/D_Documents/winfree lab/R data"
outpath <- "D:/natal/D_Documents/winfree lab/BEFresults/"

fname1 <- "blueberry_visit_table.csv"
fname2 <- "cranberry_visit_table.csv"
fname3 <- "watermelon_visit_table.csv"
fname4 <- "ca_wat_visits.csv"
fname5 <- "Wat_allsites_years_distances.csv"
fname6 <- "Site names_by_ID.csv"


## import data file as a dataframe ##
df_visits_blue <- read.csv(paste(path,fname1,sep=""),header=TRUE,stringsAsFactors=FALSE)
df_visits_cran <- read.csv(paste(path,fname2,sep=""),header=TRUE,stringsAsFactors=FALSE)
df_visits_njwat <- read.csv(paste(path,fname3,sep=""),header=TRUE,stringsAsFactors=FALSE)
df_visits_cawat <- read.csv(paste(path2,fname4,sep=""),header=TRUE,stringsAsFactors=FALSE)
# sitedistwat <- read.csv(paste(path2,fname5,sep=""),header=TRUE,stringsAsFactors=FALSE)
# sitenameswat <- read.csv(paste(path2,fname6,sep=""),header=TRUE,stringsAsFactors=FALSE)

## convert columns to correct data types ##
## blueberry data
glimpse(df_visits_blue)  # see what current data types are for each column
df_visits_blue = mutate(df_visits_blue, round=factor(round), year=as.character(year), site=factor(site))
df_visits_blue$SV_pollen_tubes <- as.numeric(df_visits_blue$SV_pollen_tubes)
df_visits_blue$SV_pollen_grains <- as.numeric(df_visits_blue$SV_pollen_grains)

# standardize species names and re-summarize across new names
df_visits_blue = mutate(df_visits_blue, gen_sp = gsub('calcarata_or_dupla','calcarata', gen_sp)) # assigned to c. calcarata
df_visits_blue = mutate(df_visits_blue, gen_sp = gsub('dup_calc','calcarata', gen_sp))
df_visits_blue = df_visits_blue %>% 
  group_by(round, year, site, gen_sp, svgroup) %>%
  summarize(visits = sum(visits), SV_pollen_tubes = mean(SV_pollen_tubes), SV_pollen_grains = mean(SV_pollen_grains))
# filter out apis mellifera
df_visits_blue = filter(df_visits_blue, gen_sp != 'Apis_mellifera')

summary(df_visits_blue)
sort(unique(df_visits_blue$gen_sp))

## cranberry data
glimpse(df_visits_cran)  # see what current data types are for each column
df_visits_cran = mutate(df_visits_cran, round=factor(round), year=as.character(year), site=factor(site))
df_visits_cran$tet_tubes <- as.numeric(df_visits_cran$tet_tubes)
df_visits_cran$tet_notubes <- as.numeric(df_visits_cran$tet_notubes)
meantubes   = mean(df_visits_cran$tet_tubes, na.rm=TRUE) ## fill in NA values with mean function for all species groups
meannotubes = mean(df_visits_cran$tet_notubes, na.rm=TRUE)
df_visits_cran$tet_tubes[is.na(df_visits_cran$tet_tubes)] <- meantubes
df_visits_cran$tet_notubes[is.na(df_visits_cran$tet_notubes)] <- meannotubes
summary(df_visits_cran)

## nj watermelon
glimpse(df_visits_njwat)
df_visits_njwat$round[df_visits_njwat$round == "NULL"] = "1"  # rename "null" rounds to "1"
df_visits_njwat$SV_pollen = as.numeric(df_visits_njwat$SV_pollen)
df_visits_njwat = mutate(df_visits_njwat, site=factor(site), year=factor(year), round=factor(round))

# standardize species names and re-summarize across new names
df_visits_njwat = mutate(df_visits_njwat, gen_sp = gsub('_seeTN','',gen_sp))
df_visits_njwat = mutate(df_visits_njwat, gen_sp = gsub('calcarata_dupla_mikmaqi','calcarata',gen_sp))
df_visits_njwat = mutate(df_visits_njwat, gen_sp = gsub('hitchensi_weemsi','weemsi',gen_sp)) # assigned to l. weemsi b/c morphogroup is small dark
df_visits_njwat = df_visits_njwat %>%
  group_by(site, year, round, genus, gen_sp, taxon_group) %>%
  summarize(visits = sum(visits), SV_pollen = mean(SV_pollen))

# fill in NA values with mean function for all species groups
meantubes  = mean(df_visits_njwat$SV_pollen, na.rm=TRUE)  
df_visits_njwat$SV_pollen[is.na(df_visits_njwat$SV_pollen)] = meantubes
df_visits_njwat = filter(df_visits_njwat, year!="0")  # remove 2 obs with no year or date listed
#filter(df_visits_njwat, site=="del")

summary(df_visits_njwat)
sort(unique(df_visits_njwat$gen_sp))

## CA watermelon ##
glimpse(df_visits_cawat)
sort(unique(df_visits_cawat$bee_sps))

# # join site names to site distances to see which names correspond to the same site
# watsites <- full_join(sitedistwat, sitenameswat, by=c("ï..SiteID" = "SiteID"))
# watsites = watsites %>% relocate(Site_common, .after = "ï..SiteID")
# watsites = watsites %>% relocate(Site_name, .after = "Site_common")
# watsites = left_join(watsites, sitenameswat, by=c("TargetSiteID" = "SiteID"))
# watsites = watsites %>% relocate(Site_common.y, .after = "TargetSiteID")
# watsites = watsites %>% relocate(Site_name.y, .after = "Site_common.y")
# 
# # see which site names are in watermelon visit table and distance table
# visitsites = sort(unique(df_visits_cawat$site))
# distsites = sort(unique(watsites$Site_common.x))
# 
# # see distribution of site distances
# hist(sitedistwat$Distance, breaks=120, xlim=c(0,20000))
# 
# # look at a subset of sites that are close to each other
# distt = 500
# watsitesnear = filter(watsites, Distance < distt)
# # based on looking at this subset, less than 500 seems reasonable cutoff

## add a column for total function (number of visits * single visit function measure)
df_visits_blue$fun = df_visits_blue$visits * df_visits_blue$SV_pollen_tubes
df_visits_cran$fun = df_visits_cran$visits * df_visits_cran$tet_tubes    # pollen tetrads with tubes
df_visits_njwat$fun = df_visits_njwat$visits * df_visits_njwat$SV_pollen
df_visits_cawat$fun = df_visits_cawat$visits * df_visits_cawat$SIGmean   # mean single visit deposition for significant visits

# combine males and females
df_visits_cawat <- df_visits_cawat %>%
  group_by(year,round,site,bee_sps) %>%
  summarise(visits = sum(visits), 
            SIGmean = mean(SIGmean),
            fun = sum(fun))
# correct column formats
df_visits_cawat = mutate(df_visits_cawat, site=factor(site), year=factor(year), round=factor(round))
df_visits_cawat = rename(df_visits_cawat, gen_sp=bee_sps)
levels(df_visits_cawat$site)

## save as R files ##
save(df_visits_blue, file=paste(rpath,'/','df_visits_blue.RData',sep=''))
save(df_visits_cran, file=paste(rpath,'/','df_visits_cran.RData',sep=''))
save(df_visits_njwat, file=paste(rpath,'/','df_visits_njwat_all.RData',sep=''))
save(df_visits_cawat, file=paste(rpath,'/','df_visits_cawat.RData',sep=''))

# filter out year 2004 for nj wat
df_visits_njwat = filter(df_visits_njwat, year != "2004")
df_visits_njwat = mutate(df_visits_njwat, site=factor(site), year=factor(year), round=factor(round))
save(df_visits_njwat, file=paste(rpath,'/','df_visits_njwat.RData',sep=''))

# filter out sites sampled only one year for ca wat
cawat_samples <- df_visits_cawat %>%
  group_by(site) %>%
  summarise(nyears=length(unique(year)),
            nrounds=length(unique(round)))

summary(df_visits_njwat)
unique(df_visits_cawat$site)
length(unique(df_visits_cawat$site))

# create tables of morphogroups and their species and pollen values
# nj watermelon
pollgroups_njwat <- df_visits_njwat %>%
  summarise(species = gen_sp,
            taxon_group = taxon_group,
            pollen = SV_pollen)
pollgroups_njwat <- distinct(pollgroups_njwat)

write.csv(pollgroups_njwat,paste(path,'morphogroups_njwat.csv',sep=''))

# nj blueberry
pollgroups_blue <- df_visits_blue %>%
  summarise(species = gen_sp,
            taxon_group = svgroup,
            pollen = SV_pollen_tubes)
pollgroups_blue <- distinct(pollgroups_blue)

write.csv(pollgroups_blue,paste(path,'morphogroups_blue.csv',sep=''))

# ca watermelon
pollgroups_cawat <- df_visits_cawat %>% ungroup() %>% select(gen_sp, SIGmean)
pollgroups_cawat <- distinct(pollgroups_cawat)

write.csv(pollgroups_cawat,paste(path,'morphogroups_cawat.csv',sep=''))

