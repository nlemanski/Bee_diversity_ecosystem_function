## load packages ##
library(plyr)
library(reshape)
library(tidyverse)

## file paths and file names ##
path <- "C:/Users/natal/OneDrive - Rutgers University/Documents/winfree lab/SQL data/"
figpath <- "C:/Users/natal/OneDrive - Rutgers University/Documents/winfree lab/figures"
rpath <- "C:/Users/natal/OneDrive - Rutgers University/Documents/winfree lab/R data/"
outpath <- "C:/Users/natal/OneDrive - Rutgers University/Documents/winfree lab/BEFresults/"

## import crop visit data as RData file ##
crops = c("blue","cran","njwat","cawat")
cropnames = c("Blueberry", "Cranberry", "NJ Watermelon", "CA Watermelon")

load(paste(rpath,'df_visits_cawat.RData',sep=''))  # load data
load(paste(rpath,'df_visits_blue.RData',sep=''))  # load data
load(paste(rpath,'df_visits_njwat.RData',sep=''))  # load data

df_species_blue = data.frame(sort(unique(df_visits_blue$gen_sp)))
df_species_njwat = data.frame(sort(unique(df_visits_njwat$gen_sp)))
df_species_cawat = data.frame(sort(unique(df_visits_cawat$gen_sp)))

save(df_species_blue, file=paste(rpath,'species_blue.Rdata',sep=''))
save(df_species_njwat, file=paste(rpath,'species_njwat.Rdata',sep=''))
save(df_species_cawat, file=paste(rpath,'species_cawat.Rdata',sep=''))


# find total abundance at each site-date
df_blue = df_visits_blue %>% group_by(round,year,site) %>% summarise(abundance = sum(visits))
df_njwat = df_visits_njwat %>% group_by(round,year,site) %>% summarise(abundance = sum(visits))
df_cawat = df_visits_cawat %>% group_by(round,year,site) %>% summarise(abundance = sum(visits))

save(df_blue, file=paste(rpath,'abundance_blue.Rdata',sep=''))
save(df_njwat, file=paste(rpath,'abundance_njwat.Rdata',sep=''))
save(df_cawat, file=paste(rpath,'abundance_cawat.Rdata',sep=''))

# find mean, min, and max abundance at each site
abund_blue = df_blue %>% 
  group_by(site) %>% 
  summarise(min.abundance = min(abundance), max.abundance = max(abundance), mean.abundance = mean(abundance), med.abundance = median(abundance))

abund_njwat = df_njwat %>% 
  group_by(site) %>% 
  summarise(min.abundance = min(abundance), max.abundance = max(abundance), mean.abundance = mean(abundance), med.abundance = median(abundance))

abund_cawat = df_cawat %>% 
  group_by(site) %>% 
  summarise(min.abundance = min(abundance), max.abundance = max(abundance), mean.abundance = mean(abundance), med.abundance = median(abundance))

save(abund_blue, file=paste(rpath,'mean_abundance_blue.Rdata',sep=''))
save(abund_njwat, file=paste(rpath,'mean_abundance_njwat.Rdata',sep=''))
save(abund_cawat, file=paste(rpath,'mean_abundance_cawat.Rdata',sep=''))
