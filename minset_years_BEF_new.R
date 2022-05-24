## Find the minimum set of species needed to meet a pollination function threshold for all years in a given timescale ##
# require threshold to be met for every round in each year
# repeat analysis for each site

## load packages ##
library(plyr)
library(reshape)
library(tidyverse)
library(gaoptim)
library(abind)

setwd("C:/Documents/Bee_diversity_ecosystem_function/R code")

## load the species required function
source(file="minfinder_function_time.R")

## file paths and file names ##
path <- "C://Documents/Bee_diversity_ecosystem_function/SQL data/"
figpath <- "C:/Documents/Bee_diversity_ecosystem_function/figures/"
rpath <- "C:/Documents/Bee_diversity_ecosystem_function/R data/"
outpath <- "C:/Documents/Bee_diversity_ecosystem_function/BEFresults/"

## import crop visit data as RData file ##
crops = c("blue","njwat","cawat")
cropnames = c("Blueberry", "NJ Watermelon", "CA Watermelon")
# make sure crop name matches data file name
crop = "cawat"  # select crop to run

# load data for crop
if (crop == "cawat") {
  load(paste(rpath,'df_visits_cawat.RData',sep=''))  # load data
  rmax = 3 # minimum number of rounds needed for inclusion
  df_visits <- df_visits_cawat  # assign crop data to df_visits
} else if (crop == "blue") {
  load(paste(rpath,'df_visits_blue.RData',sep=''))  # load data
  rmax = 3 # minimum number of rounds needed for inclusion
  df_visits <- df_visits_blue  # assign crop data to df_visits
} else if (crop == "njwat") {
  load(paste(rpath,'df_visits_njwat.RData',sep=''))  # load data
  rmax = 1 # minimum number of rounds needed for inclusion
  df_visits <- df_visits_njwat  # assign crop data to df_visits
}

# subset data to only include site-years that have same number of rounds, necessary so we have equal sampling effort
dfrmax <- df_visits %>% 
  group_by(site, year) %>%
  summarize(numr = n_distinct(round))
df_visits <- left_join(df_visits, dfrmax, by=c("site"="site", "year"="year"))
df_visits <- df_visits %>% filter(numr >= rmax)

summary(df_visits)
df_visits <- df_visits %>% mutate(year = factor(year), round = factor(round))

## list unique values of each factor ##
rounds <- sort(unique(df_visits$round))
years <- sort(unique(df_visits$year))
species <- sort(unique(df_visits$gen_sp))
sites <- sort(unique(df_visits$site))

## function threshold percentage (set here for sensitivity)
func_percent = 0.50 # value for main results
# func_percent = 0.25 # value for sensitivity analysis
# func_percent = 0.75 # value for sensitivity analysis

# function threshold value (aggregate pollen grains) from percentage, using mean across all site-year-rounds
df_totalfun = ddply(df_visits, .(round, year, site), summarize, totalfun = sum(fun))

func_level = mean(df_totalfun$totalfun, na.rm=T)*func_percent


## main program loops:
# loop through all sites. For each site-round, we will find the minimum species set at each temporal scale (number of years).
for (s in c(1:length(sites))) {
  print(sites[s])
  dfsite = df_visits[ df_visits$site == sites[s], ] # choose one site
  syears = unique(dfsite$year)
      
    # loop through number of years (temporal scale). Here h=3 would be a set of 3 years
    # for each year, need to include all rounds from that year
    for (h in c(1:length(syears))) {		
      
      spreq_list = c()							# set up an empty vector for the final species list
      print(paste("h=",h,sep=" "))		# print this to keep track of simulation progress
      
      # select the subset of h years, starting with earliest
      yearsubset = syears[1:h]
      print(yearsubset)
      
      # select the years that are in this year subset (1 to 3 years):
      dfsite_ys = dfsite %>% filter(year %in% yearsubset)
      
      # include only the first rmax rounds for each year, so that all site-years will have same number of rounds
      # syrounds = unique(dfsite_yrs$round)
      # syrounds[1:rmax]
      # dfsite_yrs = dfsite_yrs %>% filter(round %in% syrounds[1:rmax])
      dfsite_yrs = data.frame()
      for (i in 1:h) {
        dfsy = dfsite_ys %>% filter(year == yearsubset[i])
        syrounds = unique(dfsy$round)
        roundset = syrounds[1:rmax]
        dfsy = dfsy %>% filter(round %in% roundset)
        dfsite_yrs = rbind(dfsite_yrs, dfsy)
      }
      
      # combine years and rounds into year_round
      dfsite_yrs <- dfsite_yrs %>% unite("year_round",c(year,round))
      
      # convert to matrix with species as rows and year_rounds as columns
      m = daply(dfsite_yrs, .(gen_sp, year_round), function(gen_sp) gen_sp$fun, .drop_o = F, .drop_i = T)
      m[is.na(m)] = 0  # fill NA with zeros
      m = adrop(m, drop=3)  # drop 3rd array dimension to get matrix
      #print(m)
      
      # define the matrix of species and years that will go into the optimizer
      func_outputB = m
      
      
      # list the round-years that cannot meet the function threshold even with all their species.  All species present at these years will be required in the final set.
      low_sites = names(which(colSums(func_outputB)<func_level))
      
      if (length(low_sites)>0) {			# are there any low round-years? If so, all species that occur there will be pre-required.
        prereq_species = names(which(rowSums(subset(func_outputB, select=low_sites))>0))
      } else {
        prereq_species = character(0)	# if not, there will be no pre-required species.
      }
      #print(prereq_species)
      
      # list all species for this crop, in no particular order
      startvector = rownames(func_outputB)
      
      # gaoptim optimizer requires >= 3 species. If there are < 3 species at the site-date, calculate min set manually.
      if (length(startvector) < 3) {   # are there < 3 species? if so we won't run optimizer, so n and popsize are N/A
        n = NA
        popSize = NA
        
        # is there only one species present? if so, the min set is 1 regardless of total function and the required species is just the only species
        if (length(startvector) == 1) {
          new_min = 1
          spreq_list = startvector
          
        # if there are 2 species present, we need to find out if 1 or both are needed to meet function threshold
        } else {
          df_req = subset(func_outputB, rownames(func_outputB) %in% prereq_species)		# species needed because they occur at years that need all their species:
          
          # are there 2 pre-required species? If so, the min set is 2 for all years and the required species are both species
          if (length(prereq_species) > 1) {
            new_min = length(startvector)
            spreq_list = startvector
            
          # if both species aren't pre-required we need to find out if one or both are needed each year
          } else {  
            if (length(prereq_species) == 1) {  
              # if there is one pre-required species we need to evaluate it for all years
              for (z in 1:h) {
                if (sum(df_req[,z]) < func_level) { # does pre-required species fulfill function threshold for this year? If not see if the additional species does
                  if (sum(func_outputB[,z]) > func_level ) {  # if both species together fulfill function threshold in this year, add both species to required list
                    df_req = func_outputB
                  }
                }
              }
              new_min = nrow(df_req)
              spreq_list = rownames(df_req)
              
            } else {  # if there are no pre-req species, we need to see how many species are needed at each site
              # start with first permutation of 2 species
              df_req1 = func_outputB[1,,drop=F]
              for (z in 1:h) {
                if (df_req1[1,z] < func_level) {  # if the first species doesn't fulfill the function threshold for any year, both species are needed
                  df_req1 = func_outputB
                }
              }
              
              # then try second permutation of 2 species
              df_req2 = func_outputB[2,,drop=F]
              for (z in 1:h) {
                if (df_req2[1,z] < func_level) {  # if the second species doesn't fulfill the function threshold for any year, both species are needed
                  df_req2 = func_outputB
                }
              }
              
              # min set is whichever new_min is smaller
              new_min = min(nrow(df_req1),nrow(df_req2))
              if (nrow(df_req1) < nrow(df_req2)) {
                spreq_list = rownames(df_req1)
              } else {
                spreq_list = rownames(df_req2)
              }
            }
          }
        }
      } else {  # if there are > 2 species, run optimizer
        
      
      # set up the genetic algorithm optimizer (uses the gaoptim package)
      # the optimizer works on the minfinder() function, which takes a permutation of the species list and returns a number of species required to satisfy the threshold
      out = GAPerm(minfinder, length(startvector), popSize = 100)
      out$evolve(30) 							# evolve for 30 generations
      
      # continue to run the optimizer 1 generation at a time until results do not change for 30 generations
     
      n = 30
      while (   min(out$bestFit()[n:(n-29)]) !=  out$bestFit()[n] ) {
        # print(n)
        # print(out$bestFit()[n])
        n = n+1
        out$evolve(1)
      }
      
      new_min = 1/out$bestFit()[n]		# convert output of minfinder() into a number of species
      
      }  # end optimizer loop
      
      print(paste("min species", new_min))
      
      # store the list of species that were pre-required
      if (length(prereq_species)==0) {
        presets = ""
      } else {
        presets = prereq_species
      }
      # print(spreq_list)
      
      # write output to csv, one line at a time.  This way, if program is stopped before complete, the output thusfar will already be saved in the csv
      if (h == 1 & s == 1) {
        write.table(data.frame("crop" = crop, "site" = sites[s], "nyears"=h, "start_year" = yearsubset[1], "minsize" = new_min, "generations"=n, 
                               "popsize"=nrow(out$population()), "threshold" = func_level, "func_percent" = func_percent, 
                               "presets" = paste(presets, collapse=","), "run_type" = "observed_data", "minset" = paste(spreq_list, collapse=",")), 
                    file=paste(outpath,"minset_finder_results_", crop, "_observed_", func_percent*100, "_years_", rmax, "round.csv", sep=""),
                    sep=",",append=F, col.names=T, row.names=F)
      } else {
        write.table(data.frame("crop" = crop, "site" = sites[s], "nyears"=h, "start_year" = yearsubset[1], "minsize" = new_min, "generations"=n, 
                               "popsize"=nrow(out$population()), "threshold" = func_level, "func_percent" = func_percent, 
                               "presets" = paste(presets, collapse=","), "run_type" = "observed_data", "minset" = paste(spreq_list, collapse=",")), 
                    file=paste(outpath,"minset_finder_results_", crop, "_observed_", func_percent*100, "_years_", rmax, "round.csv", sep=""),
                    sep=",",append=T, col.names=F, row.names=F)
      }		
      
    } #  end h years loop (h)
  
} # end s loop


## load resulting csv files as R data frames and save for future analyses in R ##

if (crop=="blue") {
  df_years_blue  = read.csv(paste(outpath,"minset_finder_results_blue_observed_",func_percent*100,"_years_", rmax, "round.csv",sep=""),
                            header=TRUE,stringsAsFactors=F)
  df_years_blue  = mutate(df_years_blue, crop=factor(crop), site=factor(site), start_year=factor(start_year)) # convert to correct data types
  save(df_years_blue, file=paste(rpath,"df_minset_years_blue_",func_percent*100,".RData",sep='')) # save as R data for future use
} else if (crop=="njwat") {
  df_years_njwat = read.csv(paste(outpath,"minset_finder_results_njwat_observed_",func_percent*100,"_years_", rmax, "round.csv",sep=""),
                            header=TRUE,stringsAsFactors=F)
  df_years_njwat = mutate(df_years_njwat, crop=factor(crop), site=factor(site), start_year=factor(start_year)) # convert to correct data types
  save(df_years_njwat, file=paste(rpath,"df_minset_years_njwat_",rmax,"round_",func_percent*100,".RData",sep='')) # save as R data for future use
} else if (crop=="cawat") {
  df_years_cawat = read.csv(paste(outpath,"minset_finder_results_cawat_observed_",func_percent*100,"_years_", rmax, "round.csv",sep=""),
                            header=TRUE,stringsAsFactors=F)
  df_years_cawat = mutate(df_years_cawat, crop=factor(crop), site=factor(site), start_year=factor(start_year))
  save(df_years_cawat, file=paste(rpath,"df_minset_years_cawat_",func_percent*100,".RData",sep='')) # save as R data for future use
} else { print("Crop not found") }
