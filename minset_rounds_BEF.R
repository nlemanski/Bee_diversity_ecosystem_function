## Find the minimum set of species needed to meet a pollination function threshold for all rounds in a given timescale ##
# repeat analysis for each site-year
# edited NJL 5/4/2021

## load packages ##
library(plyr)
library(reshape)
library(tidyverse)
library(gaoptim)
library(abind)

setwd("~/winfree lab/R code")

## load the species required function
source(file="minfinder_function_time.R")

## file paths and file names ##
path <- "D:/natal/D_Documents/winfree lab/SQL data/"
figpath <- "D:/natal/D_Documents/winfree lab/figures"
rpath <- "D:/natal/D_Documents/winfree lab/R data/"
outpath <- "D:/natal/D_Documents/winfree lab/BEFresults/"

## import crop visit data as RData file ##
crops = c("blue","cran","njwat","cawat")
cropnames = c("Blueberry", "Cranberry", "NJ Watermelon", "CA Watermelon")
# make sure crop name matches data file name
crop = "njwat"
# fname = paste(rpath,'df_visits_',crop,'.RData',sep='')
# load(fname)
# df_visits <- df_visits_cawat

if (crop == "cawat") {
  load(paste(rpath,'df_visits_cawat.RData',sep=''))  # load data
  rmax = 1 # minimum number of rounds needed for inclusion
  df_visits <- df_visits_cawat  # assign crop data to df_visits
} else if (crop == "blue") {
  load(paste(rpath,'df_visits_blue.RData',sep=''))  # load data
  rmax = 1 # minimum number of rounds needed for inclusion
  df_visits <- df_visits_blue  # assign crop data to df_visits
} else if (crop == "njwat") {
  load(paste(rpath,'df_visits_njwat.RData',sep=''))  # load data
  rmax = 1 # minimum number of rounds needed for inclusion
  df_visits <- df_visits_njwat  # assign crop data to df_visits
}

summary(df_visits)


## list unique values of each factor ##
rounds <- sort(unique(df_visits$round))
years <- sort(unique(df_visits$year))
species <- sort(unique(df_visits$gen_sp))
sites <- sort(unique(df_visits$site))

## function threshold percentage (set here for sensitivity)
func_percent = 0.50

# function threshold value (aggregate pollen grains) from percentage, using mean across all site-year-rounds
# df_totalfun = df_visits %>%     # find the total function for each site-year-round
#                 group_by(round, year, site) %>%
#                 summarize(totalfun = sum(fun))
df_totalfun = ddply(df_visits, .(round, year, site), summarize, totalfun = sum(fun))

func_level = mean(df_totalfun$totalfun, na.rm=T)*func_percent

# rmax set to 1 because we are looking at 1 round at a time
rmax = 1

# convert dataframe to a matrix with species as rows, years as columns, and function as values
# testing with one site-round
# in full analysis, this would loop through all sites and rounds
# choose a single site-round
# dfsite1 = filter(df_visits, site == "anc") # choose one site
# srounds = unique(dfsite1$round)
# dfsiteround1 = filter(dfsite1, round == srounds[2]) # choose one round
# dfsiteround1 = mutate(dfsiteround1, year=factor(year)) # drop missing years as factors
# m1 = daply(dfsiteround1, .(gen_sp, year), function(gen_sp) gen_sp$fun, .drop_o=FALSE, .drop_i=T) # convert to matrix with species as rows and years as columns
# #m1 = m1[apply(!is.na(m1),2,any),drop=FALSE]  # drop years with no data for site s, round r
# m1[is.na(m1)] = 0  # fill NA with zeros

## main program loops:
# loop through all sites. For each site-year, we will find the minimum species set at each temporal scale (number of rounds).
for (s in c(1:length(sites))) {
  print(paste("site = ",sites[s]))
  dfsite = df_visits[ df_visits$site == sites[s], ] # choose one site
  syears = unique(dfsite$year)
  # loop through all years that were sampled at a given site
  for (y in c(1:length(syears))) {   # repeat for all years measured at site s
    print(paste("year =  ", syears[y]))
    dfsiteyear = dfsite[ dfsite$year == syears[y], ] # choose one year
    dfsiteyear = mutate(dfsiteyear, round=factor(round)) # drop missing rounds as factors
    syrounds = unique(dfsiteyear$round)
    
    m = daply(dfsiteyear, .(gen_sp, round), function(gen_sp) gen_sp$fun, .drop_o = F, .drop_i = T) # convert to matrix with species as rows and rounds as columns
    m[is.na(m)] = 0  # fill NA with zeros
    m = adrop(m, drop=3)  # drop 3rd array dimension to get matrix
    #print(m)

    # define the matrix of species and years that will go into the optimizer
    func_output = m
    
    # loop through number of rounds (temporal scale). Here h=3 would be a set of 3 sampling rounds
    for (h in c(1:length(syrounds))) {			# warning, it will take quite a while if you run everything
      
      spreq_list = c()							# set up an empty vector for the final species list
      print(paste("h=",h,sep=" "))		# print this to keep track of simulation progress
      
      # subset of h rounds, starting with earliest
      roundsubset = syrounds[1:h]
      print(roundsubset)
      
      # select the columns for the rounds that are in this round subset (1 to 9 rounds):
      func_outputB = func_output[,colnames(func_output) %in% roundsubset, drop=FALSE]
      
      # list the rounds that cannot meet the function threshold even with all their species.  All species present at these rounds will be required in the final set.
      low_sites = names(which(colSums(func_outputB)<func_level))
      
      
      if (length(low_sites)>0) {			# are there any low sites? If so, all species that occur there will be pre-required.
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
        print("< 3 species at site-date")
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
        
      print("running optimizer")
      # set up the genetic algorithm optimizer (uses the gaoptim package)
      # the optimizer works on the minfinder() function, which takes a permutation of the species list and returns a number of species required to satisfy the threshold
      out = GAPerm(minfinder, length(startvector), popSize = 100)
      out$evolve(30) 							# evolve for 30 generations
      
      # continue to run the optimizer 1 generation at a time until results do not change for 30 generations
      # while (   min(out$bestFit()[summary(out)$n:(summary(out)$n-29)]) !=  out$bestFit()[summary(out)$n] ) {
      #   out$evolve(1)
      # }  # Error: $ operator is invalid for atomic vectors
      
      n = 30
      while (   min(out$bestFit()[n:(n-29)]) !=  out$bestFit()[n] ) {
        print(n)
        print(out$bestFit()[n])
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
      print(spreq_list)
      
      # write output to csv, one line at a time.  This way, if program is stopped before complete, the output thusfar will already be saved in the csv
      if (h == 1 & s == 1 & y == 1) {
        write.table(data.frame("crop" = crop, "site" = sites[s], "year" = syears[y], "nrounds"=h, "start_round" = roundsubset[1], "minsize" = new_min, "generations"=n, "popsize"=nrow(out$population()), "threshold" = func_level, "func_percent" = func_percent, "presets" = paste(presets, collapse=","), "run_type" = "observed_data", "minset" = paste(spreq_list, collapse=",")), file=paste(outpath,"minset_finder_results_", crop, "_observed_", func_percent*100, "_rounds",".csv", sep=""),sep=",",append=F, col.names=T, row.names=F)
      } else {
        write.table(data.frame("crop" = crop, "site" = sites[s], "year" = syears[y], "nrounds"=h, "start_round" = roundsubset[1], "minsize" = new_min, "generations"=n, "popsize"=nrow(out$population()), "threshold" = func_level, "func_percent" = func_percent, "presets" = paste(presets, collapse=","), "run_type" = "observed_data", "minset" = paste(spreq_list, collapse=",")), file=paste(outpath,"minset_finder_results_", crop, "_observed_", func_percent*100, "_rounds",".csv", sep=""),sep=",",append=T, col.names=F, row.names=F)
      }		
      
    } #  end h years loop (h)
    
  } # end r loop
  
} # end s loop


## load resulting csv files as R data frames and save for future analyses in R ##
fname1 = "minset_finder_results_blue_observed_50_rounds.csv"
fname2 = "minset_finder_results_cran_observed_50_rounds.csv"
fname3 = "minset_finder_results_njwat_observed_50_rounds.csv"
fname4 = "minset_finder_results_cawat_observed_50_rounds.csv"

df_rounds_blue  = read.csv(paste(outpath,fname1,sep=""),header=TRUE,stringsAsFactors=F)
df_rounds_cran  = read.csv(paste(outpath,fname2,sep=""),header=TRUE,stringsAsFactors=F)
df_rounds_njwat = read.csv(paste(outpath,fname3,sep=""),header=TRUE,stringsAsFactors=F)
df_rounds_cawat = read.csv(paste(outpath,fname4,sep=""),header=TRUE,stringsAsFactors=F)

# view data frames
glimpse(df_rounds_blue)
glimpse(df_rounds_cran)
glimpse(df_rounds_njwat)
glimpse(df_rounds_cawat)

## convert to correct data types ##
df_rounds_blue  = mutate(df_rounds_blue, crop=factor(crop), site=factor(site), year=factor(year), start_round=factor(start_round))
df_rounds_cran  = mutate(df_rounds_cran, crop=factor(crop), site=factor(site), year=factor(year), start_round=factor(start_round))
df_rounds_njwat = mutate(df_rounds_njwat, crop=factor(crop), site=factor(site), year=factor(year), start_round=factor(start_round))
df_rounds_cawat = mutate(df_rounds_cawat, crop=factor(crop), site=factor(site), year=factor(year), start_round=factor(start_round))
summary(df_rounds_cawat)

## save as R data for future use ##
save(df_rounds_blue, file=paste(rpath,"df_minset_rounds_blue.RData",sep=''))
save(df_rounds_cran, file=paste(rpath,"df_minset_rounds_cran.RData",sep=''))
save(df_rounds_njwat, file=paste(rpath,"df_minset_rounds_njwat.RData",sep=''))
save(df_rounds_cawat, file=paste(rpath,"df_minset_rounds_cawat.RData",sep=''))

