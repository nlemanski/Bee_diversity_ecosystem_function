## null model to separate phenology effects from sampling effects in BEF over time analysis for crop data sets
# calculate size of minimum species set needed to meet 50% pollination threshold for a series of subsamples taken from the same day or from different days
# change in minimum set w/ number of samples taken from same day includes sampling effects only
# change in minimum set w/ number of samples taken from different days includes both sampling effects and phenological turnover

## load packages ##
library(plyr)
library(reshape)
library(tidyverse)
library(gaoptim)
library(abind)
library(cowplot)
library(lme4)
library(car)
library(performance)
library(effects)

setwd("C:/Documents/Bee_diversity_ecosystem_function/R code")

## load the species required function
source(file="minfinder_function_time.R")

## file paths and file names ##
path <- "C:/Documents/Bee_diversity_ecosystem_function/SQL data/"
figpath <- "C:/Documents/Bee_diversity_ecosystem_function/figures/"
rpath <- "C:/Documents/Bee_diversity_ecosystem_function/R data/"
outpath <- "C:/Documents/Bee_diversity_ecosystem_function/BEF null model results/"
altpath <- "C:/Documents/Bee_diversity_ecosystem_function/alt null model results/"

## define crops and select crop to analyze ##
crops = c("blue","njwat","cawat")
cropnames = c("Blueberry", "Eastern watermelon", "Western watermelon")
croplabs = data.frame(row.names=crops, val=cropnames)

## select crop, make sure crop name matches data file name
crop = 'blue' # 'njwat' # 'cawat' # 

## set sample parameters and set random seed ##
func_percent = 0.50 # function threshold percentage (set here for sensitivity)
nsize = 30  # set size of samples
nreps = 100 # number of reps
set.seed(42)

sampletype = 'replace'    # uncomment for sampling with replacement (bootstrapping)
# sampletype = 'noreplace' # uncomment for sampling without replacement
if (sampletype == 'replace') {repl = T} else {repl = F}

threshtype = 'same' # uncomment for same threshold between null and observed
# threshtype = ''  # uncomment for different thresholds between null and observed

# nulltype = 'global'  # uncomment for alternate null model, where samples are drawn from global pool of all days
nulltype = 'single'  # uncomment for null model where samples are drawn from a single day

## import crop visit data as RData file ##
#fname = paste(rpath,'obs_',crop,'_',x,'.Rdata',sep='')
fname = paste(rpath,'obs_',crop,'.Rdata',sep='')
load(fname)

# select crop
if (crop == 'cawat') {
  df_obs = df_obs_cawat
} else if (crop == "njwat") {
    df_obs = df_obs_njwat
} else if (crop == "blue") {
    df_obs = df_obs_blue %>% mutate(SV_pollen = SV_pollen_tubes)
} else  {
    print('crop not found')
}

# add column for site-year
df_obs = df_obs %>% mutate(siteyear = paste(site,'-',year,sep='')) # make new column for site-year

# minimum set analysis for null and observed, across different numbers of rounds ------------------------------

rmax = 1 # for within-years analysis, add 1 round at a time

## select only site-years with at least 3 or 6 rounds ##
if (crop == 'cawat') {  # for western watermelon, find site-years with >= 6 rounds
  # df_r6 = df_obs %>% filter(round == '6')
  df_r6 = df_obs %>%
    group_by(siteyear) %>% 
    summarize(nrounds = length(unique(round))) %>% 
    filter(nrounds >= 6) 
} else {               # for eastern wat and blue, find site-years with >= 3 rounds
  # df_r6 = df_obs %>% filter(round == '3') 
  df_r6 = df_obs %>%
    group_by(siteyear) %>% 
    summarize(nrounds = length(unique(round))) %>% 
    filter(nrounds >= 3) 
}
siteyears = unique(df_r6$siteyear) # list all site-years for the crop
df_obs_rounds = df_obs %>% filter(siteyear %in% siteyears)

## run null model analysis for all site-years with enough sampling rounds and save results ##
for (sy in siteyears) {
  df_site = df_obs_rounds %>% filter(siteyear == sy)  # select data from a single site-year
  sitex = df_site$site[1]
  yearx = df_site$year[1]
  
  # find min abundance, max abundance, and round with max and min abundance
  df_sum <- df_site %>% group_by(round) %>% summarise(abundance = length(gen_sp))
  minabund = min(df_sum$abundance)
  maxabund = max(df_sum$abundance)
  maxround = slice(arrange(df_sum, desc(abundance)),1)$round
  minround = slice(arrange(df_sum, abundance),1)$round
  
  # filter out rounds with abundance < nsize if sampling without replacement
  if (sampletype == 'noreplace') {
    df_site = df_site %>% filter(!round %in% (filter(df_sum, abundance < nsize)$round))
  }
  rounds = unique(df_site$round)
  
  # set number of samples equal to number of rounds included
  if (sampletype == 'noreplace') { # if sampling w/o replacement, max number of samples possible is total abundance/sample size
    nsamples = min(length(unique(df_site$round)), round(maxabund/nsize))
  } else {
    nsamples = length(unique(df_site$round)) # if bootstrapping, set number of null samples equal to max number of rounds at site-date
  }
  ## find min sets for both run types: observed (with phenology) and null (without phenology) ##
  runtypes = c("observed", "null")
  
  ## for each rep, take a random sub-sample of the data and calculate minimum set sizes for each number of rounds ##
  for (rep in 1:nreps) {
    print(paste("rep =", rep))

    for (runtype in runtypes) {
      ## choose X random samples of N individuals from different rounds (observed) or same round (null) for minimum set analysis ##
      df_sample_obs <- df_site %>% slice(0) %>% mutate(sround = round) # initialize empty data frame to hold sampled data
      if (runtype == "observed") {  # sub-sample data from each round for observed
        for (s in 1:nsamples) {
          df_sample_s <- df_site %>% filter(round == rounds[s]) %>% sample_n(size = nsize, replace = repl) %>% mutate(sround = factor(s))
          df_sample_obs <- df_sample_obs %>% bind_rows(df_sample_s)  # add sub-sampled data to df
        }
      } else if (runtype == "null") {  
        if (nulltype == "single") {  # sub-sample data from the first round for original null model
          for (s in 1:nsamples) {
            df_sample_s <- df_site %>% filter(round == rounds[1]) %>% sample_n(size = nsize, replace = repl) %>% mutate(sround = factor(s))
            df_sample_obs <- df_sample_obs %>% bind_rows(df_sample_s)  # add sub-sampled data to df
          }
        } else if (nulltype == 'global') { # sub-sample data from all rounds for global null model
          for (s in 1:nsamples) {
            df_sample_s <- df_site %>% sample_n(size = nsize, replace = repl) %>% mutate(sround = factor(s))
            df_sample_obs <- df_sample_obs %>% bind_rows(df_sample_s)  # add sub-sampled data to df
          }
        } else { print("Invalid null type") }
        
      } else {
        print("Invalid run type")
      }

      ## find number of visits and total pollination by each species at each round
      if (crop == 'cawat') {
        df_visits = df_sample_obs %>%
          group_by(sround, gen_sp, sex) %>%
          summarize(visits = length(gen_sp), SV_pollen = mean(SV_pollen)) %>%
          mutate(fun = visits*SV_pollen) %>%
          group_by(sround, gen_sp) %>%
          summarize(visits = sum(visits), fun = sum(fun))
      } else {
        df_visits = df_sample_obs %>%
          group_by(sround, gen_sp) %>%
          summarize(visits = length(gen_sp), SV_pollen = mean(SV_pollen)) %>%
          mutate(fun = visits*SV_pollen)
      }

      ## calculate the function threshold, i.e. 50% of mean function across all srounds
      if (threshtype == 'same') {  # if using same threshold, only calculate func_level for observed data, then keep it for null data
        if (runtype == 'observed') {
          meanfun = df_visits %>% group_by(sround) %>% summarize(totfun = sum(fun)) %>% summarize(meanfun = mean(totfun))
          func_level = meanfun$meanfun * func_percent }
      } else {  # if using different threshold, calculate func_level for both observed and null data
        meanfun = df_visits %>% group_by(sround) %>% summarize(totfun = sum(fun)) %>% summarize(meanfun = mean(totfun))
        func_level = meanfun$meanfun * func_percent
      }

      #print(paste("threshold:",func_level))

      ## find min set needed for 50% function for each number of srounds ##
      df_visits = mutate(df_visits, sround = factor(sround)) # drop missing rounds as factors
      srounds = unique(df_visits$sround) # list of sampled rounds

      # convert df to matrix with species as rows and rounds as columns
      m = daply(df_visits, .(gen_sp, sround), function(gen_sp) gen_sp$fun, .drop_o = F, .drop_i = T)
      m[is.na(m)] = 0  # fill NA with zeros
      m = adrop(m, drop=3)  # drop 3rd array dimension to get matrix

      # define the matrix of species and years that will go into the optimizer
      func_output = m

      # loop through number of rounds (temporal scale). Here h=3 would be a set of 3 sampling rounds
      for (h in c(1:length(srounds))) {

        spreq_list = c()							# set up an empty vector for the final species lists
        min_list = c()                # create empty vector for the minimum set size values found by the algorithm
        print(paste("h =",h,sep=" "))		# print this to keep track of simulation progress

        # subset of h rounds, starting with earliest
        roundsubset = srounds[1:h]

        # select the columns for the rounds that are in this round subset (1 to 9 rounds):
        func_outputB = func_output[,colnames(func_output) %in% roundsubset, drop=FALSE]

        # list the rounds that cannot meet the function threshold even with all their species.  All species present at these rounds will be required in the final set.
        low_sites = names(which(colSums(func_outputB)<func_level))

        if (length(low_sites)>0) {			# are there any low sites? If so, all species that occur there will be pre-required.
          prereq_species = names(which(rowSums(subset(func_outputB, select=low_sites))>0))
        } else {
          prereq_species = character(0)	# if not, there will be no pre-required species.
        }
        # print(prereq_species)

        # list all species for this crop, in no particular order
        startvector = rownames(func_outputB)

        # find total number of species observed at any round in this site-year
        rich = length(startvector)

        # gaoptim optimizer requires >= 3 species. If there are < 3 species at the site-date, calculate min set manually.
        if (length(startvector) < 3) {   # are there < 3 species? if so we won't run optimizer, so n and popsize are N/A
          n = NA
          popSize = NA
          print("< 3 species at site-date")
          # is there only one species present? if so, the min set is 1 regardless of total function and the required species is just the only species
          if (length(startvector) == 1) {
            new_min = 1
            #spreq_list = startvector
            speclist = paste(startvector, collapse = ',')

            # if there are 2 species present, we need to find out if 1 or both are needed to meet function threshold
          } else {
            df_req = subset(func_outputB, rownames(func_outputB) %in% prereq_species)		# species needed because they occur at rounds that need all their species:

            # are there 2 pre-required species? If so, the min set is 2 for all rounds and the required species are both species
            if (length(prereq_species) > 1) {
              new_min = length(startvector)
              speclist = paste(startvector, collapse=',')

              # if both species aren't pre-required we need to find out if one or both are needed each round
            } else {
              if (length(prereq_species) == 1) {
                # if there is one pre-required species we need to evaluate it for all rounds
                for (z in 1:h) {
                  if (sum(df_req[,z]) < func_level) { # does pre-required species fulfill function threshold for this year? If not see if the additional species does
                    if (sum(func_outputB[,z]) > func_level ) {  # if both species together fulfill function threshold in this year, add both species to required list
                      df_req = func_outputB
                    }
                  }
                }
                new_min = nrow(df_req)
                speclist = paste(rownames(df_req), collapse=',')

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
                  speclist = paste(rownames(df_req1), collapse = ',')
                } else {
                  speclist = paste(rownames(df_req2), collapse = ',')
                }
              }
            }
          }
        } else {  # if there are > 2 species, run optimizer

          # print("running optimizer")
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

          # find the required species list for the optimized minimum set found by genetic algorithm
          lastgen = spreq_list[(length(spreq_list)-99):length(spreq_list)] # min species list from the nth generation of genetic algorithm
          lastmin = min_list[(length(min_list)-99):length(min_list)] # min set size from the nth generation of genetic algorithm
          minsize2 = 1/max(lastmin) # the smallest min set size from nth generation, should always match new_min
          speclist = lastgen[ which.max(lastmin) ] # the required species list from the smallest min set size from nth generation
          #print(paste("min set:",speclist))

        }  # end optimizer loop

        # print(paste("min species", new_min))

        # store the list of species that were pre-required
        if (length(prereq_species)==0) {
          presets = ""
        } else {
          presets = prereq_species
        }
        #print(last(spreq_list))
        #print(1/last(min_list))

        # write output to csv, one line at a time
        if (nulltype == 'single') { 
          pth = outpath 
        } else if (nulltype == 'global') { 
          pth = altpath
        } else { 
          print('invalid null type') 
        }
        
        resname = paste(pth,"minset_subsample_results_", crop, "_", runtype, "_", func_percent*100, "_rounds_",
                        sitex,yearx,'_',"n",nsize,'_',sampletype,"_",threshtype,".csv", sep="")
        
        if (h == 1 & rep == 1) {
          write.table(data.frame("crop" = crop, "site" = sitex, "year" = yearx, "nsize" = nsize, "rep" = rep, "nrounds"= h,
                                 "start_round" = roundsubset[1], "minsize" = new_min, "richness"=rich, "generations"=n, "popsize"=nrow(out$population()),
                                 "threshold" = func_level, "func_percent" = func_percent, "presets" = paste(presets, collapse=","),
                                 "run_type" = runtype, "minset" = speclist),
                      file=resname, sep=",",append=F, col.names=T, row.names=F)
        } else {
          write.table(data.frame("crop" = crop, "site" = sitex, "year" = yearx, "nsize" = nsize, "rep" = rep, "nrounds"= h,
                                 "start_round" = roundsubset[1], "minsize" = new_min, "richness"=rich, "generations"=n, "popsize"=nrow(out$population()),
                                 "threshold" = func_level, "func_percent" = func_percent, "presets" = paste(presets, collapse=","),
                                 "run_type" = runtype, "minset" = speclist),
                      file=resname, sep=",",append=T, col.names=F, row.names=F)
        }

      } #  end h loop (number of rounds)
    } # end runtype loop
  } # end rep loop
  
  ## import results as rdata ##
  df_null = read.csv(file=paste(pth,"minset_subsample_results_", crop, "_", "null", "_", func_percent*100, "_rounds_",sitex,yearx,'_',"n",nsize,"_",sampletype,"_",threshtype,".csv", sep=""),
                     stringsAsFactors = T)
  df_obsv = read.csv(file=paste(pth,"minset_subsample_results_", crop, "_", "observed", "_", func_percent*100, "_rounds_",sitex,yearx,'_',"n",nsize,"_",sampletype,"_",threshtype,".csv", sep=""),
                     stringsAsFactors = T)
  
  #save(df_null, file=paste(pth,"minset_subsample_results_", crop, "_", "null", "_", func_percent*100, "_rounds_",sitex,yearx,'_',"n",nsize,'_x',x,"_replace.Rdata", sep=""))
  #save(df_obsv, file=paste(pth,"minset_subsample_results_", crop, "_", "observed", "_", func_percent*100, "_rounds_",sitex,yearx,'_',"n",nsize,'_x',x,"_replace.csv.Rdata", sep=""))
  
  ## plotting results ##
  # df_all_sum = df_null %>% rbind(df_obsv) %>% 
  #   group_by(nrounds, run_type) %>% 
  #   summarize(minsize = mean(minsize), threshold = mean(threshold))
  # df_all = rbind(df_null,df_obsv) # combine null and observed into one df
  # 
  # # fit a linear model and find intercept
  # lm1 = lm(minsize ~ nrounds*run_type, data=df_all)
  # 
  # # plot null vs observed
  # ggplot(data=df_all, mapping=aes(x=nrounds, y=minsize, color=run_type)) +
  #   geom_smooth() +
  #   geom_jitter(width = 0.2, height = 0.2) +
  #   geom_abline(slope=0, intercept=lm1$coefficients[1]) +
  #   # labs(x="Number of days", y="Minimum set size", color="Run type", title=paste('Simulated data - No turnover, n = ',nsize,', mean abundance =',mean(df_abund$abundance)))+ #"Site-Year",sitex,yearx)) +
  #   labs(x="Number of days", y="Minimum set size", color="Run type", title=paste(croplabs[crop,],", Site-Year = ",sitex,yearx,', n = ',nsize, ', sampled with replacement',sep='')) +
  #   #ylim(0.2,6.2) +
  #   theme_light()
  # figname = paste(figpath,'minset_subsample_results_',crop,'_rounds_',sitex,yearx,'_',"n",nsize,'_',sampletype,threshtype,'.png', sep='')
  # ggsave(filename=figname, width=7,height=5,units='in')
  # 
} # end site-year loop


# minimum set analysis for null and observed across different number of years --------

# drop data from 2004
df_obs <- df_obs %>% filter(year != '2004')

# number of rounds to include per year
if (crop == 'njwat') { rmax = 1 } else { rmax = 3 } # use only 1 round for njwat b/c there are only 3 years w/ all 3 rounds and we want 6 years

## select only site-years with at least 3 rounds, or 1 round for nj watermelon ##
df_r3 = df_obs %>% 
  group_by(siteyear) %>% 
  summarize(nrounds = length(unique(round))) %>% 
  filter(nrounds >= rmax) # find site-years with at least 3 rounds
siteyears = unique(df_r3$siteyear) # list all site-years for the crop w/ at least 3 rounds
# dfr3y3 = df_r3 %>% group_by(site) %>% summarize(nyears = length(unique(year)))
# sites = dfr3y3 %>% filter(nyears >= 3)
df_obs_years = df_obs %>% filter(siteyear %in% siteyears)

sites = unique(df_obs_years$site) # list all sites
# sites = sites[7:length(sites)]

## run null model analysis for all sites with enough sampling rounds and save results ##
for (sitex in sites) {
  df_site = df_obs_years %>% filter(site == sitex) # select single site
  print(paste('site =',sitex))
  # find min abundance, max abundance, and round with max and min abundance
  df_sum <- df_site %>% group_by(year,round) %>% summarise(abundance = length(gen_sp))
  # minabund = min(df_sum$abundance)
  # maxabund = max(df_sum$abundance)
  # maxround = slice(arrange(df_sum, desc(abundance)),1)$round
  # minround = slice(arrange(df_sum, abundance),1)$round
  
  # filter out rounds with abundance < nsize if sampling without replacement
  if (sampletype == 'noreplace') {
    df_site = df_site %>% filter(!round %in% (filter(df_sum, abundance < nsize)$round))
  }
  rounds = sort(unique(df_site$round))
  years = sort(unique(df_site$year))
  
  # set number of samples equal to number of years included
  if (sampletype == 'noreplace') {
    nsamples = min(length(unique(df_site$year)), round(maxabund/nsize))
  } else {
    nsamples = length(unique(df_site$year))
  }
  ## find min sets for both run types: observed (with phenology) and null (without phenology) ##
  runtypes = c("observed", "null")
  
  ## for each rep, take a random sub-sample of the observed data and calculate minimum set sizes for each number of years ##
  for (rep in 1:nreps) {
    for (runtype in runtypes) {
      print(paste('site =', sitex," rep =", rep, " runtype =", runtype))
      ## choose 3*X random samples of N individuals from different years (observed) or same year (null) for minimum set analysis ##
      df_sample_obs <- df_site %>% slice(0) %>% mutate(syear = factor(0)) # initialize empty data frame to hold sampled data
      if (runtype == "observed") {  # sub-sample data from 3 rounds of each year
        for (s in 1:nsamples) {
          dfsys = df_site %>% filter(year == years[s])
          yrrounds = unique(dfsys$round)
          for (r in 1:rmax) {
            df_sample_s <- df_site %>% 
                           filter(year == years[s], round == yrrounds[r]) %>% 
                           sample_n(size = nsize, replace = repl) %>% 
                           mutate(syear = factor(s))
            df_sample_obs <- df_sample_obs %>% bind_rows(df_sample_s)  # add sub-sampled data to df
          }
        }
      } else if (runtype == "null") {  # sub-sample data from 3 rounds of the first year (original null) or global pool (alt null)
        for (s in 1:nsamples) {
          dfsys = df_site %>% filter(year == years[1])
          yrrounds = unique(dfsys$round)
          for (r in 1:rmax) {
            if (nulltype == 'single') {
              df_sample_s <- df_site %>% 
                             filter(year == years[1], round == yrrounds[r]) %>% 
                             sample_n(size = nsize, replace = repl) %>% 
                             mutate(syear = factor(s))
            } else if (nulltype == 'global') {
              df_sample_s <- df_site %>% 
                             filter(round == rounds[r]) %>% 
                             sample_n(size = nsize, replace = repl) %>% 
                             mutate(syear = factor(s))
            }
            df_sample_obs <- df_sample_obs %>% bind_rows(df_sample_s)  # add sub-sampled data to df
          }
        }
      } else {
        print("Invalid run type")
      }
      
      ## find number of visits and total pollination by each species at each round-year
      if (crop == 'cawat') {
        df_visits = df_sample_obs %>%
          group_by(syear, round, gen_sp, sex) %>%
          summarize(visits = length(gen_sp), SV_pollen = mean(SV_pollen)) %>%
          mutate(fun = visits*SV_pollen) %>%
          group_by(syear, round, gen_sp) %>%
          summarize(visits = sum(visits), fun = sum(fun))
      } else {
        df_visits = df_sample_obs %>%
          group_by(syear, round, gen_sp) %>%
          summarize(visits = length(gen_sp), SV_pollen = mean(SV_pollen)) %>%
          mutate(fun = visits*SV_pollen)
      }
      
      ## calculate the function threshold, i.e. 50% of mean function across all round-years
      if (threshtype == 'same') {  # if using same threshold, only calculate func_level for observed data, then keep it for null data
        if (runtype == 'observed') {
          totfun = df_visits %>% group_by(syear, round) %>% summarize(totfun = sum(fun))
          meanfun = mean(totfun$totfun)
          func_level = meanfun * func_percent 
          }
      } else {  # if using different threshold, calculate func_level for both observed and null data
        totfun = df_visits %>% group_by(syear, round) %>% summarize(totfun = sum(fun))
        meanfun = mean(totfun$totfun)
        func_level = meanfun * func_percent
      }
      
      ## find min set needed for 50% function for all 3 rounds for each number of syears ##
      
      df_visits = mutate(df_visits, syear = factor(syear)) # drop missing years as factors
      syears = unique(df_visits$syear) # list of sampled years
      
      # # convert df to matrix with species as rows and round-years as columns
      # m = daply(df_visits, .(gen_sp, sround), function(gen_sp) gen_sp$fun, .drop_o = F, .drop_i = T)
      # m[is.na(m)] = 0  # fill NA with zeros
      # m = adrop(m, drop=3)  # drop 3rd array dimension to get matrix
      # 
      # # define the matrix of species and years that will go into the optimizer
      # func_output = m
      
      ## loop through number of years ##
      # for each year, need to include all rounds from that year
      for (h in c(1:length(syears))) {		
        
        spreq_list = c()							# set up an empty vector for the final species list
        min_list = c()                # create empty vector for the minimum set size values found by the algorithm
        print(paste("h =",h,sep=" "))		# print this to keep track of simulation progress
        
        # select the subset of h years, starting with earliest
        yearsubset = syears[1:h]
        # print(yearsubset)
        
        # select the years that are in this year subset
        dfsite_ys = df_visits %>% filter(syear %in% yearsubset)
        
        # include only the first rmax rounds for each year, so that all site-years will have same number of rounds
        # syrounds = unique(dfsite_yrs$round)
        # syrounds[1:rmax]
        # dfsite_yrs = dfsite_yrs %>% filter(round %in% syrounds[1:rmax])
        dfsite_yrs = data.frame()
        for (i in 1:h) {
          dfsy = dfsite_ys %>% filter(syear == yearsubset[i])
          syrounds = unique(dfsy$round)
          roundset = syrounds[1:rmax]
          dfsy = dfsy %>% filter(round %in% roundset)
          dfsite_yrs = rbind(dfsite_yrs, dfsy)
        }
        
        # combine syears and rounds into year_round
        dfsite_yrs <- dfsite_yrs %>% unite("year_round",c(syear,round))
        
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
        rich = length(startvector)
        
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
          
          # continue to run the optimizer 1 generation at a time until results do not change for 30 generation          
          n = 30
          while (   min(out$bestFit()[n:(n-29)]) !=  out$bestFit()[n] ) {

            n = n+1
            out$evolve(1)
          }
          
          new_min = 1/out$bestFit()[n]		# convert output of minfinder() into a number of species
          
          # find the required species list for the optimized minimum set found by genetic algorithm
          lastgen = spreq_list[(length(spreq_list)-99):length(spreq_list)] # min species list from the nth generation of genetic algorithm
          lastmin = min_list[(length(min_list)-99):length(min_list)] # min set size from the nth generation of genetic algorithm
          minsize2 = 1/max(lastmin) # the smallest min set size from nth generation, should always match new_min
          speclist = lastgen[ which.max(lastmin) ] # the required species list from the smallest min set size from nth generation

        }  # end optimizer loop
        
        # print(paste("min species", new_min))
        
        # store the list of species that were pre-required
        if (length(prereq_species)==0) {
          presets = ""
        } else {
          presets = prereq_species
        }

        # write output to csv, one line at a time.  This way, if program is stopped before complete, the output thusfar will already be saved in the csv
        if (nulltype == 'single') { pth = outpath } else { pth = altpath }
        resname = paste(pth,"minset_subsample_results_", crop, "_", runtype, "_", func_percent*100, "_years_",
                        sitex,'_',"n",nsize,'_',sampletype,"_",threshtype,".csv", sep="")
        
        if (h == 1 & rep == 1) {
          write.table(data.frame("crop" = crop, "site" = sitex, "nsize" = nsize, "rep" = rep, "nyears"=h, "start_year" = years[1],
                                 "minsize" = new_min, "richness"=rich, "generations"=n, "popsize"=nrow(out$population()),
                                 "threshold" = func_level, "func_percent" = func_percent, "presets" = paste(presets, collapse=","), 
                                 "run_type" = runtype, "minset" = speclist), 
                      file=resname, sep=",",append=F, col.names=T, row.names=F)
        } else {
          write.table(data.frame("crop" = crop, "site" = sitex, "nsize" = nsize, "rep" = rep, "nyears"=h, "start_year" = years[1],
                                 "minsize" = new_min, "richness"=rich, "generations"=n, "popsize"=nrow(out$population()), 
                                 "threshold" = func_level, "func_percent" = func_percent, "presets" = paste(presets, collapse=","), 
                                 "run_type" = runtype, "minset" = speclist), 
                      file=resname, sep=",",append=T, col.names=F, row.names=F)
        }		
        
      } #  end h years loop (h)
      
    } # end run type loop
  } # end rep loop
} # end site-year loop

# import null model results -----------------------------------------------
if (nulltype == 'global') {pth = altpath} else {pth = outpath}

## for minsize vs number of rounds, combine all site-years into a single df, then find the mean minimum set size across all reps
# make empty dataframe to hold null model results (number of rounds)
colsr = c("crop","site","year","nsize","rep","nrounds","start_round","minsize","richness","generations",
         "popsize","threshold","func_percent","presets","run_type","minset")
df_allsitesr = data.frame(matrix(nrow = 0, ncol = length(colsr)))
colnames(df_allsitesr) = colsr

# import results for all three crops and add to single df
for (crop in crops) {
  # list all site-years for the crop that were sampled >= 3 or 6 rounds
  load(paste(rpath,'obs_',crop,'.Rdata',sep=''))
  if (crop == "cawat") {
    df_res = df_obs_cawat %>% mutate(siteyear = paste(site,'-',year,sep=''))
    roundmin = 6
  } else if (crop == "njwat") {
    df_res = df_obs_njwat %>% mutate(siteyear = paste(site,'-',year,sep=''))
    roundmin = 3
  } else if (crop == "blue") {
    df_res = df_obs_blue %>% mutate(siteyear = paste(site,'-',year,sep=''), SV_pollen = SV_pollen_tubes)
    roundmin = 3
  } else  {
    print('crop not found')
  }
  df_r = df_res %>% # find site-years with at least 3 or 6 rounds
    group_by(siteyear) %>% 
    summarize(nrounds = length(unique(round))) %>% 
    filter(nrounds >= roundmin) 
  siteyears = unique(df_r$siteyear)
  
  for (sy in siteyears) {
    sitex = unlist(strsplit(sy,'-'))[1]
    yearx = unlist(strsplit(sy,'-'))[2]
    df_nullr = read.csv(file=paste(pth,"minset_subsample_results_", crop, "_", "null", "_", func_percent*100, "_rounds_",sitex,yearx,'_',"n",nsize,"_",sampletype,"_",threshtype,".csv", sep=""),
                       stringsAsFactors = T)
    df_obsvr = read.csv(file=paste(pth,"minset_subsample_results_", crop, "_", "observed", "_", func_percent*100, "_rounds_",sitex,yearx,'_',"n",nsize,"_",sampletype,"_",threshtype,".csv", sep=""),
                       stringsAsFactors = T)
    df_allr = rbind(df_nullr, df_obsvr)
    df_allsitesr = rbind(df_allsitesr, df_allr)
  }
}
# find the mean across all reps for minsize vs number of rounds
df_sitemeans = df_allsitesr %>% 
  group_by(crop, site, year, nsize, nrounds, run_type, func_percent) %>%
  summarise(minsize = mean(minsize), richness = mean(richness), threshold = mean(threshold))
df_sitemeans = df_sitemeans %>% mutate(minsizeint = round(minsize)) # nrounds results, averaged across all reps

## for minsize vs number of years, combine all sites into a single df, then find the mean minimum set size across all reps
# make empty dataframe to hold null model results (number of years)
colsy = c("crop","site","nsize","rep","nyears","start_year","minsize","richness","generations",
          "popsize","threshold","func_percent","presets","run_type","minset")
df_allsitesy = data.frame(matrix(nrow = 0, ncol = length(colsy)))
colnames(df_allsitesy) = colsy

# import results for all three crops and add to single df for minsize vs number of years
for (crop in crops) {
  # filter all site-years for the crop that were sampled >= 3 rounds
  load(paste(rpath,'obs_',crop,'.Rdata',sep=''))
  if (crop == "cawat") {
    df_res = df_obs_cawat %>% mutate(siteyear = paste(site,'-',year,sep=''))
    roundmin = 3
  } else if (crop == "njwat") {
    df_res = df_obs_njwat %>% mutate(siteyear = paste(site,'-',year,sep=''))
    roundmin = 3
  } else if (crop == "blue") {
    df_res = df_obs_blue %>% mutate(siteyear = paste(site,'-',year,sep=''), SV_pollen = SV_pollen_tubes)
    roundmin = 3
  } else  {
    print('crop not found')
  }
  df_r = df_res %>% # find site-years with at least 3 rounds
    group_by(siteyear, site) %>% 
    summarize(nrounds = length(unique(round))) %>% 
    filter(nrounds >= roundmin)
  sitelist = unique(df_r$site) # list all sites w/ >=3 rounds
  
  for (site in sitelist) {
    sitex = site
    df_nully = read.csv(file=paste(pth,"minset_subsample_results_", crop, "_", "null", "_", func_percent*100, "_years_",sitex,'_',"n",nsize,"_",sampletype,"_",threshtype,".csv", sep=""),
                        stringsAsFactors = T)
    df_obsvy = read.csv(file=paste(pth,"minset_subsample_results_", crop, "_", "observed", "_", func_percent*100, "_years_",sitex,'_',"n",nsize,"_",sampletype,"_",threshtype,".csv", sep=""),
                        stringsAsFactors = T)
    df_ally = rbind(df_nully, df_obsvy)
    df_allsitesy = rbind(df_allsitesy, df_ally)
  }
}
# find the mean across all reps for minsize vs number of years
df_sitemeansy = df_allsitesy %>% 
  group_by(crop, site, nsize, nyears, run_type, func_percent) %>%
  summarise(minsize = mean(minsize), richness = mean(richness), threshold = mean(threshold))
df_sitemeansy = df_sitemeansy %>% mutate(minsizeint = round(minsize))


# glmm analysis for minset vs number of rounds------------------------------------------------------------
## stats for minset vs nrounds results blueberry
crop='blue'
df_blue = df_allsitesr %>% filter(crop == 'blue')

hist(df_blue$minsize)
qqnorm(df_blue$minsize)

lmmblue <- lmer(minsize ~ nrounds*run_type + (1|site) + (1|year),
               data = df_blue)
hist(residuals(lmmblue))
qqnorm(residuals(lmmblue))
check_normality(lmmblue)

# glm with poisson due to non-normality
glmmblue <- glmer(minsize ~ nrounds*run_type + (1|site) + (1|year),
                  data = df_blue, family = poisson)
check_overdispersion(glmmblue)

summary(glmmblue)
Anova(glmmblue)
write.csv(summary(glmmblue)$coefficients,file=paste(pth,'glm_all_rounds_',crop,'.csv',sep=''))
write.csv(Anova(glmmblue),file=paste(pth,'anova_all_rounds_',crop,'.csv',sep=''))

effblue <- allEffects(glmmblue)[[1]] %>% as.data.frame() %>% mutate(crop = 'blue')
write.csv(effblue, file=paste(pth,'effects_all_rounds_blue.csv',sep=''))

# find percent increase in minimum set for null and observed
percincb <- effblue %>% group_by(crop, run_type) %>% summarize(percinc = (fit[5]-fit[1])/fit[1]*100)
incb <- effblue %>% group_by(crop, run_type) %>% summarize(inc = (fit[5]-fit[1]))

# find what percent of total increase is due to turnover effects (increase in observed minus increase in null over increase in observed)
propincb = (incb$inc[2]-incb$inc[1])/incb$inc[2]*100

## stats for minset vs number of rounds results nj watermelon
crop='njwat'
df_njwat = df_allsitesr %>% filter(crop == 'njwat')

hist(df_njwat$minsize)

lmmnjwat <- lmer(minsize ~ nrounds*run_type + (1|site) + (1|year),
                            data = df_njwat)
check_normality(lmmnjwat)

# run glm with poisson due to non-normality
glmmnjwat <- glmer(minsize ~ nrounds*run_type + (1|site) + (1|year),
                   data = df_njwat, family = poisson)
check_overdispersion(glmmnjwat)

# run glm with neg bin due to non-normality of residuals
# glmmnjwat <- glmer.nb(minsize ~ nrounds*run_type + (1|site), data = df_njwat)
summary(glmmnjwat)
Anova(glmmnjwat)
write.csv(summary(glmmnjwat)$coefficients,file=paste(pth,'glm_all_rounds_','njwat','.csv',sep=''))
write.csv(Anova(glmmnjwat),file=paste(pth,'anova_all_rounds_','njwat','.csv',sep=''))

# generate predicted values based on glm estimates
effnjwat <- allEffects(glmmnjwat)[[1]] %>% as.data.frame() %>% mutate(crop = 'njwat')
write.csv(effnjwat, file=paste(pth,'effects_all_rounds_njwat.csv',sep=''))

# find percent increase in minimum set for null and observed
percincn <- effnjwat %>% group_by(crop, run_type) %>% summarize(percinc = (fit[5]-fit[1])/fit[1]*100)

# find what percent of total increase is due to turnover effects (increase in observed minus increase in null over increase in observed)
incn <- effnjwat %>% group_by(crop, run_type) %>% summarize(inc = (fit[5]-fit[1]))
propincn <- (incn$inc[2]-incn$inc[1])/incn$inc[2]*100

## stats for minset vs number of rounds results california watermelon
crop = 'cawat'
df_cawat = df_allsitesr %>% filter(crop == 'cawat')

hist(df_cawat$minsize)

lmmcawat <- lmer(minsize ~ nrounds*run_type + (1|site) + (1|year),
                 data = df_cawat)
check_normality(lmmcawat)

# glmm with poisson due to non-normality
glmmcawat <- glmer(minsize ~ nrounds*run_type + (1+nrounds|site) + (1+nrounds|year),
                  data = df_cawat, family = poisson)
check_overdispersion(glmmcawat)

summary(glmmcawat)
Anova(glmmcawat)
write.csv(summary(glmmcawat)$coefficients,file=paste(pth,'glm_all_rounds_','cawat','.csv',sep=''))
write.csv(Anova(glmmcawat),file=paste(pth,'anova_all_rounds_','cawat','.csv',sep=''))

# generate predicted values based on glm estimates
effcawat <- allEffects(glmmcawat)[[1]] %>% as.data.frame() %>% mutate(crop='cawat')
write.csv(effcawat, file=paste(pth,'effects_all_rounds_cawat.csv',sep=''))

# find percent increase in minimum set for null and observed
percincc <- effcawat %>% group_by(crop, run_type) %>% summarize(percinc = (fit[5]-fit[1])/fit[1]*100)
percincc_short <- effcawat %>% group_by(crop, run_type) %>% summarize(percinc = (fit[2]-fit[1])/fit[1]*100)
incc <- effcawat %>% group_by(crop, run_type) %>% summarize(inc = (fit[5]-fit[1]))

# save percent increase in min set for all crops as a table
percincall <- rbind(percincb,percincn,percincc) %>%
              pivot_wider(names_from = crop, values_from = percinc)
write.csv(percincall, file=paste(pth,'percent_increase_all_rounds.csv',sep=''),row.names = F)

# save the increase in min set for null and obs as table
incall <- rbind(incb,incn,incc) %>%
          pivot_wider(names_from = crop, values_from = inc)
write.csv(incall, file=paste(pth,'obsvsnull_all_rounds.csv',sep=''),row.names = F)

## re-analyze ca wat with only first 3 rounds
dfcawat3 <- df_cawat %>% filter(nrounds <= 3)
hist(dfcawat3$minsize)

lmmcawat3 = lmer(minsize ~ nrounds*run_type + (1|site) + (1|year),
                data = dfcawat3)
check_normality(lmmcawat3)
hist(residuals(lmmcawat3))
effcawat3 <- allEffects(lmmcawat3)[[1]] %>% as.data.frame() %>% mutate(crop='cawat')
write.csv(effcawat3, file=paste(pth,'effects_all_rounds_cawat3.csv',sep=''))

percincc3 <- effcawat3 %>% group_by(crop, run_type) %>% summarize(percinc = (fit[5]-fit[1])/fit[1]*100, finc = fit[5]/fit[1])



# stats analysis for minset vs number of years ----------------------------
## run glm for blueberry, minset vs number of years
# subset blueberry results
dfblueyrs = filter(df_allsitesy, crop == 'blue')

# linear model on blue results
lmby = lmer(minsize ~ nyears*run_type +(1|site), data = dfblueyrs)
summary(lmby)
check_normality(lmby)
hist(residuals(lmby))

# glm with poisson due to non-normality
glmby = glmer(minsize ~ nyears*run_type +(1|site), data = dfblueyrs, family = poisson)
check_overdispersion(glmby)

# generate predicted values based on glm estimates
effblyrs <- allEffects(glmby)[[1]] %>% as.data.frame() %>% mutate(crop='blue')
write.csv(effblyrs, file=paste(pth,'effects_all_years_blue.csv',sep=''))

# find what percent of total increase is due to turnover effects (increase in observed minus increase in null over increase in observed)
incbyrs <- effblyrs %>% group_by(crop, run_type) %>% summarize(inc = (fit[5]-fit[1]))
propincbyrs = (incbyrs$inc[2]-incbyrs$inc[1])/incbyrs$inc[2]*100

## run lm for nj watermelon, minset vs number of years
# subset nj watermelon results
dfnjwatyrs = filter(df_allsitesy, crop == 'njwat')

# linear model on njwat results
lmny = lmer(minsize ~ nyears*run_type + (1|site), data = dfnjwatyrs)
summary(lmny)
check_normality(lmny)
hist(residuals(lmny))

# generalized linear model on njwat b/c non-normality detected
glmny = glmer(minsize ~ nyears*run_type+ (1|site), data = dfnjwatyrs, family = poisson)
summary(glmny)
Anova(glmny)
check_overdispersion(glmny)
hist(residuals(glmny))

# save glm for nj watermelon output
write.csv(summary(glmny)$coefficients,file=paste(pth,'glm_all_years_','njwat_1r','.csv',sep=''))
write.csv(Anova(glmny),file=paste(pth,'anova_all_years_','njwat_1r','.csv',sep=''))

# generate predicted values from lm estimates
effnwyrs <- allEffects(glmny)[[1]] %>% as.data.frame() %>% mutate(crop='njwat')
write.csv(effnwyrs, file=paste(pth,'effects_all_years_njwat_1r.csv',sep=''))

# find what percent of total increase is due to turnover effects (increase in observed minus increase in null over increase in observed)
incnyrs <- effnwyrs %>% group_by(crop, run_type) %>% summarize(inc = (fit[5]-fit[1]))
propincnyrs = (incnyrs$inc[2]-incnyrs$inc[1])/incnyrs$inc[2]*100

## run lm for ca watermelon minset vs number of years ##
# subset ca waterlon results
dfcawatyrs = filter(df_allsitesy, crop == 'cawat')

# linear model on cawat results
lmcy = lmer(minsize ~ nyears*run_type + (1|site), data = dfcawatyrs)
summary(lmcy)
check_normality(lmcy)
hist(residuals(lmcy))

# glm with poisson on cawat results because residuals are non-normal
glmcy = glmer(minsize ~ nyears*run_type + (1|site),  data = dfcawatyrs, family = poisson)
summary(glmcy)
check_overdispersion(glmcy)
hist(residuals(glmcy))

# save glm output for ca watermelon
write.csv(summary(glmcy)$coefficients, file = paste(pth,'glm_all_years_','cawat','.csv',sep=''))
write.csv(Anova(glmcy),file=paste(pth,'anova_all_years_','cawat','.csv',sep=''))

# generate predicted values from glm estimates
effcwyrs <- allEffects(glmcy)[[1]] %>% as.data.frame() %>% mutate(crop='cawat')
write.csv(effcwyrs, file=paste(pth,'effects_all_years_cawat.csv',sep=''))

# find what percent of total increase is due to turnover effects (increase in observed minus increase in null over increase in observed)
inccyrs <- effcwyrs %>% group_by(crop, run_type) %>% summarize(inc = (fit[5]-fit[1]))
propinccyrs = (inccyrs$inc[2]-inccyrs$inc[1])/inccyrs$inc[2]*100

## plot results ------------------------------------------------------------

# plot minset vs number of days for all site years on a single graph (one graph per crop)
for (crop in crops) {
  df_crop = df_sitemeans %>% filter(crop == crop)
  ggplot(data=df_crop, mapping=aes(x=nrounds, y=minsize, color=run_type)) +
    geom_smooth() +
    geom_jitter(width = 0.1, height = 0.1) +
    labs(x="Number of days", y="Minimum set size", color="Model type", title=paste(croplabs[crop,],' (sampled with replacement)',sep='')) +
    ylim(0.2,6.2) +
    scale_y_continuous(breaks=seq(0,7,1)) +
    scale_x_continuous(breaks=seq(0,9,1)) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = 'bottom')
  figname = paste(figpath,'minset_null_',crop,'.png',sep='')
  ggsave(filename=figname, width=6.5,height=5)
}

## make a multi-panel figure showing the observed and null model for all crops (min set vs number of days) ##
df_blue = df_allsitesr %>% filter(crop == 'blue')
crop='blue'
p1 <- 
  ggplot() +
  geom_jitter(data=filter(df_sitemeans,crop=='blue'), mapping=aes(x=nrounds, y=minsize, color=run_type, shape=run_type), 
              width = 0.1, height = 0.1, alpha=0.6, size=1) +
  geom_line(data=effblue, mapping=aes(x=nrounds, y=fit, color=run_type, lty=run_type)) +
  geom_ribbon(data=effblue, mapping=aes(x=nrounds, color=run_type, lty=run_type, ymin=lower,ymax=upper), alpha=0.1) +
  labs(x="Number of dates across the season", y="Number of species needed", color="Model type", lty="Model type", shape="Model type",
       title='Blueberry') +
  scale_y_continuous(breaks=seq(0,12,2)) +
  scale_x_continuous(breaks=seq(0,9,1)) +
  coord_cartesian( ylim=c(0.2,6), xlim=c(0.9, 3.1), clip='off') +
  annotate("text", x=0.75, y=6.8, label='c', size=6) +
  theme_light() +
  mtext('A', side = 3, adj = 0.05, line = -1.3) +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = 'none')

# df_njwat = df_sitemeans %>% filter(crop == 'njwat')
df_njwat = df_allsitesr %>% filter(crop == 'njwat')
crop='njwat'
p2 <- 
  ggplot() +
  geom_jitter(data=filter(df_sitemeans,crop=='njwat'), mapping=aes(x=nrounds, y=minsize, color=run_type, shape=run_type), 
              width = 0.1, height = 0.1, alpha=0.6, size=1) +
  geom_line(data=effnjwat, mapping=aes(x=nrounds,y=fit,color=run_type, lty=run_type)) +
  geom_ribbon(data=effnjwat, mapping=aes(x=nrounds, color=run_type, lty=run_type, ymin=lower,ymax=upper), alpha=0.1) +
  labs(x="Number of dates across the season", y="Number of species needed", color="Model type", lty="Model type", shape="Model type",
       title='Eastern watermelon') +
  scale_y_continuous(breaks=seq(0,12,2)) +
  scale_x_continuous(breaks=seq(0,9,1)) +
  # coord_cartesian(ylim=c(0.2,6)) +
  coord_cartesian( ylim=c(0,6), xlim=c(0.9,3.1), clip='off' ) +
  annotate("text", x=0.7, y=6.9, label='b', size=6) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = 'none')

# df_cawat = df_sitemeans %>% filter(crop == 'cawat')
df_cawat = df_allsitesr %>% filter(crop == 'cawat')
crop='cawat'
p3 <- 
  ggplot() +
  geom_jitter(data=filter(df_sitemeans,crop=='cawat'), mapping=aes(x=nrounds, y=minsize, color=run_type, shape=run_type), 
              width = 0.1, height = 0.1, alpha=0.6, size=1) +
  geom_line(data=effcawat, mapping=aes(x=nrounds, y=fit, color=run_type, lty=run_type)) +
  geom_ribbon(data=effcawat, mapping=aes(x=nrounds, color=run_type, lty=run_type, ymin=lower, ymax=upper), alpha=0.1) +
  labs(x="Number of dates across the season", y="Number of species needed", color="Model type", lty="Model type", shape="Model type",
       title='Western watermelon') +
  scale_y_continuous(breaks=seq(0,12,2)) +
  scale_x_continuous(breaks=seq(0,9,1)) +
  coord_cartesian( ylim=c(0.2,8), xlim=c(0.9,9.1), clip='off' ) +
  annotate("text", x=0.3, y=9, label='a', size=6) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = 'right')
# ggsave(filename=figname, width=8.5,height=6.5)

# Figure S1 NEE manuscript
plot_grid( p3, plot_grid(p2, p1, ncol=2), ncol=1)
figname = paste(figpath,'BEF over time final/','FigS1_nullmodel_rounds.pdf',sep='')
ggsave(filename=figname, width=8.5,height=6.5)


## make multipanel figure w/ all crops, showing min set vs number of years for null and observed ##
p2a <-
  ggplot() +
  geom_jitter(data=filter(df_sitemeansy,crop=='blue'), mapping=aes(x=nyears, y=minsize, color=run_type, shape=run_type), width = 0.1, height = 0.1, alpha=0.6) +  
  geom_ribbon(data=effblyrs, mapping=aes(x=nyears, color=run_type, lty=run_type, ymin=lower, ymax=upper), alpha=0.2) +
  geom_line(data=effblyrs, mapping=aes(x=nyears, y=fit, color=run_type, lty=run_type)) +
  labs(x="Number of years", y="Number of species needed", color="Model type", title=paste(croplabs['blue',],sep='')) +
  scale_y_continuous(breaks=seq(0,12,2)) +
  scale_x_continuous(breaks=seq(0,9,1)) +
  coord_cartesian(ylim=c(0.2,8), xlim=c(0.9,3.1), clip='off') +
  theme_light() +
  annotate("text", x=0.85, y=8.9, label='c', size=6) +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title = element_text(size=13),
        legend.position = 'none')

p2b <-
  ggplot() +
  geom_jitter(data=filter(df_sitemeansy,crop=='cawat'), mapping=aes(x=nyears, y=minsize, color=run_type,shape=run_type), 
              width = 0.1, height = 0.1, alpha=0.6) +  
  geom_ribbon(data=effcwyrs, mapping=aes(x=nyears, color=run_type, lty=run_type, ymin=lower, ymax=upper), alpha=0.2) +
  geom_line(data=effcwyrs, mapping=aes(x=nyears, y=fit, color=run_type, lty=run_type)) +
  labs(x="Number of years", y="Number of species needed", color="Model type", title=paste(croplabs['cawat',],sep='')) +
  scale_y_continuous(breaks=seq(0,12,2)) +
  scale_x_continuous(breaks=seq(0,9,1)) +
  coord_cartesian(ylim=c(0.2,6), xlim=c(0.9,3.1), clip='off') +
  theme_light() +
  annotate("text", x=0.85, y=6.7, label='b', size=6) +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title = element_text(size=13), 
        legend.position = 'none')

p2c <-
  ggplot() +
  geom_jitter(data=filter(df_sitemeansy,crop=='njwat'), mapping=aes(x=nyears, y=minsize, color=run_type, shape=run_type), width = 0.1, height = 0.01, alpha=0.6) +  
  geom_ribbon(data=effnwyrs, mapping=aes(x=nyears, color=run_type, lty=run_type, ymin=lower, ymax=upper), alpha=0.2) +
  geom_line(data=effnwyrs, mapping=aes(x=nyears, y=fit, color=run_type, lty=run_type)) +
  labs(x="Number of years", y="Number of species needed", color="Model type", lty="Model type", shape="Model type", title=paste(croplabs['njwat',],sep='')) +
  scale_y_continuous(breaks=seq(0,12,2)) +
  scale_x_continuous(breaks=seq(0,9,1)) +
  coord_cartesian(ylim=c(0.2,8.5), xlim=c(0.9,6.1), clip='off') +
  annotate("text", x=0.7, y=9.5, label='a', size=6) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title = element_text(size=13),
        legend.position = 'right')


# Figure S2 NEE manuscript
plot_grid(p2c, plot_grid(p2b, p2a, ncol=2), ncol=1)
figname = paste(figpath,'BEF over time final/','FigS2_nullmodel_years.pdf',sep='')
ggsave(filename=figname, width=8.5, height=6.5)
