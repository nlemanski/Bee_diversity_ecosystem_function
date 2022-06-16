### plot the accumulation curves showing how pollination function increases as species are added ##

## load packages ##
library(tidyverse)
library(hash)
library(RColorBrewer)
library(ggthemes)

setwd("C:/Documents/Bee_diversity_ecosystem_function/R code") # change to path where minfinder_function_time.R is located

## load the species required function
source(file="minfinder_function_time.R")

## file paths and file names ##
path <- "C:/Documents/Bee_diversity_ecosystem_function/SQL data/" # change to filepath where data is located
figpath <- "C:/Documents/Bee_diversity_ecosystem_function/figures/" # change to filepath where figures should be saved
rpath <- "C:/Documents/Bee_diversity_ecosystem_function/R data/" # change to filepath where data is located
outpath <- "C:/Documents/Bee_diversity_ecosystem_function/BEFresults/" # change to filepath where results should be saved

## color palette for plots ##
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# label crops
crops = c("blue","njwat","cawat")
cropnames = c("blue"="Blueberry","njwat"="Eastern watermelon", "cawat"="Western watermelon")
crop = 'njwat'

## import crop visit data as RData file ##
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

## list unique values of each factor ##
rounds <- sort(unique(df_visits$round))
years <- sort(unique(df_visits$year))
species <- sort(unique(df_visits$gen_sp))
sites <- sort(unique(df_visits$site))

# function threshold percentage (set here for sensitivity)
func_percent = 0.50

## find function threshold value (aggregate pollen grains) from percentage, using mean across all site-year-rounds
df_totalfun = df_visits %>%     # find the total function for each site-year-round
  group_by(round, year, site) %>%
  summarize(totalfun = sum(fun))

func_level = mean(df_totalfun$totalfun, na.rm=T)*func_percent

## order species of a single site in descending order of function and plot accumulation of pollination ##
for (s in 1:length(sites)) {
  df = df_visits %>% 
    filter(site==sites[s]) %>%
    arrange(year, round, -fun)
  df = df %>% group_by(year, round) %>% mutate(cumfun = cumsum(fun), rank=row_number())
  
  ggplot(data = df, aes(x=rank, y=cumfun, color=round)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = func_level, color='red', lty=3) +
    labs(x="Number of Species",
         y="Total Pollination",
         color="Date",
         title=paste('Cumulative pollination provided by all species (',cropnames[[crop]],', site ',sites[s],')',sep='')) +
    scale_color_viridis_d() +
    scale_x_continuous(breaks=seq(0,20,2)) +
    facet_wrap(~year) +
    theme_light()
  }

## choose a single site-year as illustrative example of high vs low abundance curves
# Figure 3
siteex = 'bri'
yearex = '2012'
crop = 'njwat'
datelabs = c('3'='High abundance','1'='Med. abundance','2'='Low abundance')

dfex = df_visits %>% 
  filter(site==siteex,year==yearex) %>%
  arrange(year, round, -fun) %>% 
  group_by(year, round) %>% 
  mutate(cumfun = cumsum(fun), rank=row_number())

dfex$round <- factor(dfex$round, levels=c('3','1','2'))

ihi = filter(dfex,round=='3',rank==2)$cumfun
imed = filter(dfex,round=='1',rank==3)$cumfun
ilow = filter(dfex,round=='2',rank==7)$cumfun

ggplot(data = dfex, aes(x=rank, y=cumfun, color=round)) +
  geom_point(size=1.5) +
  geom_line(size=0.8) +
  geom_hline(yintercept = func_level, color='red', lty=2, size=0.8) +
  geom_segment(aes(x = 2, y = 0, xend = 2, yend = ihi), color='purple4', lty=3, size=0.9) +
  geom_segment(aes(x = 3, y = 0, xend = 3, yend = imed), color='cyan4', lty=3, size=0.9) +
  geom_segment(aes(x = 7, y = 0, xend = 7, yend = ilow), color='gold2', lty=3, size=0.9) +
  labs(x="Number of Bee Species",
       y="Total Pollination Provided",
       color="Date" ) +
  scale_color_viridis_d(name="Date",labels = datelabs) +
  scale_x_continuous(breaks=seq(0,20,1)) +
  coord_cartesian(ylim = c(100,2350)) +
  theme_classic() +
  theme(legend.position = c(0.14, 0.87), 
        legend.text = element_text(size=10),
        legend.background = element_rect(fill="transparent"))
  
ggsave(filename=paste(figpath,'Fig3_accumulation_curves_',crop,'_example1.png',sep=''),width=6, height=3.6) # save as image
ggsave(filename=paste(figpath,'Fig3_accumulation_curves_',crop,'_example1.pdf',sep=''),width=6, height=3.6) # save as pdf




