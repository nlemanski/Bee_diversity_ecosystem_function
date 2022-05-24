### sensitivity analysis of pollination threshold ###
# edited by NJL 3.15.22

crop = "njwat"
funcpercs = seq(0.1,0.9,0.1)

## run minimum set analysis for different levels of function thresh --------
# across years
for (t in funcpercs){
  func_percent = t
  source("C:/Users/natal/OneDrive - Rutgers University/Documents/winfree lab/R code/minset_years_BEF_new.R")
}
# within years, across season
for (t in funcpercs){
  func_percent = t
  source("C:/Users/natal/OneDrive - Rutgers University/Documents/winfree lab/R code/minset_rounds_BEF.R")
}


## load packages ## -----------------------------------------------------
library(tidyverse)
library(lme4)
library(effects)
library(cowplot)

## file paths and file names ##
figpath <- "C:/Users/natal/OneDrive - Rutgers University/Documents/winfree lab/figures/"
rpath <- "C:/Users/natal/OneDrive - Rutgers University/Documents/winfree lab/R data/"
outpath <- "C:/Users/natal/OneDrive - Rutgers University/Documents/winfree lab/BEF sensitivity/"

## crop names and parameters ## ------------------------------------------
crops = c("blue","njwat","cawat")
cropnames = c("blue"="Blueberry","njwat"="Eastern Watermelon","cawat"="Western Watermelon")
subplot_letters <- data.frame(letter = c("a", "b", "c"), crop  = crops) # for labelling the panels in plots
# thresholds = c('25', '50', '75')
thresholds = as.character(funcpercs*100)
rmax = 1


# import results files for plotting ---------------------------------------

# import sensitivity analysis results for species needed vs number of rounds and put into a single dataframe
df_sens_rnds = data.frame()
for (i in thresholds){
  for (c in crops){
    fname = paste(rpath,"df_minset_rounds_",c,'_',i,'.RData',sep='')
    load(fname)
  }
  df_sens_rnds = rbind(df_sens_rnds, df_rounds_blue, df_rounds_njwat, df_rounds_cawat)
}
df_sens_rnds = df_sens_rnds %>% select(-minset, -presets) # drop the list of species in the minimum set
df_sens_rnds = df_sens_rnds %>% mutate(func_percentf = factor(func_percent*100)) # convert threshold to string

# import sensitivity analysis results for species needed vs number of years and put into a single dataframe
df_sens_yrs = data.frame()
for (i in thresholds) {
  for (c in c("blue","cawat")){
    fname2 = paste(rpath,"df_minset_years_",c,'_',i,".RData",sep='')
    load(fname2)
  }
  fname3 = paste(rpath,"df_minset_years_njwat_",rmax,"round_",i,".RData",sep='')
  load(fname3)
  df_sens_yrs = rbind(df_sens_yrs, df_years_blue, df_years_njwat, df_years_cawat)
}
df_sens_yrs = df_sens_yrs %>% select(-minset, -presets) # drop the list of species in the minimum set
df_sens_yrs = df_sens_yrs %>% mutate(func_percentf = factor(func_percent*100)) # convert threshold to string


# plot number of species needed vs number of rounds w/ each thresh --------
ggplot() +
  geom_smooth(data = filter(df_sens_rnds), aes(x = nrounds, y = minsize, color = func_percentf), method = glm) +
  facet_wrap(~crop, scales="free", labeller = as_labeller(cropnames)) +
  # geom_line(data=effbr1, aes(x=nrounds, y=fit), color="yellow", linetype="solid") +
  # geom_ribbon(data=effbr1, aes(x=nrounds, ymin=lower, ymax=upper), alpha= 0.3, fill="yellow") +
  scale_color_viridis_d(direction=-1) +
  labs(#title = "Sensitivity Analysis of Function Threshold",
       x= "Number of dates across the season",
       y = "Number of bee species needed",
       color = "Function Threshold (%)") +
  scale_x_continuous(breaks=seq(0,9,1)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=10), legend.position = "bottom") +
  guides(color=guide_legend(nrow=1))
ggsave(paste(figpath,"sensitivity_analysis_rounds2.png",sep=''), width=7, height=3, units="in")

# find increase in species needed w/ number of rounds at each threshold --------
slopes = data.frame()
for (c in crops) {
  for (t in thresholds) {
    df = filter(df_sens_rnds,crop==c,func_percentf==t) # %>% mutate(site=factor(site)) 
    glmr = glmer(minsize ~ nrounds + (1|site),
               data=df, family='poisson')
    slopes = rbind(slopes, data.frame(Crop = c, 
                                      Threshold = t, 
                                      Slope = summary(glmr)$coefficients[2], 
                                      Slope_SE = summary(glmr)$coefficients[4],
                                      Intercept = summary(glmr)$coefficients[1]))
  }
}

# plot the effect size (increase in species needed w/ number of rounds) for each threshold
ggplot() +
  geom_point(data=slopes, aes(x=Threshold, y=Slope, color=Crop, group=Crop)) +
  geom_line(data=slopes, aes(x=Threshold, y=Slope, color=Crop, group=Crop)) +
  geom_ribbon(data=slopes, aes(x=Threshold, ymin=(Slope-Slope_SE), ymax=(Slope+Slope_SE), color=Crop, group=Crop), alpha=0.1) +
    scale_color_manual(labels=cropnames, values=c("blue","red","green")) +
  theme_classic() +
  labs(x = "Function Threshold (%)",
       y = "Slope (Species Needed Vs. Number of Rounds)") +
  coord_cartesian(ylim = c(0,.5)) +
  theme(legend.position = "bottom")
ggsave(paste(figpath,"slopevsthreshold_rounds3.png",sep=''))

# plot number of species needed vs number of years w/ each thresh --------
ggplot() +
  geom_smooth(data = filter(df_sens_yrs), aes(x = nyears, y = minsize, color = func_percentf), method = glm) +
  facet_wrap(~crop, scales="free", labeller = as_labeller(cropnames)) +
  scale_color_viridis_d(direction=-1) +
  labs(#title = "Sensitivity Analysis of Function Threshold",
       x= "Number of years",
       y = "Number of bee species needed",
       color = "Function Threshold (%)") +
  scale_x_continuous(breaks=seq(0,6,1)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=10), legend.position = "bottom") +
  guides(color = guide_legend(nrow=1))
# tag_facet_outside(p, open = '', close = '')
ggsave(paste(figpath,"sensitivity_analysis_years3.png",sep=''), width=7, height=3, units="in")

# find increase in species needed w/ number of years at each threshold --------
slopesyr = data.frame()
for (c in crops) {
  for (t in thresholds) {
    df = filter(df_sens_yrs,crop==c,func_percentf==t) # %>% mutate(site=factor(site)) 
    glmy = glmer(minsize ~ nyears + (1|site),
                 data=df, family='poisson')
    slopesyr = rbind(slopesyr, data.frame(Crop = c, 
                                      Threshold = t, 
                                      Slope = summary(glmy)$coefficients[2],
                                      Slope_SE = summary(glmy)$coefficients[4],
                                      Intercept = summary(glmy)$coefficients[1]))
  }
}

# plot the effect size (increase in species needed w/ number of years) for each threshold
ggplot() +
  geom_point(data=slopesyr, aes(x=Threshold, y=Slope, color=Crop, group=Crop)) +
  geom_line(data=slopesyr, aes(x=Threshold, y=Slope, color=Crop, group=Crop)) +
  geom_ribbon(data=slopesyr, aes(x=Threshold, ymin=(Slope-Slope_SE), ymax=(Slope+Slope_SE), color=Crop, group=Crop), alpha=0.1) +
  scale_color_manual(labels=cropnames, values=c("blue","red","green")) +
  theme_classic() +
  labs(x = "Function Threshold (%)",
       y = "Slope (Species Needed Vs. Number of Years)") +
  coord_cartesian(ylim = c(0,0.5)) +
  theme(legend.position = "bottom")
ggsave(paste(figpath,"slopevsthreshold_years3.png",sep=''))


# effect size, number of species vs rounds --------------------------------
# glm to see if function threshold significantly affects the increase in species needed w/ number of rounds, blueberry
glmrb = glmer(minsize ~ nrounds*func_percent + (1|site), data = filter(df_sens_rnds,crop=="blue"), family = poisson)
summary(glmrb)  # there is significant interaction between function threshold and number of rounds

# find the magnitude of the increase in species needed for different thresholds, blueberry
effrb = as.data.frame(effect(term = "nrounds*func_percent", mod = glmrb, xlevels=9))
pincrb = effrb %>% group_by(func_percent) %>% 
  summarize(#finc = (last(fit)/first(fit)), 
            pinc = 100*(last(fit)-first(fit))/first(fit)) %>%
  mutate(crop='blue') %>% 
  pivot_wider(names_from = crop, values_from = pinc)

# glm to see if function threshold significantly affects the increase in species needed w/ number of rounds, nj watermelon
glmrnj = glmer(minsize ~ nrounds*func_percent + (1|site), data = filter(df_sens_rnds,crop=="njwat"), family = poisson)
summary(glmrnj)  # there is significant interaction between function threshold and number of rounds

# find the magnitude of the increase in species needed for different thresholds, nj watermelon
effrnj = as.data.frame(effect(term = "nrounds*func_percent", mod = glmrnj, xlevels=9))
pincrnj = effrnj %>% group_by(func_percent) %>% 
  summarize(#finc = (last(fit)/first(fit)), 
            pinc = 100*(last(fit)-first(fit))/first(fit)) %>%
  mutate(crop='njwat')%>% 
  pivot_wider(names_from = crop, values_from = pinc)

# glm to see if function threshold significantly affects the increase in species needed w/ number of rounds, ca watermelon
glmrca = glmer(minsize ~ nrounds*func_percent + (1|site), data = filter(df_sens_rnds,crop=="cawat"), family = poisson)
summary(glmrca)  # there is significant interaction between function threshold and number of rounds

# find the magnitude of the increase in species needed for different thresholds, nj watermelon
effrca = as.data.frame(effect(term = "nrounds*func_percent", mod = glmrca, xlevels=9))
pincrca = effrca %>% group_by(func_percent) %>% 
  summarize(#finc = (last(fit)/first(fit)), 
            pinc = 100*(last(fit)-first(fit))/first(fit)) %>%
  mutate(crop='cawat')%>% 
  pivot_wider(names_from = crop, values_from = pinc)

# make a table showing the percent increase in species needed for all rounds vs function threshold
pincr = join(pincrb, pincrnj) %>% join(pincrca) %>% mutate(func_percent = func_percent*100) %>% round(0)
write.csv(pincr, paste(outpath,"percentincrease_rounds.csv",sep=''), row.names = F)

# effect size, number of species vs years ---------------------------------
# glm to see effect of function threshold on species needed vs number of years
glmyb = glmer(minsize ~ nyears*func_percent + (1|site), data = filter(df_sens_yrs,crop=="blue"), family = poisson)
summary(glmyb)

glmynj = glmer(minsize ~ nyears*func_percent + (1|site), data = filter(df_sens_yrs,crop=="njwat"), family = poisson)
summary(glmynj)

glmyca = glmer(minsize ~ nyears*func_percent + (1|site), data = filter(df_sens_yrs,crop=="cawat"), family = poisson)
summary(glmyca)

# find the magnitude of the increase in species needed across years for different thresholds, blueberry
effyb = as.data.frame(effect(term = "nyears*func_percent", mod = glmyb, xlevels=9))
pincyb = effyb %>% group_by(func_percent) %>% 
  summarize(#finc = (last(fit)/first(fit)), 
            pinc = 100*(last(fit)-first(fit))/first(fit)) %>%
  mutate(crop='blue', func_percent = func_percent*100) %>% 
  pivot_wider(names_from = crop, values_from = pinc) %>%  round(0)

# find the magnitude of the increase in species needed across years for different thresholds, nj watermelon
effynj = as.data.frame(effect(term = "nyears*func_percent", mod = glmynj, xlevels=9))
pincynj = effynj %>% group_by(func_percent) %>% 
  summarize(#finc = (last(fit)/first(fit)), 
            pinc = 100*(last(fit)-first(fit))/first(fit)) %>%
  mutate(crop='njwat', func_percent = func_percent*100) %>% 
  pivot_wider(names_from = crop, values_from = pinc) %>%  round(0)

# find the magnitude of the increase in species needed across years for different thresholds, ca watermelon
effyca = as.data.frame(effect(term = "nyears*func_percent", mod = glmyca, xlevels=9))
pincyca = effyca %>% group_by(func_percent) %>% 
  summarize(#finc = (last(fit)/first(fit)), 
            pinc = 100*(last(fit)-first(fit))/first(fit)) %>%
  mutate(crop='cawat', func_percent = func_percent*100) %>% 
  pivot_wider(names_from = crop, values_from = pinc) %>%  round(0)

# make a table showing the percent increase in species needed for all years vs function threshold
pincy = join(pincyb, pincynj) %>% join(pincyca)
write.csv(pincy, paste(outpath,"percentincrease_years.csv",sep=''), row.names = F)
