## plot results of minimum set analysis for BEF over time ##

## load packages ##
library(tidyverse)
library(lme4)
library(effects)
library(cowplot)

## file paths and file names ##
figpath <- "C:/Documents/Bee_diversity_ecosystem_function/figures/"
rpath <- "C:/Documents/Bee_diversity_ecosystem_function/R data/"
outpath <- "C:/Documents/Bee_diversity_ecosystem_function/BEFresults/"

## functions ##
# a function to combine plots into a multi-panel plot
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## load results files ##
crops = c("blue","njwat","cawat")
cropnames = c("blue"="Blueberry","njwat"="Eastern Watermelon","cawat"="Western Watermelon")
threshold = 0.5
for (c in crops){
  fname = paste(rpath,"df_minset_years_",c,'.RData',sep='')
  fname2 = paste(rpath,"df_minset_rounds_",c,'.RData',sep='')
  load(fname)
  load(fname2)
}
rmax = 1
fname3 = paste(rpath,"df_minset_years_njwat_",rmax,"round.RData",sep='')
load(fname3)

# ## view one line for each site: minimum set vs number of rounds--------
# blueberry 
ggplot(data = df_rounds_blue, mapping = aes(x = nrounds, y = minsize, color = site)) +
  geom_jitter(width = 0.1, height = 0.1) +
  geom_smooth(mapping = aes(x = nrounds, y = minsize, color = site), method = "lm", alpha = .08) +
  labs(title = paste("Minimum Species Needed for ",threshold*100,"% of Pollination (",cropnames[1],")",sep=''),
       x= "Number of Rounds",
       y = "Number of Species",
       color= "Site") +
  coord_cartesian( ylim=c(0,10)) +
  scale_x_continuous(breaks=seq(0,3,1)) +
  scale_y_continuous(breaks=seq(0,10,2)) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = "bottom") +
  guides(col = guide_legend(nrow = 2))

figname = paste(figpath,'/minset_rounds_','blue','_',threshold*100,'.png',sep='')
ggsave(filename=figname, device="png", width=6.5,height=5)

# nj watermelon
ggplot(data = df_rounds_njwat) +
  geom_jitter(mapping = aes(x = nrounds, y = minsize, color = site), width = 0.1, height = 0.1) +
  geom_smooth(mapping = aes(x = nrounds, y = minsize, color = site), method = "lm", alpha = .08) +
  labs(title = paste("Minimum Species Needed for ",threshold*100,"% of Pollination (",cropnames[3],")",sep=''),
       x= "Number of Rounds",
       y = "Number of Species",
       color= "Site") +
  coord_cartesian(ylim=c(0,16)) +
  scale_x_continuous(breaks=seq(0,3,1)) +
  scale_y_continuous(breaks=seq(0,16,4)) +
  theme_light() +     
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = "bottom") +
  guides(col = guide_legend(nrow = 3))
  
figname = paste('minset_rounds_','njwat','_',threshold*100,'.png',sep='')
#figname = paste('minset_rounds_','njwat','_',threshold*100,'_nolegend.png',sep='')
ggsave(filename=figname, device="png", path=figpath, width=6.5,height=5)

# CA watermelon
ggplot(data = df_rounds_cawat, mapping = aes(x = nrounds, y = minsize, color = site)) +
  geom_jitter(width = 0.1, height = 0.1) +
  geom_smooth(mapping = aes(x = nrounds, y = minsize, color = site), method = "lm", alpha = .08) +
  labs(title = paste("Minimum Species Needed for ",threshold*100,"% of Pollination (",cropnames[4],")",sep=''),
       x= "Number of Rounds",
       y = "Number of Species",
       color= "Site") +
  coord_cartesian(ylim=c(0,10)) +
  scale_x_continuous(breaks=seq(0,9,1)) +
  scale_y_continuous(breaks=seq(0,10,2)) +
  theme_light() +     
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.text=element_text(size=6), 
        legend.position = "bottom") +
  guides(col = guide_legend(nrow = 4))

figname = paste('minset_rounds_','cawat','_',threshold*100,'.png',sep='')
ggsave(filename=figname, device="png", path=figpath, width=7.5,height=5)


## view as scatterplot: size of minimum set vs number of years ##
# blueberry 
ggplot(data = df_years_blue, mapping = aes(x = nyears, y = minsize, color = site)) +
  geom_jitter(width = 0.1, height = 0.1) +
  geom_smooth(mapping = aes(x = nyears, y = minsize, color = site), method = "lm", alpha = .08) +
  labs(title = paste("Minimum Species Needed for ",threshold*100,"% of Pollination (",cropnames[1],")",sep=''),
       x= "Number of Years",
       y = "Number of Species",
       color= "Site") +
  coord_cartesian( ylim=c(0,13)) +
  scale_x_continuous(breaks=seq(0,3,1)) +
  scale_y_continuous(breaks=seq(0,14,2)) +
  theme_light() +     
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = "bottom") +
  guides(col = guide_legend(nrow = 2))

figname = paste(figpath,'/minset_years_','blue','_',threshold*100,'.png',sep='')
ggsave(filename=figname, device="png", width=6.5,height=5)

# cranberry 
ggplot(data = df_years_cran, mapping = aes(x = nyears, y = minsize, color = site)) +
  geom_jitter(width = 0.1, height = 0.1) +
  geom_smooth(mapping = aes(x = nyears, y = minsize, color = site), method = "lm", alpha = .08) +
  labs(title = paste("Minimum Species Needed for ",threshold*100,"% of Pollination (",cropnames[2],")",sep=''),
       x= "Number of Years",
       y = "Number of Species",
       color= "Site") +
  coord_cartesian(xlim=c(0.7,2.2), ylim=c(0,15)) +
  scale_x_continuous(breaks=seq(0,3,1)) +
  scale_y_continuous(breaks=seq(0,16,2)) +
  theme_light() +     
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = "bottom") +
  guides(col = guide_legend(nrow = 2))

figname = paste('minset_years_','cran','_',threshold*100,'.png',sep='')
ggsave(filename=figname, device="png", path=figpath, width=8,height=5)

# nj watermelon
ggplot(data = df_years_njwat) +
  geom_jitter(mapping = aes(x = nyears, y = minsize, color = site), width = 0.1, height = 0.1) +
  geom_smooth(mapping = aes(x = nyears, y = minsize, color = site), method = "lm", alpha = .08) +
  labs(title = paste("Minimum Species Needed for ",threshold*100,"% of Pollination (",cropnames[3],")",sep=''),
       x= "Number of Years",
       y = "Number of Species",
       color= "Site") +
  coord_cartesian(ylim=c(0,20)) +
  scale_x_continuous(breaks=seq(0,6,1)) +
  scale_y_continuous(breaks=seq(0,20,4)) +
  theme_light() +     
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = "bottom") +
  guides(col = guide_legend(nrow = 3))

figname = paste('minset_years_','njwat','_',threshold*100,'_',rmax,'round.png',sep='')
#figname = paste('minset_years_','njwat','_',threshold*100,'_nolegend.png',sep='')
ggsave(filename=figname, device="png", path=figpath, width=6.5,height=5)

# CA watermelon
ggplot(data = df_years_cawat, mapping = aes(x = nyears, y = minsize, color = site)) +
  geom_jitter(width = 0.1, height = 0.1) +
  geom_smooth(mapping = aes(x = nyears, y = minsize, color = site), method = "lm", alpha = .08) +
  labs(title = paste("Minimum Species Needed for ",threshold*100,"% of Pollination (",cropnames[4],")",sep=''),
       x= "Number of Years",
       y = "Number of Species",
       color= "Site") +
  coord_cartesian(ylim=c(0,10)) +
  scale_x_continuous(breaks=seq(0,9,1)) +
  scale_y_continuous(breaks=seq(0,10,2)) +
  theme_light() +     
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.text=element_text(size=6), 
        legend.position = "bottom") +
  guides(col = guide_legend(nrow = 4))

figname = paste('minset_years_','cawat','_',threshold*100,'.png',sep='')
ggsave(filename=figname, device="png", path=figpath, width=6.5,height=5)



# manuscript figures 1 and 2 ----------------------------------------------

### view scatterplot with a single trend line based on glmm estimates ###
## min set vs number of rounds (within-years) with site-years as reps


glmmblueround <- glmer(minsize ~ nrounds + (1|site),
                       data=df_rounds_blue, family='poisson')
effblue = as.data.frame(effect(term="nrounds", mod=glmmblueround))

p1a <-
  ggplot() +
  geom_jitter(data = df_rounds_blue, mapping = aes(x = nrounds, y = minsize), width = 0.1, height = 0.1, size=1.3) +
  geom_line(data=effblue, aes(x=nrounds, y=fit), color="black", linetype="solid") +
  geom_ribbon(data=effblue, aes(x=nrounds, ymin=lower, ymax=upper), alpha= 0.3, fill="black") +
  # geom_abline(slope=slb, intercept=intb) +
  #labs(title = paste("Minimum Species Needed for ",threshold*100,"% of Pollination (",cropnames[1],")",sep=''),
  labs(title = "Blueberry",
       x= "Number of Dates Across the Season",
       y = "Number of Species") +
       #color= "Site") +
  coord_cartesian( ylim=c(0,8), xlim=c(0.9, 3.1), clip='off') +
  annotate("text", x=0.7, y=9.2, label='c', size=6) +
  scale_x_continuous(breaks=seq(0,3,1)) +
  scale_y_continuous(breaks=seq(0,10,2)) +
  theme_classic() +     
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = "bottom") +
  guides(col = guide_legend(nrow = 2))

figname = paste(figpath,'/minset_rounds_','blue','_',threshold*100,'_glmmline.png',sep='')
ggsave(filename=figname, device="png", width=6.5,height=5)

# nj watermelon glmm
glmmnjwatround <- glmer.nb(minsize ~ nrounds + (1|site),  # neg bin was used because of overdispersion
                           data = df_rounds_njwat)
effnjwat = as.data.frame(effect(term="nrounds", mod=glmmnjwatround))

p1b <- 
  ggplot() +
  geom_jitter(data = df_rounds_njwat, mapping = aes(x = nrounds, y = (minsize)), width = 0.1, height = 0.1, size=1.3) +
  geom_line(data=effnjwat, aes(x=nrounds, y=fit), color="black", linetype="solid") +
  geom_ribbon(data=effnjwat, aes(x=nrounds, ymin=lower, ymax=upper), alpha= 0.3, fill="black") +
  #labs(title = paste("Minimum Species Needed for ",threshold*100,"% of Pollination (",cropnames[3],")",sep=''),
  labs(title = "Eastern watermelon",
       x= "Number of Dates Across the Season",
       y = "Number of Species") +
       #color= "Site") +
  coord_cartesian( ylim=c(0,16), xlim=c(0.9,3.1), clip='off' ) +
  annotate("text", x=0.68, y=18.5, label='b', size=6) +
  scale_x_continuous(breaks=seq(0,3,1)) +
  scale_y_continuous(breaks=seq(0,16,4)) +
  theme_classic() +     
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = "bottom") +
  guides(col = guide_legend(nrow = 3))

figname = paste(figpath,'/minset_rounds_','njwat','_',threshold*100,'_glmmline.png',sep='')
ggsave(filename=figname, device="png", width=6.5,height=5)

# ca watermelon glmm
glmmcawatround <-  glmer(minsize ~ nrounds + (1|site),
                         data=df_rounds_cawat, family='poisson')

glmcawatround <-  glm(minsize ~ nrounds,
                         data=df_rounds_cawat, family='poisson')

effcawat = as.data.frame(effect(term="nrounds", mod=glmcawatround))

p1c <- 
  ggplot() +
  geom_jitter(data = df_rounds_cawat, mapping = aes(x = nrounds, y = minsize), width = 0.1, height = 0.1, size=1.3) +
  geom_line(data=effcawat, aes(x=nrounds, y=fit), color="black", linetype="solid") +
  geom_ribbon(data=effcawat, aes(x=nrounds, ymin=lower, ymax=upper), alpha= 0.3, fill="black") +
 # labs(title = paste("Minimum Species Needed for ",threshold*100,"% of Pollination (",cropnames[4],")",sep=''),
  labs(title = "Western watermelon",
       x= "Number of Dates Across the Season",
       y = "Number of Species") +
       #color= "Site") +
  coord_cartesian( ylim=c(0,10), xlim=c(0.9,9.1), clip='off' ) +
  annotate("text", x=0.3, y=11.3, label='a', size=6) +
  scale_x_continuous(breaks=seq(0,9,1)) +
  scale_y_continuous(breaks=seq(0,10,2)) +
  theme_classic() +     
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = "bottom", legend.text=element_text(size=6), ) +
  guides(col = guide_legend(nrow = 4))

figname = paste(figpath,'/minset_rounds_','cawat','_',threshold*100,'_',rmax,'round_glmline.png',sep='')
ggsave(filename=figname, device="png", width=7.5,height=5)


### make multi-panel figure for min species across days ##
# plot min species vs Number of Dates Across the Season
plot_grid(p1c, plot_grid(p1b, p1a, ncol=2), ncol=1)
figname = paste(figpath,'/','Fig1_minset_rounds_',threshold*100,'new.png',sep='')
ggsave(filename=figname, device='png', width=8.5, height=6.5)

### scatterplot with a single trend line based on glmm estimates ###
## min set vs number of years (across-years) with all rounds included per year

# blueberry glmm
glmmblueyear <- glmer(minsize ~ nyears + (1|site),
                       data=df_years_blue, family='poisson')

effblue = as.data.frame(effect(term="nyears", mod=glmmblueyear))

p2c <-
  ggplot() +
  geom_jitter(data = df_years_blue, mapping = aes(x = nyears, y = minsize), width = 0.1, height = 0.1, size=1.3) +
  geom_line(data=effblue, aes(x=nyears, y=fit), color="black", linetype="solid") +
  geom_ribbon(data=effblue, aes(x=nyears, ymin=lower, ymax=upper), alpha= 0.3, fill="black") +
  labs(title = "Blueberry",
       x= "Number of Years",
       y = "Number of Species") +
       #color= "Site") +
  coord_cartesian( ylim=c(0,11), xlim=c(0.9,3.1), clip='off') +
  annotate("text", x=0.7, y=12.5, label='c', size=6) +
  scale_x_continuous(breaks=seq(0,3,1)) +
  scale_y_continuous(breaks=seq(0,12,2)) +
  theme_classic() +     
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = "bottom") +
  guides(col = guide_legend(nrow = 2))

figname = paste(figpath,'/minset_years_','blue','_',threshold*100,'_glmmline.png',sep='')
ggsave(filename=figname, device="png", width=6.5,height=5)

# nj watermelon glmm
glmmfullnjwatyear <-  glmer(minsize ~ nyears + (1 + nyears | site),
                        data=df_years_njwat, family = "poisson")

effnjwat = as.data.frame(effect(term="nyears", mod=glmmfullnjwatyear))

p2a <-
  ggplot() +
  geom_jitter(data = df_years_njwat, mapping = aes(x = nyears, y = (minsize)), width = 0.1, height = 0.1, size=1.3) +
  geom_line(data=effnjwat, aes(x=nyears, y=fit), color="black", linetype="solid") +
  geom_ribbon(data=effnjwat, aes(x=nyears, ymin=lower, ymax=upper), alpha= 0.3, fill="black") +
  labs(title = "Eastern watermelon",
       x= "Number of Years",
       y = "Number of Species") +
       #color= "Site") +
  coord_cartesian( ylim=c(0,16), xlim=c(0.8,6.2), clip='off' ) + 
  annotate("text", x=0.4, y=18, label='a', size=6) +
  scale_x_continuous(breaks=seq(0,6,1)) +   
  scale_y_continuous(breaks=seq(0,20,4)) +
  theme_classic() +     
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = "bottom") +
  guides(col = guide_legend(nrow = 3))

figname = paste(figpath,'/minset_years_','njwat','_',threshold*100,'_',rmax,'round_glmmlinefull.png',sep='')
ggsave(filename=figname, device="png", width=6.5,height=5)

# ca watermelon glmm
glmmcawatyear <-  glmmTMB(minsize ~ nyears + (1|site),
                        data=df_years_cawat, family = "poisson")
glmcawatyear <- glm(minsize ~ nyears,
                data=df_years_cawat, family = "poisson")

effcawat = as.data.frame(effect(term="nyears", mod=glmmcawatyear)) # changed to mixed model w/ random intercept

p2b <-
  ggplot() +
  geom_jitter(data = df_years_cawat, mapping = aes(x = nyears, y = minsize), width = 0.1, height = 0.1, size=1.3) +
  geom_line(data=effcawat, aes(x=nyears, y=fit), color="black", linetype="solid") +
  geom_ribbon(data=effcawat, aes(x=nyears, ymin=lower, ymax=upper), alpha= 0.3, fill="black") +
  labs(title = "Western watermelon",
       x= "Number of Years",
       y = "Number of Species") +
       #color= "Site") +
  coord_cartesian( ylim=c(0,8), xlim=c(0.9,3.1), clip='off' ) +
  annotate("text", x=0.7, y=9.3, label='b', size=6) +
  scale_x_continuous(breaks=seq(0,9,1)) +
  scale_y_continuous(breaks=seq(0,8,2)) +
  theme_classic() +     
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = "bottom", legend.text=element_text(size=6), ) +
  guides(col = guide_legend(nrow = 4))

figname = paste(figpath,'/minset_years_','cawat','_',threshold*100,'_glmline.png',sep='')
ggsave(filename=figname, device="png", width=6.5,height=5)

### make multi-panel figure for min species across years ##
# plot min species vs number of years
plot_grid(p2a, plot_grid(p2b, p2c, ncol=2), ncol=1)
figname = paste(figpath,'/','Fig2_minset_years_',threshold*100,'new.png',sep='')
ggsave(filename=figname, device='png', width=8.5, height=6.5)


### make plot of pre-required species ---------------------------------------
# add column for number of pre-required species to within-year results
df_rounds_blue <- df_rounds_blue %>% mutate(prereq = 0)
for (i in 1:nrow(df_rounds_blue)) {
  prereqn = length(unlist(strsplit(df_rounds_blue$presets[i],',')))
  df_rounds_blue$prereq[i] = prereqn
}
df_rounds_njwat <- df_rounds_njwat %>% mutate(prereq = 0)
for (i in 1:nrow(df_rounds_njwat)) {
  prereqn = length(unlist(strsplit(df_rounds_njwat$presets[i],',')))
  df_rounds_njwat$prereq[i] = prereqn
}
df_rounds_cawat <- df_rounds_cawat %>% mutate(prereq = 0)
for (i in 1:nrow(df_rounds_cawat)) {
  prereqn = length(unlist(strsplit(df_rounds_cawat$presets[i],',')))
  df_rounds_cawat$prereq[i] = prereqn
}

# add column for number of pre-required species to across-year results
df_years_blue <- df_years_blue %>% mutate(prereq = 0)
for (i in 1:nrow(df_years_blue)) {
  prereqn = length(unlist(strsplit(df_years_blue$presets[i],',')))
  df_years_blue$prereq[i] = prereqn
}
df_years_njwat <- df_years_njwat %>% mutate(prereq = 0)
for (i in 1:nrow(df_years_njwat)) {
  prereqn = length(unlist(strsplit(df_years_njwat$presets[i],',')))
  df_years_njwat$prereq[i] = prereqn
}
df_years_cawat <- df_years_cawat %>% mutate(prereq = 0)
for (i in 1:nrow(df_years_cawat)) {
  prereqn = length(unlist(strsplit(df_years_cawat$presets[i],',')))
  df_years_cawat$prereq[i] = prereqn
}

# plot minimum set vs rounds along with prereq species vs rounds
A <-
  ggplot() +
  geom_smooth(data=df_rounds_cawat, aes(nrounds,minsize,color='blue'), method='loess') +
  geom_smooth(data=df_rounds_cawat, aes(nrounds,prereq,color='red'), method='loess') +
  labs(title = 'Western watermelon',
       x= "Number of Sampling Days",
       y = "Number of Species") +
  scale_color_manual(name='',values=c('blue'='blue', 'red'='red'), labels=c('Minimum set','Prerequired species')) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = "bottom", legend.text=element_text(size=11))

B <-
  ggplot() +
  geom_smooth(data=df_rounds_njwat, aes(nrounds,minsize,color='blue'), method='loess') +
  geom_smooth(data=df_rounds_njwat, aes(nrounds,prereq,color='red'), method='loess') +
  labs(title = 'Eastern watermelon',
       x= "Number of Sampling Days",
       y = "Number of Species") +
  scale_color_manual(name='',values=c('blue'='blue', 'red'='red'), labels=c('Minimum set','Prerequired species')) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = "bottom", legend.text=element_text(size=11))

C <- 
  ggplot() +
  geom_smooth(data=df_rounds_blue, aes(nrounds,minsize,color='blue'), method='loess') +
  geom_smooth(data=df_rounds_blue, aes(nrounds,prereq,color='red'), method='loess') +
  labs(title = 'Blueberry',
       x= "Number of Sampling Days",
       y = "Number of Species") +
  scale_color_manual(name='',values=c('blue'='blue', 'red'='red'), labels=c('Minimum set','Prerequired species')) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = "bottom", legend.text=element_text(size=11))

plot_grid( A, plot_grid(B, C, ncol=2), ncol=1)
figname = paste(figpath,'/','Prereqspec_rounds','.png',sep='')
ggsave(filename=figname, device='png', width=8.5, height=6.5)

# plot minimum set vs years along with prereq species vs years
B <-
  ggplot() +
  geom_smooth(data=df_years_cawat, aes(nyears,minsize,color='blue'), method='loess') +
  geom_smooth(data=df_years_cawat, aes(nyears,prereq,color='red'), method='loess') +
  labs(title = 'Western watermelon',
       x= "Number of Years",
       y = "Number of Species") +
  scale_color_manual(name='',values=c('blue'='blue', 'red'='red'), labels=c('Minimum set','Prerequired species')) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = "bottom", legend.text=element_text(size=11))

A <-
  ggplot() +
  geom_smooth(data=df_years_njwat, aes(nyears,minsize,color='blue'), method='loess') +
  geom_smooth(data=df_years_njwat, aes(nyears,prereq,color='red'), method='loess') +
  labs(title = 'Eastern watermelon',
       x= "Number of Years",
       y = "Number of Species") +
  scale_color_manual(name='',values=c('blue'='blue', 'red'='red'), labels=c('Minimum set','Prerequired species')) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = "bottom", legend.text=element_text(size=11))

C <- 
  ggplot() +
  geom_smooth(data=df_years_blue, aes(nyears,minsize,color='blue'), method='loess') +
  geom_smooth(data=df_years_blue, aes(nyears,prereq,color='red'), method='loess') +
  # geom_jitter(data=df_years_blue, aes(nyears,minsize,color='blue'), width=0.1, height=0.1) +
  # geom_jitter(data=df_years_blue, aes(nyears,prereq,color='red'), width=0.1, height=0.1) +
  labs(title = 'Blueberry',
       x= "Number of Years",
       y = "Number of Species") +
  scale_color_manual(name='',values=c('blue'='blue', 'red'='red'), labels=c('Minimum set','Prerequired species')) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = "bottom", legend.text=element_text(size=11))

plot_grid( A, plot_grid(B, C, ncol=2), ncol=1)
figname = paste(figpath,'/','Prereqspec_years','.png',sep='')
ggsave(filename=figname, device='png', width=8.5, height=6.5)


### make plots of species needed vs total abundance --------------------
# make df for all 3 crops
df_funall = data.frame() # total function by site-date for just 1st rounds
df_funallrnds = data.frame() # total function by site-date for all rounds

for (crop in crops) {
  # import visit data for crop
  fname = paste(rpath,'df_visits_',crop,'.RData',sep='')
  load(fname)
  if (crop == "cawat") {
    df_visits <- df_visits_cawat  # assign crop data to df_visits
    df_rounds <- df_rounds_cawat
  } else if (crop == "blue") {
    df_visits <- df_visits_blue  
    df_rounds <- df_rounds_blue
  } else if (crop == "njwat") {
    df_visits <- df_visits_njwat  
    df_rounds <- df_rounds_njwat
  }
  
  # find total function and richness at each site date
  df_fun = df_visits %>%
    group_by(site,year,round) %>%
    summarize(visits = sum(visits), fun=sum(fun), rich=length(unique(gen_sp)))
  
  # add df with min set to df with richness and total function for all rounds
  df_funrnd = left_join(df_fun,df_rounds)
  
  # add crop result to df with all other crops
  df_funallrnds <-rbind(df_funallrnds,df_funrnd)
  
  # filter out just the first round for each site-year
  df_fun1 = df_fun %>%
    group_by(site,year) %>%
    summarize(round = first(round), visits=first(visits), fun=first(fun), rich=first(rich))
  
  # for minimum set output, filter out just values for one round
  df_minset1 = df_rounds %>% filter(nrounds==1)
  
  # add min set to df with richness and total function for 1 round
  df_fun1 = left_join(df_fun1, df_minset1, by = c("site","year"))
  
  # add single round crop results to df w/ other crops
  df_funall <- rbind(df_funall,df_fun1)
}

## plots for single crops ##
crop = 'cawat'
# plot site-date total function vs number of species needed to meet threshold
ggplot(data=df_fun1, aes(x=visits, y=minsize)) +
  geom_point() +
  geom_smooth(method='loess') +
  labs(x="Total pollination provided by all species",
       y="Number of species needed to meet threshold") +
  theme_light()

ggplot(data=df_fun1, aes(x=fun, y=(minsize/rich*100))) +
  geom_point() +
  geom_smooth(method='loess') +
  labs(x="Total pollination provided by all species",
       y="Percent of species needed to meet threshold",
       title=paste("Effect of total function on species needed (",cropnames[[crop]],")",sep='')) +
  coord_cartesian(ylim=c(0,100)) +
  # scale_y_log10(breaks=seq(0,100,20)) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5))
# ggsave(filename=paste(figpath,'totalfunvspecneed_',crop,'.png'))

ggplot(data=df_fun1, aes(x=visits, y=(minsize/rich*100))) +
  geom_point() +
  geom_smooth(method='loess') +
  labs(x="Total abundance of all species",
       y="Percent of species needed to meet threshold",
       title=paste("Effect of bee abundance on species needed (",cropnames[[crop]],")",sep='')) +
  coord_cartesian(ylim=c(0,100)) +
  # scale_y_log10(breaks=seq(0,100,20)) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5))
# ggsave(filename=paste(figpath,'abundvspecneed_',crop,'.png'))

## plot all three crops w/ multipanels
ggplot(data=df_funall, aes(x=fun, y=(minsize/rich*100))) +
  geom_point() +
  geom_smooth(method='loess') +
  geom_vline(data=df_funall, aes(xintercept=threshold), color='red', lty=3, size=1) +
  facet_wrap(~ crop, scales='free',labeller=labeller(crop = cropnames)) +
  labs(x="Total pollination provided by all species",
       y="Percent of species needed to meet threshold",
       title=paste("Effect of total function on proportion of species needed to meet threshold",sep='')) +
  coord_cartesian(ylim=c(0,100)) +
  # scale_y_log10(breaks=seq(0,100,20)) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename=paste(figpath,'totalfunvspecneed_all.png'),width=8.5, height=4)

ggplot(data=df_funall, aes(x=visits, y=(minsize/rich*100))) +
  geom_point() +
  geom_smooth(method='loess') +
  facet_wrap(~ crop, scales='free',labeller=labeller(crop = cropnames)) +
  labs(x="Total abundance of all species",
       y="Percent of species needed to meet threshold",
       title=paste("Effect of bee abundance on proportion of species needed to meet threshold",sep='')) +
  coord_cartesian(ylim=c(0,100)) +
  # scale_y_log10(breaks=seq(0,100,20)) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename=paste(figpath,'abundvspecneed_all.png'),width=8.5, height=4)

## report mean number and proportion of species needed to meet function at a single site
df_funall <- df_funall %>% mutate(propspec = minsize/rich)
df_funall %>% group_by(crop) %>% summarize(specneed = mean(minsize), propneed = mean(propspec))

## divide dates into high and low function
dfmeds <- df_funall %>% group_by(crop) %>% summarize(medfun = median(fun))
df_funall <- left_join(df_funall,dfmeds)
df_funall <- df_funall %>% mutate(abundance = if_else(fun>medfun,'high','low'))
# divide dates into needing all species or not
df_funall <- df_funall %>% mutate(allspec = if_else(minsize==rich,'Yes','No'))

# find how many high vs low sites needed all species
allspectable <- df_funall %>% 
  group_by(crop,abundance,allspec) %>% 
  summarize(num = length(fun)) %>% 
  pivot_wider(names_from = allspec, values_from = num)
allspectable[is.na(allspectable)] <- 0
allspectable <- allspectable %>% mutate(propall = 100*all /(all+notall))
allspectable <- select(allspectable, crop, abundance, propall) %>% 
  pivot_wider(names_from=abundance, values_from=propall)
write.table(allspectable,paste(outpath,'allspectable.csv',sep=''),sep=',',row.names=F)

df_allspec <- df_funall %>% 
  group_by(crop,abundance,allspec) %>% 
  summarize(num = length(fun))

ggplot(data=df_allspec) +
  geom_bar(stat='identity',aes(x=abundance, y=num,fill=allspec)) +
  facet_wrap(~crop, labeller = labeller(crop = cropnames)) +
  labs(x = 'Total bee abundance',
       y = 'Number of site-dates',
       fill = 'All species needed') +
  theme_clean() +
  theme(legend.background = element_rect(color = NA))
ggsave(paste(figpath,'abund_v_specneed.png',sep=''))
