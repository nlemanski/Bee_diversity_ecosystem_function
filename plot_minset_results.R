## plot results of minimum set analysis for BEF over time ##
# edited NJL, 4/26/21

## load packages ##
library(tidyverse)
library(lme4)
library(effects)

## file paths and file names ##
figpath <- "D:/natal/D_Documents/winfree lab/figures"
rpath <- "D:/natal/D_Documents/winfree lab/R data/"
outpath <- "D:/natal/D_Documents/winfree lab/BEFresults/"

## load results files ##
crops = c("blue","cran","njwat","cawat")
cropnames = c("Blueberry", "Cranberry", "Eastern Watermelon", "Western Watermelon")
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

## view as scatterplot: size of minimum set vs number of rounds ##
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

# cranberry 
ggplot(data = df_rounds_cran, mapping = aes(x = nrounds, y = minsize, color = site)) +
  geom_jitter(width = 0.1, height = 0.1) +
  geom_smooth(mapping = aes(x = nrounds, y = minsize, color = site), method = "lm", alpha = .08) +
  labs(title = paste("Minimum Species Needed for ",threshold*100,"% of Pollination (",cropnames[2],")",sep=''),
       x= "Number of Rounds",
       y = "Number of Species",
       color= "Site") +
  coord_cartesian(xlim=c(0.7,2.2), ylim=c(0,15)) +
  scale_x_continuous(breaks=seq(0,3,1)) +
  scale_y_continuous(breaks=seq(0,16,2)) +
  theme_light() +     
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = "bottom") +
  guides(col = guide_legend(nrow = 2))

figname = paste('minset_rounds_','cran','_',threshold*100,'.png',sep='')
ggsave(filename=figname, device="png", path=figpath, width=6.5,height=5)

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

### view scatterplot with a single trend line based on glmm estimates ###
## min set vs number of rounds (within-years) with site-years as reps

# blueberry glmm
glmmblueround <- glmer(minsize ~ nrounds + (1|site),
                       data=df_rounds_blue, family='poisson')
effblue = as.data.frame(effect(term="nrounds", mod=glmmblueround))

ggplot() +
  geom_jitter(data = df_rounds_blue, mapping = aes(x = nrounds, y = minsize, color = site), width = 0.1, height = 0.1) +
  geom_line(data=effblue, aes(x=nrounds, y=fit), color="black", linetype="dashed") +
  geom_ribbon(data=effblue, aes(x=nrounds, ymin=lower, ymax=upper), alpha= 0.3, fill="black") +
  # geom_abline(slope=slb, intercept=intb) +
  labs(title = paste("Minimum Species Needed for ",threshold*100,"% of Pollination (",cropnames[1],")",sep=''),
       x= "Number of Rounds",
       y = "Number of Species",
       color= "Site") +
  coord_cartesian( ylim=c(0,8)) +
  scale_x_continuous(breaks=seq(0,3,1)) +
  scale_y_continuous(breaks=seq(0,10,2)) +
  theme_light() +     
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = "bottom") +
  guides(col = guide_legend(nrow = 2))

figname = paste(figpath,'/minset_rounds_','blue','_',threshold*100,'_glmmline.png',sep='')
ggsave(filename=figname, device="png", width=6.5,height=5)

# nj watermelon glmm
glmmnjwatround <- glmer.nb(minsize ~ nrounds + (1|site),  # neg bin was used because of overdispersion
                           data = df_rounds_njwat)

effnjwat = as.data.frame(effect(term="nrounds", mod=glmmnjwatround))

ggplot() +
  geom_jitter(data = df_rounds_njwat, mapping = aes(x = nrounds, y = (minsize), color = site), width = 0.1, height = 0.1) +
  geom_line(data=effnjwat, aes(x=nrounds, y=fit), color="black", linetype="dashed") +
  geom_ribbon(data=effnjwat, aes(x=nrounds, ymin=lower, ymax=upper), alpha= 0.3, fill="black") +
  labs(title = paste("Minimum Species Needed for ",threshold*100,"% of Pollination (",cropnames[3],")",sep=''),
       x= "Number of Rounds",
       y = "Number of Species",
       color= "Site") +
  coord_cartesian( ylim=c(0,16) ) +
  scale_x_continuous(breaks=seq(0,3,1)) +
  scale_y_continuous(breaks=seq(0,16,4)) +
  theme_light() +     
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = "bottom") +
  guides(col = guide_legend(nrow = 3))

figname = paste(figpath,'/minset_rounds_','njwat','_',threshold*100,'_glmmline.png',sep='')
ggsave(filename=figname, device="png", width=6.5,height=5)

# ca watermelon glmm
glmmcawatround <-  glmer(minsize ~ nrounds + (1|site),
                         data=df_rounds_cawat, family='poisson')

effcawat = as.data.frame(effect(term="nrounds", mod=glmmcawatround))

ggplot() +
  geom_jitter(data = df_rounds_cawat, mapping = aes(x = nrounds, y = minsize, color = site), width = 0.1, height = 0.1) +
  geom_line(data=effcawat, aes(x=nrounds, y=fit), color="black", linetype="dashed") +
  geom_ribbon(data=effcawat, aes(x=nrounds, ymin=lower, ymax=upper), alpha= 0.3, fill="black") +
  labs(title = paste("Minimum Species Needed for ",threshold*100,"% of Pollination (",cropnames[4],")",sep=''),
       x= "Number of Rounds",
       y = "Number of Species",
       color= "Site") +
  coord_cartesian( ylim=c(0,10) ) +
  scale_x_continuous(breaks=seq(0,9,1)) +
  scale_y_continuous(breaks=seq(0,10,2)) +
  theme_light() +     
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = "bottom", legend.text=element_text(size=6), ) +
  guides(col = guide_legend(nrow = 4))

figname = paste(figpath,'/minset_rounds_','cawat','_',threshold*100,'_',rmax,'round_glmmline.png',sep='')
ggsave(filename=figname, device="png", width=7.5,height=5)

### view scatterplot with a single trend line based on glmm estimates ###
## min set vs number of years (across-years) with all rounds included per year

# blueberry glmm
glmmblueyear <- glmer(minsize ~ nyears + (1|site),
                       data=df_years_blue, family='poisson')

effblue = as.data.frame(effect(term="nyears", mod=glmmblueyear))

ggplot() +
  geom_jitter(data = df_years_blue, mapping = aes(x = nyears, y = minsize, color = site), width = 0.1, height = 0.1) +
  geom_line(data=effblue, aes(x=nyears, y=fit), color="black", linetype="dashed") +
  geom_ribbon(data=effblue, aes(x=nyears, ymin=lower, ymax=upper), alpha= 0.3, fill="black") +
  labs(title = paste("Minimum Species Needed for ",threshold*100,"% of Pollination (",cropnames[1],")",sep=''),
       x= "Number of Years",
       y = "Number of Species",
       color= "Site") +
  coord_cartesian( ylim=c(0,11)) +
  scale_x_continuous(breaks=seq(0,3,1)) +
  scale_y_continuous(breaks=seq(0,12,2)) +
  theme_light() +     
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = "bottom") +
  guides(col = guide_legend(nrow = 2))

figname = paste(figpath,'/minset_years_','blue','_',threshold*100,'_glmmline.png',sep='')
ggsave(filename=figname, device="png", width=6.5,height=5)

# nj watermelon glmm
glmmnjwatyear <-  glmer(minsize ~ nyears + (1|site),
                        data=df_years_njwat, family = "poisson")

effnjwat = as.data.frame(effect(term="nyears", mod=glmmnjwatyear))

ggplot() +
  geom_jitter(data = df_years_njwat, mapping = aes(x = nyears, y = (minsize), color = site), width = 0.1, height = 0.1) +
  geom_line(data=effnjwat, aes(x=nyears, y=fit), color="black", linetype="dashed") +
  geom_ribbon(data=effnjwat, aes(x=nyears, ymin=lower, ymax=upper), alpha= 0.3, fill="black") +
  labs(title = paste("Minimum Species Needed for ",threshold*100,"% of Pollination (",cropnames[3],")",sep=''),
       x= "Number of Years",
       y = "Number of Species",
       color= "Site") +
  coord_cartesian( ylim=c(0,20) ) + 
  scale_x_continuous(breaks=seq(0,6,1)) +   
  scale_y_continuous(breaks=seq(0,20,4)) +
  theme_light() +     
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = "bottom") +
  guides(col = guide_legend(nrow = 3))

figname = paste(figpath,'/minset_years_','njwat','_',threshold*100,'_',rmax,'round_glmmline.png',sep='')
ggsave(filename=figname, device="png", width=6.5,height=5)

# ca watermelon glmm
glmmcawatyear <-  glmer(minsize ~ nyears + (1|site),
                        data=df_years_cawat, family = "poisson")

effcawat = as.data.frame(effect(term="nyears", mod=glmmcawatyear))

ggplot() +
  geom_jitter(data = df_years_cawat, mapping = aes(x = nyears, y = minsize, color = site), width = 0.1, height = 0.1) +
  geom_line(data=effcawat, aes(x=nyears, y=fit), color="black", linetype="dashed") +
  geom_ribbon(data=effcawat, aes(x=nyears, ymin=lower, ymax=upper), alpha= 0.3, fill="black") +
  labs(title = paste("Minimum Species Needed for ",threshold*100,"% of Pollination (",cropnames[4],")",sep=''),
       x= "Number of Years",
       y = "Number of Species",
       color= "Site") +
  coord_cartesian( ylim=c(0,8) ) +
  scale_x_continuous(breaks=seq(0,9,1)) +
  scale_y_continuous(breaks=seq(0,8,2)) +
  theme_light() +     
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size=13), legend.position = "bottom", legend.text=element_text(size=6), ) +
  guides(col = guide_legend(nrow = 4))

figname = paste(figpath,'/minset_years_','cawat','_',threshold*100,'_glmmline.png',sep='')
ggsave(filename=figname, device="png", width=6.5,height=5)

