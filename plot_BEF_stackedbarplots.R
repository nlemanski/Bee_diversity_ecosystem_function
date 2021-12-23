## make stacked barplots showing function contributions by each species at each site-date ##

## load packages ##
library(plyr)
library(tidyverse)
library(ggplot2)
library(tidytext)
library(forcats)
library(RColorBrewer)

## file paths ##
figpath <- "C:/Documents/winfree lab/figures"
rpath <- "C:/Documents/winfree lab/R data/"

## import crop visit data as RData files ##
crops = c("blue","njwat","cawat")
cropnames = c("Blueberry", "NJ Watermelon", "CA Watermelon")

for (c in crops){
  fname = paste(rpath,"df_visits_",c,".RData",sep='')
  load(fname)
}

## function threshold percentage
func_percent = 0.50

## make stacked bar plots showing total function contributed by each species at each site-year ##
## blueberry
# arrange data by date
df_date_blue = df_visits_blue %>% 
  unite("date",c(year,round))

# create a variable to order bars by descending height for each site and join it to dataframe
df_allsp = df_date_blue %>%
            group_by(date, site) %>%
            summarize(totalfun = sum(fun)) %>%
            arrange(site, -totalfun)
df_allsp = df_allsp %>% 
            group_by(site) %>%
            mutate(order = row_number())
df_date_blue = df_date_blue %>% left_join(df_allsp, by = c("date"="date", "site"="site"))
df_date_blue = df_date_blue %>% 
                arrange(site, order, -fun)

# function threshold value (aggregate pollen grains) from percentage, using mean across all site-year-rounds
df_totalfun = ddply(df_visits_blue, .(round, year, site), summarize, totalfun = sum(fun))
func_level_blue = mean(df_totalfun$totalfun, na.rm=T)*func_percent

#create custom color scale
nspec = length(unique(df_date_blue$gen_sp))
# spec = sort(unique(df_date_blue$gen_sp))
# mycolors = rainbow(nspec)
# names(mycolors) = sort(unique(df_date_blue$gen_sp))
mycolors = sample(colors(nspec)[1:nspec])
colors = c("green4", "green3", "mediumpurple", "mediumpurple2", "mediumpurple4", "red", "red3",
               "orangered", "orange" , "sienna", "sienna3", "darkorange", "orange3", "chocolate1",
               "chocolate", "grey20" ,"black", "aquamarine1", "aquamarine3","brown4", "coral4",
               "orangered3", "blueviolet" , "aquamarine" , "tomato", "mediumblue","burlywood2", "brown1",
               "blue4", "steelblue" , "royalblue" , "blue" , "steelblue2","cadetblue4", "skyblue3",
               "skyblue1", "skyblue4", "hotpink4", "powderblue", "lavender", "gold3", "brown3",
               "gold", "deeppink", "goldenrod1", "brown", "slateblue1", "goldenrod", "hotpink",
               "cadetblue1", "orchid")

# stacked barplots of functional contributions of species at each site-date
ggplot(data = df_date_blue, aes(fill = fct_reorder(gen_sp, fun, mean), y = fun, x = order)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(~site, scales = "free") +
  geom_hline(yintercept = func_level_blue, color = "red") +
  theme(legend.position = "bottom", axis.text.x = element_blank()) +
  labs(x = "Date", y = "Function Contribution by Species", fill = "Species") +
  # scale_fill_manual(name=spec, values=mycolors)
  scale_fill_manual(values = colors)

figname = "stackedbar_blue.png"
ggsave(filename=figname, device="png", path=figpath, width=10, height=7.5)

# stacked barplots with species stacked alphabetically
nspec = length(unique(df_date_blue$gen_sp))
ggplot(data = df_date_blue, aes(fill = gen_sp, y = fun, x = order)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(~site, scales = "free") +
  geom_hline(yintercept = func_level_blue, color = "red") +
  theme(legend.position = "bottom", axis.text.x = element_blank()) +
  labs(x = "Date", y = "Function Contribution by Species", fill = "Species") +
  scale_fill_manual(values=colors)

figname = "stackedbar_blue_v2.png"
ggsave(filename=figname, device="png", path=figpath, width=10, height=7.5)

## nj watermelon

# arrange data by date
df_date_njwat = df_visits_njwat %>% 
  unite("date",c(year,round))

# create a variable to order bars by descending height for each site and join it to dataframe
df_allsp = df_date_njwat %>%
  group_by(date, site) %>%
  summarize(totalfun = sum(fun)) %>%
  arrange(site, -totalfun)
df_allsp = df_allsp %>% 
  group_by(site) %>%
  mutate(order = row_number())
df_date_njwat = df_date_njwat %>% left_join(df_allsp, by = c("date"="date", "site"="site"))
df_date_njwat = df_date_njwat %>% 
  arrange(site, order, -fun)

# function threshold value (aggregate pollen grains) from percentage, using mean across all site-year-rounds
df_totalfun = ddply(df_visits_njwat, .(round, year, site), summarize, totalfun = sum(fun))
func_level_njwat = mean(df_totalfun$totalfun, na.rm=T)*func_percent

#create custom color scale
nspec = length(unique(df_date_njwat$gen_sp))
mycolors = sample(colors(nspec)[1:nspec])

# stacked barplots of functional contributions of species at each site-date
ggplot(data = df_date_njwat, aes(fill = fct_reorder(gen_sp, fun, mean), y = fun, x = order)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(~site, scales = "free") +
  geom_hline(yintercept = func_level_njwat, color = "red") +
  theme(legend.position = "bottom", legend.text = element_text(size=7), legend.title = element_text(size=9),
        axis.text.x = element_blank()) +
  labs(x = "Date", y = "Function Contribution by Species", fill = "Species") +
  scale_fill_manual(values=mycolors)

figname = "stackedbar_njwat.png"
ggsave(filename=figname, device="png", path=figpath, width=10, height=8.5)

# stacked barplots with species stacked alphabetically
ggplot(data = df_date_njwat, aes(fill = gen_sp, y = fun, x = order)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(~site, scales = "free") +
  geom_hline(yintercept = func_level_njwat, color = "red") +
  theme(legend.position = "bottom", legend.text = element_text(size=7), legend.title = element_text(size=9),
        axis.text.x = element_blank()) +
  labs(x = "Date", y = "Function Contribution by Species", fill = "Species")

figname = "stackedbar_njwat_v2.png"
ggsave(filename=figname, device="png", path=figpath)


## western watermelon

# arrange data by date
df_date_cawat = df_visits_cawat %>% 
  unite("date",c(year,round))

# create a variable to order bars by descending height for each site and join it to dataframe
df_allsp = df_date_cawat %>%
  group_by(date, site) %>%
  summarize(totalfun = sum(fun)) %>%
  arrange(site, -totalfun)
df_allsp = df_allsp %>% 
  group_by(site) %>%
  mutate(order = row_number())
df_date_cawat = df_date_cawat %>% left_join(df_allsp, by = c("date"="date", "site"="site"))
df_date_cawat = df_date_cawat %>% 
  arrange(site, order, -fun)

# function threshold value (aggregate pollen grains) from percentage, using mean across all site-year-rounds
df_totalfun = ddply(df_visits_cawat, .(round, year, site), summarize, totalfun = sum(fun))
func_level_cawat = mean(df_totalfun$totalfun, na.rm=T)*func_percent

# stacked barplots of functional contributions of species at each site-date
nspec = length(unique(df_date_cawat$gen_sp))
mycolors = sample(rainbow(nspec)[1:nspec])
ggplot(data = df_date_cawat, aes(fill = fct_reorder(gen_sp, fun, mean), y = fun, x = order)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(~site, scales = "free") +
  geom_hline(yintercept = func_level_cawat, color = "red") +
  theme(legend.position = "bottom", legend.text = element_text(size=7), legend.title = element_text(size=9),
        axis.text.x = element_blank()) +
  labs(x = "Date", y = "Function Contribution by Species", fill = "Species") +
  scale_fill_manual(values=mycolors)

figname = "stackedbar_cawat.png"
ggsave(filename=figname, device="png", path=figpath, width=10, height=8.5)

# stacked barplots with species stacked alphabetically
ggplot(data = df_date_cawat, aes(fill = gen_sp, y = fun, x = order)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(~site, scales = "free") +
  geom_hline(yintercept = func_level_cawat, color = "red") +
  theme(legend.position = "bottom", legend.text = element_text(size=7), legend.title = element_text(size=9),
        axis.text.x = element_blank()) +
  labs(x = "Date", y = "Function Contribution by Species", fill = "Species")

figname = "stackedbar_cawat_v2.png"
ggsave(filename=figname, device="png", path=figpath)
