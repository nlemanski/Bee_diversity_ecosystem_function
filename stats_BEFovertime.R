### statistical analysis of BEF over time analysis results, in all crop systems ###

## load packages ##
library(ggplot2)
library(lme4)
library(car)
library(MASS)
library(MuMIn)
library(fitdistrplus)
library(ggfortify) # for autoplot
library(broom) 	# for augment
library(performance)
library(effects)

## file paths ##
figpath <- "C:/Documents/Bee_diversity_ecosystem_function/figures/"
rpath <- "C:/Documents/Bee_diversity_ecosystem_function/R data/"
respath <- "C:/Documents/Bee_diversity_ecosystem_function/BEFresults/"

## load data ##
threshold = 0.5
crops = c("blue","njwat","cawat")

for (c in crops){
  fname1 = paste(rpath,"df_minset_years_",c,'.RData',sep='')
  fname2 = paste(rpath,"df_minset_rounds_",c,'.RData',sep='')
  load(fname1)
  load(fname2)
}

# load across-years data for nj wat with 1 and 3 rounds included
fname3 = paste(rpath,"df_minset_years_njwat_1round",'.RData',sep='')
load(fname3)
df_years_njwat_1r = df_years_njwat
fname4 = paste(rpath,"df_minset_years_njwat_3round",'.RData',sep='')
load(fname4)
df_years_njwat_3r = df_years_njwat

# stats analysis of minimum species vs number of years --------
## visualize blueberry data to see if linear model will fit ##
hist(df_years_blue$minsize)  # data is right tailed
hist(log(df_years_blue$minsize))
bluenorm = fitdistr(df_years_blue$minsize, densfun="normal")
curve(dnorm(x, bluenorm$estimate[1], bluenorm$estimate[2]), col="red", lwd=2, add=T)

# fit glmer with poisson family, random intercept for site
glmmblue = glmer(minsize ~ nyears + (1|site), 
                   data=df_years_blue, family = "poisson")
summary(glmmblue)
Anova(glmmblue)

# look at distribution of residuals
qqnorm(residuals(glmmblue))
abline(0,1)

# check for model assumptions
check_overdispersion(glmmblue) # no overdispersion detected
check_singularity(glmmblue) # no singularity
check_heteroscedasticity(glmmblue) # no heteroscedasticity detected
check_model(glmmblue)

# fit glmm with random slopes and intercepts
glmmblue_full <- glmer(minsize ~ nyears + (1 + nyears | site), 
                         data=df_years_blue, family = "poisson")
isSingular(glmmblue_full)  # singular fit
summary(glmmblue_full)

## compare AIC of full vs intercept only mixed models
aicblue = AIC(glmmblue, glmmblue_full)
aicblue = aicblue[order(aicblue$AIC),]
aicblue  # intercept only model has lower AIC
compareblue = compare_performance(glmmblue, glmmblue_full, rank=T)
write.csv(compareblue, paste(respath,"modelcompare_years_blue.csv",sep=''))
AICtable <- data.frame(aicblue) %>% mutate(mod = row.names(aicblue), crop = 'blue', timescale = 'across')

# save output of glmm w/ lowest AIC
fnamecoeff = paste(respath,"glmm_coeff_","blue","years",".csv",sep="")
fnameanova = paste(respath,"glmm_anova_","blue","years",".csv",sep="")
write.csv(summary(glmmblue)$coefficients,fnamecoeff)
write.csv(Anova(glmmblue),fnameanova)
r2(glmmblue)

effblue = as.data.frame(effect(term="nyears", mod=glmmblue))
write.csv(effblue,paste(respath,'blue_years_effects.csv',sep=''))

# find factor increase and percent increase in species needed w/ number of years over whole time period
pincyb = effblue %>% summarize(finc = last(fit)/first(fit), pinc = 100*(last(fit)-first(fit))/first(fit))

# put results from best model in data frame dfr2
dfr2 = data.frame(crop='blue', timescale = 'across', model='int_only', r2m=r2(glmmblue)[2], r2c=r2(glmmblue)[1],
                  slope = summary(glmmblue)$coefficients[2], slopese = summary(glmmblue)$coefficients[4], 
                  int=summary(glmmblue)$coefficients[1], intse = summary(glmmblue)$coefficients[3],
                  percinc = pincyb$pinc)

## fit a glmm for NJ watermelon  with number of years as x var and number of species in minimum set as y var ##
# stats with 1 sampling round per year
# visualize data
hist(df_years_njwat_1r$minsize) # data is sort of bimodal
hist(log(df_years_njwat_1r$minsize))

# glmm with site as random intercept, using data for nj watermelon w/ only 1 sampling round per year
glmmnjwat <- glmer(minsize ~ nyears + (1|site),
                   data=df_years_njwat_1r, family = "poisson")
summary(glmmnjwat)

# check for model assumptions
check_singularity(glmmnjwat)
check_heteroscedasticity(glmmnjwat)  # no heterosc. detected, p = 0.208
check_overdispersion(glmmnjwat) # no overdisp detected, disp ratio = 0.516, pearson's chi sq = 56.226, p = 1
check_model(glmmnjwat)

# nj watermelon model with site as random slope and intercept
glmmnjwatfull <- glmer(minsize ~ nyears + (1 + nyears | site),
                   data=df_years_njwat_1r, family = "poisson")
summary(glmmnjwatfull)

check_singularity(glmmnjwatfull) # fit not singular
check_overdispersion(glmmnjwatfull) # no overdispersion
check_heteroscedasticity(glmmnjwatfull) # not heterosced.

# compare models for nj wat data
aicnjwat = AIC(glmmnjwat, glmmnjwatfull)
aicnjwat = aicnjwat[order(aicnjwat$AIC),]
comparenjwat = compare_performance(glmmnjwat, glmmnjwatfull, rank=T)  
comparenjwat # lowest AIC is for glmm w/ random slope and intercept
write.csv(comparenjwat, paste(respath,"modelcompare_years_njwat_1round.csv",sep=''), row.names = F)
# save the aic comparison to a table for Table S6
AICtable <- rbind(AICtable, data.frame(aicnjwat) %>% mutate(mod = row.names(aicnjwat), crop = 'njwat', timescale = 'across'))

# run anova for lowest AIC model and save output to file
# glmm with random slope and intercepts
fnamecoeff = paste(respath,"glmmfull_coeff_","njwat","years_1round",".csv",sep="")
fnameanova = paste(respath,"glmmfull_anova_","njwat","years_1round",".csv",sep="")
write.csv(summary(glmmnjwatfull)$coefficients,fnamecoeff)
write.csv(Anova(glmmnjwatfull),fnameanova)
r2(glmmnjwatfull)

# save estimated effect sizes for lowest aic model to file
effnjwat = as.data.frame(effect(term="nyears", mod=glmmnjwatfull))
write.csv(effnjwat, paste(respath,"njwat_years_effects.csv",sep=''))

# find factor increase and percent increase in species needed w/ number of years over whole time period
pincynj = effnjwat %>% summarize(finc = last(fit)/first(fit), pinc = 100*(last(fit)-first(fit))/first(fit))

# add r2 and effect sizes to dataframe
dfr2 = rbind(dfr2, 
             data.frame(crop='njwat', timescale = 'across', model='full', r2m=r2(glmmnjwatfull)[2], r2c=r2(glmmnjwatfull)[1],
                        slope = summary(glmmnjwatfull)$coefficients[2], slopese = summary(glmmnjwatfull)$coefficients[4], 
                        int=summary(glmmnjwatfull)$coefficients[1], intse = summary(glmmnjwatfull)$coefficients[3],
                        percinc = pincynj$pinc )) %>% remove_rownames()

## fit a glmm for CA watermelon  with number of years as x var and number of species in minimum set as y var ##
hist(df_years_cawat$minsize, breaks=7) # left skewed
hist(log(df_years_cawat$minsize))

# glmm with site as random slope and intercept, poisson family with log link
glmmcawat_full <- glmmTMB(minsize ~ nyears + (1 + nyears | site),
                        data=df_years_cawat, family = "poisson")

# glmm with site as random slope, poisson family with log link
glmmcawat <- glmmTMB(minsize ~ nyears + (1|site),
                   data=df_years_cawat, family = "poisson")

# check for model assumptions
check_singularity(glmmcawat) # singular fit
check_heteroscedasticity(glmmcawat) # no heterosc, detected, p = 0.642
check_overdispersion(glmmcawat) # no overdispersion detected, p=1
check_model(glmmcawat)

# compare models
comparecawat = compare_performance(glmmcawat, glmmcawat_full, rank=T)  # best fit for glmm w/ random effect
comparecawat  # random intercept only mixed model has lower AIC than full model
aiccawat = AIC(glmmcawat, glmmcawat_full)
write.csv(comparecawat,paste(respath,"modelcompare_years_cawat.csv",sep=''), row.names = F)
# save the aic comparison to a table for Table S6
AICtable <- rbind(AICtable, data.frame(aiccawat) %>% mutate(mod = row.names(aiccawat), crop = 'cawat', timescale = 'across'))

# perform anova and save glmm results for lowest aic model to file
fnamecoeff = paste(respath,"glmm_coeff_","cawat","years",".csv",sep="")
fnameanova = paste(respath,"glmm_anova_","cawat","years",".csv",sep="")
summary(glmmcawat)
write.csv(summary(glmmcawat)$coefficients$cond,fnamecoeff)
write.csv(Anova(glmmcawat),fnameanova)
r2(glmmcawat, tolerance = 1e-11)

# save estimated effect sizes to file
effcawat = as.data.frame(effect(term="nyears", mod=glmmcawat))
write.csv(effcawat, paste(respath,"cawat_years_effects.csv",sep=''))

# find factor increase and percent increase in species needed w/ number of years over whole time period
pincyca = effcawat %>% summarize(finc = last(fit)/first(fit), pinc = 100*(last(fit)-first(fit))/first(fit))

# add r2 and effect sizes to dataframe
dfr2 = rbind(dfr2, 
             data.frame(crop='cawat', timescale = 'across', model='int_only', r2m=r2(glmmcawat)[2], r2c=r2(glmmcawat)[1],
                        slope = summary(glmmcawat)$coefficients$cond[2], slopese = summary(glmmcawat)$coefficients$cond[4], 
                        int=summary(glmmcawat)$coefficients$cond[1], intse = summary(glmmcawat)$coefficients$cond[3],
                        percinc = pincyca$pinc )) %>% remove_rownames()


### stats for change in min species with number of rounds ### ---------------

## fit a linear mixed model for blueberry with number of rounds as x var and number of species in minimum set as y var
# view histogram of blueberry data
hist(df_rounds_blue$minsize,breaks=6)
hist(log(df_rounds_blue$minsize),breaks=6)

# glmm with site as random intercept
glmmblueround <- glmer(minsize ~ nrounds + (1|site),
                data=df_rounds_blue, family='poisson')
summary(glmmblueround)

# visualize model fit
qqnorm(residuals(glmmblueround))
hist(residuals(glmmblueround))
check_model(glmmblueround)
check_singularity(glmmblueround)
isSingular(glmmblueround)
check_heteroscedasticity(glmmblueround) # no heterosced., p = 0.096
check_overdispersion(glmmblueround) # no overdispersion, p = 1

# glmm with site as random slope and intercept
glmmblueround_full <- glmer(minsize ~ nrounds + (1 + nrounds |site),
                       data=df_rounds_blue, family='poisson')
check_singularity(glmmblueround_full)
isSingular(glmmblueround_full)
check_heteroscedasticity(glmmblueround_full)
check_overdispersion(glmmblueround_full)

# compare models
aicbluernd = AIC(glmmblueround, glmmblueround_full)
compareblueround = compare_performance(glmmblueround, glmmblueround_full, rank=T) # model with random effect has better score
plot(compareblueround)
compareblueround # AIC value lower for random intercept only glmm
write.csv(compareblueround,paste(respath,"modelcompare_rounds_blue.csv",sep=''), row.names = F)
AICtable <- rbind(AICtable, data.frame(aicbluernd) %>% mutate(mod = row.names(aicbluernd), crop = 'blue', timescale = 'within'))

# save output of lowest aic model to file
summary(glmmblueround)
fnamecoeff = paste(respath,"glmm_coeff_","blue","rounds",".csv",sep="")
fnameanova = paste(respath,"glmm_anova_","blue","rounds",".csv",sep="")
write.csv(summary(glmmblueround)$coefficients,fnamecoeff)
write.csv(Anova(glmmblueround),fnameanova)
r2(glmmblueround)
dfr2rd = data.frame(crop='blue', model='int_only', r2m=r2(glmmblueround)[2], row.names=NULL)

effblueround = as.data.frame(effect(term="nrounds", mod=glmmblueround))
write.csv(effblueround, paste(respath,"blue_rounds_effects.csv",sep=''))

percincblueround = (effblueround$fit[5]-effblueround$fit[1])/effblueround$fit[1]*100

# add r2 and effect sizes for blueberry within years model to dataframe
dfr2 = rbind(dfr2, 
             data.frame(crop='blue', timescale = 'within', model='int_only', r2m=r2(glmmblueround)[2], r2c=r2(glmmblueround)[1],
                        slope = summary(glmmblueround)$coefficients[2], slopese = summary(glmmblueround)$coefficients[4], 
                        int=summary(glmmblueround)$coefficients[1], intse = summary(glmmblueround)$coefficients[3],
                        percinc = percincblueround)) %>% remove_rownames()


## fit a glmm for NJ watermelon with number of rounds as x var and number of species in minimum set as y var
# view nj watermelon data
hist(df_rounds_njwat$minsize)
hist(log(df_rounds_njwat$minsize))

# fit glmm with random intercept for effect of site
glmmnjwatround <- glmer(minsize ~ nrounds + (1|site), 
                        data = df_rounds_njwat, family = 'poisson')

# view model fit
hist(residuals(glmmnjwatround))
check_model(glmmnjwatround)
check_singularity(glmmnjwatround)
check_heteroscedasticity(glmmnjwatround) # no heterosced., p = 0.145
check_overdispersion(glmmnjwatround) # overdispersion detected

# fit full glmm with random slope intercept for effect of site
glmmnjwatroundfull <- glmer(minsize ~ nrounds + (1 + nrounds |site),
                        data = df_rounds_njwat, family = 'poisson')
check_overdispersion(glmmnjwatroundfull) # overdispersion detected 

# fit glmm with negative binomial to account for overdispersion
glmmnjwatround2 <- glmer.nb(minsize ~ nrounds + (1|site),
                         data = df_rounds_njwat)

# view model fit
hist(residuals(glmmnjwatround2))
check_model(glmmnjwatround2)
check_singularity(glmmnjwatround2)

# fit glmm full model with negative binomial to account for overdispersion
glmmnjwatround2full <- glmer.nb(minsize ~ nrounds + (1 + nrounds |site),
                            data = df_rounds_njwat)
check_model(glmmnjwatround2full)
isSingular(glmmnjwatround2full) # singular fit

# compare models
comparenjwatround = compare_performance(glmmnjwatround, glmmnjwatroundfull, glmmnjwatround2, glmmnjwatround2full,
                    rank=T) 
comparenjwatround 
aicnjwatrnd = AIC(glmmnjwatround, glmmnjwatroundfull, glmmnjwatround2, glmmnjwatround2full)
aicnjwatrnd = aicnjwatrnd[order(aicnjwatrnd$AIC),] # poisson model with random intercept only has lowest aic
write.csv(comparenjwatround,paste(respath,"modelcompare_rounds_njwat.csv",sep=''), row.names = F)
AICtable <- rbind(AICtable, data.frame(aicnjwatrnd) %>% mutate(mod = row.names(aicnjwatrnd), crop = 'nwjat', timescale = 'within'))

# save results from negative binomial due to overdispersion
fnamecoeff = paste(respath,"glmm_coeff_negbin_","njwat","rounds",".csv",sep="")
fnameanova = paste(respath,"glmm_anova_negbin_","njwat","rounds",".csv",sep="")
write.csv(summary(glmmnjwatround2)$coefficients,fnamecoeff)
write.csv(Anova(glmmnjwatround2),fnameanova)
r2(glmmnjwatround2)


effnjwatrounds = as.data.frame(effect(term = "nrounds", mod = glmmnjwatround)) # predicted values from glmm w/ poisson
effnjwatrounds2 = as.data.frame(effect(term = "nrounds", mod = glmmnjwatround2)) # predicted values from neg bin 
write.csv(effnjwatrounds2, paste(respath, "njwat_effects_rounds.csv",sep=''))

percincnjwatround = (effnjwatrounds2$fit[5]-effnjwatrounds2$fit[1])/effnjwatrounds2$fit[1]*100 # est 76% increase over all 3 rounds for nj watermelon

# add r2 value and effect sizes for nj watermelon within years to dataframe
dfr2 = rbind(dfr2, 
             data.frame(crop='njwat', timescale = 'within', model='int_only', r2m=r2(glmmnjwatround2)[2], r2c=r2(glmmnjwatround2)[1],
                        slope = summary(glmmnjwatround2)$coefficients[2], slopese = summary(glmmnjwatround2)$coefficients[4], 
                        int=summary(glmmnjwatround2)$coefficients[1], intse = summary(glmmnjwatround2)$coefficients[3],
                        percinc = percincnjwatround )) %>% remove_rownames()

## fit a glmm for CA watermelon with number of rounds as x var and number of species in minimum set as y var
# view data
hist(df_rounds_cawat$minsize,breaks=8)
hist(log(df_rounds_cawat$minsize),breaks=6)

# site as random effect w/ random slope and intercept (full model)
glmmcawatroundfull <- glmer(minsize ~ nrounds + (1 + nrounds | site),
                        data=df_rounds_cawat, family='poisson')  # singular fit
summary(glmmcawatroundfull)

# site as random effect, intercept only
glmmcawatround <- glmer(minsize ~ nrounds + (1|site),
                data=df_rounds_cawat, family='poisson')
summary(glmmcawatround)

# view model fit
hist(residuals(glmmcawatround))
check_model(glmmcawatround)
check_singularity(glmmcawatround)
check_heteroscedasticity(glmmcawatround)  # hetero detected
check_overdispersion(glmmcawatround) # no overdispersion, dispersion ratio =  0.477, Chi-Squared = 93, p-value = 1

# compare model fit
comparecawatround = compare_performance(glmmcawatround, glmmcawatroundfull, rank=T)
plot(comparecawatround)
comparecawatround  # model w/ only random intercept has lowest aic
aiccawatrnd = AIC(glmmcawatround, glmmcawatroundfull)
aiccawatrnd
write.csv(comparecawatround,paste(respath,"modelcompare_rounds_cawat.csv",sep=''), row.names = F)
AICtable <- rbind(AICtable, data.frame(aiccawatrnd) %>% mutate(mod = row.names(aiccawatrnd), crop = 'cawat', timescale = 'within'))

# anova and save output of model with lowest aic
summary(glmmcawatround)
fnamecoeff = paste(respath,"glmm_coeff_","cawat","rounds",".csv",sep="")
fnameanova = paste(respath,"glmm_anova_","cawat","rounds",".csv",sep="")
write.csv(summary(glmmcawatround)$coefficients,fnamecoeff)
write.csv(Anova(glmmcawatround),fnameanova)
r2(glmmcawatround) # get r^2 values
# dfr2rd = rbind(dfr2rd, data.frame(crop='cawat', model='int_only', r2m=r2(glmmcawatround)[2], row.names=NULL))

# tab_model(glmmcawatround, file=paste(figpath,'glm_table_rnds_cawat')) # create a table of model output

effcawatround = as.data.frame(effect(term="nrounds", mod=glmmcawatround))
write.csv(effcawatround, paste(respath,"cawat_rounds_effects.csv",sep=''))
percinccawatround = (effcawatround$fit[5]-effcawatround$fit[1])/effcawatround$fit[1]*100 # 

# add r2 value and effect sizes for ca watermelon within years to dataframe
dfr2 = rbind(dfr2, 
             data.frame(crop='cawat', timescale = 'within', model='int_only', r2m=r2(glmmcawatround)[2], r2c=r2(glmmcawatround)[1],
                        slope = summary(glmmcawatround)$coefficients[2], slopese = summary(glmmcawatround)$coefficients[4], 
                        int=summary(glmmcawatround)$coefficients[1], intse = summary(glmmcawatround)$coefficients[3],
                        percinc = percinccawatround )) %>% remove_rownames()

## save the dataframe w/ effect sizes and r2 values and df with aic values
write.csv(dfr2, paste(respath, "effects_r2_all.csv", sep=''))
write.csv(AICtable, paste(respath, "aic_table_all.csv", sep=''))

### report number of species needed to meet threshold for single site-date --------
df_rounds_njwat = df_rounds_njwat[,c(1:6,8:14)]
df_rounds = rbind(df_rounds_blue, df_rounds_njwat, df_rounds_cawat)
minset1r <- df_rounds %>% 
  filter(nrounds==1) %>%
  group_by(crop) %>% 
  summarise(spec_mean = mean(minsize), spec_sd = sd(minsize))
write.csv(minset1r, paste(respath,'minset_1sitedate.csv',sep=''))
