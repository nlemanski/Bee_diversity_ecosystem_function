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
# library(lmerTest)

## file paths ##
figpath <- "C:/Documents/winfree lab/figures/"
rpath <- "C:/Documents/winfree lab/R data/"
respath <- "C:/Documents/winfree lab/BEFresults/"

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

## functions ##
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

### stats analysis of minimum species needed for 50% pollination  vs. number of years in analysis ###

## visualize blueberry data to see if linear model will fit ##
hist(df_years_blue$minsize)  # data is right tailed
hist(log(df_years_blue$minsize))
bluenorm = fitdistr(df_years_blue$minsize, densfun="normal")
curve(dnorm(x, bluenorm$estimate[1], bluenorm$estimate[2]), col="red", lwd=2, add=T)

## fit lmm for blueberry with number of years as x var and number of species in minimum set as y var ##
lmer_blue <- lmer(minsize ~ nyears + (1|site),
                data=df_years_blue)
summary(lmer_blue)
Anova(lmer_blue)

qqnorm(residuals(lmer_blue))
abline(0,1)

# plot fitted values vs residuals
plot(fitted(lmer_blue), residuals(lmer_blue), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(2), residuals(lmer_blue)))

# check model assumptions
check_singularity(lmer_blue) # no singularity
check_heteroscedasticity(lmer_blue) # no heteroscedasticity detected
check_model(lmer_blue) # residuals look bimodal
plot(lmer_blue)

## fit a linear model with no random effects
lmblue2 <- lm(minsize ~ nyears, data=df_years_blue)
check_model(lmblue2)
autoplot(lmblue2)

## fit a glmm with a poisson distribution ##
# poisson makes sense because data is based on counts and has a lower bound of 1
glm_blue = glm(minsize ~ nyears, data=df_years_blue, family = "poisson")
summary(glm_blue)
Anova(glm_blue)

# look at distribution of residuals
qqnorm(residuals(glm_blue))
abline(0,1)
autoplot(glm_blue)
hist(glm_blue$resid)

# fit glmer with poisson family, random effect of site
glmer_blue = glmer(minsize ~ nyears + (1|site), 
                   data=df_years_blue, family = "poisson")
summary(glmer_blue)
Anova(glmer_blue)

# look at distribution of residuals
qqnorm(residuals(glmer_blue))
abline(0,1)
plot(glmer_blue)

# check for model assumptions
check_overdispersion(glmer_blue) # no overdispersion detected
check_singularity(glmer_blue) # no singularity
check_heteroscedasticity(glmer_blue) # no heteroscedasticity detected
check_model(glmer_blue)

# fit glmm with random slopes and intercepts
glmer_blue_full <- glmer(minsize ~ nyears + (1 + nyears | site), 
                         data=df_years_blue, family = "poisson")
isSingular(glmer_blue_full)  # singular fit

summary(glmer_blue_full)

## compare AIC of different models
aicblue = AIC(lmblue2, lmer_blue, glm_blue, glmer_blue, glmer_blue_full)
aicblue = aicblue[order(aicblue$AIC),]
# linear model with site as random effect has lowest AIC
compareblue = compare_performance(lmblue2, lmer_blue, glm_blue, glmer_blue, glmer_blue_full, rank=T)
write.csv(compareblue, paste(respath,"modelcompare_years_blue.csv",sep=''))

# save output of lmm
fnamecoeff = paste(respath,"lmm_coeff_","blue","years",".csv",sep="")
fnameanova = paste(respath,"lmm_anova_","blue","years",".csv",sep="")
write.csv(summary(lmer_blue)$coefficients,fnamecoeff)
write.csv(Anova(lmer_blue),fnameanova)

# save output of glmm
fnamecoeff = paste(respath,"glmm_coeff_","blue","years",".csv",sep="")
fnameanova = paste(respath,"glmm_anova_","blue","years",".csv",sep="")
write.csv(summary(glmer_blue)$coefficients,fnamecoeff)
write.csv(Anova(glmer_blue),fnameanova)

effblue = as.data.frame(effect(term="nyears", mod=glmer_blue))
write.csv(effblue,paste(respath,'blue_years_effects.csv',sep=''))


## fit a lmm for NJ watermelon  with number of years as x var and number of species in minimum set as y var ##
# stats with 1 sampling round per year
# visualize data
hist(df_years_njwat_1r$minsize) # data is sort of bimodal
hist(log(df_years_njwat_1r$minsize))

# lmm with site as random effect
lmmnjwat <- lmer(minsize ~ nyears + (1|site),
                 data=df_years_njwat_1r)
summary(lmmnjwat)
plot(lmmnjwat)
qqnorm(residuals(lmmnjwat))
abline(0,1)

# check for model assumptions
check_singularity(lmmnjwat)
check_heteroscedasticity(lmmnjwat)  # no heterosc. detected, p = 0.894
check_model(lmmnjwat) # residuals look super normal

# lm
lmnjwat <- lm(minsize ~ nyears,
              data=df_years_njwat_1r)
summary(lmnjwat)
plot(lmnjwat)
qqnorm(residuals(lmnjwat))
abline(0,1)

# check for model assumptions
check_singularity(lmnjwat)
check_heteroscedasticity(lmnjwat)  # no heterosc. detected, p = 0.431
check_model(lmnjwat)

# glmm with site as random effect
glmmnjwat <- glmer(minsize ~ nyears + (1|site),
                   data=df_years_njwat_1r, family = "poisson")
summary(glmmnjwat)
plot(glmmnjwat)

# check for model assumptions
check_singularity(glmmnjwat)
check_heteroscedasticity(glmmnjwat)  # no heterosc. detected, p = 0.208
check_overdispersion(glmmnjwat) # no overdisp detected, disp ratio = 0.516, pearson's chi sq = 56.226, p = 1
check_model(glmmnjwat)

# glm with no random effect
glmnjwat <- glm(minsize ~ nyears,
                   data=df_years_njwat_1r, family = "poisson")
summary(glmnjwat)

# check for model assumptions
check_singularity(glmnjwat)
check_heteroscedasticity(glmnjwat)  # no heterosc. detected, p = 0.179
check_overdispersion(glmnjwat) # overdisp detected for glm, p < 0.001
check_model(glmnjwat)

# check model with random slope and intercept
glmmnjwatfull <- glmer(minsize ~ nyears + (1 + nyears | site),
                   data=df_years_njwat_1r, family = "poisson")
summary(glmmnjwatfull)

check_singularity(glmmnjwatfull) # fit not singular
check_overdispersion(glmmnjwatfull) # no overdispersion
check_heteroscedasticity(glmmnjwatfull) # not heterosced.

# compare models for nj wat data
comparenjwat = compare_performance(lmnjwat,lmmnjwat,glmnjwat,glmmnjwat, glmmnjwatfull, rank=T)  
# highest score is for glmm w/ random effects
write.csv(comparenjwat, paste(respath,"modelcompare_years_njwat_1round.csv",sep=''), row.names = F)

# run anova and save output to file
# glmm with random intercepts
fnamecoeff = paste(respath,"glmm_coeff_","njwat","years_1round",".csv",sep="")
fnameanova = paste(respath,"glmm_anova_","njwat","years_1round",".csv",sep="")
write.csv(summary(glmmnjwat)$coefficients,fnamecoeff)
write.csv(Anova(glmmnjwat),fnameanova)

# glmm with random slope and intercepts
fnamecoeff = paste(respath,"glmmfull_coeff_","njwat","years_1round",".csv",sep="")
fnameanova = paste(respath,"glmmfull_anova_","njwat","years_1round",".csv",sep="")
write.csv(summary(glmmnjwatfull)$coefficients,fnamecoeff)
write.csv(Anova(glmmnjwatfull),fnameanova)

effnjwat = as.data.frame(effect(term="nyears", mod=glmmnjwatfull))
write.csv(effnjwat, paste(respath,"njwat_years_effects.csv",sep=''))

# lmm with random intercepts
fnamecoeff = paste(respath,"lmm_coeff_","njwat","years_1round",".csv",sep="")
fnameanova = paste(respath,"lmm_anova_","njwat","years_1round",".csv",sep="")
write.csv(summary(lmmnjwat)$coefficients,fnamecoeff)
write.csv(Anova(lmmnjwat),fnameanova)

## fit a lmm for CA watermelon  with number of years as x var and number of species in minimum set as y var ##
hist(df_years_cawat$minsize, breaks=7) # left skewed
hist(log(df_years_cawat$minsize))

# lmm with site as random effect
lmmcawat <- lmer(minsize ~ nyears + (1|site),
                data=df_years_cawat)
summary(lmmcawat)
plot(lmmcawat)
qqnorm(residuals(lmmcawat))
abline(0,1)
hist(residuals(lmmcawat))

# check for model assumptions
check_singularity(lmmcawat)
check_heteroscedasticity(lmmcawat)  # no heterosc. detected, p = 0.798
check_model(lmmcawat)  # residuals look normal

# lm with no random effect
lmcawat <- lm(minsize ~ nyears,
                data=df_years_cawat)
summary(lmcawat)
plot(lmcawat)
qqnorm(residuals(lmcawat))
abline(0,1)

# check for model assumptions
check_singularity(lmcawat)
check_heteroscedasticity(lmcawat)  # no heterosc. detected, p = 0.232
check_model(lmcawat)  # residuals look normal

# glmm with site as random slope and intercept, poisson family with log link
glmmcawat_full <- glmer(minsize ~ nyears + (1 + nyears | site),
                   data=df_years_cawat, family = "poisson")

# glmm with site as random effect, poisson family with log link
glmmcawat <- glmer(minsize ~ nyears + (1|site),
                  data=df_years_cawat, family = "poisson")
summary(glmmcawat)

# check for model assumptions
check_singularity(glmmcawat) # singular fit
check_heteroscedasticity(glmmcawat) # no heterosc, detected, p = 0.642
check_overdispersion(glmmcawat) # no overdispersion detected, p=1
check_model(glmmcawat)

# glm with no random effect
glmcawat <- glm(minsize ~ nyears,
                   data=df_years_cawat, family = "poisson")
summary(glmcawat)

# check for model assumptions
check_singularity(glmcawat)
check_heteroscedasticity(glmcawat) # no heterosc, detected, p = 0.65
check_overdispersion(glmcawat) # no overdispersion detected, p = 1
check_model(glmcawat)

# compare models
comparecawat = compare_performance(lmcawat, lmmcawat, glmcawat, glmmcawat, rank=T)  # best fit for lmm w/ random effect
plot(comparecawat)
write.csv(comparecawat,paste(respath,"modelcompare_years_cawat.csv",sep=''), row.names = F)

# perform anova and save glmm results to file
# save glm output
fnamecoeff = paste(respath,"glm_coeff_","cawat","years",".csv",sep="")
fnameanova = paste(respath,"glm_anova_","cawat","years",".csv",sep="")
write.csv(summary(glmcawat)$coefficients,fnamecoeff)
write.csv(Anova(glmcawat),fnameanova)

effcawat = as.data.frame(effect(term="nyears", mod=glmcawat))
write.csv(effcawat, paste(respath,"cawat_years_effects.csv",sep=''))

# save glmm output
fnamecoeff = paste(respath,"glmm_coeff_","cawat","years",".csv",sep="")
fnameanova = paste(respath,"glmm_anova_","cawat","years",".csv",sep="")
write.csv(summary(glmmcawat)$coefficients,fnamecoeff)
write.csv(Anova(glmmcawat),fnameanova)

# save lmm output
fnamecoeff = paste(respath,"lmm_coeff_","cawat","years",".csv",sep="")
fnameanova = paste(respath,"lmm_anova_","cawat","years",".csv",sep="")
write.csv(summary(lmmcawat)$coefficients,fnamecoeff)
write.csv(Anova(lmmcawat),fnameanova)

### stats for change in min species with number of rounds ###

## fit a linear mixed model for blueberry with number of rounds as x var and number of species in minimum set as y var
# view histogram of blueberry data
hist(df_rounds_blue$minsize,breaks=6)
hist(log(df_rounds_blue$minsize),breaks=6)

# glmm with site as random effect
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

# glm with no random effect
glmblueround <- glm(minsize ~ nrounds,
                       data=df_rounds_blue, family='poisson')
summary(glmblueround)

# visualize model fit
qqnorm(residuals(glmblueround))
hist(residuals(glmblueround))
check_model(glmblueround)
check_singularity(glmblueround)
check_heteroscedasticity(glmblueround) # no heterosced., p = 0.739
check_overdispersion(glmblueround) # no overdispersion, p = 1

# glmm with site as random slope and intercept
glmmblueround_full <- glmer(minsize ~ nrounds + (1 + nrounds |site),
                       data=df_rounds_blue, family='poisson')
check_singularity(glmmblueround_full)
isSingular(glmmblueround_full)
check_heteroscedasticity(glmmblueround_full)
check_overdispersion(glmmblueround_full)

# lm w/ no random effects
lmblueround <- lm(minsize ~ nrounds, data=df_rounds_blue)
check_model(lmblueround)
check_heteroscedasticity(lmblueround) # heterosced detected, p=0.015
check_homogeneity(lmblueround)

# lmm w/ random intercept
lmmblueround <- lmer(minsize ~ nrounds + (1|site),
                      data=df_rounds_blue)
summary(lmmblueround)
#check model
check_model(lmmblueround)
check_heteroscedasticity(lmmblueround) 
check_homogeneity(lmmblueround)

# full lmm w/ random slope and intercept
lmmblueround_full <- lmer(minsize ~ nrounds + (1 + nrounds |site),
                     data=df_rounds_blue)
summary(lmmblueround_full)

# check model assumptions
isSingular(lmmblueround_full)  # singular fit
check_model(lmmblueround_full)
check_heteroscedasticity(lmmblueround_full) 
check_homogeneity(lmmblueround_full)

# compare models
compareblueround = compare_performance(lmblueround, lmmblueround, lmmblueround_full, 
                                       glmblueround, glmmblueround, glmmblueround_full, rank=T) # model with random effect has better score
plot(compareblueround)
write.csv(compareblueround,paste(respath,"modelcompare_rounds_blue.csv",sep=''), row.names = F)

# do chi sqr and save output
Anova(glmmblueround)
fnamecoeff = paste(respath,"glmm_coeff_","blue","rounds",".csv",sep="")
fnameanova = paste(respath,"glmm_anova_","blue","rounds",".csv",sep="")
write.csv(summary(glmmblueround)$coefficients,fnamecoeff)
write.csv(Anova(glmmblueround),fnameanova)

effblueround = as.data.frame(effect(term="nrounds", mod=glmmblueround))
write.csv(effblueround, paste(respath,"blue_rounds_effects.csv",sep=''))
percincblueround = (effblueround$fit[5]-effblueround$fit[1])/effblueround$fit[1]*100

Anova(lmmblueround)
fnamecoeff = paste(respath,"lmm_coeff_","blue","rounds",".csv",sep="")
fnameanova = paste(respath,"lmm_anova_","blue","rounds",".csv",sep="")
write.csv(summary(lmmblueround)$coefficients,fnamecoeff)
write.csv(Anova(lmmblueround),fnameanova)

Anova(lmmblueround_full)
fnamecoeff = paste(respath,"lmmfull_coeff_","blue","rounds",".csv",sep="")
fnameanova = paste(respath,"lmmfull_anova_","blue","rounds",".csv",sep="")
write.csv(summary(lmmblueround_full)$coefficients,fnamecoeff)
write.csv(Anova(lmmblueround_full),fnameanova)


## fit a glmm for NJ watermelon with number of rounds as x var and number of species in minimum set as y var
# view watermelon data
hist(df_rounds_njwat$minsize)
hist(log(df_rounds_njwat$minsize))

# fit glmm with random effect of site
glmmnjwatround <- glmer(minsize ~ nrounds + (1|site),
                        data = df_rounds_njwat, family = 'poisson')
summary(glmmnjwatround)

# view model fit
hist(residuals(glmmnjwatround))
check_model(glmmnjwatround)
check_singularity(glmmnjwatround)
check_heteroscedasticity(glmmnjwatround) # no heterosced., p = 0.145
check_overdispersion(glmmnjwatround) # overdispersion detected, p < 0.001

# fit glm with no random effect
glmnjwatround <-  glm(minsize ~ nrounds,
                        data = df_rounds_njwat, family = 'poisson')
summary(glmnjwatround)

check_model(glmnjwatround)
check_singularity(glmnjwatround)
check_heteroscedasticity(glmnjwatround)
check_overdispersion(glmnjwatround) # overdisp detected

# fit glmm with negative binomial to account for overdispersion
glmmnjwatround2 <- glmer.nb(minsize ~ nrounds + (1|site),
                         data = df_rounds_njwat)
summary(glmmnjwatround2)

# view model fit
hist(residuals(glmmnjwatround2))
check_model(glmmnjwatround2)
check_singularity(glmmnjwatround2)
check_heteroscedasticity(glmmnjwatround2) # no heterosced., p = 0.410

# fit glmm full model with negative binomial to account for overdispersion
glmmnjwatround2full <- glmer.nb(minsize ~ nrounds + (1 + nrounds |site),
                            data = df_rounds_njwat)
summary(glmmnjwatround2full)
check_model(glmmnjwatround2full)
isSingular(glmmnjwatround2full) # singular fit

# fit lm w/o random effect
lmnjwatround <-  lm(minsize ~ nrounds,
                      data = df_rounds_njwat)
summary(lmnjwatround)

check_model(lmnjwatround)
check_heteroscedasticity(lmnjwatround) # hetersced detected
check_homogeneity(lmnjwatround) # non homogeneous variance detected

# fit lmm w/ random int
lmmnjwatround <-  lmer(minsize ~ nrounds + (1|site),
                    data = df_rounds_njwat)
check_model(lmmnjwatround)

lmmnjwatroundfull <- lmer(minsize ~ nrounds + (1 + nrounds | site),
                          data = df_rounds_njwat)  # singular fit

# compare models
comparenjwatround = compare_performance(glmnjwatround, glmmnjwatround, glmmnjwatround2, glmmnjwatround2full,
                    lmnjwatround, lmmnjwatround, lmmnjwatroundfull, rank=T) # model with random effect and neg binom does best

write.csv(comparenjwatround,paste(respath,"modelcompare_rounds_njwat.csv",sep=''), row.names = F)

# save results from negative binomial due to overdispersion
Anova(glmmnjwatround2)
fnamecoeff = paste(respath,"glmm_coeff_negbin_","njwat","rounds",".csv",sep="")
fnameanova = paste(respath,"glmm_anova_negbin_","njwat","rounds",".csv",sep="")
write.csv(summary(glmmnjwatround2)$coefficients,fnamecoeff)
write.csv(Anova(glmmnjwatround2),fnameanova)

effnjwatrounds = as.data.frame(effect(term = "nrounds", mod = glmmnjwatround2))
write.csv(effnjwatrounds, paste(respath, "njwat_effects_rounds.csv",sep=''))
percincnjwatround = (effnjwatrounds$fit[5]-effnjwatrounds$fit[1])/effnjwatrounds$fit[1]*100


## fit a lmm for CA watermelon with number of rounds as x var and number of species in minimum set as y var
# view data
hist(df_rounds_cawat$minsize,breaks=8)
hist(log(df_rounds_cawat$minsize),breaks=6)

# site as random effect (full model)
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

# glm w/ no random effect
glmcawatround <- glm(minsize ~ nrounds,
                        data=df_rounds_cawat, family='poisson')
check_heteroscedasticity(glmcawatround) # heterosced detected
check_overdispersion(glmcawatround) # no overdispersion
check_model(glmcawatround)

# lmm w/ random slope and intercept
lmmcawatroundfull <- lmer(minsize ~ nrounds + (1 + nrounds | site),
                            data=df_rounds_cawat) # singular fit
check_model(lmmcawatroundfull)
check_heteroscedasticity(lmmcawatroundfull)  # heterosced detected
check_homogeneity(lmmcawatroundfull) # non-homogeneity detected

# lmm w/ random intercept
lmmcawatround <- lmer(minsize ~ nrounds + (1 | site),
                          data=df_rounds_cawat)
check_model(lmmcawatround)

# lm w/ no random effect
lmcawatround <- lm(minsize ~ nrounds,
                      data=df_rounds_cawat)
check_model(lmcawatround)

# compare model fit
comparecawatround = compare_performance(lmcawatround, lmmcawatround, lmmcawatroundfull,
                                  glmcawatround, glmmcawatround, glmmcawatroundfull, rank=T)
plot(comparecawatround)
write.csv(comparecawatround,paste(respath,"modelcompare_rounds_cawat.csv",sep=''), row.names = F)

# anova and save output
Anova(glmmcawatround)
fnamecoeff = paste(respath,"glmm_coeff_","cawat","rounds",".csv",sep="")
fnameanova = paste(respath,"glmm_anova_","cawat","rounds",".csv",sep="")
write.csv(summary(glmmcawatround)$coefficients,fnamecoeff)
write.csv(Anova(glmmcawatround),fnameanova)

# glm w/ no random effects has lowest AIC
Anova(glmcawatround)
fnamecoeff = paste(respath,"glm_coeff_","cawat","rounds",".csv",sep="")
fnameanova = paste(respath,"glm_anova_","cawat","rounds",".csv",sep="")
write.csv(summary(glmcawatround)$coefficients,fnamecoeff)
write.csv(Anova(glmcawatround),fnameanova)

effcawatround = as.data.frame(effect(term="nrounds", mod=glmmcawatround))
write.csv(effcawatround, paste(respath,"cawat_rounds_effects.csv",sep=''))


### report number of species needed to meet threshold for single site-date --------
df_rounds = rbind(df_rounds_blue, df_rounds_njwat, df_rounds_cawat)
minset1r <- df_rounds %>% 
  filter(nrounds==1) %>%
  group_by(crop) %>% 
  summarise(spec_mean = mean(minsize), spec_sd = sd(minsize))
write.csv(minset1r, paste(respath,'minset_1sitedate.csv',sep=''))

