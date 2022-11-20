library(data.table)
library(ggplot2)
library(ez)
library(pwr)
library(effsize)
library(lme4)
library(lmerTest)
library(nlme)
library(viridis)
library(sjstats)
library(sjPlot)
library(clipr)

# here is data.table, dplyr will be appeared later (or not)

# locale

Sys.setlocale(category = 'LC_ALL', locale = "RU_ru")
Sys.getlocale()

# color

viridis(2)
scale_color_viridis()

#----cathodal dt--------------------------------------------------------

cathodal <-
  fread("tDCS_cathodal.csv",
        sep = ";",
        header = TRUE)

# View(cathodal)

#dt with rejected only
cathodal_rejected <- cathodal[Choice == "Rejected", ]

#dt with selected only
cathodal_selected <- cathodal[Choice == "Selected", ]


#dt with single difference insteas of repeated measured Evaluation

cathodal_dif <-
  cathodal[Rating == "R1", ]  #just for creating dataframe with certain number of rows


cathodal_dif[, 'Evaluation'] <-
  cathodal[Rating == "R2", Evaluation] -
  cathodal[Rating == "R1", Evaluation]
colnames(cathodal_dif)[colnames(cathodal_dif) == "Evaluation"] <-
  "Difference"
cathodal_dif[, 'Rating'] <- NULL
View(cathodal_dif)


#----anodal dt--------------------------------------------------------

anodal <-
  fread("tDCS_anodal2.csv",
        sep = ";",
        header = TRUE)

# View(anodal)

#dt with rejected only
anodal_rejected <- anodal[Choice == "Rejected", ]

#dt with selected only
anodal_selected <- anodal[Choice == "Selected", ]


#dt with single difference insteas of repeated measured Evaluation

anodal_dif <-
  anodal[Rating == "R1", ]  #just for creating dataframe with certain number of rows


anodal_dif[, 'Evaluation'] <-  anodal[Rating == "R2", Evaluation] -
  anodal[Rating == "R1", Evaluation]
colnames(anodal_dif)[colnames(anodal_dif) == "Evaluation"] <-
  "Difference"
anodal_dif[, 'Rating'] <- NULL
View(anodal_dif)

#Post-Ex is processed separately, remove it
# anodal_dif_short <- subset(anodal_dif, Type != "Post-Ex")

length_cat <- length(unique(cathodal_dif$Sub, incomparables = FALSE))
length_cat
length_an <- length(unique(anodal_dif$Sub, incomparables = FALSE))
length_an
#----   function for extracting t values and DFs (not mine)-----------------

pairwise.t.test.with.t.and.df <- function (x, g, p.adjust.method = p.adjust.methods, pool.sd = !paired, 
                                           paired = FALSE, alternative = c("two.sided", "less", "greater"), 
                                           ...) 
{
  if (paired & pool.sd) 
    stop("pooling of SD is incompatible with paired tests")
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
  g <- factor(g)
  p.adjust.method <- match.arg(p.adjust.method)
  alternative <- match.arg(alternative)
  if (pool.sd) {
    METHOD <- "t tests with pooled SD"
    xbar <- tapply(x, g, mean, na.rm = TRUE)
    s <- tapply(x, g, sd, na.rm = TRUE)
    n <- tapply(!is.na(x), g, sum)
    degf <- n - 1
    total.degf <- sum(degf)
    pooled.sd <- sqrt(sum(s^2 * degf)/total.degf)
    compare.levels <- function(i, j) {
      dif <- xbar[i] - xbar[j]
      se.dif <- pooled.sd * sqrt(1/n[i] + 1/n[j])
      t.val <- dif/se.dif
      if (alternative == "two.sided") 
        2 * pt(-abs(t.val), total.degf)
      else pt(t.val, total.degf, lower.tail = (alternative == 
                                                 "less"))
    }
    compare.levels.t <- function(i, j) {
      dif <- xbar[i] - xbar[j]
      se.dif <- pooled.sd * sqrt(1/n[i] + 1/n[j])
      t.val = dif/se.dif 
      t.val
    }       
  }
  else {
    METHOD <- if (paired) 
      "paired t tests"
    else "t tests with non-pooled SD"
    compare.levels <- function(i, j) {
      xi <- x[as.integer(g) == i]
      xj <- x[as.integer(g) == j]
      t.test(xi, xj, paired = paired, alternative = alternative, 
             ...)$p.value
    }
    compare.levels.t <- function(i, j) {
      xi <- x[as.integer(g) == i]
      xj <- x[as.integer(g) == j]
      t.test(xi, xj, paired = paired, alternative = alternative, 
             ...)$statistic
    }
    compare.levels.df <- function(i, j) {
      xi <- x[as.integer(g) == i]
      xj <- x[as.integer(g) == j]
      t.test(xi, xj, paired = paired, alternative = alternative, 
             ...)$parameter
    }
  }
  PVAL <- pairwise.table(compare.levels, levels(g), p.adjust.method)
  TVAL <- pairwise.table.t(compare.levels.t, levels(g), p.adjust.method)
  if (pool.sd) 
    DF <- total.degf
  else
    DF <- pairwise.table.t(compare.levels.df, levels(g), p.adjust.method)           
  ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL, 
              p.adjust.method = p.adjust.method, t.value = TVAL, dfs = DF)
  class(ans) <- "pairwise.htest"
  ans
}


pairwise.table.t <- function (compare.levels.t, level.names, p.adjust.method) 
{
  ix <- setNames(seq_along(level.names), level.names)
  pp <- outer(ix[-1L], ix[-length(ix)], function(ivec, jvec) sapply(seq_along(ivec), 
                                                                    function(k) {
                                                                      i <- ivec[k]
                                                                      j <- jvec[k]
                                                                      if (i > j)
                                                                        compare.levels.t(i, j)               
                                                                      else NA
                                                                    }))
  pp[lower.tri(pp, TRUE)] <- pp[lower.tri(pp, TRUE)]
  pp
}








#----   1.1.1 ✓ paired t test for 6 hypothesis: main and 4 controls  -----------------

# 1) difficult rejected cathodal отличаются от difficult rejected sham -- main
# 2)difficult rejected less  than difficult selected -- нельзя объединять группы с tDCS и sham в один сэмпл
# 3)difficult rejected  отличаются от easy rejected 
# 4)difficult rejected  отличаются от computer rejected 
# 5)difficult rejected  отличаются от post-ex rejected ????
# 6) there is no difference in each condition: difficult rejected, easy rejected... 


# cathodal ----

alpha_corrected <- 0.05/4
# 0.01/4
# 0.001/4


# 1) for full, main:
# difficult rejected cathodal differ from difficult rejected sham (less) 

cathodal_dif_tDCS_rej <- cathodal_dif[Type == "Difficult" & Stimulation == "tDCS_cat" & Choice == 'Rejected', Difference]
cathodal_dif_sham_rej <- cathodal_dif[Type == "Difficult" & Stimulation == "sham" & Choice == 'Rejected', Difference]



t.test(x = cathodal_dif_sham_rej,  
       y = cathodal_dif_tDCS_rej,  
       paired = T, var.equal = FALSE,
       conf.level = 1-alpha_corrected, alternative = "less")

cohen.d( cathodal_dif_tDCS_rej, cathodal_dif_sham_rej, paired = T, hedges.correction = F)
cohen.d( cathodal_dif_tDCS_rej, cathodal_dif_sham_rej, paired = T, hedges.correction = T)


#2) for all where difficult: 
# rejected less than selected in only tDCS trials

cathodal_difficult_rej <- cathodal_dif[Type == "Difficult" & Stimulation == "tDCS_cat" & Choice == 'Rejected', Difference]
cathodal_difficult_sel <- cathodal_dif[Type == "Difficult" & Stimulation == "tDCS_cat" & Choice == 'Selected', Difference]

ttest_difficult <- t.test(x = cathodal_difficult_rej,  
                          y = cathodal_difficult_sel,  
                          paired = T, var.equal = FALSE,
                          conf.level = 1-alpha_corrected, alternative = "less")


ttest_difficult

length(cathodal_difficult_rej)
length(cathodal_difficult_sel)

View(cathodal_dif)


cohen.d( cathodal_difficult_sel, cathodal_difficult_rej, paired = T, hedges.correction = F)

#3) for all where rejected only tDCS: 
# difficult отличаются от easy 
# difficult отличаются от computer 
# difficult отличаются от post-ex ?


cathodal_easy_rej <- cathodal_dif[Type == "Easy" & Stimulation == "tDCS_cat" & Choice == 'Rejected', Difference]
cathodal_comp_rej <- cathodal_dif[Type == "Computer" & Stimulation == "tDCS_cat" & Choice == 'Rejected', Difference]
cathodal_postex_rej <- cathodal_dif[Type == "Post-Ex"  & Stimulation == "tDCS_cat" & Choice == 'Rejected', Difference]

t.test(x = cathodal_difficult_rej,  
       y = cathodal_easy_rej,  
       paired = T, var.equal = FALSE,
       conf.level = 1-alpha_corrected, alternative = "less")


cohen.d( cathodal_easy_rej, cathodal_difficult_rej, paired = T, hedges.correction = F)

t.test(x = cathodal_difficult_rej,  
       y = cathodal_comp_rej,  
       paired = T, var.equal = FALSE,
       conf.level = 1-alpha_corrected, alternative = "less")

cohen.d( cathodal_comp_rej, cathodal_difficult_rej, paired = T, hedges.correction = F)


t.test(x = cathodal_difficult_rej,  
       y = cathodal_postex_rej,  
       paired = T, var.equal = FALSE,
       conf.level = 1-alpha_corrected, alternative = "less")

cohen.d( cathodal_postex_rej, cathodal_difficult_rej, paired = T, hedges.correction = F)


#pairwise t tests separated and only for searching for differece with cathodal stimulation ?
#does not make sence: we are not interested in precesely cathodal stimulation in all particular
#cases here. So, I don't include this part (see appendix).


cathodal_dif_cropped_for_ttest_interaction <- cathodal_dif


cathodal_dif_cropped_for_ttest_interaction[, 'Stimulation_Type_Choice'] <-
  paste0(cathodal_dif_cropped_for_ttest_interaction[, Stimulation], '_',
         cathodal_dif_cropped_for_ttest_interaction[, Type], '_',
         cathodal_dif_cropped_for_ttest_interaction[, Choice])

View(cathodal_dif_cropped_for_ttest_interaction)


pairwTTest_fdr_dif <-
  pairwise.t.test.with.t.and.df(
    x = cathodal_dif_cropped_for_ttest_interaction$Difference,
    g = cathodal_dif_cropped_for_ttest_interaction$Stimulation_Type_Choice,
    p.adjust.method = 'fdr',
    paired = T
  )

View(pairwTTest_fdr_dif$p.value)

round(pairwTTest_fdr_dif$p.value,2)

# write.table(round(pairwTTest_fdr_dif$p.value,2),"clipboard",sep="\t",col.names=NA)

clipr::write_clip(round(pairwTTest_fdr_dif$p.value,2))

#  anodal -----

# 1) for full, main:
# difficult rejected cathodal differ from difficult rejected sham (less)

anodal_dif_tDCS_rej <- anodal_dif[Type == "Difficult" & Stimulation == "tDCS" 
                                  & Choice == 'Rejected', Difference]
anodal_dif_sham_rej <- anodal_dif[Type == "Difficult" & Stimulation == "sham" 
                                  & Choice == 'Rejected', Difference]

t.test(x = anodal_dif_tDCS_rej,  
       y = anodal_dif_sham_rej,  
       paired = T, var.equal = FALSE,
       conf.level = 1-alpha_corrected, alternative = "less")

alpha_corrected


#2) for all where difficult: 
# rejected less than selected 

View(anodal_dif)

anodal_difficult_rej <- anodal_dif[Type == "Difficult" & Stimulation == "tDCS" & Choice == 'Rejected', Difference]
anodal_difficult_sel <- anodal_dif[Type == "Difficult" & Stimulation == "tDCS" & Choice == 'Selected', Difference]

t.test(x = anodal_difficult_rej,  
       y = anodal_difficult_sel,  
       paired = T, var.equal = FALSE,
       conf.level = 1-alpha_corrected, alternative = "less")

cohen.d( anodal_difficult_rej, anodal_difficult_sel, paired = T, hedges.correction = F)

#3) for all where rejected: 
# difficult отличаются от easy 
# difficult отличаются от computer 
# difficult отличаются от post-ex ?

anodal_difficult_rej <- anodal_dif[Type == "Difficult" & Stimulation == "tDCS" & Choice == 'Rejected', Difference]
anodal_easy_rej <- anodal_dif[Type == "Easy" & Stimulation == "tDCS" & Choice == 'Rejected', Difference]
anodal_comp_rej <- anodal_dif[Type == "Computer" & Stimulation == "tDCS" & Choice == 'Rejected', Difference]
anodal_postex_rej <- anodal_dif[Type == "Post-Ex" & Stimulation == "tDCS" & Choice == 'Rejected', Difference]



t.test(x = anodal_difficult_rej,  
       y = anodal_easy_rej,  
       paired = T, var.equal = FALSE,
       conf.level = 1-alpha_corrected, alternative = "less")


cohen.d( anodal_difficult_rej, anodal_easy_rej, paired = T, hedges.correction = F)

t.test(x = anodal_difficult_rej,  
       y = anodal_comp_rej,  
       paired = T, var.equal = FALSE,
       conf.level = 1-alpha_corrected, alternative = "less")

cohen.d( anodal_difficult_rej, anodal_comp_rej, paired = T, hedges.correction = F)


t.test(x = anodal_difficult_rej,  
       y = anodal_postex_rej,  
       paired = T, var.equal = FALSE,
       conf.level = 1-alpha_corrected, alternative = "less")

cohen.d( anodal_difficult_rej, anodal_postex_rej, paired = T, hedges.correction = T)


anodal_dif_cropped_for_ttest_interaction <- anodal_dif


anodal_dif_cropped_for_ttest_interaction[, 'Stimulation_Type_Choice'] <-
  paste0(anodal_dif_cropped_for_ttest_interaction[, Stimulation], '_',
         anodal_dif_cropped_for_ttest_interaction[, Type], '_',
         anodal_dif_cropped_for_ttest_interaction[, Choice])

View(anodal_dif_cropped_for_ttest_interaction)


an_pairwTTest_fdr_dif <-
  pairwise.t.test.with.t.and.df(
    x = anodal_dif_cropped_for_ttest_interaction$Difference,
    g = anodal_dif_cropped_for_ttest_interaction$Stimulation_Type_Choice,
    p.adjust.method = 'fdr',
    paired = T
  )

View(an_pairwTTest_fdr_dif$p.value)

round(an_pairwTTest_fdr_dif$p.value,2)

clipr::write_clip(round(an_pairwTTest_fdr_dif$p.value,2))

#----   1.1.4 ✓  LME  cathodal -------------------------------------------------------

#Cathodal full data

cathodal_wide <- dcast(cathodal, Sub + Stimulation + Type + Choice ~ Rating)
cathodal_wide <- as.data.table(cathodal_wide)

#unsuccessful trials for covariate design------

lmerka1 <- lmer(R2 ~ R1 + Stimulation*Type*Choice + (1 + Stimulation|Sub), 
                data = cathodal_wide, REML = F)
anova(lmerka1)
summary(lmerka1)

lmerka5 <- lmer(R2 ~ Stimulation*Type*Choice + (R1 + Stimulation|Sub), 
                data = cathodal_wide, REML = F)
anova(lmerka5)
summary(lmerka5)

lmerka4 <- lmer(R2 ~ Stimulation*Type*Choice + (R1 + R1|Sub), 
                data = cathodal_wide, REML = F)
anova(lmerka4)
summary(lmerka4)

lmerka2 <- lmer(R2 ~ Stimulation*Type*Choice + (R1 + R1|Sub) + (Stimulation|Sub), 
                data = cathodal_wide, REML = F)
anova(lmerka2)
summary(lmerka2)

# lmerka3 <- lmer(R2 ~ R1 + Stimulation*Type*Choice + (1 + 1|Sub), 
#                 data = cathodal_wide, REML = F)


# lmerka_gls <- gls(R2 ~ R1 + Stimulation*Type*Choice, data = cathodal_wide, 
#                   random =  ~ (1 + Stimulation|Sub), 
#                   method = 'ML')

lmerka_lme <- lme(R2 ~ R1 + Stimulation*Type*Choice, data = cathodal_wide,
                  random =  ~ 1 + Stimulation|Sub,
                  method = 'ML')

anova(lmerka_lme )
summary(lmerka_lme )

pairwise.t.test()

View(cathodal_dif)


summary(lmerka1)
summary(lmerka2)
summary(lmerka_dif)
anova(lmerka2)
anova(lmerka_dif)
View(cathodal_wide)

ggplot(fortify(lmerka1))

#successful lme model on the full data ------

#if we are comparing models with same fixed factors and different random factors,
#we should use REML (restricted max likelihood)

#https://stats.stackexchange.com/questions/13166/rs-lmer-cheat-sheet

#model with only random intercept for sub only (has huge Type 1 error rate)
##why this does not work??? -- because of brackets???
#https://stats.stackexchange.com/questions/58745/using-lmer-for-repeated-measures-linear-mixed-effect-model

lmerka_random_int_sub <- lmer(Difference ~ Stimulation*Type*Choice + (1|Sub), 
                              data = cathodal_dif, REML = T)
lmerka_random_int_sub 
anova(lmerka_random_int_sub)
summary(lmerka_random_int_sub)
isSingular(lmerka_random_int_sub)

#model with random intercept and slope for sub only
lmerka_random_int_slope_sub <- lmer(Difference ~ Stimulation*Type*Choice + (1 + 1|Sub), 
                                    data = cathodal_dif, REML = T)
lmerka_random_int_slope_sub 
anova(lmerka_random_int_slope_sub)
summary(lmerka_random_int_slope_sub)
isSingular(lmerka_random_int_slope_sub)


#chosed model

#model with correlated random intercept and ransom slope by stimulation 
# Conditional Akaike information criterion:   386.79
lmerka_random_int_slope_stim_sub <- lmer(Difference ~ 
                                           Stimulation*Type*Choice + (1 + Stimulation|Sub), 
                                         data = cathodal_dif, REML = T)
lmerka_random_int_slope_stim_sub 
anova_lmerka_random_int_slope_stim_sub <- anova(lmerka_random_int_slope_stim_sub)
#why do we need anova(lm)
#https://stats.stackexchange.com/questions/115304/interpreting-output-from-anova-when-using-lm-as-input
anova_lmerka_random_int_slope_stim_sub
summary(lmerka_random_int_slope_stim_sub)
isSingular(lmerka_random_int_slope_stim_sub)

tab_model(lmerka_random_int_slope_stim_sub)
# sjt.lmer(lmerka_random_int_slope_stim_sub)
clipr::write_clip(tab_model(lmerka_random_int_slope_stim_sub))



eta_sq(anova_lmerka_random_int_slope_stim_sub, partial = T,ci.lvl = NULL)



#model with only random slope for stimulation / sub 
# Conditional Akaike information criterion:   386.79
lmerka_random_slope_stim_sub <- lmer(Difference ~ Stimulation*Type*Choice + (Stimulation|Sub), 
                                     data = cathodal_dif, REML = T)
lmerka_random_slope_stim_sub 
anova(lmerka_random_slope_stim_sub)
summary(lmerka_random_slope_stim_sub)
isSingular(lmerka_random_slope_stim_sub)




#model with  ransom slope by stimulationn, type and choice 
lmerka_random_int_slope_all_slopes_sub <- lmer(Difference ~ 
                                                 Stimulation*Type*Choice + (Stimulation|Sub) +
                                                 (Type|Sub) + (Choice|Sub), 
                                               data = cathodal_dif, REML = T)

lmerka_random_int_slope_all_slopes_sub 
anova(lmerka_random_int_slope_all_slopes_sub)
summary(lmerka_random_int_slope_all_slopes_sub)
isSingular(lmerka_random_int_slope_all_slopes_sub)

#model with correlated random intercept and ransom slope by stimulationn, type and choice 
lmerka_random_int_slope_all_sub <- 
  lmer(Difference ~ Stimulation*Type*Choice + (1 + Stimulation|Sub) +
         (1 + Type|Sub) + (1 + Choice|Sub), data = cathodal_dif, REML = T)

lmerka_random_int_slope_all_sub 
anova(lmerka_random_int_slope_all_sub)
summary(lmerka_random_int_slope_all_sub)
isSingular(lmerka_random_int_slope_all_sub)
allFit(lmerka_random_int_slope_all_sub)

# lmerka_random_slope_all_in_one <- lmer(Difference ~ Stimulation*Type*Choice + 
#                                          (1 + Stimulation*Type*Choice|Sub), 
#                                      data = cathodal_dif, REML = T)

lmerka_random_slope_all_in_one <- lmer(Difference ~ Stimulation*Type*Choice + 
                                         (1 + Stimulation + Type + Choice|Sub), 
                                       data = cathodal_dif, REML = T)
isSingular(lmerka_random_slope_all_in_one)
anova(lmerka_random_slope_all_in_one)
summary(lmerka_random_slope_all_in_one)
isSingular(lmerka_random_slope_all_in_one)

allFit(lmerka_random_slope_all_in_one)
optimizeLmer(Difference ~ Stimulation*Type*Choice + 
               (1 + Stimulation + Type + Choice|Sub), 
             data = cathodal_dif, REML = T, optimizer = 'bobyqa')



View(cathodal_dif)

model.matrix(lmerka_random_int_slope_stim_sub)

#--Cathodal model selection----------

#by AIC and LogLik?

anova(lmerka_random_int_sub, lmerka_random_int_slope_sub, lmerka_random_slope_stim_sub, 
      lmerka_random_int_slope_stim_sub,lmerka_random_int_slope_all_slopes_sub,
      lmerka_random_int_slope_all_sub, lmerka_random_slope_all_in_one)
#this comparisson are going through REML=F (ML instead of REML), which is not good for 
#comparisons of model with random factros only. For random factors model comparison better to use 
#REML

# 
# contourplot(Difference ~ Stimulation*Type*Choice | factor(Sub), data=cathodal_dif)
# 
# xyplot(Difference ~ Stimulation*Type | Stimulation*Type, data=cathodal_dif)

# fft <- lme4::fortify(lmerka_dif)

cAIC(lmerka_random_int_sub)
cAIC(lmerka_random_int_slope_sub)
cAIC(lmerka_random_slope_stim_sub)
cAIC(lmerka_random_int_slope_stim_sub)
cAIC(lmerka_random_int_slope_all_slopes_sub)
cAIC(lmerka_random_int_slope_all_sub)
cAIC(lmerka_random_slope_all_in_one)

lmerControl(lmerka_random_int_slope_all_slopes_sub)


#--1.1.4 ✓  LME  anodal ------------------


#successful anodal lme model on the full data------------------------------------------------------

#if we are comparing models with same fixed factors and different random factors,
#we should use REML (restricted max likelihood)

#https://stats.stackexchange.com/questions/13166/rs-lmer-cheat-sheet



#model with only random intercept for sub only (has huge Type 1 error rate)
##why this does not work??? -- because of breakets???
#https://stats.stackexchange.com/questions/58745/using-lmer-for-repeated-measures-linear-mixed-effect-model

lmerka_random_int_sub_an <- lmer(Difference ~ Stimulation*Type*Choice + (1|Sub), 
                                 data = anodal_dif, REML = T)
lmerka_random_int_sub_an 
anova(lmerka_random_int_sub_an)
summary(lmerka_random_int_sub_an)
isSingular(lmerka_random_int_sub_an)

#model with random intercept and slope for sub only
lmerka_random_int_slope_sub_an <- lmer(Difference ~ Stimulation*Type*Choice + (1 + 1|Sub), 
                                       data = anodal_dif, REML = T)
lmerka_random_int_slope_sub_an 
anova(lmerka_random_int_slope_sub_an)
summary(lmerka_random_int_slope_sub_an)
isSingular(lmerka_random_int_slope_sub_an)


#anova chosed model

#model with correlated random intercept and ransom slope by stimulation 
# Conditional Akaike information criterion:   212.05
lmerka_random_int_slope_stim_sub_an <- lmer(Difference ~ 
                                              Stimulation*Type*Choice + (1 + Stimulation|Sub), 
                                            data = anodal_dif, REML = T)
lmerka_random_int_slope_stim_sub_an 
anova_lmerka_random_int_slope_stim_sub_an <- anova(lmerka_random_int_slope_stim_sub_an)
#why do we need anova(lm)
#https://stats.stackexchange.com/questions/115304/interpreting-output-from-anova-when-using-lm-as-input
anova_lmerka_random_int_slope_stim_sub_an
summary(lmerka_random_int_slope_stim_sub_an)
isSingular(lmerka_random_int_slope_stim_sub_an)


eta_sq(anova_lmerka_random_int_slope_stim_sub_an, partial = T,ci.lvl = NULL)


tab_model(lmerka_random_int_slope_stim_sub_an)
# sjt.lmer(lmerka_random_int_slope_stim_sub)
# clipr::write_clip(tab_model(lmerka_random_int_slope_stim_sub_an))

#model with only random slope for stimulation / sub 
# Conditional Akaike information criterion:   212.05
lmerka_random_slope_stim_sub_an <- lmer(Difference ~ Stimulation*Type*Choice + (Stimulation|Sub), 
                                        data = anodal_dif, REML = T)
lmerka_random_slope_stim_sub_an 
anova(lmerka_random_slope_stim_sub_an)
summary(lmerka_random_slope_stim_sub_an)
isSingular(lmerka_random_slope_stim_sub_an)




#model with  ransom slope by stimulationn, type and choice 
lmerka_random_int_slope_all_slopes_sub_an <- lmer(Difference ~ 
                                                    Stimulation*Type*Choice + (Stimulation|Sub) +
                                                    (Type|Sub) + (Choice|Sub), 
                                                  data = anodal_dif, REML = T)

lmerka_random_int_slope_all_slopes_sub_an 
anova(lmerka_random_int_slope_all_slopes_sub_an)
summary(lmerka_random_int_slope_all_slopes_sub_an)
isSingular(lmerka_random_int_slope_all_slopes_sub_an)

#model with correlated random intercept and ransom slope by stimulationn, type and choice 
lmerka_random_int_slope_all_sub_an <- 
  lmer(Difference ~ Stimulation*Type*Choice + (1 + Stimulation|Sub) +
         (1 + Type|Sub) + (1 + Choice|Sub), data = anodal_dif, REML = T)

lmerka_random_int_slope_all_sub_an 
anova(lmerka_random_int_slope_all_sub_an)
summary(lmerka_random_int_slope_all_sub_an)
isSingular(lmerka_random_int_slope_all_sub_an)
allFit(lmerka_random_int_slope_all_sub_an)



lmerka_random_slope_all_in_one_an <- lmer(Difference ~ Stimulation*Type*Choice + 
                                            (1 + Stimulation + Type + Choice|Sub), 
                                          data = anodal_dif, REML = T)
isSingular(lmerka_random_slope_all_in_one_an)
anova(lmerka_random_slope_all_in_one_an)
summary(lmerka_random_slope_all_in_one_an)
isSingular(lmerka_random_slope_all_in_one_an)



cAIC(lmerka_random_int_sub_an)
cAIC(lmerka_random_int_slope_sub_an)
cAIC(lmerka_random_slope_stim_sub_an)
cAIC(lmerka_random_int_slope_stim_sub_an)
cAIC(lmerka_random_int_slope_all_slopes_sub_an)
cAIC(lmerka_random_int_slope_all_sub_an)
cAIC(lmerka_random_slope_all_in_one_an)


#----  ✓ 1.5 plots for cathodal part------------------------------------------


#plot for full anova for visualization searate factors

ggplot(data = cathodal_dif, aes(x=Stimulation,y = Difference, fill = Stimulation)) +
  geom_boxplot(outlier.size=-1) +
  theme_minimal() +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  scale_y_continuous(breaks = seq(-2.5, 1.5, 0.5), limits=c(-2.5, 1.5))

ggplot(data = cathodal_dif, aes(x=Type,y = Difference, fill = Type)) +
  geom_boxplot(outlier.size=-1) +
  theme_minimal() +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  scale_y_continuous(breaks = seq(-2.5, 1.5, 0.25), limits=c(-2.5, 1.5))

# cathodal_dif_not_postex <- subset(cathodal_dif, Type != "Post-Ex")


# for reordering boxes
# cathodal_dif_plot <- cathodal_dif_not_postex
cathodal_dif_plot <- cathodal_dif

cathodal_dif_plot[Type == 'Computer', 'Type'] <- 
  paste0('F', cathodal_dif_plot[Type == 'Computer', Type]) 

cathodal_dif_plot[Choice == 'Selected', 'Type'] <- 
  paste0('S_', cathodal_dif_plot[Choice == 'Selected', Type])  
cathodal_dif_plot[Choice == 'Rejected', 'Type'] <- 
  paste0('R_', cathodal_dif_plot[Choice == 'Rejected', Type]) 


ggplot(data = cathodal_dif_plot, aes(x=Type,y = Difference, fill = Stimulation)) +
  geom_boxplot(outlier.size=-1) +
  scale_y_continuous(breaks = seq(-2.5, 1.5, 0.25), limits=c(-2.5, 1.5))


ggplot(data = cathodal_dif_plot, aes(x=Type,y = Difference, fill = Type)) +
  geom_boxplot(outlier.size=-1) +
  scale_y_continuous(breaks = seq(-2.5, 1.5, 0.25), limits=c(-2.5, 1.5)) +
  theme_minimal() 

viridis(4)

# plot for 2x2 ANOVA only for difficult choices

# View(cathodal_dif[Type == "Difficult",])

ggplot(data = cathodal_dif[Type == "Difficult",], aes(x=Choice,y = Difference, fill = Stimulation)) +
  geom_boxplot(outlier.size=-1) +
  theme_minimal() +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  scale_y_continuous(breaks = seq(-2.5, 1.5, 0.25), limits=c(-2.5, 1.5))

ggplot(data = cathodal_dif, aes(x=Choice,y = Difference, fill = Choice)) +
  geom_boxplot(outlier.size=-1) +
  theme_minimal() +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  scale_y_continuous(breaks = seq(-2.5, 1.5, 0.5), limits=c(-2.5, 1.5))


#----  ✓ 1.5 plots for anodal part------------------------------------------

# different stimulation
ggplot(data = anodal_dif, aes(x = Stimulation, y = Difference, fill = Stimulation)) +
  geom_boxplot(outlier.size = -1) +
  theme_minimal() +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  scale_y_continuous(breaks = seq(-2.5, 1.5, 0.5),
                     limits = c(-2.5, 1.5))

# different type of choice
ggplot(data = anodal_dif, aes(x = Type, y = Difference, fill = Type)) +
  geom_boxplot(outlier.size = -1) +
  theme_minimal() +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  scale_y_continuous(breaks = seq(-2.5, 1.5, 0.25),
                     limits = c(-2.5, 1.5))


# plot only for difficult choices


ggplot(data = anodal_dif[Type == "Difficult",], aes(x=Choice,y = Difference, fill = Stimulation)) +
  geom_boxplot(outlier.size=-1) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(-2.5, 1.5, 0.25), limits=c(-2.5, 1.5))



# final plot with all conditions
anodal_dif_plot <- anodal_dif

anodal_dif_plot[Type == 'Computer', 'Type'] <-
  paste0('F', anodal_dif_plot[Type == 'Computer', Type])

anodal_dif_plot[Choice == 'Selected', 'Type'] <-
  paste0('S_', anodal_dif_plot[Choice == 'Selected', Type])
anodal_dif_plot[Choice == 'Rejected', 'Type'] <-
  paste0('R_', anodal_dif_plot[Choice == 'Rejected', Type])


ggplot(data = anodal_dif_plot, aes(x = Type, y = Difference, fill = Stimulation)) +
  geom_boxplot(outlier.size = -1) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(-2.5, 1.5, 0.25),
                     limits = c(-2.5, 1.5))


viridis(4)

#----  ✓ 1.5 plots for summary of two experiments------------------------------------------

cathodal_dif_temp <- cathodal_dif[order(Sub, Stimulation),]
cathodal_dif_temp[Stimulation == 'sham', 'Stimulation'] <-
  paste0(cathodal_dif_temp[Stimulation == 'sham', Stimulation], sep='_', 'cat')


anodal_dif_temp <- anodal_dif[order(Sub, Stimulation),]
# anodal_dif_temp[, 'Sub'] <- paste0(anodal_dif_temp[, Sub], sep='_', 'an')
anodal_dif_temp[, 'Stimulation'] <-
  paste0(anodal_dif_temp[, Stimulation], sep='_', 'an')

View(cathodal_dif_temp)
View(anodal_dif_temp)

# colnames(cathodal_dif_temp)
# names(cathodal_dif_temp)[names(cathodal_dif_temp) == 'Stimulation']

cat_an_dif_temp <- merge(cathodal_dif_temp, anodal_dif_temp, 
                         all.x = T, all.y = T)
View(cat_an_dif_temp)

cat_an_plot <- cat_an_dif_temp[Type == "Difficult" | Type == "Easy", ]
# View(cat_an_dif[Choice == 'Rejected' & Stimulation != 'sham_an' |  Stimulation != 'sham_cat')
# 
View(cat_an_plot)

cat_an_plot[Choice == 'Selected', 'Type'] <-
  paste0('S_', cat_an_plot[Choice == 'Selected', Type])
cat_an_plot[Choice == 'Rejected', 'Type'] <-
  paste0('R_', cat_an_plot[Choice == 'Rejected', Type])

cat_an_plot[Stimulation == 'sham_cat', 'Stimulation'] <-
  paste0('1_', cat_an_plot[Stimulation == 'sham_cat', Stimulation])
cat_an_plot[Stimulation == 'tDCS_cat', 'Stimulation'] <-
  paste0('2_', cat_an_plot[Stimulation == 'tDCS_cat', Stimulation])
cat_an_plot[Stimulation == 'sham_an', 'Stimulation'] <-
  paste0('3_', cat_an_plot[Stimulation == 'sham_an', Stimulation])
cat_an_plot[Stimulation == 'tDCS_an', 'Stimulation'] <-
  paste0('4_', cat_an_plot[Stimulation == 'tDCS_an', Stimulation])


cat_an_plot2 <- cat_an_plot[Type != 'S_Easy', ]


ggplot(data = cat_an_plot2, aes(x = Type, y = Difference, fill = Stimulation)) +
  geom_boxplot(outlier.size = -1) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(-2.5, 1.5, 0.25),
                     limits = c(-2.5, 1.5))


#---    t tests for difference in difficult rejected for checking--------------

anodal_dif_sham_rej <- anodal_dif[Type == "Difficult" & Stimulation == "sham" & Choice == 'Rejected', Difference]
anodal_dif_tDCS_rej <- anodal_dif[Type == "Difficult" & Stimulation == "tDCS" & Choice == 'Rejected', Difference]

cohensD <- cohensD(cathodal_dif_sham_rej, cathodal_dif_tDCS_rej)
cohensD(cathodal_dif_sham_rej, cathodal_dif_tDCS_rej)
sdpool <- sqrt( ( sd(cathodal_dif_sham_rej)^2 + sd(cathodal_dif_tDCS_rej)^2 ) /2 )
sqrt( ( sd(cathodal_dif_sham_rej)^2 + sd(cathodal_dif_tDCS_rej)^2 ) /2 )


t.test(x = anodal_dif_tDCS_rej,  
       y = anodal_dif_sham_rej,  
       paired = T, var.equal = FALSE,
       conf.level = 0.95, alternative = "less")

#----   t tests for comparisson two experiments------------------

t.test(x = cat_an_plot[Type == "R_Difficult" & Stimulation == "1_sham_cat", Difference],  
       y = cat_an_plot[Type == "R_Difficult" & Stimulation == "3_sham_an", Difference],  
       paired = F, var.equal = FALSE,
       conf.level = 0.95, alternative = "two.sided")


t.test(x = cat_an_plot[Type == "R_Difficult" & Stimulation == "2_tDCS_cat", Difference],  
       y = cat_an_plot[Type == "R_Difficult" & Stimulation == "4_tDCS_an", Difference],  
       paired = F, var.equal = FALSE,
       conf.level = 0.95, alternative = "two.sided")


t.test(x = cat_an_plot[Type == "R_Difficult" & Stimulation == "1_sham_cat", Difference],  
       y = cat_an_plot[Type == "R_Difficult" & Stimulation == "2_tDCS_cat", Difference],  
       paired = T, var.equal = FALSE,
       conf.level = 0.95, alternative = "two.sided")
