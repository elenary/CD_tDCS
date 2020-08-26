# install.packages(c('data.table', 'ez'))
# install.packages("installr")
install.packages('lsr')
install.packages('effsize')
install.packages("tidyr")
install.packages('plotly')
install.packages('TeachingDemos')
install.packages("pwr")
install.packages("psychReport")
install.packages("viridis")
install.packages('HH')
install.packages('rstatix')
install.packages('nlme')
install.packages('lme4')
install.packages('emmeans')
install.packages('contrast')
install.packages('stats')
install.packages('RColorBrewer')
install.packages('lmerTest')
install.packages('nlme')

# install.packages("rlang")

library(data.table)
library(ez)
library(pwr)
library(effsize)
library(tidyr)
library(lsr)
library(plotly)
# library(TeachingDemos)
library(psychReport)
library(rColorBrewer, lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
# library(reshape2)
library(viridis)
library('HH')
library(reshape2)
library(rstatix)
library(nlme)
library(emmeans)
library(contrast)
library(stats)
library(lme4)
library(lmerTest)
library(nlme)

# here is data.table, dplyr ver will appear later

# locale

Sys.setlocale(category = 'LC_ALL', locale = "RU_ru")
Sys.getlocale()

# color

display.brewer.all()
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







#----   function for gain score ANOVA output (by me) -----------

gs_ANOVA_full_and_posthoc <- function(dataset, dv_numbers, PostHocFactor1, PostHocFactor2) 
  {
  #ANOVA
  if (dv_numbers == 3) {
    gs_ANOVA_full <- ezANOVA(data = dataset, dv = Difference, wid = Sub,
            within = .(Stimulation, Type, Choice), return_aov = T, type = 3)
  } 
  
  if (dv_numbers == 2) {
    gs_ANOVA_full <- ezANOVA(data = dataset, dv = Difference, wid = Sub,
                             within = .(Stimulation, Choice), return_aov = T, type = 3)
  }
  #(effect size = ges)
  
  #residuals
  #(extracting residuals for check it on normality)
  
    gs_ANOVA_resid <-
    proj(gs_ANOVA_full$aov)[[3]][, "Residuals"]
    attr(gs_ANOVA_resid, "ATT") <- NULL
  
  #areResidualsNormal
  gs_ANOVA_AreResidNormal <- shapiro.test(gs_ANOVA_resid)
  
  #(residHist -- better do after function call hist(a$residuals))

  
  # #posthoc for interaction
  # dataset <- cathodal_dif
  # dataset_interaction <- as.data.table(dataset)
  # PostHocFactor1 <- 'Type'
  # PostHocFactor2 <- 'Choice'
  # dataset_interaction[, 'PostHocFactor1','_','PostHocFactor2'] <-
  #   paste0(dataset_interaction[, PostHocFactor1], '_', dataset_interaction[, PostHocFactor1])
  # 
  # 
  
  
  #result
  list(gs_ANOVA = gs_ANOVA_full,
       residuals = gs_ANOVA_resid,
       areResidualsNormal = gs_ANOVA_AreResidNormal)
      # newDTwithInteraction = dataset_interaction)
  
}
#----   0 ANALYSIS PLANS AND INTERATIONS     --------------------------

# Iterations:
# 1. - Only Anova separately for different conditions (easy, difficult, computer, post-ex)
# 2. - gain score / repeated measurment Anova for whole dataset
# 3  + ANCOVA for whole dataset + ANCOVA for different conditions
# 4. + LME for whole dataset for differences + LME for differen conditions 
#+ pairwise t test for seprate hypothesis

#----   1 CATHODAL ANALYSIS PART   --------------------------
#---- 1.1.1 ✖ Cathodal 2x4x2 gain score ANOVA for difference    ---------

gs_ANOVA_cathodal <- gs_ANOVA_full_and_posthoc(cathodal_dif, 3, Type, Choice)
gs_ANOVA_cathodal$gs_ANOVA
gs_ANOVA_cathodal$gs_ANOVA$ANOVA
gs_ANOVA_cathodal$residuals
hist(gs_ANOVA_cathodal$residuals) 
# View(cathodal_dif)



#-----  1.1.2 ✓ cathodal 2x4x2 ANCOVA (main)---------------------------

cathodal_wide <- dcast(cathodal, Sub + Stimulation + Type + Choice ~ Rating)
cathodal_wide <- as.data.table(cathodal_wide)

# ezANOVA(data = cathodal, dv = Evaluation, wid = Sub,
#         within = .(Stimulation, Type, Choice, Rating), 
#         return_aov = T, type = 3)

ezANOVA(data = cathodal_wide, dv = R2, wid = Sub,
        within = .(Stimulation, Type, Choice), 
        within_covariates = R1, return_aov = T, type = 3)

head(cathodal)


#----   post hoc for cathodal significant interaction Type:Choice    ------------------------------------

cathodal_dif_interaction <- cathodal_dif
cathodal_dif_interaction[, 'Type_Choice'] <-
  paste0(cathodal_dif_interaction[, Type], '_', cathodal_dif_interaction[, Choice])

View(cathodal_dif_interaction)

# pwtt_interaction_fdr_gs <- pairwise.t.test(
#   cathodal_dif_interaction$Difference,
#   g = cathodal_dif_interaction$Type_Choice,
#   p.adjust.method = 'fdr',
#   paired = T
# )
# 
# pwtt_interaction_fdr_gs


pairwTTest_fdr_gs_withtvalues <-
  pairwise.t.test.with.t.and.df(
    x = cathodal_dif_interaction$Difference,
    g = cathodal_dif_interaction$Type_Choice,
    p.adjust.method = 'fdr',
    paired = T
  ) 

pairwTTest_fdr_gs_withtvalues
pairwTTest_fdr_gs_withtvalues$p.value
View(as.data.table(pairwTTest_fdr_gs_withtvalues$p.value, keep.rownames = T))
# str(pairwTTest_fdr_gs_withtvalues)
pairwTTest_fdr_gs_withtvalues
pairwTTest_fdr_gs_withtvalues$t.value
View(as.data.table(pairwTTest_fdr_gs_withtvalues$t.value, keep.rownames = T))
pairwTTest_fdr_gs_withtvalues$dfs
View(as.data.table(pairwTTest_fdr_gs_withtvalues$dfs, keep.rownames = T))


#  post hoc interaction Stimulation:Type:Choice (it's not significant but interesting to see)

cathodal_dif_interaction2 <- cathodal_dif
cathodal_dif_interaction2[, 'Stimulation_Type_Choice'] <-
  paste0(cathodal_dif_interaction2[, Stimulation], '_', cathodal_dif_interaction2[, Type], '_', cathodal_dif_interaction2[, Choice])

View(cathodal_dif_interaction2)

pairwTTest_fdr_gs_withtvalues2 <-
  pairwise.t.test.with.t.and.df(
    x = cathodal_dif_interaction2$Difference,
    g = cathodal_dif_interaction2$Stimulation_Type_Choice,
    p.adjust.method = 'fdr',
    paired = T
  ) 

pairwTTest_fdr_gs_withtvalues2$p.value

#----   1.1.3 ✖ TRY cathodal 2x4x2 ANCOVA with planned contrasts --------------------
View(cathodal)
View(cathodal_dif)

View(cathodal_wide)

cathodal_wide_factors <- cathodal_wide

levels(as.factor(cathodal_wide$Stimulation))
levels(as.factor(cathodal_wide$Type))
levels(as.factor(cathodal_wide$Choice))

cathodal_wide_factors$Stimulation <- as.factor(cathodal_wide$Stimulation)
cathodal_wide_factors$Type <- as.factor(cathodal_wide$Type)
cathodal_wide_factors$Choice <- as.factor(cathodal_wide$Choice)

head(cathodal_wide_factors)
# View(cathodal_wide_factors)

dummies_Stimulation <- c(rep(-1, 8), rep(1,8), rep(-1, 8), rep(1,8))
dummies_Type <- c(rep(-1/3, 2), rep(1, 2), rep(-1/3, 2), rep(-1/3,2), rep(-1/3, 2), rep(1, 2), rep(-1/3, 2), rep(-1/3,2))
dummies_Choice <- c(rep(c(1,-1),8))

contrasts_Stimulation <- c(-1,1)
contrasts_Type <- c(-1/3,1,-1/3, -1/3)
contrasts_Choice <- c(1,-1)
# contrasts_Type



contrasts(cathodal_wide_factors$Stimulation) <- dummies_Stimulation
contrasts(cathodal_wide_factors$Type) <- dummies_Type
contrasts(cathodal_wide_factors$Choice) <- dummies_Choice

attributes()

# options( contrasts = c('contr.sum', 'contr.poly'))

View(cathodal_wide_factors)

anovka <- ezANOVA(data = cathodal_wide_factors, dv = R2, wid = Sub,
        within = .(Stimulation, Type, Choice), 
        within_covariates = R1, type = 3, return_aov = T)

anovka

#here i realized that actial class of anovka$aov is [1] "aovlist" "listof" , NOT aov!!
#and that ruined all the analysis

# anovka$aov
# summary(anovka$aov)
# str(anovka$aov)


#and switching to aov or lim functions.....
aovka <- aov(R2 ~ Stimulation * Type * Choice + R1 + Error(Sub/(Stimulation*Type*Choice)), 
             data = cathodal_wide_factors)
class(aovka)

aovka2 <- aov(R2 ~ Stimulation * Type * Choice + R1 + Error(Sub/(Stimulation*Type*Choice)),
              data = cathodal_wide_factors)




summary(aovka)
summary(anovka$aov)
anovka

anova(aovka)
class(aovka2)


lmodelka <- lm(R2 ~ Stimulation : Type : Choice + R1, data = cathodal_wide_factors)

summary(lmodelka)
summary.lm(lmodelka)
anova.lm(lmodelka)

emmeans(aovka,  ~ Stimulation * Type * Choice) 
emmeans(aovka2,  ~ Stimulation * Type * Choice)

emmeans(lmodelka,  ~ Stimulation * Type * Choice)

help("models", package = "emmeans")



head(cathodal)
head(cathodal_wide)
View(cathodal)
View(cathodal_wide)

#----   1.1.4 ✓ cathodal LME  -----------------

# modelka <- lmer(R2~R1+Stimulation+Type+Choice, data = cathodal_wide)
# summary(lm(R2~R1+Stimulation+Type+Choice, data = cathodal_wide))
# summary.lm(lm(R2~R1+Stimulation+Type+Choice, data = cathodal_wide))

ezANOVA(data = cathodal_wide, dv = R2, wid = Sub,
        within = .(Stimulation, Type, Choice), 
        within_covariates = R1, return_aov = T, type = 3)


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

#worst by AIC, BIC and LogLik
lmerka_dif <- lmer(Difference ~ Stimulation*Type*Choice + (1 + Stimulation|Sub), 
                   data = cathodal_dif, REML = F)
anova(lmerka_dif)
summary(lmerka_dif)

model.matrix(lmerka_dif)



# lmerka_dif2 <- lmer(Difference ~ Stimulation*Type*Choice + (1 + Stimulation|Sub), 
#                     data = cathodal_dif, REML = F)

lmerka_dif2 <- lmer(Difference ~ Stimulation*Type*Choice + (1 + Stimulation|Sub), 
                   data = cathodal_dif, REML = T)

anova(lmerka_dif2)
summary(lmerka_dif)



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


# model on the subset of data (only data involved in key comparisonsn)


# 1) for full, main:
# difficult rejected cathodal отличаются от difficult rejected sham
#2) for all where difficult: (control these is keeping for cathodal only and for sham only)
# rejected отличаются от selected 
#3) for all where rejected: (control these is keeping for cathodal only and for sham only)
# difficult отличаются от easy 
# difficult отличаются от computer 
# difficult отличаются от post-ex ?

#comparissons
# 1) difficult rejected cathodal отличаются от difficult rejected sham
# 2)difficult rejected cathodal отличаются от difficult selected cathodal
# 3)difficult rejected cathodal отличаются от easy rejected cathodal
# 4)difficult rejected cathodal отличаются от computer rejected cathodal
# 5)difficult rejected cathodal отличаются от post-ex rejected cathodal


# cathodal_wide_cropped <- cathodal_wide[Stimulation == 'tDCS_cat' & Choice == "Rejected" | 
#                                      Stimulation == 'tDCS_cat' & Type == 'Difficult' & Choice == "Selected" | 
#                                      Stimulation == 'sham' & Type == 'Difficult' & Choice == "Rejected" ,]

cathodal_dif_cropped <- cathodal_dif[Choice == "Rejected" | 
                                       Type == 'Difficult' ,]

View(cathodal_dif_cropped)

# lmerka_cropped <- lmer(R2 ~ R1 + Stimulation*Type*Choice + (1 + Stimulation|Sub), 
#                  data = cathodal_dif_cropped, REML = T)


lmerka_dif_cropped <- lmer(Difference ~ Stimulation*Type*Choice + (1 + Stimulation|Sub), 
                           data = cathodal_dif_cropped, REML = F)

anova(lmerka_dif_cropped)
anova(lmerka_cropped, type = 1) #whithout error but still no calculations for interaction of factors

summary(lmerka_cropped)
anova(lmerka_cropped)
head(cathodal_wide_cropped)
View(cathodal_dif_cropped)

#unworking models

# lmer(R2 ~ R1 + Stimulation*Type*Choice + 
#                  (1|Sub) + (1|Stimulation:Sub) + (1|Type:Sub) + (1|Choice:Sub), 
#                data=cathodal_wide)

# lmer(R2 ~ R1 + Stimulation*Type*Choice + 
#        (R1|Sub) + (R1|Stimulation:Sub) + (R1|Type:Sub) + (R1|Choice:Sub), 
#      data=cathodal_wide)


# lmerka2 <- lmer(R2 ~ Stimulation*Type*Choice + 
#                  (1|Sub) + (1|Stimulation:Sub) + (1|Type:Sub) + (1|Choice:Sub), 
#                data=cathodal_wide)



#----   1.1.5 ✓ paired t test for 6 hypothesis  -----------------



# separeted t test for checking 



# 1) for full, main:
# difficult rejected cathodal отличаются от difficult rejected sham


#2) for all where difficult: (control these is keeping for cathodal only and for sham only)
# rejected отличаются от selected 
#3) for all where rejected: (control these is keeping for cathodal only and for sham only)
# difficult отличаются от easy 
# difficult отличаются от computer 
# difficult отличаются от post-ex ?


# 1) difficult rejected cathodal отличаются от difficult rejected sham
# 2)difficult rejected cathodal отличаются от difficult selected cathodal
# 3)difficult rejected cathodal отличаются от easy rejected cathodal
# 4)difficult rejected cathodal отличаются от computer rejected cathodal
# 5)difficult rejected cathodal отличаются от post-ex rejected cathodal

cathodal_dif_cropped_for_ttest <- cathodal_dif[Stimulation == 'tDCS_cat' & Choice == "Rejected" |
                                                 Stimulation == 'tDCS_cat' & Type == 'Difficult' & Choice == "Selected" |
                                                 Stimulation == 'sham' & Type == 'Difficult' & Choice == "Rejected" ,]


cathodal_dif_cropped_for_ttest_interaction <- cathodal_dif_cropped_for_ttest
cathodal_dif_cropped_for_ttest_interaction[, 'Stimulation_Type_Choice'] <-
  paste0(cathodal_dif_cropped_for_ttest_interaction[, Stimulation], '_', 
         cathodal_dif_cropped_for_ttest_interaction[, Type], '_', 
         cathodal_dif_cropped_for_ttest_interaction[, Choice])

View(cathodal_dif_cropped_for_ttest_interaction)





#pairwise t test for cathodal simulation mostly

pairwTTest_fdr_dif_cropped <-
  pairwise.t.test.with.t.and.df(
    x = cathodal_dif_cropped_interaction$Difference,
    g = cathodal_dif_cropped_interaction$Stimulation_Type_Choice,
    p.adjust.method = 'fdr',
    paired = T
  ) 

View(pairwTTest_fdr_dif_cropped$p.value)


pairwTTest_bh_dif_cropped <-
  pairwise.t.test.with.t.and.df(
    x = cathodal_dif_cropped_interaction$Difference,
    g = cathodal_dif_cropped_interaction$Stimulation_Type_Choice,
    p.adjust.method = 'BH',
    paired = T
  ) 

pairwTTest_bh_dif_cropped$p.value


pairwTTest_fdr_dif_cropped <-
  pairwise.t.test.with.t.and.df(
    x = cathodal_dif_cropped_interaction$Difference,
    g = cathodal_dif_cropped_interaction$Stimulation_Type_Choice,
    p.adjust.method = '',
    paired = T
  ) 

pairwTTest_fdr_dif_cropped$p.value


#----   1.2 t tests for difference in difficult rejected for checking--------

#t test  for cathodal, cohens'D and stat power

cathodal_dif_sham_rej <- cathodal_dif[Type == "Difficult" & Stimulation == "sham" & Choice == 'Rejected', Difference]
cathodal_dif_tDCS_rej <- cathodal_dif[Type == "Difficult" & Stimulation == "tDCS_cat" & Choice == 'Rejected', Difference]

cathodal_dif_tDCS_sel <- cathodal_dif[Type == "Difficult" & Stimulation == "tDCS_cat" & Choice == 'Selected', Difference]
cathodal_easy_tDCS_rej <- cathodal_dif[Type == "Easy" & Stimulation == "tDCS_cat" & Choice == 'Rejected', Difference]
cathodal_comp_tDCS_rej <- cathodal_dif[Type == "Computer" & Stimulation == "tDCS_cat" & Choice == 'Rejected', Difference]
cathodal_postex_tDCS_rej <- cathodal_dif[Type == "Post-Ex" & Stimulation == "tDCS_cat" & Choice == 'Rejected', Difference]




# 
# cohensD <- cohensD(cathodal_dif_sham_rej, cathodal_dif_tDCS_rej)
# cohensD(cathodal_dif_sham_rej, cathodal_dif_tDCS_rej)
# sdpool <- sqrt( ( sd(cathodal_dif_sham_rej)^2 + sd(cathodal_dif_tDCS_rej)^2 ) /2 )
# sqrt( ( sd(cathodal_dif_sham_rej)^2 + sd(cathodal_dif_tDCS_rej)^2 ) /2 )

cohen.d(cathodal_dif_sham_rej, cathodal_dif_tDCS_rej, paired = T, hedges.correction = T)
cohen.d(cathodal_dif_sham_rej, cathodal_dif_tDCS_rej, paired = T, hedges.correction = F)


citation("effsize")

t.test(x = cathodal_dif_sham_rej,  
       y = cathodal_dif_tDCS_rej,  
       paired = T, var.equal = FALSE,
       conf.level = 0.95, alternative = "less")

t.test(x = cathodal_dif_tDCS_sel,  
       y = cathodal_dif_tDCS_rej,  
       paired = T, var.equal = FALSE,
       conf.level = 0.95, alternative = "two.sided")

t.test(x = cathodal_dif_tDCS_sel,  
       y = cathodal_dif_tDCS_rej,  
       paired = T, var.equal = FALSE,
       conf.level = 0.95, alternative = "two.sided")


main_test <- t.test(x = cathodal_dif_sham_rej,  
       y = cathodal_dif_tDCS_rej,  
       paired = T, var.equal = FALSE,
       conf.level = 0.95, alternative = "two.sided")

main_test$statistic

t = seq(-5, 5, 0.01)
sim_df <- data.frame(t = t, cdf = dt(t, 16))
ggplot(sim_df, aes(x = t, y = cdf))+
  geom_line()+
  geom_vline(xintercept = main_test$statistic, colour = "#FF2277", size = 3)


Pvalue.norm.sim() 



#stat power 

pwr.t.test(n = 17, d=0.29, sig.level = 0.05)

#----   1.3.1 cathodal 2x2 gain score ANOVA by different choice factor---------------------------

# for checking is gain score ANOVA differ from the initial RM Marco by factors

ANOVA_cathodal_gs_dif <- ezANOVA(data = cathodal_dif[Type == 'Difficult', ], dv = Difference, wid = Sub,
                                 within = .(Stimulation, Choice), return_aov = T, type = 3)
ANOVA_cathodal_gs_dif$ANOVA

#interaction for only difficult choices is significant

ANOVA_cathodal_gs_easy <- ezANOVA(data = cathodal_dif[Type == 'Easy', ], dv = Difference, wid = Sub,
                                 within = .(Stimulation, Choice), return_aov = T, type = 3)
ANOVA_cathodal_gs_easy$ANOVA

ANOVA_cathodal_gs_computer <- ezANOVA(data = cathodal_dif[Type == 'Computer', ], dv = Difference, wid = Sub,
                                  within = .(Stimulation, Choice), return_aov = T, type = 3)
ANOVA_cathodal_gs_computer$ANOVA


ANOVA_cathodal_gs_PostEx <- ezANOVA(data = cathodal_dif[Type == 'Post-Ex', ], dv = Difference, wid = Sub,
                                      within = .(Stimulation, Choice), return_aov = T, type = 3)
ANOVA_cathodal_gs_PostEx$ANOVA


#----   post hoc Stimulation:Choice FOR GAIN SCORE ANOVA by different choice type-----------------

#difficult choices (target)

cathodal_dif_interaction3 <- cathodal_dif[Type == 'Difficult', ]
cathodal_dif_interaction3[, 'Stim_Choice'] <-
  paste0(cathodal_dif_interaction3[, Stimulation], '_', cathodal_dif_interaction3[, Choice])
View(cathodal_dif)
View(cathodal_dif_interaction3)

# pwtt_interaction_fdr_gs_diff <- pairwise.t.test(
#   cathodal_dif_interaction3$Difference,
#   g = cathodal_dif_interaction3$Stim_Choice,
#   p.adjust.method = 'fdr',
#   paired = T
# )
# pwtt_interaction_fdr_gs_diff
# View(as.data.table(pwtt_interaction_fdr_gs[3], keep.rownames=T))



pairwTTest_fdr_gs_difficult <-
  pairwise.t.test.with.t.and.df(
    x = cathodal_dif_interaction3$Difference,
    g = cathodal_dif_interaction3$Stim_Choice,
    p.adjust.method = 'fdr',
    paired = T
  ) 


pairwTTest_fdr_gs_difficult$p.value
View(pairwTTest_fdr_gs_difficult$p.value)
pairwTTest_fdr_gs_difficult$t.value
View(pairwTTest_fdr_gs_difficult$t.value)
pairwTTest_fdr_gs_difficult$dfs



#easy choices

cathodal_dif_interaction4 <- cathodal_dif[Type == 'Easy', ]
cathodal_dif_interaction4[, 'Stim_Choice'] <-
  paste0(cathodal_dif_interaction4[, Stimulation], '_', cathodal_dif_interaction4[, Choice])
View(cathodal_dif_interaction4)


pairwTTest_fdr_gs_easy<-
  pairwise.t.test.with.t.and.df(
    x = cathodal_dif_interaction4$Difference,
    g = cathodal_dif_interaction4$Stim_Choice,
    p.adjust.method = 'fdr',
    paired = T
  ) 


pairwTTest_fdr_gs_easy$p.value
View(pairwTTest_fdr_gs_easy$p.value)
pairwTTest_fdr_gs_easy$t.value
View(pairwTTest_fdr_gs_easy$t.value)
pairwTTest_fdr_gs_easy$dfs



#----   1.3.2 cathodal ANCOVA only difficult  ---------------------------

cathodal_difficult <- cathodal[Type == 'Difficult', ]

cathodal_difficult_wide <- dcast(cathodal_difficult, 
                                 Sub + Stimulation + Type + Choice ~ Rating)

ezANOVA(data = cathodal_difficult_wide, dv = R2, wid = Sub,
        within = .(Stimulation, Choice), 
        within_covariates = R1, return_aov = T)


#----   1.3.3 cathodal ANCOVA only rejected ---------------------------

cathodal_rejected_wide <- dcast(cathodal_rejected, 
                                Sub + Stimulation + Type + Choice ~ Rating)


cathodal_difficult_rejected <- cathodal_rejected[Type == 'Difficult', ]
cathodal_difficult_rejected_wide <- dcast(cathodal_difficult_rejected, 
                                          Sub + Stimulation + Type + Choice ~ Rating)

# View(cathodal_difficult_rejected_wide)
head(cathodal_difficult_rejected_wide)

ezANOVA(data = cathodal_rejected_wide, dv = R2, wid = Sub,
        within = .(Stimulation, Type),
        within_covariates = R1, return_aov = T)

cathodal_rejected_difference <- cathodal_dif[Choice == 'Rejected',]
# View(cathodal_rejected_difference)

ezANOVA(data = cathodal_rejected_difference, dv = Difference, wid = Sub,
        within = .(Stimulation, Type), 
        return_aov = T, type = 3)






#----   1.4 plots for cathodal part------------------------------------------


#--plot for full anova for visualization searate factors

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


# ggplot(data = cogdis_dif_plot, aes(x=Type,y = Difference, fill = Stimulation)) +
#   geom_bar(data = cogdis_dif_plot, position=position_dodge(), stat="identity") +
#   scale_y_continuous(breaks = seq(-2,0.75,0.25)) +
# 
#   geom_errorbar(aes(ymin=soa0-SEm0, ymax=soa0+SEm0),
#                 size=.3,    # Thinner lines
#                   width=.2,
#                 position=position_dodge(.9)) +
# 
#   scale_x_discrete(limits=c("S_Difficult", "S_Easy", "S_Computer",
#                             "R_Difficult", "R_Easy", "R_Computer"))











ggplot(data = cathodal_dif, aes(x=Choice,y = Difference, fill = Choice)) +
  geom_boxplot(outlier.size=-1) +
  theme_minimal() +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  scale_y_continuous(breaks = seq(-2.5, 1.5, 0.5), limits=c(-2.5, 1.5))






#----   ANODAL GAIN SCORE 2x4x2 ANOVA FOR DIFFERENCE--------------------------------------


gs_ANOVA_anodal <- gs_ANOVA_full_and_posthoc(anodal_dif, 3, Type, Choice)
gs_ANOVA_anodal$gs_ANOVA
gs_ANOVA_anodal$gs_ANOVA$ANOVA
gs_ANOVA_anodal$residuals
gs_ANOVA_anodal$areResidualsNormal 
hist(gs_ANOVA_anodal$residuals) 


#----   POST HOC for anodal interaction Type:Choice    ------------------------------------

anodal_dif_interaction <- anodal_dif
anodal_dif_interaction[, 'Type_Choice'] <-
  paste0(anodal_dif_interaction[, Type], '_', anodal_dif_interaction[, Choice])

View(anodal_dif_interaction)

# pwtt_interaction_fdr_gs <- pairwise.t.test(
#   cathodal_dif_interaction$Difference,
#   g = cathodal_dif_interaction$Type_Choice,
#   p.adjust.method = 'fdr',
#   paired = T
# )
# 
# pwtt_interaction_fdr_gs


pairwTTest_fdr_gs_anodal <-
  pairwise.t.test.with.t.and.df(
    x = anodal_dif_interaction$Difference,
    g = anodal_dif_interaction$Type_Choice,
    p.adjust.method = 'fdr',
    paired = T
  ) 

pairwTTest_fdr_gs_anodal
pairwTTest_fdr_gs_anodal$p.value
View(pairwTTest_fdr_gs_anodal$p.value)
# str(pairwTTest_fdr_gs_withtvalues)
pairwTTest_fdr_gs_anodal$t.value
View(pairwTTest_fdr_gs_anodal$t.value)
pairwTTest_fdr_gs_anodal$dfs



#  post hoc interaction Stimulation:Type:Choice (it's not significant but interesting)

anodal_dif_interaction_an2 <- anodal_dif
anodal_dif_interaction_an2[, 'Stimulation_Type_Choice'] <-
  paste0(anodal_dif_interaction_an2[, Stimulation], '_', anodal_dif_interaction_an2[, Type], '_', anodal_dif_interaction_an2[, Choice])

View(anodal_dif_interaction_an2)


#----   ANODAL 2x2 GAIN SCORE ANOVA by different choice factor---------------------------

# for checking is gain score ANOVA differ from the initial RM Marco by factors

ANOVA_anodal_gs_dif <- ezANOVA(data = anodal_dif[Type == 'Difficult', ], dv = Difference, wid = Sub,
                                 within = .(Stimulation, Choice), return_aov = T, type = 3)
ANOVA_anodal_gs_dif$ANOVA
ANOVA_cathodal_gs_dif$ANOVA

#interaction for only difficult choices is significant

ANOVA_anodal_gs_easy <- ezANOVA(data = anodal_dif[Type == 'Easy', ], dv = Difference, wid = Sub,
                                  within = .(Stimulation, Choice), return_aov = T, type = 3)
ANOVA_anodal_gs_easy$ANOVA

ANOVA_anodal_gs_computer <- ezANOVA(data = anodal_dif[Type == 'Computer', ], dv = Difference, wid = Sub,
                                      within = .(Stimulation, Choice), return_aov = T, type = 3)
ANOVA_anodal_gs_computer$ANOVA


ANOVA_anodal_gs_PostEx <- ezANOVA(data = anodal_dif[Type == 'Post-Ex', ], dv = Difference, wid = Sub,
                                    within = .(Stimulation, Choice), return_aov = T, type = 3)
ANOVA_anodal_gs_PostEx$ANOVA


#----   POST HOC Choice FOR GAIN SCORE ANOVA by different choice type-----------------

#difficult choices (target)


pairwTTest_fdr_gs_anodal_difficult <-
  pairwise.t.test.with.t.and.df(
    x = anodal_dif[Type == 'Difficult', ]$Difference,
    g = anodal_dif[Type == 'Difficult', ]$Choice,
    p.adjust.method = 'fdr',
    paired = T
  ) 


pairwTTest_fdr_gs_anodal_difficult$p.value
View(pairwTTest_fdr_gs_anodal_difficult$p.value)
pairwTTest_fdr_gs_anodal_difficult$t.value
View(pairwTTest_fdr_gs_anodal_difficult$t.value)
pairwTTest_fdr_gs_anodal_difficult$dfs

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



#----   plots anodal------------------

#--plot for full anova for visualization searate factors

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


# plot for 2x2 ANOVA only for difficult choices


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


#----   plots for summary of two experiments------------------


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



#-------___________________________-------------------------------
#----   GENERAL COMMENTS  AND EXPLORATORY PART  -------------------------------------------

#----   CATHODAL REP MEASURMENT ANOVA   -------------------------------------------

View(cathodal)
ANOVA_cathodal_rm <- ezANOVA(data = cathodal, dv = Evaluation, wid = Sub,
                             within = .(Stimulation, Type, Choice, Rating), return_aov = T, type = 3)
ANOVA_cathodal_rm
summary(ANOVA_cathodal_rm$aov)


cathodal_interaction_rm <- cathodal
cathodal_interaction_rm[, 'Stimulation_Type_Choice_Rating'] <-
  paste0(cathodal_interaction_rm[, Stimulation], '_', cathodal_interaction_rm[, Type], '_', cathodal_interaction_rm[, Choice], '_', cathodal_interaction_rm[, Rating])

View(cathodal_interaction_rm)

pwtt_interaction_rm <- pairwise.t.test(
  cathodal_interaction_rm$Evaluation,
  g = cathodal_interaction_rm$Stimulation_Type_Choice_Rating,
  p.adjust.method = 'bonf',
  paired = T
)

cathodal_interaction_rm$Stimulation_Type_Choice_Rating

# View(pwtt_interaction2)
pwtt_interaction_rm 
str(pwtt_interaction_rm )

pwtt_interaction_rm <- as.data.table(pwtt_interaction_rm[3], keep.rownames=T)
View(pwtt_interaction_rm)




#----   ANODAL ANCOVA ---------------------------

anodal_wide <- dcast(anodal, Sub + Stimulation + Type + Choice ~ Rating)
View(anodal_wide)

anodal_difficult <- cathodal[Type == 'Difficult', ]

anodal_difficult_wide <- dcast(anodal_difficult, 
                               Sub + Stimulation + Type + Choice ~ Rating)

ezANOVA(data = anodal, dv = Evaluation, wid = Sub,
        within = .(Stimulation, Type, Choice, Rating), 
        return_aov = T, type = 3)


ezANOVA(data = anodal_wide, dv = R2, wid = Sub,
        within = .(Stimulation, Type, Choice), 
        within_covariates = R1, return_aov = T)


ezANOVA(data = anodal_difficult_wide, dv = R2, wid = Sub,
        within = .(Stimulation, Choice), 
        within_covariates = R1, return_aov = T)

#----   CATHODAL REP MEASURMENT ANOVA by different choice factor (as Marco did)-----------------------------------------

cathodal_difficult <- cathodal[Type == 'Difficult', ]

anova_rm_diff <- ezANOVA(data = cathodal[Type == 'Difficult', ], dv = Evaluation, wid = Sub,
        within = .(Stimulation, Choice, Rating), return_aov = T, type = 3)

anova_rm_diff$ANOVA

ezANOVA(data = cathodal[Type == 'Easy', ], dv = Evaluation, wid = Sub,
        within = .(Stimulation, Choice, Rating), return_aov = T, type = 3)

ezANOVA(data = cathodal[Type == 'Computer', ], dv = Evaluation, wid = Sub,
        within = .(Stimulation, Choice, Rating), return_aov = T, type = 3)

ezANOVA(data = cathodal[Type == 'Post-Ex', ], dv = Evaluation, wid = Sub,
        within = .(Stimulation, Choice, Rating), return_aov = T, type = 3)

ANOVA_cathodal_difficult_rm <- ezANOVA(data = cathodal_difficult, dv = Evaluation, wid = Sub,
                                       within = .(Stimulation, Choice, Rating), return_aov = T, type = 3)

ezANOVA(data = cathodal, dv = Evaluation, wid = Sub,
        within = .(Stimulation, Choice, Rating, Type), return_aov = T, type = 3)
ANOVA_cathodal_difficult_rm
summary(ANOVA_cathodal_rm$aov)


cathodal_interaction_rm_diff <- cathodal[Type == 'Difficult', ]
cathodal_interaction_rm_diff[, 'Stimulation_Choice_Rating'] <-
  paste0(cathodal_interaction_rm_diff[, Stimulation], '_', 
         cathodal_interaction_rm_diff[, Choice], '_', 
         cathodal_interaction_rm_diff[, Rating])

View(cathodal_interaction_rm_diff)

pairwTTest_fdr_rm_diff <-
  pairwise.t.test.with.t.and.df(
    x = cathodal_interaction_rm_diff$Evaluation,
    g = cathodal_interaction_rm_diff$Stimulation_Choice_Rating,
    p.adjust.method = 'fdr',
    paired = T
  ) 

View(pairwTTest_fdr_rm_diff$p.value)
pairwTTest_fdr_gs_withtvalues$p.value




#----   looking at differnt corrections-------------------------------

#bonferonni correction

pwtt_interaction2 <- pairwise.t.test(
  cathodal_dif_interaction2$Difference,
  g = cathodal_dif_interaction2$Stimulation_Type_Choice,
  p.adjust.method = 'bonf',
  paired = T
)


#pool.sd = False by default

# View(pwtt_interaction2)
pwtt_interaction2
str(pwtt_interaction2)

pwtt_interaction2_pv <- as.data.table(pwtt_interaction2[3], keep.rownames=T)
View(pwtt_interaction2_pv)

#holm correction

pwtt_interaction2_holm <- pairwise.t.test(
  cathodal_dif_interaction2$Difference,
  g = cathodal_dif_interaction2$Stimulation_Type_Choice,
  p.adjust.method = 'holm',
  paired = T
)

pwtt_interaction2_pv_holm <- as.data.table(pwtt_interaction2_holm[3], keep.rownames=T)
View(pwtt_interaction2_pv_holm)

#WITHOUT  correction

pwtt_interaction2_none <- pairwise.t.test(
  cathodal_dif_interaction2$Difference,
  g = cathodal_dif_interaction2$Stimulation_Type_Choice,
  p.adjust.method = 'none',
  paired = T
)

pwtt_interaction2_pv_none <- as.data.table(pwtt_interaction2_none[3], keep.rownames=T)
View(pwtt_interaction2_pv_none)


# pwtt_interaction2_pv_none[10,4]

pwtt_interaction2_pv_none_tvalues <- pairwise.t.test.with.t.and.df(x = cathodal_dif_interaction2$Difference, 
                                                                   g = cathodal_dif_interaction2$Stimulation_Type_Choice, 
                                                                   p.adjust.method = 'none', paired = T) 
str(pwtt_interaction2_pv_none_tvalues)
View(pwtt_interaction2_pv_none_tvalues$t.value)




#----   and other corrections--------------
# pwtt_interaction2_hochberg <- pairwise.t.test(
#   cathodal_dif_interaction2$Difference,
#   g = cathodal_dif_interaction2$Stimulation_Type_Choice,
#   p.adjust.method = 'hochberg',
#   paired = T
# )
# 
# pwtt_interaction2_pv_hochberg <- as.data.table(pwtt_interaction2_hochberg[3], keep.rownames=T)
# View(pwtt_interaction2_pv_hochberg)


# pwtt_interaction2_BH <- pairwise.t.test(
#   cathodal_dif_interaction2$Difference,
#   g = cathodal_dif_interaction2$Stimulation_Type_Choice,
#   p.adjust.method = 'BH',
#   paired = T
# )
# 
# pwtt_interaction2_pv_BH <- as.data.table(pwtt_interaction2_BH[3], keep.rownames=T)
# View(pwtt_interaction2_pv_BH)
# 
# 
# 
# pwtt_interaction2_BY <- pairwise.t.test(
#   cathodal_dif_interaction2$Difference,
#   g = cathodal_dif_interaction2$Stimulation_Type_Choice,
#   p.adjust.method = 'BY',
#   paired = T
# )
# 
# pwtt_interaction2_pv_BY <- as.data.table(pwtt_interaction2_BY[3], keep.rownames=T)
# View(pwtt_interaction2_pv_BY)
#----   trying one more long format for post hoc pairwise t tests-------

# cathodal_dif_verylong <-
#   melt(
#     cathodal_dif,
#     id.vars = c("Sub", "Difference"),
#     measure.vars = c("Stimulation", "Type", "Choice")
#   )
# View(cathodal_dif_verylong)
# 
# 
# pairwise.t.test(Cathodal_dif_long[, Difference],
#                 g = Cathodal_dif_long[, value],
#                 p.adjust.method = 'bonf') 


#-------is there difference between ANOVA type 3 and 2? No (there)---------

ANOVA_cathodal_gs <- ezANOVA(data = cathodal_dif, dv = Difference, wid = Sub,
                             within = .(Stimulation, Type, Choice), return_aov = T, type = 2)
ANOVA_cathodal_gs
summary(ANOVA_cathodal_gs$aov)

gs_ANOVA_full_and_posthoc()