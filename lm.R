library('cAIC4')
library('tidyverse')
# different tryes of different slopes and intercepts

cathodal_dif %>% 
  lmer(Difference ~ Stimulation*Type*Choice + (1|Sub), 
       data = ., REML = T) ->lmer_random_int_only2 # ok


cathodal_dif %>% 
  lmer(Difference ~ Stimulation*Type*Choice + (1 + 1|Sub), 
       data = ., REML = T) ->lmer_random_int_only # ok


cathodal_dif %>% 
  lmer(Difference ~ Stimulation*Type*Choice + (1||Sub), 
       data = ., REML = T) -> lmer_random_int_only_uncor #error

cathodal_dif %>% 
  lmer(Difference ~ Stimulation*Type*Choice + (Stimulation|Sub), 
       data = ., REML = T) -> lmer_random_slope_stim #ok, chosen


cathodal_dif %>% 
  filter(Type !="Post-Ex") %>% 
  lmer(Difference ~ Stimulation*Type*Choice + (Stimulation|Sub), 
       data = ., REML = T) -> lmer_random_slope_stim_without_postex #ok, chosen


cathodal_dif %>% 
  filter(Type !="Post-Ex") %>% 
  filter(Choice =="Selected") %>% 
  lmer(Difference ~ Stimulation*Type + (Stimulation|Sub), 
       data = ., REML = T) -> lmer_random_slope_stim_without_postex_selected
anova(lmer_random_slope_stim_without_postex_selected)

cathodal_dif %>% 
  filter(Type !="Post-Ex") %>% 
  filter(Type !="Computer") %>% 
  lmer(Difference ~ Stimulation*Type*Choice + (Stimulation|Sub), 
       data = ., REML = T) -> lmer_random_slope_stim_without_postex_compt #singular
anova(lmer_random_slope_stim_without_postex_compt)

cathodal_dif %>% 
  lmer(Difference ~ Stimulation*Type*Choice + (Stimulation||Sub), 
       data = ., REML = T) -> lmer_random_slope_stim_uncor #ok


cathodal_dif %>% 
  lmer(Difference ~ Stimulation*Type*Choice + (1 + Stimulation|Sub), 
       data = ., REML = T) -> lmer_random_int_slope_stim #ok, chosen

cathodal_dif %>% 
  lmer(Difference ~ Stimulation*Type*Choice + (1 + Stimulation||Sub), 
       data = ., REML = T) -> lmer_random_int_slope_stim_uncor #ok

cathodal_dif %>% 
  lmer(Difference ~ Stimulation*Type*Choice + (1 + Stimulation*Type*Choice|Sub), 
       data = ., REML = T) -> lmer_random_slopes_all #the random-effects parameters unidentifiable

cathodal_dif %>% 
  lmer(Difference ~ Stimulation*Type*Choice + (1 + Stimulation*Type*Choice||Sub), 
       data = ., REML = T) -> lmer_random_slopes_all_uncor #the random-effects parameters unidentifiable

cathodal_dif %>% 
  lmer(Difference ~ Stimulation*Type*Choice + (Stimulation*Type*Choice||Sub), 
       data = ., REML = T)  #the random-effects parameters unidentifiable

cathodal_dif %>% 
  lmer(Difference ~ Stimulation*Type*Choice + (Stimulation + Type + Choice|Sub), 
       data = ., REML = T) -> lmer_random_slopes1_all #singular

cathodal_dif %>% 
  lmer(Difference ~ Stimulation*Type*Choice + (1+ Stimulation + Type + Choice|Sub), 
       data = ., REML = T) -> lmer_random_int_slopes1_all #singular

cathodal_dif %>% 
  lmer(Difference ~ Stimulation*Type*Choice + (Stimulation + Type + Choice|Sub), 
       data = ., REML = T) -> lmer_random_int_all #singular

cathodal_dif %>% 
  lmer(Difference ~ Stimulation*Type*Choice + (Stimulation + Type + Choice||Sub), 
       data = ., REML = T) -> lmer_random_int_all_uncor #singular

cathodal_dif %>% 
  lmer(Difference ~ Stimulation*Type*Choice + ((Stimulation|Sub) + (Type|Sub) + (Choice|Sub)), 
       data = ., REML = T) -> lmer_random_slopes_all_separet #singular

anova(lmer_random_slopes_all_separet)

anova(lmer_random_int_only2,
      lmer_random_int_only, 
      lmer_random_slope_stim,
      lmer_random_slope_stim_uncor,
      lmer_random_int_slope_stim,
      lmer_random_int_slope_stim_uncor,
      # lmer_random_int_only_uncor,
      # lmer_random_slopes_all,
      # lmer_random_slopes_all_uncor,
      # lmer_random_int_all,
      # lmer_random_slopes_all_separet, 
      test = "LRT")

lmer_random_slope_stim_without_postex
summary(lmer_random_slope_stim_without_postex)
anova(lmer_random_slope_stim_without_postex)

# ---c AIC ----

cAIC(lmer_random_int_only2)
cAIC(lmer_random_int_only)
cAIC(lmer_random_slope_stim)
cAIC(lmer_random_slope_stim_uncor)
cAIC(lmer_random_int_slope_stim)
cAIC(lmer_random_int_slope_stim_uncor)
cAIC(lmer_random_slope_stim_without_postex) #final


cathodal_dif %>% 
  filter(Type !="Post-Ex") %>% 
  lmer(Difference ~ Stimulation*Type*Choice + (Stimulation|Sub), 
       data = ., REML = T) -> lmer_random_slope_stim_without_postex #ok, chosen

options(scipen = 999)

r.squaredGLMM(lmer_random_slope_stim_without_postex)
lmer_random_slope_stim_without_postex
anova(lmer_random_slope_stim_without_postex)
effectsize::eta_squared(lmer_random_slope_stim_without_postex, partial = TRUE)

summary(lmer_random_slope_stim_without_postex)

# only rejected

cathodal_dif %>% 
  filter(Choice == "Rejected") %>% 
  lmer(Difference ~ Stimulation*Type + (1 + Stimulation|Sub), 
       data = ., REML = T) ->lmer6 #singular
isSingular(lmer6)

cAIC(lmer6)

cathodal_dif %>% 
  filter(Choice == "Rejected") %>% 
  lmer(Difference ~ Stimulation*Type + (1 + Stimulation||Sub), 
       data = ., REML = T) ->lmer13 #singular
isSingular(lmer13)

cAIC(lmer13)


cathodal_dif %>% 
  filter(Choice == "Rejected") %>% 
  lmer(Difference ~ Stimulation*Type + (1 + Stimulation + Stimulation|Sub), 
       data = ., REML = T) ->lmer11 #singular

cathodal_dif %>% 
  filter(Choice == "Rejected") %>% 
  lmer(Difference ~ Stimulation*Type + (1|Sub), 
       data = ., REML = T) ->lmer10 #ok


anova(lmer10)
anova(lmer6, lmer13, lmer10, lmer11)

effectsize::eta_squared(lmer10, partial = TRUE)
r.squaredGLMM(lmer10)


cAIC(lmer10)


cathodal_dif %>% 
  filter(Choice == "Rejected") %>% 
  lmer(Difference ~ Stimulation*Type + (Stimulation*Type|Sub), 
       data = ., REML = T) -> lmer15 #nottttt
anova(lmer12)

#covariate

cathodal_wide %>% 
  lmer(R2 ~ R1 + Stimulation*Type*Choice + (1 + Stimulation|Sub), data = ., REML = T) %>% 
  anova() # yes but no



anova(lmer0, lmer0_2, lmer7, lmer8, lmer9)
anova(lmer6, lmer13, lmer11, lmer10, lmer12)
# new variables -- Pradigm and Difficulty

View(cathodal_dif)

cathodal_dif %>% 
  mutate(Paradigm = ifelse(Type == "Difficult", "experimental", "control")) -> 
  cathodal_dif_exp

cathodal_dif %>% 
  mutate(Paradigm = ifelse(Type == "Difficult" | Type == "Easy", "experimental","control")) %>% 
  mutate(Difficulty = ifelse(Type == "Difficult" | Type == "Easy", Type, "Difficult")) -> 
  cathodal_dif_exp2

View(cathodal_dif_exp2)


lmer2 <- lmer(Difference ~ Stimulation*Paradigm*Choice + (1 + Stimulation|Sub), 
              data = cathodal_dif_exp, REML = T)

anova(lmer2)
car::Anova(lmer2)

cathodal_dif_exp %>% 
  filter(Choice == "Rejected") %>% 
  lmer(Difference ~ Stimulation*Paradigm + (1 + Stimulation|Sub), 
       data = ., REML = T) -> lmer3

anova(lmer3)
car::Anova(lmer3)

cathodal_dif_exp2 %>% 
  lmer(Difference ~ Stimulation*Choice*Paradigm*Difficulty + (1 + Stimulation|Sub), 
       data = ., REML = T) ->lmer4

anova(lmer4)
car::Anova(lmer4)

cathodal_dif_exp2 %>% 
  filter(Choice == "Rejected") %>% 
  lmer(Difference ~ Stimulation*Paradigm*Difficulty + (1 + Stimulation|Sub), 
       data = ., REML = T) ->lmer5

anova(lmer5)
car::Anova(lmer5)


#--------unsuccessful trials for covariate design-------

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


# model for anodal

anodal_dif %>% 
  filter(Type !="Post-Ex") %>% 
  lmer(Difference ~ Stimulation*Type*Choice + (Stimulation|Sub), 
       data = ., REML = T) -> lmer_random_slope_stim_anodal #ok, chosen
  
anova(lmer_random_slope_stim_anodal)
summary(lmer_random_slope_stim_anodal)

cAIC(lmer_random_slope_stim_anodal) #final

r.squaredGLMM(lmer_random_slope_stim_anodal)
effectsize::eta_squared(lmer_random_slope_stim_anodal, partial = TRUE)




