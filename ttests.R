library("tidyverse")
library(data.table)

options(scipen = 999)

# 1.1 cathodal ------------------------------------------------------------

alpha_corrected <- 0.05/4
alpha_corrected

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


#2) for all where rejected difficult: 

# rejected in difficult less than selected in difficult

cathodal_difficult_rej <- cathodal_dif[Type == "Difficult" & Choice == 'Rejected', Difference]
cathodal_difficult_sel <- cathodal_dif[Type == "Difficult" & Choice == 'Selected', Difference]

t.test(x = cathodal_difficult_rej,  
        y = cathodal_difficult_sel,  
        paired = T, var.equal = FALSE,
        conf.level = 1-alpha_corrected, alternative = "less")


cathodal_difficult_rej
cathodal_difficult_sel

cohen.d( cathodal_difficult_sel, cathodal_difficult_rej, paired = T, hedges.correction = F)


# rejected in difficult less than rejected in easy

cathodal_easy_rej <- cathodal_dif[Type == "Easy" & Stimulation == "tDCS_cat" & Choice == 'Rejected', Difference]
cathodal_comp_rej <- cathodal_dif[Type == "Computer" & Stimulation == "tDCS_cat" & Choice == 'Rejected', Difference]


t.test(x = cathodal_difficult_rej,  
       y = cathodal_easy_rej,  
       paired = T, var.equal = FALSE,
       conf.level = 1-alpha_corrected, alternative = "less")

cohen.d( cathodal_easy_rej, cathodal_difficult_rej, paired = T, hedges.correction = F)


# rejected in difficult less than rejected in computer

t.test(x = cathodal_difficult_rej,  
       y = cathodal_comp_rej,  
       paired = T, var.equal = FALSE,
       conf.level = 1-alpha_corrected, alternative = "less")

cohen.d( cathodal_comp_rej, cathodal_difficult_rej, paired = T, hedges.correction = F)


# 1.2 anodal --------------------------------------------------------------


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


#2) for all where difficult: 
# rejected in difficult less than selected in difficult

anodal_difficult_rej <- anodal_dif[Type == "Difficult" & Choice == 'Rejected', Difference]
anodal_difficult_sel <- anodal_dif[Type == "Difficult" & Choice == 'Selected', Difference]

t.test(x = anodal_difficult_rej,  
       y = anodal_difficult_sel,  
       paired = T, var.equal = FALSE,
       conf.level = 1-alpha_corrected, alternative = "less")

# options(scipen = 999)

cohen.d( anodal_difficult_rej, anodal_difficult_sel, paired = T, hedges.correction = F)

# rejected in difficult less than rejected in easy

anodal_easy_rej <- anodal_dif[Type == "Easy" & Choice == 'Rejected', Difference]
anodal_comp_rej <- anodal_dif[Type == "Computer" & Choice == 'Rejected', Difference]
# anodal_postex_rej <- anodal_dif[Type == "Post-Ex" & Choice == 'Rejected', Difference]


t.test(x = anodal_difficult_rej,  
       y = anodal_easy_rej,  
       paired = T, var.equal = FALSE,
       conf.level = 1-alpha_corrected, alternative = "less")

cohen.d(anodal_difficult_rej, anodal_easy_rej, paired = T, hedges.correction = F)

# rejected in difficult less than rejected in computer

t.test(x = anodal_difficult_rej,  
       y = anodal_comp_rej,  
       paired = T, var.equal = FALSE,
       conf.level = 1-alpha_corrected, alternative = "less")

cohen.d( anodal_difficult_rej, anodal_comp_rej, paired = T, hedges.correction = F)


# 2. ANOVA to probe RRC vs RCR -----------------------------------------------


cathodal_dif %>% 
  mutate(Paradigm = ifelse(Type == "Difficult" | 
                             Type == "Easy", "RCR","RRC")) -> cathodal_dif_exp3

cathodal_dif_exp3[, 'Choice_Paradigm'] <- paste0(cathodal_dif_exp3[, Choice], '_',
                                                 cathodal_dif_exp3[, Paradigm])

anodal_dif %>% 
  mutate(Paradigm = ifelse(Type == "Difficult" | 
                             Type == "Easy", "RCR","RRC")) -> anodal_dif_exp3

anodal_dif_exp3[, 'Choice_Paradigm'] <- paste0(anodal_dif_exp3[, Choice], '_',
                                               anodal_dif_exp3[, Paradigm])

  
library(ez)
library(rstatix)
detach(package:rstatix)
library(effectsize)

anova_paradigm_cathodal <- ezANOVA(data = cathodal_dif_exp3, dv = Difference, wid = Sub, 
                          within = c(Choice, Paradigm))  
anova_paradigm_cathodal


pairwise.t.test.with.t.and.df(cathodal_dif_exp3$Difference, cathodal_dif_exp3$Choice_Paradigm,
                              p.adjust.method = 'fdr', paired = T)

rstatix::pairwise_t_test(cathodal_dif_exp3, Difference ~ Choice_Paradigm, paired = TRUE)


anova_paradigm_anodal <- ezANOVA(data = anodal_dif_exp3, dv = Difference, wid = Sub, 
                          within = c(Choice, Paradigm))  
anova_paradigm_anodal


rstatix::pairwise_t_test(anodal_dif_exp3, Difference ~ Choice_Paradigm, paired = TRUE)


