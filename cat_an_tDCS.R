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

# locale

Sys.setlocale(category = 'LC_ALL', locale = "RU_ru")
Sys.getlocale()

# color

display.brewer.all()
viridis(2)
scale_color_viridis()


#----cathodal df--------------------------------------------------------

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

#Post-Ex is processed separatly, remove it
# cathodal_dif_short <- subset(cathodal_dif, Type != "Post-Ex")

#----anodal df--------------------------------------------------------

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

#Post-Ex is processed separatly, remove it
# anodal_dif_short <- subset(anodal_dif, Type != "Post-Ex")

length_cat <- length(unique(cathodal_dif$Sub, incomparables = FALSE))
length_an <- length(unique(anodal_dif$Sub, incomparables = FALSE))

#----   function for extracting t values and DFs-----------------

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








#----   CATHODAL GAIN SCORE ANOVA for difference    -----main(1)-----

#ANOVA for cathodal

ANOVA_cathodal_gs <- ezANOVA(data = cathodal_dif, dv = Difference, wid = Sub,
                             within = .(Stimulation, Type, Choice), return_aov = T, type = 3)
ANOVA_cathodal_gs
summary(ANOVA_cathodal_gs$aov)

# effct size
ANOVA_cathodal_gs

#check residuals on normality

ANOVA_cathodal_gs_resid <-
  proj(ANOVA_cathodal_gs$aov)[[3]][, "Residuals"]
attr(ANOVA_cathodal_gs_resid, "ATT") <- NULL


ANOVA_cathodal_gs_resid
shapiro.test(ANOVA_cathodal_gs_resid)
hist(ANOVA_cathodal_gs_resid)

#----   without post-ex----
# ANOVA_cathodal_gs2 <- ezANOVA(data = cathodal_dif_not_postex, dv = Difference, wid = Sub,
#                              within = .(Stimulation, Type, Choice), return_aov = T, type = 3)
# ANOVA_cathodal_gs2
# summary(ANOVA_cathodal_gs2$aov)


#is there difference between type 3 and 2? No

# ANOVA_cathodal_gs <- ezANOVA(data = cathodal_dif, dv = Difference, wid = Sub,
#                              within = .(Stimulation, Type, Choice), return_aov = T, type = 2)
# ANOVA_cathodal_gs
# summary(ANOVA_cathodal_gs$aov)



#----   POST HOC FOR GS CATHODAL ANOVA    --------------------------------

#----   trying one more long format----

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



#----   post hoc interaction Type:Choice (it's significant)   ------------------------------------main(2)-----

cathodal_dif_interaction <- cathodal_dif
cathodal_dif_interaction[, 'Type_Choice'] <-
  paste0(cathodal_dif_interaction[, Type], '_', cathodal_dif_interaction[, Choice])

View(cathodal_dif_interaction)

pwtt_interaction_fdr_gs <- pairwise.t.test(
  cathodal_dif_interaction$Difference,
  g = cathodal_dif_interaction$Type_Choice,
  p.adjust.method = 'fdr',
  paired = T
)

View(as.data.table(pwtt_interaction_fdr_gs[3], keep.rownames=T))



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
pairwTTest_fdr_gs_withtvalues$dfs



#----   post hoc interaction Stimulation:Type:Choice (it's not significant but interesting)--------main(3.1)---------


cathodal_dif_interaction2 <- cathodal_dif
cathodal_dif_interaction2[, 'Stimulation_Type_Choice'] <-
  paste0(cathodal_dif_interaction2[, Stimulation], '_', cathodal_dif_interaction2[, Type], '_', cathodal_dif_interaction2[, Choice])

View(cathodal_dif_interaction2)

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

#--   false discovery rate (fdr) correction-------------------------------------------------------------main(3.2)------

pwtt_interaction2_fdr <- pairwise.t.test(
  cathodal_dif_interaction2$Difference,
  g = cathodal_dif_interaction2$Stimulation_Type_Choice,
  p.adjust.method = 'fdr',
  paired = T
)

pwtt_interaction2_pv_fdr <- as.data.table(pwtt_interaction2_fdr[3], keep.rownames=T)
View(pwtt_interaction2_pv_fdr)

pwtt_interaction2_pv_fdr_tvalues <- pairwise.t.test.with.t.and.df(x = cathodal_dif_interaction2$Difference, 
                                                                   g = cathodal_dif_interaction2$Stimulation_Type_Choice, 
                                                                   p.adjust.method = 'fdr', paired = T) 
str(pwtt_interaction2_pv_fdr_tvalues)
View(pwtt_interaction2_pv_fdr_tvalues$t.value)


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


#----   t tests for difference in difficult rejected for checking--------

#t test  for cathodal

cathodal_dif_sham_rej <- cathodal_dif[Type == "Difficult" & Stimulation == "sham" & Choice == 'Rejected', Difference]
cathodal_dif_tDCS_rej <- cathodal_dif[Type == "Difficult" & Stimulation == "tDCS_cat" & Choice == 'Rejected', Difference]



cohensD <- cohensD(cathodal_dif_sham_rej, cathodal_dif_tDCS_rej)
cohensD(cathodal_dif_sham_rej, cathodal_dif_tDCS_rej)
sdpool <- sqrt( ( sd(cathodal_dif_sham_rej)^2 + sd(cathodal_dif_tDCS_rej)^2 ) /2 )
sqrt( ( sd(cathodal_dif_sham_rej)^2 + sd(cathodal_dif_tDCS_rej)^2 ) /2 )


t.test(x = cathodal_dif_sham_rej,  
       y = cathodal_dif_tDCS_rej,  
       paired = T, var.equal = FALSE,
       conf.level = 0.95, alternative = "less")

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

#PLOT WITH p VALUES AND NUMBER OF SUBJECTS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#don't remember why



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



#----   CATHODAL REP MEASURMENT ANOVA by different choice factor (as Marco did)-----------------------------------------

cathodal_difficult <- cathodal[Type == 'Difficult', ]

ezANOVA(data = cathodal[Type == 'Difficult', ], dv = Evaluation, wid = Sub,
        within = .(Stimulation, Choice, Rating), return_aov = T, type = 3)

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


cathodal_interaction_rm <- cathodal
cathodal_interaction_rm[, 'Stimulation_Type_Choice_Rating'] <-
  paste0(cathodal_interaction_rm[, Stimulation], '_', cathodal_interaction_rm[, Type], '_', cathodal_interaction_rm[, Choice], '_', cathodal_interaction_rm[, Rating])






#----   CATHODAL GAIN SCORE ANOVA by different choice factor---------------------------

# for checking is gain score ANOVA differ from the initial RM Marco by factors

ANOVA_cathodal_gs_dif <- ezANOVA(data = cathodal_dif[Type == 'Difficult', ], dv = Difference, wid = Sub,
                                 within = .(Stimulation, Choice), return_aov = T, type = 3)
ANOVA_cathodal_gs_dif

#interaction for only difficult choices is significant

ANOVA_cathodal_gs_easy <- ezANOVA(data = cathodal_dif[Type == 'Easy', ], dv = Difference, wid = Sub,
                                 within = .(Stimulation, Choice), return_aov = T, type = 3)
ANOVA_cathodal_gs_easy

ANOVA_cathodal_gs_computer <- ezANOVA(data = cathodal_dif[Type == 'Computer', ], dv = Difference, wid = Sub,
                                  within = .(Stimulation, Choice), return_aov = T, type = 3)
ANOVA_cathodal_gs_computer


ANOVA_cathodal_gs_PostEx <- ezANOVA(data = cathodal_dif[Type == 'Post-Ex', ], dv = Difference, wid = Sub,
                                      within = .(Stimulation, Choice), return_aov = T, type = 3)
ANOVA_cathodal_gs_PostEx



#----   plots for cathodal------------------------------------------main(4)-------


#--plot for anova for visualization searate factors

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



#----
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



#----   plots anodal------------------

anodal_dif_plot <- anodal_dif

anodal_dif_plot[Type == 'Computer', 'Type'] <- 
  paste0('F', anodal_dif_plot[Type == 'Computer', Type]) 

anodal_dif_plot[Choice == 'Selected', 'Type'] <- 
  paste0('S_', anodal_dif_plot[Choice == 'Selected', Type])  
anodal_dif_plot[Choice == 'Rejected', 'Type'] <- 
  paste0('R_', anodal_dif_plot[Choice == 'Rejected', Type]) 


ggplot(data = anodal_dif_plot, aes(x=Type,y = Difference, fill = Stimulation)) +
  geom_boxplot() +
  theme_minimal() +
  scale_y_continuous()
