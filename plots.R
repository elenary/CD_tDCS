# main plot for cathodal

cathodal_dif_plot %>% 
  filter(Type != "S_Post-Ex") %>%
  filter(Type != "R_Post-Ex") %>%
  # filter(!(Type == "S_Easy" & Choice == "Selected")) %>% 
  # filter(!(Type == "S_FComputer" & Choice == "Selected")) %>% 
  mutate(Stimulation = ifelse(Stimulation == "sham", "sham", "cathodal tDCS")) %>% 
  ggplot(aes(x=Type,y = Difference, fill = Stimulation)) +
  geom_boxplot(outlier.size=-1) +
  scale_y_continuous(breaks = seq(-2.5, 1.5, 0.25), limits=c(-2.5, 1.5)) +
  theme_minimal() +
  scale_fill_manual(values = c("lightseagreen", "gray80")) +
  # scale_x_discrete(labels = function(x) str_wrap(c("self-difficult", "self-easy", 
  #                                                  "computer", "self-difficult"), 
  #                                                width = 100)) +
  scale_x_discrete(labels = function(x) str_wrap(c("self-difficult", "self-easy",
                                                   "computer", "self-difficult", "self-easy",
                                                   "computer"),
                                                 width = 100)) +
  ggtitle("Influence of cathodal tDCS on preference changes") +
  xlab("Choice type") + 
  ylab("Preference changes (points)") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, margin = margin(t = 0, r = 00, b = 15, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 00, b = 0, l = 0)),
        axis.text = element_text(size = 12))

  


cathodal_dif_plot %>% 
  filter(Type != "S_Post-Ex") %>% 
  filter(Type != "R_Post-Ex") %>% 
  filter(!(Type == "S_Easy" & Choice == "Selected")) %>% 
  filter(!(Type == "S_FComputer" & Choice == "Selected")) %>% 
  ggplot(aes(x=Type,y = Difference, color = Stimulation)) +
  stat_summary(fun = mean, geom = 'point', 
               position = position_dodge(0.5)) +
  stat_summary(fun.data = mean_se, geom = 'errorbar', 
               position = position_dodge(0.5), width = 0.4) +
  scale_y_continuous(breaks = seq(-2.5, 1.5, 0.25), limits=c(-2.5, 1.5)) +
  theme_minimal() +
  scale_color_manual(values = c("gray80", "lightseagreen"))

# subjects 1-9

cathodal_dif_plot %>% 
  filter(Type != "S_Post-Ex") %>% 
  filter(Type != "R_Post-Ex") %>% 
  # filter(!(Type == "S_Easy" & Choice == "Selected")) %>% 
  # filter(!(Type == "S_FComputer" & Choice == "Selected")) %>% 
  filter(Sub<=9) %>% 
  ggplot(aes(x=Type,y = Difference, color = Stimulation)) +
  stat_summary(fun = mean, geom = 'point', shape=19, size = 3,
               position = position_dodge(0.5)) +
  scale_y_continuous(breaks = seq(-2, 1, 0.25), limits=c(-2, 1)) +
  theme_minimal() +
  scale_color_manual(values = c("gray", "lightseagreen")) +
  facet_wrap(~Sub, ncol = 3)

# subjects 9-17

cathodal_dif_plot %>% 
  filter(Type != "S_Post-Ex") %>% 
  filter(Type != "R_Post-Ex") %>% 
  # filter(!(Type == "S_Easy" & Choice == "Selected")) %>% 
  # filter(!(Type == "S_FComputer" & Choice == "Selected")) %>% 
  filter(Sub > 9) %>% 
  ggplot(aes(x=Type,y = Difference, color = Stimulation)) +
  stat_summary(fun = mean, geom = 'point', shape=19, size = 3,
               position = position_dodge(0.5)) +
  scale_y_continuous(breaks = seq(-2, 1, 0.25), limits=c(-2, 1)) +
  theme_minimal() +
  scale_color_manual(values = c("gray", "lightseagreen")) +
  facet_wrap(~Sub, ncol = 3)

# final plot with all conditions

anodal_dif_plot <- anodal_dif

anodal_dif_plot[Type == 'Computer', 'Type'] <-
  paste0('F', anodal_dif_plot[Type == 'Computer', Type])

anodal_dif_plot[Choice == 'Selected', 'Type'] <-
  paste0('S_', anodal_dif_plot[Choice == 'Selected', Type])
anodal_dif_plot[Choice == 'Rejected', 'Type'] <-
  paste0('R_', anodal_dif_plot[Choice == 'Rejected', Type])

# spaghetti plot cathodal

cathodal_dif_plot %>% 
  filter(Type != "S_Post-Ex") %>% 
  filter(Type != "R_Post-Ex") %>% 
  ggplot(aes(x=Stimulation,y = Difference, group = as.factor(Sub))) +
  geom_line() +
  scale_y_continuous(breaks = seq(-2.5, 1.5, 0.25), limits=c(-2.5, 1.5)) +
  theme_minimal()

# main plot for anodal

anodal_dif_plot %>% 
  filter(Type != "S_Post-Ex") %>% 
  filter(Type != "R_Post-Ex") %>% 
  filter(!(Type == "S_Easy" & Choice == "Selected")) %>% 
  filter(!(Type == "S_FComputer" & Choice == "Selected")) %>% 
  mutate(Stimulation = ifelse(Stimulation == "sham", "sham", "anodal tDCS")) %>% 
  ggplot(aes(x=Type,y = Difference, fill = Stimulation)) +
  geom_boxplot(outlier.size=-1) +
  scale_y_continuous(breaks = seq(-2.5, 1.5, 0.25), limits=c(-2.5, 1.5)) +
  theme_minimal() +
  scale_fill_manual(values = c("indianred3", "gray80")) +
  scale_x_discrete(labels = function(x) str_wrap(c("self-difficult", "self-easy", 
                                                   "computer", "self-difficult"), 
                                                 width = 100)) +
  ggtitle("Influence of anodal tDCS on preference changes") +
  xlab("Choice type") + 
  ylab("Preference changes (points)") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, margin = margin(t = 0, r = 00, b = 15, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 00, b = 0, l = 0)),
        axis.text = element_text(size = 12))

anodal_dif_plot %>% 
  filter(Type != "S_Post-Ex") %>% 
  filter(Type != "R_Post-Ex") %>% 
  filter(!(Type == "S_Easy" & Choice == "Selected")) %>% 
  filter(!(Type == "S_FComputer" & Choice == "Selected")) %>% 
  ggplot(aes(x=Type,y = Difference, fill = Stimulation)) +
  geom_boxplot(outlier.size=-1) +
  scale_y_continuous(breaks = seq(-2.5, 1.5, 0.25), limits=c(-2.5, 1.5)) +
  theme_minimal() +
  scale_fill_manual(values = c("gray80", "indianred3")) 

# spaghetti plot anodal

anodal_dif_plot %>% 
  filter(Type != "S_Post-Ex") %>% 
  filter(Type != "R_Post-Ex") %>% 
  ggplot(aes(x=Stimulation,y = Difference, group = as.factor(Sub))) +
  geom_line() +
  scale_y_continuous(breaks = seq(-2.5, 1.5, 0.25), limits=c(-2.5, 1.5)) +
  theme_minimal()

# subjects 1-9

View(anodal_dif_plot)
anodal_dif_plot %>% 
  filter(Type != "S_Post-Ex") %>% 
  filter(Type != "R_Post-Ex") %>% 
  filter(!(Type == "S_Easy" & Choice == "Selected")) %>% 
  filter(!(Type == "S_FComputer" & Choice == "Selected")) %>% 
  filter(Sub<=12) %>% 
  ggplot(aes(x=Type,y = Difference, color = Stimulation)) +
  stat_summary(fun = mean, geom = 'point', shape=19, size = 3,
               position = position_dodge(0.5)) +
  scale_y_continuous(breaks = seq(-2, 1, 0.25), limits=c(-2, 1)) +
  theme_minimal() +
  scale_color_manual(values = c("gray", "indianred3")) +
  facet_wrap(~Sub, ncol = 3)


# subjects 1-9

View(anodal_dif_plot)
anodal_dif_plot %>% 
  filter(Type != "S_Post-Ex") %>% 
  filter(Type != "R_Post-Ex") %>% 
  filter(!(Type == "S_Easy" & Choice == "Selected")) %>% 
  filter(!(Type == "S_FComputer" & Choice == "Selected")) %>% 
  filter(Sub>12) %>% 
  ggplot(aes(x=Type,y = Difference, color = Stimulation)) +
  stat_summary(fun = mean, geom = 'point', shape=19, size = 3,
               position = position_dodge(0.5)) +
  scale_y_continuous(breaks = seq(-2, 1, 0.25), limits=c(-2, 1)) +
  theme_minimal() +
  scale_color_manual(values = c("gray", "indianred3")) +
  facet_wrap(~Sub, ncol = 3)
