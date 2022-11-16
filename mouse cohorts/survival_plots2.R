#overall survival#####################################


leg_order<-
  read_csv("survival_compiled.csv")%>%
  count(genes_ko)%>%
  select(1)%>%
  pull(genes_ko)

library(survminer)
library(survival)

survfit(Surv(time = age_death_capped)~genes_ko, data = read_csv("survival_compiled.csv"))%>%
  
  ggsurvplot(xlim = c(0, 150),
             ylim = c(0, 1.02),
             size =3,
             alpha = .9,
             break.x.by = 25,
             break.y.by = .25,
             
             axes.offset = FALSE,
             legend = "right",
             ggtheme = theme_classic(),
             palette = custom_lz,
             xlab = "Time Post Injection (Days)",
             legend.title = "Genes KO",
             legend.lab = leg_order
  )


#################################
#now filtered by tumor survival. 




#survival time after day tumor was noticed. (Aka tumor-free survival)
#this plot is not necessarily useful, since there was no intervention,
#also because most tumors were not visible prior to sac for end date. 


library(ggsci)
library(survival)
library(survminer)


leg_order<-
  read_csv("survival_compiled.csv")%>%
  count(genes_ko)%>%
  select(1)%>%
  pull(genes_ko)

d<-
survfit(Surv(time = elapsed_from_tumor)~genes_ko, 
        data = read_csv("survival_compiled.csv")%>%
          mutate(elapsed_from_tumor = as.numeric(elapsed_from_tumor)))#%>%
          #filter(tumor_noticed != is.na(tumor_noticed)))%>%
  
ggsurvplot(d, xlim = c(0, 60),
             ylim = c(0, 1.02),
             size =3,
             alpha = .9,
             break.x.by = 10,
             break.y.by = .25,
             axes.offset = FALSE,
             legend = "right",
             ggtheme = theme_classic(),
              palette = custom_lz,
             xlab = "Time After Tumor Observed (Days)",
             legend.title = "Genes KO")  


custom_lz

color_hex

#####################################
library(viridis)

tumor_x_order<-
  read_csv("compiled_cohorts3.csv")%>%
  group_by(genes_ko)%>%
  count(had_tumor)%>%
  drop_na()%>%
  mutate(had_tumor = as.logical(had_tumor))%>%
  mutate(per = n/sum(n)*100)%>%
  ungroup()%>%
  group_by(had_tumor)%>%
  arrange(per)%>%
  filter(had_tumor == FALSE)%>%
  arrange(per)%>%
  pull(genes_ko)


tumor_color_order<-
  read_csv("compiled_cohorts3.csv")%>%
  group_by(genes_ko)%>%
  count(had_tumor)%>%
  drop_na()%>%
  mutate(had_tumor = as.logical(had_tumor))%>%
  mutate(per = n/sum(n)*100)%>%
  filter(had_tumor == FALSE)%>%
  inner_join(read_csv("cohort_colors.csv"))%>%
  pull(rgb)


read_csv("compiled_cohorts3.csv")%>%
  group_by(genes_ko)%>%
  count(had_tumor)%>%
  drop_na()%>%
  mutate(had_tumor = as.logical(had_tumor))%>%
  mutate(per = n/sum(n)*100)%>%
  filter(had_tumor == FALSE)%>%
  
  
  ggplot(aes((100-per), factor(genes_ko, level = rev(tumor_x_order)), color = genes_ko))+
  geom_segment(aes(x = -0, 
                   xend = 100-per, 
                   y = factor(genes_ko, level = rev(tumor_x_order)), 
                   yend = factor(genes_ko, level = rev(tumor_x_order))), 
               size = 13)+
  scale_color_manual(values = color_hex)+
  geom_segment(aes(x = 100-per,
                   xend = 100,
                   y = factor(genes_ko, level = rev(tumor_x_order)),
                   yend = factor(genes_ko, level = rev(tumor_x_order))),
               color = "grey", size = 3, alpha = 1)+
  geom_segment(aes(x = -1,
                   xend = -.7,
                   y = factor(genes_ko, level = rev(tumor_x_order)),
                   yend = factor(genes_ko, level = rev(tumor_x_order))),
               size = 13)+
  

  #geom_point(size = 5, color = "black", alpha = 1)+
  
xlim(0,100)+
  theme_classic()+
  # theme(panel.grid.major.x = element_blank(),
  #   panel.border = element_blank(),
  #   axis.ticks.x = element_blank(),
  #   )+
  scale_x_continuous(expand=(c(0,1)))



##########################
