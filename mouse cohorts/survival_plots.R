
#####overall survival



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
             size =2,
             alpha = .75,
             break.x.by = 25,
             break.y.by = .25,
             axes.offset = FALSE,
             legend = "right",
             ggtheme = theme_classic(),
             palette = c("locuszoom"),
             xlab = "Time Post Injection (Days)",
             legend.title = "Genes KO",
             legend.lab = leg_order
  )


#################################
#now filtered by tumor survival. 




#survival time after day tumor was noticed. (Aka tumor-free survival)
#this plot is not necessarily useful, since there was no intervention,
#also because most tumors were not visible prior to sac for end date. 


leg_order<-
  read_csv("survival_compiled.csv")%>%
  count(genes_ko)%>%
  select(1)%>%
  pull(genes_ko)



survfit(Surv(time = elapsed_from_tumor)~genes_ko, 
        data = read_csv("survival_compiled.csv")%>%
          filter(tumor_noticed != is.na(tumor_noticed)))%>%
  
  ggsurvplot(xlim = c(0, 150),
             ylim = c(0, 1.02),
             size =2,
             alpha = .75,
             break.x.by = 25,
             break.y.by = .25,
             axes.offset = FALSE,
             legend = "right",
             ggtheme = theme_classic(),
             palette = c("locuszoom"),
             xlab = "Time Post Injection (Days)",
             legend.title = "Genes KO",
             #legend.lab = leg_order
  )  
