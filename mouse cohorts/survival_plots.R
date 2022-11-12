
##############plotting survival curve with full cohort data. 
library(tidyverse)
library(janitor)
#1st thing, creating a clean dataset that can be easily used for plotting. 
#pull data from compiled_cohorts2.csv, which has all cleaned and unified data.
read_csv("compiled_cohorts2.csv", show_col_types = FALSE)%>%
  
  #selecting variables that are needed for survival curve plots, 
  #plus identifiers if needed to link back to main data. 
  select_at(vars(mouse_num,
                 dob,
                 virus,
                 injection_date,
                 death_date,
                 tumor_noticed,
                 behavior_noticed,
                 strain,
                 src,
                 age_death_capped,
                 elapsed_from_tumor,
                 elapsed_from_behavior,
                 exclude, 
                 genes_ko))%>%
  
  group_by(mouse_num)%>%
  arrange(desc(death_date))%>%
  slice(1)%>%
  ungroup()%>%
  
  filter(is.na(exclude))%>%
  select(-exclude)%>%
  relocate(virus, .after = strain)%>%
  
  write_csv("survival_compiled.csv")

##############################################################
#trying to create a yes/no of had tumor or not. 
#select.csv is a list of any variables that could help tell me if there was a tumor present. 

#extract list of variable names from select.csv, to use for selection from main dataset. 
cols<-
  read_csv("select.csv", show_col_types = FALSE)%>%
  clean_names()%>%
  filter(x1 == "1")%>%
  select(x3)%>%
  pull()




tt<-
read_csv("compiled_cohorts2.csv", show_col_types = FALSE)%>%
  select((cols))%>%
  mutate(had_tumor = 0)%>%
  relocate(had_tumor, .before = tumor_id)%>%
  as.data.frame()
 
#transpose and keep columns as first row. Will use to check how to use each variable.
setNames(data.frame(t(tt[,-1])), tt[,1])%>%

write.csv("tumor_test.csv")

###################################################################
real1<-
read_csv("compiled_cohorts2.csv", show_col_types = FALSE)%>%
  #subset db
  select((cols))%>%
  #create had_tumor col.
  mutate(had_tumor = NA)%>%
  relocate(had_tumor)%>%
  mutate(had_tumor = case_when(is.na(tumor_id)==FALSE|
                                 is.na(tumor)==FALSE|
                                 is.na(tumor_noticed)==FALSE ~ 1))%>%
  
  #case_when has issues with complex self-referential when nested in one
  #mutate command. Separate out if there are multiple complex and/or statements. 
  
  mutate(had_tumor = case_when(is.na(tumor_locations)==FALSE&
                                 grepl("Pre-tumorigenic", tumor_locations)==FALSE&
                                 grepl("no obs", tumor_locations)==FALSE ~ 1))%>%
  
  mutate(pretumorigenic = if_else(grepl("Pre-tumorigenic", tumor_locations), 1, NA_real_))
  
  real1%>%
    mutate(had_tumor = case_when(pretumorigenic=="1" ~ "0", 
                                 TRUE ~ as.character(had_tumor)))%>%
  
  view()
  








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
