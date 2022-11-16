
#####dataset cleaning: 1st, make accessible dataset with all relevant information,
#pulled from main dataset. 
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
#############################################################
#trying to create a yes/no of had tumor or not. 
#select.csv is a list of any variables that could help tell me if there was a tumor present. 

#extract list of variable names from select.csv, to use for selection from main dataset. 
cols<-
  read_csv("select.csv", show_col_types = FALSE)%>%
  clean_names()%>%
  filter(x1 == "1")%>%
  select(x3)%>%
  pull()



#use variable names for below transposition of dataframe. 

tt<-
read_csv("compiled_cohorts2.csv", show_col_types = FALSE)%>%
  select((cols))%>%
  mutate(had_tumor = 0)%>%
  relocate(had_tumor, .before = tumor_id)%>%
  as.data.frame()
 
#transpose and keep columns as first row. Will use to check how to use each variable.
setNames(data.frame(t(tt[,-1])), tt[,1])%>%

write.csv("tumor_test.csv")


#used output of above to look through data variables. 
#made notes on conditions to look for for tumor determination. 




#########

#continue: looking to create a variable for whether mice had any tumors, at any time. 
#reading in cleaned cohorts data, selecting for any variables that can be informative. 

#greyed out were used to test dataset conditions. 
#once code was verified on subset data, easier to visualise, 
#code then used to perform operations on full compiled cohorts data. 

read_csv("compiled_cohorts2.csv", show_col_types = FALSE)%>%
  #subset db
  #select(c((cols), exclude))%>%
  #create had_tumor col.
  mutate(had_tumor = NA)%>%
  relocate(had_tumor)%>%
  
  #some low-hanging fruit. 
  mutate(had_tumor = case_when(is.na(tumor_id)==FALSE|
                                 is.na(tumor)==FALSE|
                                 is.na(tumor_noticed)==FALSE ~ "1"))%>%
  
  #case_when is having issues with complex self-referential when nested in one
  #mutate command. Separate out if there are multiple complex and/or statements. 
  
  #first "complex" conditional. Entry can exist but may also say that there is no tumor. 
  #sorts for any entries that contain values, but NOT specific phrases denoting no tumors. 
  mutate(had_tumor = case_when(is.na(tumor_locations)==FALSE&
                                 grepl("Pre-tumorigenic", tumor_locations)==FALSE&
                                 grepl("no obs", tumor_locations)==FALSE ~ "1",
                                TRUE ~ as.character(had_tumor)))%>%
  
  #make new column that indicates if mouse was pretumorigenic. 
  #may be useful later, if I want to combine both tumor+pretumorigenic groups. 
  mutate(pretumorigenic = if_else(grepl("Pre-tumorigenic", tumor_locations), 1, NA_real_))%>%
  

  #####have to use as.character(x) in order for the false condition to ignore NAs
  #mark 0 for tumor if was pretumorigenic. 
  #leave remaining alone. 
  mutate(had_tumor = case_when(pretumorigenic=="1" ~ "0", 
                                 TRUE ~ as.character(had_tumor)))%>%
  
  mutate(had_tumor = case_when(notes == "No tumor. Sac'd for end date." ~ "0",
                               TRUE ~ as.character(had_tumor)))%>%
  
  
  mutate(had_tumor = case_when(tumor_locations == "no obs" ~ "0",
                               TRUE ~ as.character(had_tumor)))%>%
  
  mutate(had_tumor = case_when(notes == "no evidence of tumor" ~ "0",
                               TRUE ~ as.character(had_tumor)))%>%
  
  mutate(had_tumor = case_when(notes == "no evidence of tumor"|
                                notes == "slightly hunched back but no evidence of tumor"|
                                notes == "ruffle fur noticed, but no obs. Tumor. found dead"|
                                notes == "Lost weight, extreme hunchback, head twitching, slow moving. Soft skull and brain. No obs. Tumor, no other CNS tumor found"|
                                notes == "Lost weight, mild hunchback, slow moving. Soft skull and brain. No obs. Tumor, no other CNS tumor found. Enlarged optic nerve" ~ "0",
                               
                               TRUE ~ as.character(had_tumor)))%>%
  
  mutate(had_tumor = case_when(necropsy_behavior_comments == "end of experiment, looked normal" ~ "0",
                               TRUE ~ as.character(had_tumor)))%>%

  mutate(had_tumor = case_when(proof_of_mass == "yes" ~ "1",
                               proof_of_mass == "no" ~"0",
                               TRUE ~ as.character(had_tumor)))%>%
  
  
  mutate(had_tumor = case_when(is.na(tumor_noticed) == FALSE ~ "1",
                               TRUE ~ as.character(had_tumor)))%>%
  
  write_csv("compiled_cohorts3.csv")
  
  #relocate(exclude)%>%
  #filter(is.na(exclude))%>%
  
  #filter(is.na(had_tumor))%>%
  #select(genes_ko)%>%
  #count(genes_ko)%>%
  

##################################################
#genes_ko'd groups
#"heatmap" of genes ko's for easy visualisation of groups. 

read_csv("compiled_cohorts3.csv")%>%
  group_by(genes_ko)%>%
  slice(1)%>%
  select(genes_ko, nf1_ko, pten_ko, ink_ko, atrx_ko)


kos<-
read_csv("kos.csv", show_col_types = FALSE)%>%
  drop_na()%>%
  row_to_names(1)%>%clean_names()%>%
  bind_cols(read_csv("compiled_cohorts3.csv", show_col_types = FALSE)%>%
              group_by(genes_ko)%>%
              slice(1)%>%
              select(genes_ko, nf1_ko, pten_ko, ink_ko, atrx_ko))%>%
  select(1:5)%>%
  relocate(genes_ko)%>%
  mutate_at(vars(2:5), as_factor)

kos%>%
  pivot_longer(2:5, names_to = "gene", values_to = "ko")%>%


ggplot(aes(x = factor(gene, levels = c("nf1", "pten", "ink", "atrx")), genes_ko, fill = ko))+
  geom_point(shape = 22, size = 8)+
  facet_grid(.~factor(gene, levels = c("nf1", "pten", "ink", "atrx")), scales = "free")+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.text.x.bottom = element_blank(),
        legend.position = "none",
        panel.spacing = unit(-2, "lines"),
        axis.line = element_blank(),
        panel.grid = element_blank())+
  scale_fill_manual(values = c("#00FF81", "#ffffff"))
  






###############################
#######################################
#creating common color palette for all plots. 
library(paletteer)
library(ggsci)

lz<-
paletteer_d("ggsci::default_locuszoom", n=6)
lz
rev(lz)

custom_lz<-
c("#9632B8FF",  "#357EBDFF", "#5CB85CFF",  "#EEA236FF", "#D43F3AFF", "#46B8DAFF")


leg_order%>%
  as_tibble()%>%
  mutate(rgb = custom_lz)%>%
  mutate(nickname = c("purple", "dk blue", "green", "yellow", "red", "light blue"))%>%
  rename(genes_ko = value)%>%
  write_csv("cohort_colors.csv")


df_col<-
  read_csv("cohort_colors.csv")%>%
          select(rgb, genes_ko)

color_hex<-
df_col%>%
  pull(rgb)

color_genes<-
  df_col%>%
  pull(genes_ko)

names(color_hex) = color_genes
color_hex





