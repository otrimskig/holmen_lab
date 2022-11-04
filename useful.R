
################          USEFUL FUNCTIONS      ##############################
#hit control + shift + C to comment-out entire highlighted section of code (and un-comment it)
################################################################################

# #standard csv read-in (remove 1st column for file with row numbers)
read.csv("ki67_joined.csv")%>%
  as_tibble()%>%
  select(-1)




#standard csv write with no rownames. 
write.csv("filename.csv",
          rownames = FALSE)


#############################################################
######df manipulation

#move row to column names, uses janitor package. takes number for row number.
#clean_names, also from janitor. makes titles R-friendly.
df%>%
  row_to_names(1)%>%clean_names()



###########################################################
#####plotting


#re-order x-axis:
            #vector with desired order of x-axis
x_order <- c("pre_tumorigenic", "glioma","spinal_cord", "head", "nerve")
              #by setting "level" with above vector within factor() function in x aesthetic. 
ggplot(aes(x = factor(tumor_location, level = x_order), 
            y = this))
  




ggplot(aes(timepoint, avg_conc, color = group)) +
  
  geom_point(shape = 1, size = 2, stroke = 1.5, alpha = .6) +
  geom_errorbar(aes(ymin = avg_conc + sd, ymax = avg_conc - sd), width = 0.2, alpha = .2)+ #position = position_dodge(width = 1)) +
  
  #geom_ribbon(aes(fill = group, ymin = avg_conc + sd, ymax = avg_conc - sd), linetype = 0, alpha = .2) + 
  #geom_line(alpha = .5) +
  theme_classic()


#shape 1 allows border. Shape 16 can be fill. 
#use position = position_dodge(width = .9)) to dodge points.


mutate(`injection_date` = as.Date(`injection_date`, format = "%m/%d/%Y"))

  