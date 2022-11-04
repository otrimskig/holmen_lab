
##############################################

#main

home<-"X:/Holmen Lab/Lab Personnel/Garrett/R/"
setwd(home)



library(tidyverse)
library(ggplot2)
library(janitor)
library(fuzzyjoin)
library(tools)
library(readr)
library(fs)
library(stringr)
library(googlesheets4)

##############################################


##colors and other packages for plots


library(viridis)
library(ggsci)
library(RColorBrewer)
library(survival)



#############################################


install.packages("clipr")
install.packages("googlesheets4")




install.packages("tidyverse")
library(tidyverse)

install.packages("ggplot2")
library(ggplot2)

install.packages("janitor")
library(janitor)

install.packages("survival")
library(survival)

install.packages("RColorBrewer")
library(RColorBrewer)

install.packages("viridis")
library(viridis)

install.packages("ggsci")
library(ggsci)

install.packages("fuzzyjoin")
library(fuzzyjoin)

library(tools)
library(readr)

install.packages("fs")
library(fs)







coalesce_join <- function(x, y, 
                          by = NULL, suffix = c(".x", ".y"), 
                          join = dplyr::full_join, ...) {
  joined <- join(x, y, by = by, suffix = suffix, ...)
  # names of desired output
  cols <- union(names(x), names(y))
  
  to_coalesce <- names(joined)[!names(joined) %in% cols]
  suffix_used <- suffix[ifelse(endsWith(to_coalesce, suffix[1]), 1, 2)]
  # remove suffixes and deduplicate
  to_coalesce <- unique(substr(
    to_coalesce, 
    1, 
    nchar(to_coalesce) - nchar(suffix_used)
  ))
  
  coalesced <- purrr::map_dfc(to_coalesce, ~dplyr::coalesce(
    joined[[paste0(.x, suffix[1])]], 
    joined[[paste0(.x, suffix[2])]]
  ))
  names(coalesced) <- to_coalesce
  
  dplyr::bind_cols(joined, coalesced)[cols]
}



