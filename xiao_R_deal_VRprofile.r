## this is Xiao's first R script
## this script is used to visualize the VRprofile output data

#___1___this function get the input data virulence disturibtion
virulence_dist <- function(exact_directory, tag = NULL)
{
  library(readxl)
  VRprofile <- read_excel(exact_directory, col_names = FALSE, na = "-", skip = 2)
  library(dplyr)
  names(VRprofile) <- c( "1", "2", "gene_name", "length_aa",
                          "VF_hit", "6", "AR_hit", "8", 
                          "T3SE_hit", "10", "T4SE_hit", "12",
                          "T6SE_hit",  "14", "T7SE_hit", "16",
                          "prophage_hit",  "18", "ICE_hit", "20",
                          "T3SS_hit",  "22", "T4SS_hit", "24",
                          "T6SS_hit",  "26", "T7SS_hit", "28",
                          "Integron_hit",  "30", "IS_hit", "32",
                          "PAI_hit", "34", "ARI_hit", "36"
                         )
  
  wash_VRprofile <- select(VRprofile, c("gene_name", "length_aa",
                                        "VF_hit", "AR_hit", 
                                        "T3SE_hit", "T4SE_hit",
                                        "T6SE_hit", "T7SE_hit",
                                        "prophage_hit", "ICE_hit",
                                        "T3SS_hit", "T4SS_hit",
                                        "T6SS_hit", "T7SS_hit",
                                        "Integron_hit","IS_hit",
                                        "PAI_hit", "ARI_hit"
                                        )) 
  
  library(tidyr)
  #NA_VRprofile <- gather(wash_VRprofile, key = virulence_parts, value = hit_name, -gene_name, -length_aa)
  #NA_VRprofile <- NA_VRprofile[is.na(NA_VRprofile$hit_name),]
  #NA_VRprofile ["data_tag"]<- tag
  wash_VRprofile <- gather(wash_VRprofile, key = virulence_parts, value = hit_name, -gene_name, -length_aa, na.rm = TRUE)
  wash_VRprofile["data_tag"]  <- tag
  
  library(ggplot2)
  graph <- ggplot(wash_VRprofile, aes(virulence_parts, fill = virulence_parts)) +
    geom_bar() + theme(axis.text.x = element_text(angle = 45))
  ggsave("virulence_dist.pdf", dpi = 600)
  
 
  
  #View(wash_VRprofile)
  #View(VRprofile)
  return(wash_VRprofile)
}


#___2___
virulence_dist_group <- function(inputdata, baseline = NULL)
{
  library(ggplot2)
  inputdata_text <- table(inputdata$data_tag,inputdata$virulence_parts)
  View(inputdata_text)
  write.csv(inputdata_text, file = "virulence_dist_group.csv")
  
  graph <- ggplot(inputdata, aes(virulence_parts, fill = data_tag)) +
    geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 45)) +
    labs(y = "distribution rate")
    
  if (!is.null(baseline)) {
    graph <- graph + geom_hline(yintercept= baseline, linetype="dashed", color = "red")
  }
  ggsave("virulence_dist_group_fill.pdf", dpi = 600)
  
  graph <- ggplot(inputdata, aes(virulence_parts, fill = data_tag)) +
    geom_bar(position = "dodge") + theme(axis.text.x = element_text(angle = 45))
    
  ggsave("virulence_dist_group_dodge.pdf", dpi = 600)
  
  }
  
  
  #__3__ return the na part of dataframe
virulence_na_dataframe <- function(a1, a2, tag = NULL)
{
  a1["data_tag"] <- tag
  a2["data_tag"] <- tag
  a1["length_aa"] <- 999
  a2["length_aa"] <- 999
  a1["hit_name"] <- "hit"
  a2["hit_name"] <- "hit"
  
  #View(a1)
  washed_set <- setdiff(a1,a2)
  return(washed_set)
}
  
  
  
  
  
  
  
