library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# plot the base figures (output plot_base_stack, plot_base_fill)
Interpro_all <- read.delim("./Interpro_all.Rtab.combine", header = FALSE)
names(Interpro_all) <- c("Tag_name", "DB", "ex_confirmed", "reference")
#plot_base_stack <- ggplot(Interpro_all, aes(DB, fill=ex_confirmed))+geom_bar(position = "stack")
#plot_base_fill <- ggplot(Interpro_all, aes(DB, fill=ex_confirmed))+geom_bar(position = "fill")

# plot annotation modify
#plot_base_stack <- plot_base_stack+theme(axis.text.x = element_text(angle = 45))
#plot_base_fill <- plot_base_fill+theme(axis.text.x = element_text(angle = 45))

# analysis the raw data
#Interpro_all_summarized <- summarise(group_by(Interpro_all, DB, ex_confirmed), number= n(), percent=n()/nrow(DB))
Interpro_all_summarized_1 <- Interpro_all %>% group_by(DB, ex_confirmed, reference) %>% summarise( number =n())
#Interpro_all_summarized_stat <- summarise(group_by(Interpro_all_summarized_1, ex_confirmed), median_num=median(number))
#median_num_Yes = as.numeric(filter(Interpro_all_summarized_stat, ex_confirmed == "Yes")[1,2])
#median_num_No = as.numeric(filter(Interpro_all_summarized_stat, ex_confirmed == "No")[1,2])
Interpro_all_summarized_2 <- spread(Interpro_all_summarized_1, key = ex_confirmed, value = number)
Interpro_all_summarized_2 <- replace_na(Interpro_all_summarized_2, list(No=0, Yes=0))
Interpro_all_summarized_2 <- mutate(Interpro_all_summarized_2, ex_confirmed_rate = Yes/(Yes+No))
Interpro_all_summarized_2$reference <- factor(Interpro_all_summarized_2$reference, 
                                              levels = c('Clostridium','S.aureus','C.jejuni','Shigella','E.coli','Salmonella','Non-Salmonella','SetA'))
#median_num_Yes_rate = median(Interpro_all_summarized_2$ex_confirmed_rate)
# plot the base figures using analyzed data
#plot_base_stack <- ggplot(Interpro_all_summarized_1, aes(x= DB, y= number, fill= ex_confirmed, label = number))+geom_bar(stat = "identity")+
#  geom_text(size = 3, position = position_stack(vjust = 0.5))+
#  theme(axis.text.x = element_text(angle = 45))

plot_base_stack <- ggplot(Interpro_all_summarized_2, aes(x=DB, y= ex_confirmed_rate, fill=reference)) +
  geom_bar(stat = "identity")+
  #geom_text(aes(label = paste("number",as.character(Yes))), size = 3, position = position_stack(vjust = 0.5))+
  geom_text(aes(label = paste(label_percent()(ex_confirmed_rate),"(",as.character(Yes),")")), size = 3, position = position_stack(vjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45))

plot_base_fill <- ggplot(Interpro_all_summarized_2, aes(x=DB, y= ex_confirmed_rate, fill=reference)) +
  geom_bar(stat = "identity", position = "fill")+
  #geom_text(aes(label = paste("number",as.character(Yes))), size = 3, position = position_stack(vjust = 0.5))+
  #geom_text(aes(label = paste(label_percent()(ex_confirmed_rate),"(",as.character(Yes),")")), size = 3, position = position_fill(vjust = 0.65))+
  geom_text(aes(label = paste(label_percent()(ex_confirmed_rate))), size = 3, position = position_fill(vjust = 0.65))+
  geom_point(position = position_fill(vjust = 0.3), aes(size = Yes))+
  #scale_size_continuous(breaks = breaks_extended(n=9))+
  theme(axis.text.x = element_text(angle = 45))

plot_base_dodge <- ggplot(Interpro_all_summarized_2, aes(x=DB, y= ex_confirmed_rate, fill=reference)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5))+
  #geom_text(aes(label = paste("number",as.character(Yes))), size = 3, position = "dodge"+
  #geom_text(aes(label = paste(label_percent()(ex_confirmed_rate),"(",as.character(Yes),")")), size = 3, position = position_dodge(width = 0.9))+
  theme(axis.text.x = element_text(angle = 45))

plot_base_line <- ggplot(Interpro_all_summarized_2, aes(x=DB, y= ex_confirmed_rate, colour=reference, group=reference)) +
  geom_line(stat = "identity",)+
  theme(axis.text.x = element_text(angle = 45))

plot_base_fill_num <- ggplot(Interpro_all_summarized_2, aes(x=DB, y= Yes, fill=reference)) +
  geom_bar(stat = "identity", position = "fill")+
  #geom_text(aes(label = paste("number",as.character(Yes))), size = 3, position = position_stack(vjust = 0.5))+
  #geom_text(aes(label = paste(label_percent()(ex_confirmed_rate),"(",as.character(Yes),")")), size = 3, position = position_fill(vjust = 0.65))+
  geom_text(aes(label = paste(as.character(Yes))), size = 3, position = position_fill(vjust = 0.65))+
  geom_point(position = position_fill(vjust = 0.3), aes(size = ex_confirmed_rate))+
  #scale_size_continuous(breaks = breaks_extended(n=9))+
  theme(axis.text.x = element_text(angle = 45))

# add annotation lines
#plot_base_stack <- plot_base_stack+geom_hline(yintercept = median_num_Yes, linetype="dashed", color = "blue", size=0.4)
#plot_base_fill <- plot_base_fill+geom_hline(yintercept = median_num_Yes_rate, linetype="dashed", color = "blue", size=0.4)
#plot_base_fill <- plot_base_fill + geom_hline(yintercept = 167/4548, linetype="dashed", color = "red", size=0.4, show.legend = TRUE)

###manually
#plot_base_stack + coord_cartesian(ylim = c(0,100))
plot_base_stack <- plot_base_stack + scale_x_discrete(limits = c("GENE3D", "SUPERFAMILY", "PFAM", "SMART", "TIGRFAM", "PIRSF", "SFLD", "HAMAP", "PROSITE_PROFILES", "CDD", "PRINTS", "PROSITE_PATTERNS", "MOBIDB_LITE", "COILS"))
plot_base_fill <- plot_base_fill + scale_x_discrete(limits = c("GENE3D", "SUPERFAMILY", "PFAM", "SMART", "TIGRFAM", "PIRSF", "SFLD", "HAMAP", "PROSITE_PROFILES", "CDD", "PRINTS", "PROSITE_PATTERNS", "MOBIDB_LITE", "COILS"))
plot_base_fill_num <- plot_base_fill_num + scale_x_discrete(limits = c("GENE3D", "SUPERFAMILY", "PFAM", "SMART", "TIGRFAM", "PIRSF", "SFLD", "HAMAP", "PROSITE_PROFILES", "CDD", "PRINTS", "PROSITE_PATTERNS", "MOBIDB_LITE", "COILS"))

### post marker modifiecation
plot_base_fill <- plot_base_fill + labs(y= "Rate distribution of experimentally verified virulence factor", x = "Seed_database name")
plot_base_fill <- plot_base_fill + scale_fill_discrete("Reference_proteins name") + scale_size_continuous("Experimentally verified \nvirulence factor number", breaks = breaks_extended(n=9))

plot_base_fill_num <- plot_base_fill_num + labs(y= "Number distribution of experimentally verified virulence factor", x = "Seed_database name")
plot_base_fill_num <- plot_base_fill_num + scale_fill_discrete("Reference_proteins name") + scale_size_continuous("Experimentally verified \nvirulence factor rate", breaks = breaks_extended(n=9))
