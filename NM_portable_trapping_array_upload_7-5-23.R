# 'A portable and wind resistant drift fence array for arid environments'
#submitted to Herpetology Notes 7/8/23

#BT Camper

library('dplyr') #version 1.0.10
library('tidyr') #version 1.2.1
library('ggplot2') #version 3.3.6
library('vegan') #version 2.6-2

#set working directory
setwd("insert working directory here")

#read in 2022 full data
reptiles_2022_full <- read.csv("NM_reptiles_2022_full_data.csv", header=TRUE)


#---

### ABUNDANCE AND RICHNESS BY SITE

#generates species richness by site
richness_by_site <- aggregate(reptiles_2022_full, Species~Site, function(x) length(unique(x)))

#generates abundance by site
abundance_by_site <- aggregate(reptiles_2022_full, Species~Site, length)


#generate Appendix 1
abundance <- data.frame(reptiles_2022_full %>% group_by(Site, Species) %>% dplyr::count(Site))
  #pivot abundance data to look like Appendix 1
abund <- pivot_wider(abundance, names_from=Site, values_from=n) #appendix 1 data

#---------------------

### SAMPLING EFFORT DATA

#Sampling effort in Days
effort <- reptiles_2022_full %>% group_by(Site,Date) %>% summarise(Freq=n())
effort_days <- effort %>% group_by(Site) %>% count()

#Sampling effort in trap Checks
effort1 <- reptiles_2022_full %>% group_by(Site,Date,Time) %>% summarise(Freq=n())
effort_checks <- effort1 %>% group_by(Site, Time) %>% count()


#-----------------------
#-----------------------


### COMPARING DRIFT FENCE TYPES


#-----------------------

### CREATING A COMMUNITY MATRIX

#read in comparison data (silt fence [2021] vs. novel drift fence [2022])
SNE_2021_2022 <- read.csv("NM_reptiles_2021-2022_sub_data.csv", header=TRUE)

data_table <- data.frame(table(SNE_2021_2022$Species, SNE_2021_2022$Date, SNE_2021_2022$year)) #creating frequency table by Species/Date/Year
data_matrix <- pivot_wider(data_table, names_from=Var1, values_from=Freq) #pull frequency table into a matrix

#renaming columns
names(data_matrix)[names(data_matrix) == "Var2"] <- "Date"
names(data_matrix)[names(data_matrix) == "Var3"] <- "Year"

#calculating abundance and richness
data_matrix$abundance <- rowSums(data_matrix[,sapply(data_matrix, is.integer)]) #calculating abundance per replicate
data_matrix <- data_matrix[data_matrix$abundance!=0,] #removing rows that have no reptiles
data_matrix$richness <- apply(data_matrix[names(data_matrix) %in% unique(SNE_2021_2022$Species)],1,function(x) sum(x > 0)) #calculating richness per replicate

#controlling for amount of traps and fencing
data_matrix$abundance_trap <- if_else(data_matrix$Year == "2021", data_matrix$abundance/8, data_matrix$abundance/6) #controlling for number of traps per drift fence type
data_matrix$abundance_fence <- if_else(data_matrix$Year == "2021", data_matrix$abundance/36.576, data_matrix$abundance/24) #controlling for number of traps per drift fence type
data_matrix$richness_trap <- if_else(data_matrix$Year == "2021", data_matrix$richness/8, data_matrix$richness/6) #controlling for amount of fencing per drift fence type
data_matrix$richness_fence <- if_else(data_matrix$Year == "2021", data_matrix$richness/36.576, data_matrix$richness/24) #controlling for amount of fencing per drift fence type

#-----------------------

### COMPARISONS OF ABUNDANCE AND RICHNESS

#Kruskal-Wallis Test

#abundance normal
kruskal.test(abundance ~ Year, data = data_matrix)
#abundance adjusted #traps
kruskal.test(abundance_trap ~ Year, data = data_matrix)
#abundance adjusted #fence
kruskal.test(abundance_fence ~ Year, data = data_matrix)

#richness normal
kruskal.test(richness ~ Year, data = data_matrix)
#richness adjusted #traps
kruskal.test(richness_trap ~ Year, data = data_matrix)
#richness adjusted #fence
kruskal.test(richness_fence ~ Year, data = data_matrix)


#Box-and-whisker plots

#creating ggplot-friendly dataframe
#filtering dataframe by abundance/richness type, adding 'fence_type' column
norm <- data.frame(data_matrix[c("Year", "abundance", "richness")], info_type = rep("Normal", nrow(data_matrix)))
trap <- data.frame(data_matrix[c("Year", "abundance_trap", "richness_trap")], info_type = rep("#Traps", nrow(data_matrix)))
fence <- data.frame(data_matrix[c("Year", "abundance_fence", "richness_fence")], info_type = rep("Fencing (m)", nrow(data_matrix)))

#renaming columns to match each other (for rbind function)
names(trap) <- names(norm)
names(fence) <- names(norm)

#creating new dataframe
plot_fence <- rbind(norm, trap, fence)
#reworking year column to fence_type column
plot_fence$Year <- if_else(plot_fence$Year == "2021", "Silt Fence", "Novel Array")
plot_fence <- rename(plot_fence, "fence_type" = "Year") #renaming Year column

#abundance figure 3a
tiff("../figures/abundance_box&whisker_figure_3a_v2.tiff", units="in", width=7.921687, height=5, res=600)

ggplot(data=plot_fence, aes(x=fence_type, y=abundance, colour=info_type)) + geom_boxplot(position = position_dodge(width = 0.9)) +
  labs(x="Abundance", y="Fence Type", colour="Data Type") + theme_classic(22)

dev.off()

#species richness figure 3b
tiff("../figures/richness_box&whisker_figure_3b_v2.tiff", units="in", width=7.921687, height=5, res=600)

ggplot(data=plot_fence, aes(x=fence_type, y=richness, colour=info_type)) + geom_boxplot(position = position_dodge(width = 0.9)) + 
  labs(x="Species Richness", y="Fence Type", colour="Data Type") + theme_classic(22)

dev.off()

#-----------------------

### ORDINATION SPACE ANALYSES

#rescaling raw abundance data to proportion data, each row = 1 (100%)
prop_matrix<-cbind((data_matrix[names(data_matrix) %in% unique(SNE_2021_2022$Species)]/rowSums(data_matrix[names(data_matrix) %in% unique(SNE_2021_2022$Species)]))*100,
                   data_matrix[c("Date", "Year")])

#PERMANOVA test for bray-curtis and jaccard disimilarity
adonis2(prop_matrix[,names(prop_matrix) %in% unique(SNE_2021_2022$Species)]~prop_matrix$Year, method ="bray", perm=999) #bray-curtis
adonis2(prop_matrix[,names(prop_matrix) %in% unique(SNE_2021_2022$Species)]~prop_matrix$Year, method ="jaccard", perm=999) #jaccard

