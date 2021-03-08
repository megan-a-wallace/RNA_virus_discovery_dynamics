library(ggplot2)
library(ggthemes)
library(plyr)
library(dplyr)
library(tidyr)
library(gridExtra)
library(scales)

##importing coverage data across new virus genomes from SeptOct2016 data
#setwd("C:/Users/s1667991/Dropbox/PhD - 1st Year/Sequencing/05_19_totalRNASeq/genome_presentation")
# Reads mapped across new virus genomes - very-sensitive mapping ----------

##importing read depth data across new viruses generated using --very-sensitive read mapping to bowtie database, SeptOct2016 data
SeptOct2016_new_virus_readdepth_s <- read.table("SeptOct2016_11_19_new_viruses_5_datasets_ReadDepth_sensitive.tsv", sep = "")
cols<-c("virus","position","depth")
colnames(SeptOct2016_new_virus_readdepth_s)<-cols
SeptOct2016_new_virus_readdepth_s$position<-as.numeric(as.character(SeptOct2016_new_virus_readdepth_s$position))
SeptOct2016_new_virus_readdepth_s$depth<-as.numeric(as.character(SeptOct2016_new_virus_readdepth_s$depth))
##adding column which specifys the dataset
SeptOct2016_new_virus_readdepth_s$dataset="SeptOct2016"

Dec2016June2017_new_virus_readdepth_s <- read.table("Dec2016June2017_11_19_new_viruses_5_datasets_ReadDepth_sensitive.tsv", sep = "")
cols<-c("virus","position","depth")
colnames(Dec2016June2017_new_virus_readdepth_s)<-cols
Dec2016June2017_new_virus_readdepth_s$position<-as.numeric(as.character(Dec2016June2017_new_virus_readdepth_s$position))
Dec2016June2017_new_virus_readdepth_s$depth<-as.numeric(as.character(Dec2016June2017_new_virus_readdepth_s$depth))
##adding column which specifys the dataset
Dec2016June2017_new_virus_readdepth_s$dataset="Dec2016June2017"

DrosMS_JulOct17_new_virus_readdepth_s <- read.table("JulOct17_11_19_new_viruses_5_datasets_ReadDepth_sensitive.tsv", sep = "")
cols<-c("virus","position","depth")
colnames(DrosMS_JulOct17_new_virus_readdepth_s)<-cols
DrosMS_JulOct17_new_virus_readdepth_s$position<-as.numeric(as.character(DrosMS_JulOct17_new_virus_readdepth_s$position))
DrosMS_JulOct17_new_virus_readdepth_s$depth<-as.numeric(as.character(DrosMS_JulOct17_new_virus_readdepth_s$depth))
##adding column which specifys the dataset
DrosMS_JulOct17_new_virus_readdepth_s$dataset="DrosMS_JulOct17"

DrosMS_AprJun18_new_virus_readdepth_s <- read.table("AprJun18_11_19_new_viruses_5_datasets_ReadDepth_sensitive.tsv", sep = "")
cols<-c("virus","position","depth")
colnames(DrosMS_AprJun18_new_virus_readdepth_s)<-cols
DrosMS_AprJun18_new_virus_readdepth_s$position<-as.numeric(as.character(DrosMS_AprJun18_new_virus_readdepth_s$position))
DrosMS_AprJun18_new_virus_readdepth_s$depth<-as.numeric(as.character(DrosMS_AprJun18_new_virus_readdepth_s$depth))
##adding column which specifys the dataset
DrosMS_AprJun18_new_virus_readdepth_s$dataset="DrosMS_AprJun18"

DrosMS_JulOct18_new_virus_readdepth_s <- read.table("JulOct18_11_19_new_viruses_5_datasets_ReadDepth_sensitive.tsv", sep = "")
cols<-c("virus","position","depth")
colnames(DrosMS_JulOct18_new_virus_readdepth_s)<-cols
DrosMS_JulOct18_new_virus_readdepth_s$position<-as.numeric(as.character(DrosMS_JulOct18_new_virus_readdepth_s$position))
DrosMS_JulOct18_new_virus_readdepth_s$depth<-as.numeric(as.character(DrosMS_JulOct18_new_virus_readdepth_s$depth))
##adding column which specifys the dataset
DrosMS_JulOct18_new_virus_readdepth_s$dataset="DrosMS_JulOct18"

# Plotting coverage for new viruses across all datasets - senstitive --------

###Now we have all of the datasets we want to bind the dataframes together and then access each virus (across all datasets) separately
rbind.data.frame(SeptOct2016_new_virus_readdepth_s,Dec2016June2017_new_virus_readdepth_s,DrosMS_JulOct17_new_virus_readdepth_s,DrosMS_AprJun18_new_virus_readdepth_s,DrosMS_JulOct18_new_virus_readdepth_s)->all_datasets_read_depth_s

##Summing across the positions to create a single coverage plot per virus
coverage_all_datasets_by_position_s <- all_datasets_read_depth_s %>%
  group_by(virus,position) %>%
  summarise(depth=sum(depth))

glimpse(coverage_all_datasets_by_position_s)
class(coverage_all_datasets_by_position_s$depth)

##creating loop to plot coverage for each virus/segment in the table
levels(coverage_all_datasets_by_position_s$virus)<-c("Sunshine_virus","Burdiehouse_burn_virus","Cockenzie_virus","Crammond_virus","Dalkeith_virus","Dansoman_like_v_s2","Dansoman_like_v_s1","Phasmav_L","Hubei_odonate_like_s3","Hubei_odonate_like_s2","Hubei_odonate_like_s5","Hubei_odonate_like_s1","Hubei_odonate_like_s4","Hubei_toti_42.4","Hubei_toti_42.7","Inveresk_virus","Lasswade_virus","Loreto_like_virus","Negev_like_virus","Phasmav_M","Phasmav_S","Goose_discist_picornav","Sighthill_L","Sighthill_M","Tranent_L","Tranent_M","Tranent_S","Vogrie_1","Vogrie_2","Vogrie_5","Narna_v")
#converting coverage_all_datasets_by_position to a proper data frame
df_coverage_all_datasets_by_position_s<- data.frame(coverage_all_datasets_by_position_s$virus,coverage_all_datasets_by_position_s$position,coverage_all_datasets_by_position_s$depth)
colnames(df_coverage_all_datasets_by_position_s)=cols
str(df_coverage_all_datasets_by_position_s)

#input data for loop
x<-list(levels(df_coverage_all_datasets_by_position_s$virus))
data <-df_coverage_all_datasets_by_position_s

#function for transforming non-zero depth vals to log scale
non_zero_trans <- function(x)
{
  x <- log10(x+1)
  return(x)
}

#clearing plot list
pltList <- list()

for (i in 1:length(x[[1]])){
  ##filtering data to subject virus
  dat<-filter(data,virus==x[[1]][i])
  
  dat$depth<-non_zero_trans(dat$depth)
  
  ##setting limits
  xlim<-11000
  
  ##plotting
  pltList[[i]]<- print(ggplot(dat)+
                         geom_area(aes(y=depth, x = position), fill ="#009E73", show.legend = FALSE, stat="identity")+
                         theme(axis.text = element_text(size = rel(.6), margin = margin(t=.5,r=.5,unit = "cm")),
                               panel.background = element_blank(),
                               plot.subtitle = element_text(size = rel(1),margin = margin(b=.3)),
                               axis.ticks.length = unit(0.01,"cm"),
                               panel.border = element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.line = element_line(colour = "black", size = .3),
                               plot.margin = unit(c(.1,.3,.1,.3), "cm"))+
                         scale_x_continuous(limits = c(0,xlim), expand = expansion(mult = c(.01,.01)))+
                         scale_y_continuous(labels = math_format(10^.x))+
                         xlab("")+
                         ylab("")+
                         ggtitle("")+
                         labs(subtitle = paste(gsub("_", " ",x[[1]][i]))))
  print(paste(x[[1]][i],"coverage plot printed", sep=" "))
} 

#printing the coverage plots for new viruses to multiple pages 
ml <- marrangeGrob(pltList, nrow=4, ncol=2)
ggsave(plot = ml, file = "multipage_newvirus_readmapping_alldatasets_sensitive_MW.pdf")

