#'---
#'title: Sampling location map of spatial and temporal sampling 
#'author: Megan Wallace
#'date: 25th July 2017
#'---

setwd("C:/Users/s1667991/Dropbox/PhD - 1st Year/thesis/ch2_metagenomic_virus_discovery")

#'Loading required mapping packages into R

require(ggplot2)
require(scales)
require(ggmap) 
require(dplyr)
require(tidyr)

# For google map, you have to give the center of the window you are looking at.
# Possibility for the map type argument: terrain / satellite / roadmap / hybrid

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Load the sampling site data
sampling_site_data<-read.csv(file = "sampling_site_coords.csv", header = TRUE)

# get the map info
#I registered this API before 
#register_google(key = "AIzaSyCqBPPYeRoLIrdzbWwMZghNdYsO2yis2lI", write = TRUE)
edin_map <- get_googlemap("Musselburgh, Scotland", zoom = 10.5, maptype = "satellite")

# Plot it
ggmap(edin_map) + 
        geom_point(data=sampling_site_data, aes(x=lon, y=lat),size=7, color=cbPalette[5], alpha=0.7) +
        ylim(55.845,56.025) + xlim(-3.35,-2.82) + coord_map() +
        xlab("Longitude") + ylab("Lattitude") +
        theme(axis.title.x = element_text(vjust = -0.35, size = 14,
                                          colour = "black"),
              axis.title.y = element_text(vjust = 1.2,
                                          size = 14,
                                          colour = "black"),
              axis.text = element_text(size = 12,
                                       colour = "black")) 

##Now making it a bit more fancy, with cirles sized by how many flies I collected there
read.csv(file = "../../Data/Fly Collections/Spatial_temporal.prev/data_processing_Aug17_Oct18/virus_presence_absence_PCRs/Indiv_virus_PCR_results.csv", header = TRUE, )->virus_pa

virus_pa$site<-gsub(" ","",virus_pa$site) #there was aome trailing whitespace in SI

#Summarizing the number of flies collected at each site
flies_per_site <- virus_pa %>% 
        group_by(site) %>%
        summarise(no_flies = sum(no_flies))
colnames(flies_per_site)<-c("code","no_flies")

#combining the BB and BR sites, because they're about 50m from one another
flies_per_site[4, 2] <- flies_per_site[4, 2] + flies_per_site[1, 2]
flies_per_site<-flies_per_site[2:21,]

##Combining this w the sampling site coords
coords_flies<-merge(sampling_site_data,flies_per_site, by = "code")

# Create breaks for the color scale
mybreaks <- c(10, 25, 150, 500)
Darren_palette<-colorRampPalette(c("firebrick3","darkgoldenrod1","steelblue","white"))(10)
Darren_scale<-colorRampPalette(c("white","steelblue","darkgoldenrod1","firebrick3"))

# Reorder data to show biggest cities on top
coords_flies <- coords_flies %>%
        arrange(no_flies) %>%
        mutate(code=factor(code, unique(code))) 

# Build the map
ggmap(edin_map) + 
        geom_point(data=coords_flies,aes(x=lon, y=lat, size=no_flies, color=no_flies, alpha=no_flies), shape=20, stroke=FALSE) +
        scale_size_continuous(name="No. of flies collected", trans="log", range=c(5,20), breaks=mybreaks) +
        scale_alpha_continuous(name="No. of flies collected", trans="log", range=c(0.5, .9), breaks=mybreaks) +
        scale_color_gradientn(colours = Darren_scale(11)[-(11)], trans="log", breaks=mybreaks, name="No. of flies collected") +
        ylim(55.845,56.025) + xlim(-3.35,-2.82) + coord_map() +
        xlab("Longitude") + ylab("Lattitude") +
        theme(axis.title.x = element_text(vjust = -0.35, size = 14,
                                          colour = "black"),
              axis.title.y = element_text(vjust = 1.2,
                                          size = 14,
                                          colour = "black"),
              axis.text = element_text(size = 12,
                                       colour = "black")) +
        guides(colour = guide_legend()) +
        theme(text = element_text(color = "#22211d"),
              legend.background = element_rect(fill = "white", color = NA))

##And labeling w how many times I went to each site...eg. sampling effort
#Summarizing the number of flies collected at each site

##importing the effort data, how many times I visited a site regardless of whether I collected flies there
read.csv(file = "../../Data/Fly Collections/Spatial_temporal.prev/sampling_effort.csv", header = TRUE, )->effort_dat

effort_dat$month_year<-as.factor(effort_dat$month_year)

effort_per_site <- effort_dat %>% 
        group_by(code) %>%
        summarise(no_efforts = length(unique(month_year)))
colnames(effort_per_site)<-c("code","effort")

##Combining this w the sampling site coords
coords_flies<-merge(sampling_site_data,flies_per_site,effort_per_site, by = "code")
coords_flies<-Reduce(merge, list(sampling_site_data,flies_per_site,effort_per_site))

coords_flies <- coords_flies %>%
        arrange(no_flies) %>%
        mutate(code=factor(code, unique(code))) 

#############
## PANEL A ##
#############

#longer colour-blind pallette
pal <- c("#004949","#009292","#ff6db6","#ffb6db",
         "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
         "#920000","#924900","#db6d00","#24ff24","#ffff6d","#000000")

# Build the map
ggmap(edin_map) + 
        geom_point(data=coords_flies,aes(x=lon, y=lat, size=no_flies, color=no_flies, alpha=no_flies), shape=20, stroke=FALSE) +
        geom_text(data=coords_flies, aes(x=lon, y=lat, label=effort), size=5) +
        scale_size_continuous(name="No. of flies collected", trans="log", range=c(5,20), breaks=mybreaks) +
        scale_alpha_continuous(name="No. of flies collected", trans="log", range=c(0.5, .9), breaks=mybreaks) +
        scale_color_gradientn(colours = Darren_scale(11)[-(11)], trans="log", breaks=mybreaks, name="No. of flies collected") +
        ylim(55.845,56.025) + xlim(-3.35,-2.82) + coord_map() +
        xlab("Longitude") + ylab("Lattitude") +
        theme(axis.title.x = element_text(vjust = -0.35, size = 14,
                                          colour = "black"),
              axis.title.y = element_text(vjust = 1.2,
                                          size = 14,
                                          colour = "black"),
              axis.text = element_text(size = 12,
                                       colour = "black")) +
        guides(colour = guide_legend()) +
        theme(text = element_text(color = "#22211d"),
              legend.background = element_rect(fill = "white", color = NA),
              legend.key = element_blank(),
              legend.text = element_text(size = 12,colour = "black"),
              legend.title = element_text(size = 14,colour = "black"),
              legend.box.just = "left")
##I shifted the label off the smallest circle in inkscape

##############################################################
## PANEL B : Stacked boxplot of species collected over time ##
##############################################################
virus_pa$date<-as.Date(virus_pa$date, format = "%d/%m/%Y", ordered = TRUE)

virus_pa$species<-factor(virus_pa$species,levels = c("imm","sub","obs","sil","mel","fun","pha","tri","hyd","hel","def","chy","cam","vir","bus"),labels = c("D.immigrans", "D.subobscura","D.obscura","D.subsilvestris","D.melanogaster","D.funebris","D.phalerata","D.tristis","D.hydei","D.helvetica","S.deflexa","C.costata","H.cameraria","Virilis sp.","D.buskii"))

#combining the BB and BR sites, because they're about 50m from one another
flies_per_site[4, 2] <- flies_per_site[4, 2] + flies_per_site[1, 2]
flies_per_site<-flies_per_site[2:21,]

#reducing virus_pa to only the cols needed for stacked bar plot
temp_fly_stack<-select(virus_pa,date,species,no_flies)
##because the x intervals overlap for the two pilot collections, combining them for the purpose of this plot - as the stacks colours don't work otherwise
temp_fly_stack$date<-as.Date(gsub("2016-09-23","2016-09-26",temp_fly_stack$date),format = "%Y-%m-%d", ordered = TRUE)


#stacked bar plot
stack_plot<-ggplot(temp_fly_stack, aes(x = date, y = no_flies, fill = species)) + 
        geom_bar(stat = "identity", width = 13)+
        scale_x_date(date_labels ="%b-%y",date_breaks='1 month')+
        xlab("Collection Date")+
        ylab("No. of Flies Collected")+
        scale_fill_manual(values = pal)+
        theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size=14),axis.text.y = element_text(), axis.text.x = element_text(size=10,angle=45,hjust=1),plot.title = element_text(hjust = 0.5))+
        labs(fill = "Species")+
        annotate("rect", xmin = as.Date("2016-09-01"), xmax = as.Date("2016-11-30"), ymin = -30, ymax = -10,
                 alpha = .5,fill = "#009E73")+
        annotate("rect", xmin = as.Date("2016-12-02"), xmax = as.Date("2017-06-30"), ymin = -30, ymax = -10,
                 alpha = .5,fill = "#F0E442")+
        annotate("rect", xmin = as.Date("2017-07-02"), xmax = as.Date("2017-10-31"), ymin = -30, ymax = -10,
                 alpha = .5,fill = "#009E73")+
        annotate("rect", xmin = as.Date("2017-11-02"), xmax = as.Date("2018-06-30"), ymin = -30, ymax = -10,
                 alpha = .5,fill = "#F0E442")+
        annotate("rect", xmin = as.Date("2018-07-02"), xmax = as.Date("2018-10-31"), ymin = -30, ymax = -10,
                 alpha = .5,fill = "#009E73")+
        annotate("text", x = as.Date("2016-10-15"), y = -20, label = "1")+
        annotate("text", x = as.Date("2017-03-15"), y = -20, label = "2")+
        annotate("text", x = as.Date("2017-09-01"), y = -20, label = "3")+
        annotate("text", x = as.Date("2018-03-01"), y = -20, label = "4")+
        annotate("text", x = as.Date("2018-09-01"), y = -20, label = "5")

