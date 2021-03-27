###################################
## Known virus host range matrix ## 
###################################

##Feb 2021
##Megan A. Wallace

##packages
#install.packages('pheatmap')
library(pheatmap)
library(dplyr)
library(grid)

##setting working directory
setwd("C:/Users/s1667991/Dropbox/PhD - 1st Year/Data/Fly Collections/Spatial_temporal.prev/data_processing_Aug17_Oct18/virus_presence_absence_PCRs")

##importing data
read.csv(file = "all_viruses_host_range_matrix_by_yr_binary_minimal.csv", 
         header = TRUE, row.names = "species_collection", stringsAsFactors = FALSE, strip.white = TRUE)->host_range_dat
host_range_dat$species<-factor(host_range_dat$species, levels = c("imm","mel","fun","pha","obs","sub","sil","tri","hyd","hel","chy","cam","vir","def","bus"))
host_range_dat$time<-factor(host_range_dat$time, levels = c("2016","2017","2018","all"))
host_range_dat<-host_range_dat %>% arrange(species,time)

##turning the binary cols to numeric in a loop
for (i in 3:57) {
  host_range_dat[,i]<-as.numeric(host_range_dat[,i])
}

#converting to a numeric matrix, and transposing
host_range_mat<-t(as.matrix(host_range_dat[,4:57]))

#Setting up the annotations for collection
Time_col_annotations <- data.frame(Collections = host_range_dat$time)
row.names(Time_col_annotations) <- colnames(host_range_mat)

Time_col_annotations$Collections<-ifelse(Time_col_annotations$Collections == "2016","2016",
                                              ifelse(Time_col_annotations$Collections == "2017", "2017",
                                                     ifelse(Time_col_annotations$Collections == "2018", "2018",
                                                            ifelse(Time_col_annotations$Collections == "all", "Whole dataset", no = "missing")
                                                            )
                                                     )
                                  )


ann_colours<-list(Collections = c("2016"="#E69F00", "2017"="#56B4E9","2018"="#F0E442", "Whole dataset"="#0072B2"), Order = c("Mononegavirales"="#004949","Bunyavirales"="#009292","Picornavirales"="#ff6db6","Jingchuvirales"="#ffb6db","Martellivirales"="#490092", "Ghabrivirales"="#006ddb","Wolframvirales"="#b66dff","Reovirales"="#6db6ff","Nodamuvirales"="#b6dbff","Sobelivirales"="#920000","Permutotetraviridae"="#924900","Durnavirales"="#db6d00","Birnaviridae"="#24ff24","Amarillovirales"="#ffff6d","Piccovirales"="#000000","Nudiviridae"="#1B9E77","Flaviviridae"="#D95F02",), Genome = c("ssRNA (+)"="#7570B3","ssRNA (-)"="#E7298A","dsRNA"="#66A61E","dsDNA"="#E6AB02","ssDNA"="#A6761D"))

virus_labels<-c("Inveresk virus","Sunshine virus","Sighthill virus","Vogrie virus","Tranent virus","Lasswade virus","Cockenzie DCV variant","Burdiehouse burn virus","Crammond virus","Dalkeith virus","River Almond virus","Hillwood park virus","Craighall virus","Inverleith virus","Midmar virus","Glencorse burn virus","North Esk virus","Gosford virus","Dimm Sigma virus","Dimm Nora virus","Dmel Nora virus","Dsub Nora virus","Prestney Burn virus","Muthill virus","Grom virus","Motts mill virus","Galbut virus","Chaq (satellite of Galbut)", "Dmel Sigma virus","Drosophila A virus","Drosophila C virus","American Noda virus","Drosophila X virus","Charvil virus","Twyford virus","Thika virus","Killifi virus","Brandeis virus","Berkeley virus","Vesanto virus","Kallithea virus","Cherry Gardens virus","Pow Burn virus","Craigies hill virus","Craigmillar park virus","La Jolla virus","Corseley virus","Lye green virus","Kinkell virus","Buckhurst virus","Hermitage virus","Withyham virus","Larkfield virus","Tartou virus")

virus_order_row_annotations <- data.frame(Order = c("Mononegavirales","Bunyavirales","Bunyavirales","Reovirales","Bunyavirales","Mononegavirales","Picornavirales","Jingchuvirales","Martellivirales","Jingchuvirales","Martellivirales","Martellivirales","Ghabrivirales","Ghabrivirales","Nodamuvirales","Reovirales","Bunyavirales","Wolframvirales","Mononegavirales","Picornavirales","Picornavirales","Picornavirales","Sobelivirales","Martellivirales","Sobelivirales","Sobelivirales","Durnavirales","Durnavirales","Mononegavirales","Permutotetraviridae","Picornavirales","Nodamuvirales","Birnaviridae","Amarillovirales","Picornavirales","Picornavirales","Picornavirales","Martellivirales","Picornavirales","Piccovirales","Nudiviridae","Mononegavirales","Picornavirales","Nodamuvirales","Nodamuvirales","Picornavirales","Nodamuvirales","Mononegavirales","Picornavirales","Martellivirales","Flaviviridae","Mononegavirales","Ghabrivirales","Picornavirales"), Genome = c("ssRNA (-)","ssRNA (-)","ssRNA (-)","dsRNA","ssRNA (-)","ssRNA (-)","ssRNA (+)","ssRNA (-)","ssRNA (+)","ssRNA (-)","ssRNA (+)","ssRNA (+)","dsRNA","dsRNA","ssRNA (+)","dsRNA","ssRNA (-)","ssRNA (+)","ssRNA (-)","ssRNA (+)","ssRNA (+)","ssRNA (+)","ssRNA (+)","ssRNA (+)","ssRNA (+)","ssRNA (+)","dsRNA","dsRNA","ssRNA (-)","ssRNA (+)","ssRNA (+)","ssRNA (+)","dsRNA","ssRNA (+)","ssRNA (+)","ssRNA (+)","ssRNA (+)","ssRNA (+)","ssRNA (+)","ssDNA","dsDNA","ssRNA (-)","ssRNA (+)","ssRNA (+)","ssRNA (+)","ssRNA (+)","ssRNA (+)","ssRNA (-)","ssRNA (+)","ssRNA (+)","ssRNA (+)","ssRNA (-)","dsRNA","ssRNA (+)"))
row.names(virus_order_row_annotations) <- rownames(host_range_mat)

pheatmap(host_range_mat, scale = "none",
         color = c("white","#009E73"), legend = FALSE, cluster_cols = FALSE,gaps_col = c(3,4,5,6,9,12,15,16,17,18,19,20,21,22), annotation_col = Time_col_annotations, annotation_row = virus_order_row_annotations, annotation_colors = ann_colours, labels_col = c("","D.immigrans","","D.melanogaster","D.funebris","D.phalerata","","D.obscura","","","D.subobscura","","","D.subsilvestris","","D.tristis","D.hydei","D.helvetica","C.costata","H.cameraria","virilis sp.","S.deflexa","D.busckii"), angle_col = "45", fontsize_col = 13, labels_row = virus_labels, fontsize = 12, fontsize_row = 13, border_color = "grey60",cellwidth = 17, cellheight = 14)
