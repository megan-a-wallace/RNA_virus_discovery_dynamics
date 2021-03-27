#################################################
## Heatmep of newly-described virus host range ##
#################################################

##Dec 2020
##Megan A. Wallace

##packages
#install.packages('pheatmap')
library(pheatmap)
library(dplyr)
library(grid)

##setting working directory
setwd("C:/Users/s1667991/Dropbox/PhD - 1st Year/Data/Fly Collections/Spatial_temporal.prev/data_processing_Aug17_Oct18/virus_presence_absence_PCRs")

##importing data
read.csv(file = "novel_viruses_host_range_matrix_binary_minimal.csv", 
         header = TRUE, row.names = "species_collection", stringsAsFactors = FALSE, strip.white = TRUE)->host_range_dat
host_range_dat$species<-factor(host_range_dat$species, levels = c("imm","mel","fun","pha","obs","sub","sil","tri","hyd","hel","chy","cam","vir","def","bus"))
host_range_dat$collection<-factor(host_range_dat$collection, levels = c("SO16","DJ","JO17","AJ18","JO18"))
host_range_dat<-host_range_dat %>% arrange(species,collection)

##turning the binary cols to numeric in a loop
for (i in 3:32) {
  host_range_dat[,i]<-as.numeric(host_range_dat[,i])
}

##making a version of the data where the segmented viruses are combined into a single results column
#combined
host_range_dat_comb_seg<-host_range_dat

#tranent
host_range_dat_comb_seg$tranent<-rowSums(host_range_dat_comb_seg[,grep("tranent_", colnames(host_range_dat_comb_seg))])
host_range_dat_comb_seg$tranent[host_range_dat_comb_seg$tranent > 0] <- 1
#sighthill
host_range_dat_comb_seg$sighthill<-rowSums(host_range_dat_comb_seg[,grep("sighthill_", colnames(host_range_dat_comb_seg))])
host_range_dat_comb_seg$sighthill[host_range_dat_comb_seg$sighthill > 0] <- 1
#vogrie
host_range_dat_comb_seg$vogrie<-rowSums(host_range_dat_comb_seg[,grep("vogrie_", colnames(host_range_dat_comb_seg))])
host_range_dat_comb_seg$vogrie[host_range_dat_comb_seg$vogrie > 0] <- 1
#north_esk
host_range_dat_comb_seg$north_esk<-rowSums(host_range_dat_comb_seg[,grep("north_esk_", colnames(host_range_dat_comb_seg))])
host_range_dat_comb_seg$north_esk[host_range_dat_comb_seg$north_esk > 0] <- 1
#midmar
host_range_dat_comb_seg$midmar<-rowSums(host_range_dat_comb_seg[,grep("midmar_", colnames(host_range_dat_comb_seg))])
host_range_dat_comb_seg$midmar[host_range_dat_comb_seg$midmar > 0] <- 1
#glencorse_burn
host_range_dat_comb_seg$glencorse_burn<-rowSums(host_range_dat_comb_seg[,grep("glencorse_", colnames(host_range_dat_comb_seg))])
host_range_dat_comb_seg$glencorse_burn[host_range_dat_comb_seg$glencorse_burn > 0] <- 1

#converting to a numeric matrix, and transposing
host_range_mat<-t(as.matrix(host_range_dat[,3:32]))
host_range_mat_comb_seg<-t(as.matrix(host_range_dat_comb_seg[,c(3,4,13:21,32,34:39)]))

#Setting up the annotations for collection
collection_col_annotations <- data.frame(Collection = host_range_dat$collection)
row.names(collection_col_annotations) <- colnames(host_range_mat)

collection_col_annotations$Collection<-ifelse(collection_col_annotations$Collection == "SO16","Sep - Oct '16",
                                              ifelse(collection_col_annotations$Collection == "DJ", "Dec '16 - Jun '17",
                                                     ifelse(collection_col_annotations$Collection == "JO17", "Jul - Oct '17",
                                                            ifelse(collection_col_annotations$Collection == "AJ18", "Apr - Jun '18",
                                                                   ifelse(collection_col_annotations$Collection == "JO18", "Jul - Oct '18", no = "missing")
                                                                   )
                                                            )
                                                     )
                                              )

collection_col_annotations_comb_seg <- data.frame(Collection = host_range_dat_comb_seg$collection)
row.names(collection_col_annotations_comb_seg) <- colnames(host_range_mat_comb_seg)

collection_col_annotations_comb_seg$Collection<-ifelse(collection_col_annotations_comb_seg$Collection == "SO16","Sep - Oct '16",
                                              ifelse(collection_col_annotations_comb_seg$Collection == "DJ", "Dec '16 - Jun '17",
                                                     ifelse(collection_col_annotations_comb_seg$Collection == "JO17", "Jul - Oct '17",
                                                            ifelse(collection_col_annotations_comb_seg$Collection == "AJ18", "Apr - Jun '18",
                                                                   ifelse(collection_col_annotations_comb_seg$Collection == "JO18", "Jul - Oct '18", no = "missing")
                                                            )
                                                     )
                                              )
)

#Setting up the annotations for virus group, orders taken from 2019 release of ICTV
virus_order_row_annotations <- data.frame(Order = c("Mononegavirales","Bunyavirales","Bunyavirales","Bunyavirales","Reovirales","Reovirales","Reovirales","Bunyavirales","Bunyavirales","Bunyavirales","Mononegavirales","Picornavirales","Jingchuvirales","Martellivirales","Jingchuvirales","Martellivirales","Martellivirales","Ghabrivirales","Ghabrivirales","Nodamuvirales","Nodamuvirales","Reovirales","Reovirales","Reovirales","Reovirales","Reovirales","Bunyavirales","Bunyavirales","Bunyavirales","Wolframvirales"),Genome = c("ssRNA (-)","ssRNA (-)","ssRNA (-)","ssRNA (-)","dsRNA","dsRNA","dsRNA","ssRNA (-)","ssRNA (-)","ssRNA (-)","ssRNA (-)","ssRNA (+)","ssRNA (-)","ssRNA (+)","ssRNA (-)","ssRNA (+)","ssRNA (+)","dsRNA","dsRNA","ssRNA (+)","ssRNA (+)","dsRNA","dsRNA","dsRNA","dsRNA","dsRNA","ssRNA (-)","ssRNA (-)","ssRNA (-)","ssRNA (+)"))
row.names(virus_order_row_annotations) <- rownames(host_range_mat)

#Setting up the annotations for virus group, combined segments plot
virus_order_row_annotations_comb_seg <- data.frame(Order = c("Mononegavirales","Bunyavirales","Mononegavirales","Picornavirales","Jingchuvirales","Martellivirales","Jingchuvirales","Martellivirales","Martellivirales","Ghabrivirales","Ghabrivirales","Wolframvirales","Bunyavirales","Bunyavirales","Reovirales","Bunyavirales","Nodamuvirales","Reovirales"), Genome = c("ssRNA (-)","ssRNA (-)","ssRNA (-)","ssRNA (+)","ssRNA (-)","ssRNA (+)","ssRNA (-)","ssRNA (+)","ssRNA (+)","dsRNA","dsRNA","ssRNA (+)","ssRNA (-)","ssRNA (-)","dsRNA","ssRNA (-)","ssRNA (+)","dsRNA"))
row.names(virus_order_row_annotations_comb_seg) <- rownames(host_range_mat_comb_seg)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ann_colours<-list(Collection = c("Sep - Oct '16"="#E69F00", "Dec '16 - Jun '17"="#56B4E9","Jul - Oct '17"="#F0E442", "Apr - Jun '18"="#0072B2","Jul - Oct '18"="#CC79A7"), Order = c("Mononegavirales"="#004949","Bunyavirales"="#009292","Picornavirales"="#ff6db6","Jingchuvirales"="#ffb6db","Martellivirales"="#490092", "Ghabrivirales"="#006ddb","Wolframvirales"="#b66dff","Reovirales"="#6db6ff","Nodamuvirales"="#b6dbff"), Genome = c("ssRNA (+)"="#920000","ssRNA (-)"="#924900","dsRNA"="#db6d00"))

#longer colour-blind pallette
pal <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
         "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
         "#920000","#924900","#db6d00","#24ff24","#ffff6d")

virus_labels<-c("Inveresk virus","Sunshine virus","Lasswade virus","Cockenzie virus","Burdiehouse burn virus","Crammond virus","Dalkeith virus","River almond virus","Hillwood park virus","Craighall virus","Inverleith virus","Gosford virus","Tranent virus","Sighthill virus","Vogrie virus","North esk virus","Midmar virus","Glencorse burn virus")

virus_seg_labels<-c("Inveresk virus","Sunshine virus","Sighthill virus (L)","Sighthill virus (M)","Vogrie virus (1)","Vogrie virus (5)","Vogrie virus (2)","Tranent virus (L)","Tranent virus (M)","Tranent virus (S)","Lasswade virus","Cockenzie virus","Burdiehouse burn virus","Crammond virus","Dalkeith virus","River almond virus","Hillwood park virus","Craighall virus","Inverleith virus","Midmar virus (1)","Midmar virus (2)","Glencorse burn virus (1)","Glencorse burn virus (2)","Glencorse burn virus (3)","Glencorse burn virus (4)","Glencorse burn virus (5)","North esk virus (L)","North esk virus (M)","North esk virus (S)","Gosford virus")

#With all segments separate
pheatmap(host_range_mat,scale = "none", 
         color = c("white","#009E73"), legend = FALSE, cluster_cols = FALSE, gaps_col = c(5,9,12,17,22,27,32,35,37,38,40,41,42,44), annotation_col = collection_col_annotations,annotation_row = virus_order_row_annotations, annotation_colors = ann_colours,labels_col = c("","","D.immigrans","","","","D.melanogaster","","","","D.funebris","","","","D.phalerata","","","","","D.obscura","","","","","D.subobscura","","","","","D.subsilvestris","","","","D.tristis","","D.hydei","","D.helvetica","C.costata","","H.cameraria","virilis sp.","S.deflexa","","D.busckii"), angle_col = "45", fontsize_col = 12, labels_row = virus_seg_labels, fontsize = 10, fontsize_row = 12, border_color = "grey60",cellwidth = 13, cellheight = 17)

downViewport("col_annotation.3-3-3-3")
grid.text(c("No. of flies",1333,161,28,24,197,286,164,20,3,3,2,1,1,3,1), x = c(-0.05, 0.051, 0.15, 0.226, 0.314, 0.422, 0.528, 0.635, 0.723, 0.780, 0.815, 0.855, 0.89, 0.916, 0.951, 0.988), y = rep.int(1.3, times = 16), gp = gpar(fontsize=11), vjust = 0.2, hjust = 0.5)

#With segments combined into viruses
pheatmap(host_range_mat_comb_seg, scale = "none",
         color = c("white","#009E73"), legend = FALSE, cluster_cols = FALSE,gaps_col = c(5,9,12,17,22,27,32,35,37,38,40,41,42,44),annotation_col = collection_col_annotations_comb_seg,annotation_row = virus_order_row_annotations_comb_seg, annotation_colors = ann_colours, labels_col = c("","","D.immigrans","","","","D.melanogaster","","","","D.funebris","","","","D.phalerata","","","","","D.obscura","","","","","D.subobscura","","","","","D.subsilvestris","","","","D.tristis","","D.hydei","","D.helvetica","C.costata","","H.cameraria","virilis sp.","S.deflexa","","D.busckii"), angle_col = "45", fontsize_col = 12, labels_row = virus_labels, fontsize = 11, fontsize_row = 12, border_color = "grey60",cellwidth = 13, cellheight = 17)

downViewport("col_annotation.3-3-3-3")
grid.text(c("No. of flies",1333,161,28,24,197,286,164,20,3,3,2,1,1,3,1), x = c(-0.05, 0.051, 0.15, 0.226, 0.314, 0.422, 0.528, 0.635, 0.723, 0.780, 0.815, 0.855, 0.89, 0.916, 0.951, 0.988), y = rep.int(1.3, times = 16), gp = gpar(fontsize=11), vjust = 0.2, hjust = 0.5)
