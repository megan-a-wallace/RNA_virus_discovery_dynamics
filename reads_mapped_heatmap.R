#######################################################
## Matrix of relative virus reads mapped per dataset ##
#######################################################

##Mar 2021
##Megan A. Wallace

##packages
#install.packages('pheatmap')
library(pheatmap)
library(dplyr)
library(grid)

##setting working directory
setwd("C:/Users/s1667991/Dropbox/PhD - 1st Year/Sequencing/05_19_totalRNASeq/TotalRNAseqdata_Edinburghgenomics_09_19/")

##importing data
read.csv(file = "relative_virus_mapped_reads.csv", 
         header = TRUE, row.names = "virus", stringsAsFactors = FALSE, strip.white = TRUE)->reads_mapped_dat

#converting to a numeric matrix, and transposing
reads_mapped_mat<-as.matrix(reads_mapped_dat)

col_pallette<-colorRampPalette(c("white","white","azure","steelblue","darkgoldenrod1","firebrick3"))(20)

###Trying with a log scale on the matrix
log_mat<-reads_mapped_mat
log_mat[log_mat>0] <- log10(log_mat[log_mat>0])
log_mat[log_mat==0] <- -7

pheatmap(log_mat, cluster_cols = FALSE)->hm_1

pheatmap(log_mat,cluster_cols = FALSE, gaps_col = c(1,2,3,4),
         labels_col = c("Sept - Oct '16", "Dec '16 - Jun '17","Jul - Oct '17", "Apr - Jun '18","Jul - Oct '18"),
         angle_col = 45, fontsize = 16, 
         cellwidth = 45, cellheight = 15,
         color = col_pallette, border_color = NA, 
         labels_row = gsub(pattern = "_",replacement = " ", x = hm_1$tree_row$labels), 
         legend_breaks = c(-5,-3,-1,1,2))->hm_2

pdf(file = "../../../thesis/ch2_metagenomic_virus_discovery/reads_mapped_heatmap.pdf",height = 14,width = 9)
hm_2
dev.off()

