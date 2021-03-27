##########################################################################
### Annotating Phylogenies of Sanger sequencing from spat-temp analyses ##
##########################################################################

##Megan A. Wallace
##March 2021


##installing the latest Bioconductor version
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")

#'installing ggtree
BiocManager::install("ggtree")
install.packages("treeio")
install.packages("phytools")
install.packages("rlang")
#'loading ggtree into session 
library(ggplot2)
library(ggtree)
library(treeio)
library(scales)
library(rlang)
library(grid)
library(ape)
library(stringr)
library(phytools)

##ImmSV linked tree L and N gene
Imm_SV_tree<-read.beast(file = "Imm_SV/Imm_SV_LN_ST_exponpop.MCC.trees.nexus")

Imm_SV_tree_info<-fortify(Imm_SV_tree)

#pallette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#annotation fields
get.fields(Imm_SV_tree)

#plotting tree with filled circles for posterior values >=0.6
cols <- c("imm" = "#E69F00", "mel" = "#56B4E9")
Imm_SV_tree_plot<-ggtree(Imm_SV_tree, ladderize = TRUE, mrsd = "2018-08-21", size = 1.1, aes(color=host_species)) +
  scale_colour_manual(values = cols,name = "Host Species", na.translate = F, breaks = c("imm","mel"),labels = c("Dimm","Dmel")) +
  geom_point2(aes(subset=!is.na(posterior) & posterior > 0.59),color="grey20", size=3.5) +
  theme_tree2(legend.position=c(.1, .8),legend.text=element_text(size=19),legend.title=element_text(size=20))+geom_tippoint(size = 4)+
  theme(axis.line.x = element_line(size = 1.3), axis.text.x = element_text(size=20,color = "black"), axis.ticks.length.x = unit(0.4,"lines"), plot.margin = unit(c(1,1,1,1),"cm")) + scale_x_continuous(breaks = c(2015,2017.5,2020),minor_breaks = seq(2015,2020,2.5),limits = c(2014,2020))

pdf(file = "../../../thesis/ch2_metagenomic_virus_discovery/Imm_SV_tree_plot.pdf",height = 6,width = 13)
Imm_SV_tree_plot
dev.off()

##Imm Nora
Imm_nora_tree<-read.beast(file = "Imm_Nora/trees/Imm_nora_SS_wo_outgroup_codonaln_exponpop_GR1_combined.MCC.trees")

Imm_nora_tree_info<-fortify(Imm_nora_tree)

#pallette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#annotation fields
get.fields(Imm_nora_tree)

#plotting tree with filled circles for posterior values >=0.6
cols <- c("imm" = "#E69F00", "mel" = "#56B4E9", "obs" = "#009E73", "pha" = "#F0E442", "sub" = "#0072B2", "sil" = "#D55E00")
Imm_nora_tree_plot<-ggtree(Imm_nora_tree, ladderize = TRUE, mrsd = "2018-08-21", size = 1.1, aes(color=species)) +
  scale_colour_manual(values = cols,name = "Host Species", na.translate = F, breaks = c("imm","mel","obs","pha","sub","sil"),labels = c("Dimm","Dmel","Dobs","Dpha","Dsub","Dsus")) +
  geom_point2(aes(subset=!is.na(posterior) & posterior > 0.59),color="grey20", size=3.5) +
  theme_tree2(legend.position=c(.1, .8),legend.text=element_text(size=19),legend.title=element_text(size=20))+geom_tippoint(size = 4)+
  theme(axis.line.x = element_line(size = 1.3), axis.text.x = element_text(size=20,color = "black"), axis.ticks.length.x = unit(0.4,"lines"), plot.margin = unit(c(1,1,1,1),"cm")) + scale_x_continuous(breaks = c(2010,2015,2020),minor_breaks = seq(2010,2020,2.5),limits = c(2010,2020))

pdf(file = "../../../thesis/ch2_metagenomic_virus_discovery/Imm_nora_tree_plot.pdf",height = 6,width = 13)
Imm_nora_tree_plot
dev.off()

##Prestney burn
PB_tree<-read.beast(file = "Prestney_burn/Prestney_burn_con_16_seqs_exponpop.mcc.trees.nexus.txt")

PB_tree_info<-fortify(PB_tree)

##tidying the labels on the tree, removing underscores
#removing accession numbers and underscores - this is ugly code - fix!
PB_tree@phylo$tip.label<-gsub("[0-9]{4}_[0-9]{2}_","",PB_tree@phylo$tip.label)
PB_tree@phylo$tip.label<-gsub("_PB*","",PB_tree@phylo$tip.label)

#pallette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#annotation fields
get.fields(PB_tree)

#plotting tree with filled circles for posterior values >=0.6

cols <- c("hyd" = "#E69F00", "mel" = "#56B4E9", "obs" = "#009E73", "tri" = "#F0E442", "sub" = "#0072B2", "sil" = "#D55E00")
PB_tree_plot<-ggtree(PB_tree, ladderize = TRUE, mrsd = "2018-09-30", size = 1.1, aes(color=host_species)) +
  scale_colour_manual(values = cols,name = "Host Species", na.translate = F, breaks = c("hyd","mel","obs","tri","sub","sil"),labels = c("Dhyd","Dmel","Dobs","Dtri","Dsub","Dsus")) +
  geom_point2(aes(subset=!is.na(posterior) & posterior > 0.59),color="grey20", size=3.5) +
  theme_tree2(legend.position=c(.1, .8),legend.text=element_text(size=19),legend.title=element_text(size=20))+geom_tippoint(size = 4)+
  theme(axis.line.x = element_line(size = 1.3), axis.text.x = element_text(size=20,color = "black"), axis.ticks.length.x = unit(0.4,"lines"), plot.margin = unit(c(1,1,1,1),"cm")) + scale_x_continuous(breaks = c(2010,2015,2020),minor_breaks = seq(2010,2020,2.5),limits = c(2007,2020))

pdf(file = "../../../thesis/ch2_metagenomic_virus_discovery/PB_tree_plot.pdf",height = 8,width = 13)
PB_tree_plot
dev.off()

##Grom 
grom_tree<-read.beast(file = "Grom/Grom_P1_30_con_seqs_codon_exponpop_combined.MCC.trees.nexus")
#grom_tree<-midpoint.root(grom_tree)

grom_tree_info<-fortify(grom_tree)

##tidying the labels on the tree, removing underscores
#removing accession numbers and underscores - this is ugly code - fix!
grom_tree@phylo$tip.label<-gsub("[0-9]{4}_[0-9]{2}_","",grom_tree@phylo$tip.label)
grom_tree@phylo$tip.label<-gsub("_GR*","",grom_tree@phylo$tip.label)
grom_tree@phylo$tip.label<-gsub("_R","",grom_tree@phylo$tip.label)
grom_tree@phylo$tip.label<-gsub("_F","",grom_tree@phylo$tip.label)
grom_tree@phylo$tip.label<-gsub("_[0-9]{2}[A-Z]{1}","",grom_tree@phylo$tip.label)
grom_tree@phylo$tip.label<-gsub("_[0-9]{1}","",grom_tree@phylo$tip.label)

#pallette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#annotation fields
get.fields(grom_tree)

#plotting tree with filled circles for posterior values >=0.6

cols <- c("imm" = "#E69F00", "mel" = "#56B4E9", "obs" = "#009E73", "tri" = "#F0E442", "sub" = "#0072B2", "sil" = "#D55E00")
grom_tree_plot<-ggtree(grom_tree, ladderize = TRUE, mrsd = "2018-10-20", size = 1.1, aes(color=species)) +
  scale_colour_manual(values = cols,name = "Host Species", na.translate = F, breaks = c("imm","mel","obs","tri","sub","sil"),labels = c("Dimm","Dmel","Dobs","Dtri","Dsub","Dsus")) +
  geom_point2(aes(subset=!is.na(posterior) & posterior > 0.59),color="grey20", size=3.5) +
  theme_tree2(legend.position=c(.1, .8),legend.text=element_text(size=19),legend.title=element_text(size=20)) + 
  theme(axis.line.x = element_line(size = 1.3), axis.text.x = element_text(size=20,color = "black"), axis.ticks.length.x = unit(0.4,"lines"), plot.margin = unit(c(1,1,1,1),"cm")) + scale_x_continuous(breaks = c(2010,2015,2020),minor_breaks = seq(2010,2020,2.5),limits = c(2010,2020))+geom_tippoint(size = 4)

pdf(file = "../../../thesis/ch2_metagenomic_virus_discovery/grom_tree_plot.pdf",height = 8,width = 13)
grom_tree_plot
dev.off()

