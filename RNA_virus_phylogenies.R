#########################
## New RNA virus trees ##
#########################

##Using ggtree to make nicer looking phylogenies of new RNA viruses. 

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

##################################
## Negative sense ssRNA viruses ##
##################################

##Burdiehouse burn chuvirus
Burdiehouse_burn_tree<-read.tree(file = "Burdiehouse_burn_virus/02_21/Burdiehouse_burn_virus_5BLASTp.aln.fas.treefile")
Burdiehouse_burn_tree<-midpoint.root(Burdiehouse_burn_tree)

Burdiehouse_burn_tree_info<-fortify(Burdiehouse_burn_tree)

##tidying the labels on the tree, removing underscores
#removing accession numbers and underscores - this is ugly code - fix!
Burdiehouse_burn_tree$tip.label<-gsub("[A-Z]{3}[0-9]{5}","",Burdiehouse_burn_tree$tip.label)
Burdiehouse_burn_tree$tip.label<-gsub("[A-Z]{2}_[0-9]{9}","",Burdiehouse_burn_tree$tip.label)
Burdiehouse_burn_tree$tip.label<-gsub("\\.[1-9]{1}","",Burdiehouse_burn_tree$tip.label)
Burdiehouse_burn_tree$tip.label<-gsub("_"," ",Burdiehouse_burn_tree$tip.label)
Burdiehouse_burn_tree$tip.label<-gsub("^ ","",Burdiehouse_burn_tree$tip.label)

#generating a vector of colours to colour the new virus in red
Burdiehouse_cols<-ifelse(str_detect(Burdiehouse_burn_tree$tip.label,"Burdiehouse"),"#56B4E9",
                         ifelse(str_detect(Burdiehouse_burn_tree$tip.label,"polygyrus"),"#F0E442","black")) #to get label order to label with black and red colours
Burdiehouse_nodes<-Burdiehouse_burn_tree$node.label[str_detect(Burdiehouse_burn_tree$node.label,"[0-9]{2}")] #make subset of labels containing the bootstrap values > 59
Burdiehouse_nodes<-Burdiehouse_burn_tree$node.label %in% Burdiehouse_nodes[as.numeric(Burdiehouse_nodes)>59]

#plotting the tree
Burdiehouse_burn_tree_plot<-ggtree(Burdiehouse_burn_tree, ladderize = TRUE, size = 2.3) + 
                            geom_tiplab(size=12.75, color=Burdiehouse_cols, hjust = -0.02)+
                            geom_nodepoint(aes(subset=Burdiehouse_nodes), color="#0072B2", size=7.5)+
                            geom_treescale(color = "grey", width = 0.5, x=0, y=8.5, linesize = 4, fontsize = 0)+
                            xlim(0, 4)

pdf(file = "../../../../thesis/ch2_metagenomic_virus_discovery/Burdiehouse_burn_tree_plot.pdf",height = 6,width = 15)
Burdiehouse_burn_tree_plot
dev.off()

##North esk phasmavirus
North_esk_tree<-read.tree(file = "Phasmavirus/02_21/North_esk_virus_5BLASTp.aln.fas.treefile")
North_esk_tree<-midpoint.root(North_esk_tree)

North_esk_tree_info<-fortify(North_esk_tree)

##tidying the labels on the tree, removing underscores
#removing accession numbers and underscores - this is ugly code - fix!
North_esk_tree$tip.label<-gsub("[A-Z]{3}[0-9]{5}","",North_esk_tree$tip.label)
North_esk_tree$tip.label<-gsub("[A-Z]{2}_[0-9]{9}","",North_esk_tree$tip.label)
North_esk_tree$tip.label<-gsub("\\.[1-9]{1}","",North_esk_tree$tip.label)
North_esk_tree$tip.label<-gsub("_"," ",North_esk_tree$tip.label)
North_esk_tree$tip.label<-gsub("^ ","",North_esk_tree$tip.label)

#generating a vector of colours to colour the new virus in red
North_esk_cols<-ifelse(str_detect(North_esk_tree$tip.label,"North esk"),"#56B4E9","black") 
North_esk_nodes<-North_esk_tree$node.label[str_detect(North_esk_tree$node.label,"[0-9]{2}")] #make subset of labels containing the bootstrap values > 59
North_esk_nodes<-North_esk_tree$node.label %in% North_esk_nodes[as.numeric(North_esk_nodes)>59]

#plotting the tree
North_esk_tree_plot<-ggtree(North_esk_tree, ladderize = TRUE, size = 2.3) + 
  geom_tiplab(size=12.75, color=North_esk_cols, hjust = -0.02)+
  geom_nodepoint(aes(subset=North_esk_nodes), color="#0072B2", size=7.5)+
  geom_treescale(color = "grey", width = 0.5, x=0, y=8.5, linesize = 4, fontsize = 0)+
  xlim(0, 10)

pdf(file = "../../../../thesis/ch2_metagenomic_virus_discovery/North_esk_tree_plot.pdf",height = 7,width = 15)
North_esk_tree_plot
dev.off()

##Phleboviruses - Sighthill + Tranent
Phlebo_tree<-read.tree(file = "Phleboviruses/02_21/Phleboviruses_5BLASTp.aln.fas.treefile")
Phlebo_tree<-midpoint.root(Phlebo_tree)
Phlebo_tree_info<-fortify(Phlebo_tree)

##tidying the labels on the tree, removing underscores
#removing accession numbers and underscores - this is ugly code - fix!
Phlebo_tree$tip.label<-gsub("[A-Z]{3}[0-9]{5}","",Phlebo_tree$tip.label)
Phlebo_tree$tip.label<-gsub("[A-Z]{2}_[0-9]{9}","",Phlebo_tree$tip.label)
Phlebo_tree$tip.label<-gsub("\\.[1-9]{1}","",Phlebo_tree$tip.label)
Phlebo_tree$tip.label<-gsub("_"," ",Phlebo_tree$tip.label)
Phlebo_tree$tip.label<-gsub("^ ","",Phlebo_tree$tip.label)

#generating a vector of colours to colour the new virus in red
Phlebo_cols<-ifelse(str_detect(Phlebo_tree$tip.label,"Tranent|Sighthill"),"#56B4E9","black") 
Phlebo_nodes<-Phlebo_tree$node.label[str_detect(Phlebo_tree$node.label,"[0-9]{2}")] #make subset of labels containing the bootstrap values > 59
Phlebo_nodes<-Phlebo_tree$node.label %in% Phlebo_nodes[as.numeric(Phlebo_nodes)>59]

#plotting the tree
Phlebo_tree_plot<-ggtree(Phlebo_tree, ladderize = TRUE, size = 2.3) + 
  geom_tiplab(size=12.75, color=Phlebo_cols, hjust = -0.02)+
  geom_nodepoint(aes(subset=Phlebo_nodes), color="#0072B2", size=7.5)+
  geom_treescale(color = "grey", width = 0.5, x=0, y=8.5, linesize = 4, fontsize = 0)+
  xlim(0, 9)

pdf(file = "../../../../thesis/ch2_metagenomic_virus_discovery/Phlebo_tree_plot.pdf",height = 10,width = 18)
Phlebo_tree_plot
dev.off()

##Sunshine
Sunshine_tree<-read.tree(file = "Sunshine_L4_buny_virus/02_21/Sunshine_virus_5BLASTp.aln.fas.treefile")
Sunshine_tree<-midpoint.root(Sunshine_tree)
Sunshine_tree_info<-fortify(Sunshine_tree)

##tidying the labels on the tree, removing underscores
#removing accession numbers and underscores - this is ugly code - fix!
Sunshine_tree$tip.label<-gsub("[A-Z]{3}[0-9]{5}","",Sunshine_tree$tip.label)
Sunshine_tree$tip.label<-gsub("[A-Z]{2}_[0-9]{9}","",Sunshine_tree$tip.label)
Sunshine_tree$tip.label<-gsub("\\.[1-9]{1}","",Sunshine_tree$tip.label)
Sunshine_tree$tip.label<-gsub("_"," ",Sunshine_tree$tip.label)
Sunshine_tree$tip.label<-gsub("^ ","",Sunshine_tree$tip.label)

#generating a vector of colours to colour the new virus in red
Sunshine_cols<-ifelse(str_detect(Sunshine_tree$tip.label,"Sunshine"),"#56B4E9",
                      ifelse(str_detect(Sunshine_tree$tip.label,"Edhazardia aedis"),"#F0E442","black"))

Sunshine_nodes<-Sunshine_tree$node.label[str_detect(Sunshine_tree$node.label,"[0-9]{2}")] #make subset of labels containing the bootstrap values > 59
Sunshine_nodes<-Sunshine_tree$node.label %in% Sunshine_nodes[as.numeric(Sunshine_nodes)>59]

#plotting the tree
Sunshine_tree_plot<-ggtree(Sunshine_tree, ladderize = TRUE, size = 2.3) + 
  geom_tiplab(size=12.75, color=Sunshine_cols, hjust = -0.02)+
  geom_nodepoint(aes(subset=Sunshine_nodes), color="#0072B2", size=7.5)+
  geom_treescale(color = "grey", width = 0.5, x=0, y=6.5, linesize = 4, fontsize = 0)+
  xlim(0, 7)

pdf(file = "../../../../thesis/ch2_metagenomic_virus_discovery/Sunshine_tree_plot.pdf",height = 5,width = 12)
Sunshine_tree_plot
dev.off()

##Inveresk
Inveresk_tree<-read.tree(file = "Inveresk_virus_updated/02_21/Inveresk_Nyamivirus_5BLASTp.aln.fas.treefile")
Inveresk_tree<-midpoint.root(Inveresk_tree)
Inveresk_tree_info<-fortify(Inveresk_tree)

##tidying the labels on the tree, removing underscores
#removing accession numbers and underscores - this is ugly code - fix!
Inveresk_tree$tip.label<-gsub("[A-Z]{3}[0-9]{5}","",Inveresk_tree$tip.label)
Inveresk_tree$tip.label<-gsub("[A-Z]{2}_[0-9]{9}","",Inveresk_tree$tip.label)
Inveresk_tree$tip.label<-gsub("\\.[1-9]{1}","",Inveresk_tree$tip.label)
Inveresk_tree$tip.label<-gsub("_"," ",Inveresk_tree$tip.label)
Inveresk_tree$tip.label<-gsub("^ ","",Inveresk_tree$tip.label)

#generating a vector of colours to colour the new virus in red
Inveresk_cols<-ifelse(str_detect(Inveresk_tree$tip.label,"Inveresk"),"#56B4E9","black")

Inveresk_nodes<-Inveresk_tree$node.label[str_detect(Inveresk_tree$node.label,"[0-9]{2}")] #make subset of labels containing the bootstrap values > 59
Inveresk_nodes<-Inveresk_tree$node.label %in% Inveresk_nodes[as.numeric(Inveresk_nodes)>59]

#plotting the tree
Inveresk_tree_plot<-ggtree(Inveresk_tree, ladderize = TRUE, size = 2.3) + 
  geom_tiplab(size=12.75, color=Inveresk_cols, hjust = -0.02)+
  geom_nodepoint(aes(subset=Inveresk_nodes), color="#0072B2", size=7.5)+
  geom_treescale(color = "grey", width = 0.5, x=0, y=11.5, linesize = 4, fontsize = 0)+
  xlim(0, 12)

pdf(file = "../../../../thesis/ch2_metagenomic_virus_discovery/Inveresk_tree_plot.pdf",height = 8,width = 13)
Inveresk_tree_plot
dev.off()

##Lasswade
Lasswade_tree<-read.tree(file = "Lasswade_virus/02_21/Lasswade_virus_5BLASTp.aln.fas.treefile")
Lasswade_tree<-midpoint.root(Lasswade_tree)
Lasswade_tree_info<-fortify(Lasswade_tree)

##tidying the labels on the tree, removing underscores
#removing accession numbers and underscores - this is ugly code - fix!
Lasswade_tree$tip.label<-gsub("[A-Z]{3}[0-9]{5}","",Lasswade_tree$tip.label)
Lasswade_tree$tip.label<-gsub("[A-Z]{2}_[0-9]{9}","",Lasswade_tree$tip.label)
Lasswade_tree$tip.label<-gsub("\\.[1-9]{1}","",Lasswade_tree$tip.label)
Lasswade_tree$tip.label<-gsub("_"," ",Lasswade_tree$tip.label)
Lasswade_tree$tip.label<-gsub("^ ","",Lasswade_tree$tip.label)

#generating a vector of colours to colour the new virus in red
Lasswade_cols<-ifelse(str_detect(Lasswade_tree$tip.label,"Lasswade"),"#56B4E9", 
                      ifelse(str_detect(Lasswade_tree$tip.label,"Drosophila"),"#009E73","black")
                      )

Lasswade_nodes<-Lasswade_tree$node.label[str_detect(Lasswade_tree$node.label,"[0-9]{2}")] #make subset of labels containing the bootstrap values > 59
Lasswade_nodes<-Lasswade_tree$node.label %in% Lasswade_nodes[as.numeric(Lasswade_nodes)>59]

#plotting the tree
Lasswade_tree_plot<-ggtree(Lasswade_tree, ladderize = TRUE, size = 2.3) + 
  geom_tiplab(size=12.75, color=Lasswade_cols, hjust = -0.02)+
  geom_nodepoint(aes(subset=Lasswade_nodes), color="#0072B2", size=7.5)+
  geom_treescale(color = "grey", width = 0.5, x=0, y=9, linesize = 4, fontsize = 0)+
  xlim(0, 12)

pdf(file = "../../../../thesis/ch2_metagenomic_virus_discovery/Lasswade_tree_plot.pdf",height = 7.5,width = 15)
Lasswade_tree_plot
dev.off()

##Dalkeith
Dalkeith_tree<-read.tree(file = "Dalkeith_virus/02_21/Dalkeith_virus_5BLASTp.aln.fas.treefile")
Dalkeith_tree<-midpoint.root(Dalkeith_tree)
Dalkeith_tree_info<-fortify(Dalkeith_tree)

##tidying the labels on the tree, removing underscores
#removing accession numbers and underscores - this is ugly code - fix!
Dalkeith_tree$tip.label<-gsub("[A-Z]{3}[0-9]{5}","",Dalkeith_tree$tip.label)
Dalkeith_tree$tip.label<-gsub("[A-Z]{2}_[0-9]{9}","",Dalkeith_tree$tip.label)
Dalkeith_tree$tip.label<-gsub("\\.[1-9]{1}","",Dalkeith_tree$tip.label)
Dalkeith_tree$tip.label<-gsub("_"," ",Dalkeith_tree$tip.label)
Dalkeith_tree$tip.label<-gsub("^ ","",Dalkeith_tree$tip.label)

#generating a vector of colours to colour the new virus in red
Dalkeith_cols<-ifelse(str_detect(Dalkeith_tree$tip.label,"Dalkeith"),"#56B4E9","black")

Dalkeith_nodes<-Dalkeith_tree$node.label[str_detect(Dalkeith_tree$node.label,"[0-9]{2}")] #make subset of labels containing the bootstrap values > 59
Dalkeith_nodes<-Dalkeith_tree$node.label %in% Dalkeith_nodes[as.numeric(Dalkeith_nodes)>59]

#plotting the tree
Dalkeith_tree_plot<-ggtree(Dalkeith_tree, ladderize = TRUE, size = 2.3) + 
  geom_tiplab(size=12.75, color=Dalkeith_cols, hjust = -0.02)+
  geom_nodepoint(aes(subset=Dalkeith_nodes), color="#0072B2", size=7.5)+
  geom_treescale(color = "grey", width = 0.5, x=0, y=6.5, linesize = 4, fontsize = 0)+
  xlim(0, 19)

pdf(file = "../../../../thesis/ch2_metagenomic_virus_discovery/Dalkeith_tree_plot.pdf",height = 6,width = 13)
Dalkeith_tree_plot
dev.off()

##################################
## Positive sense ssRNA viruses ##
##################################

##Cockenzie
Cockenzie_tree<-read.tree(file = "Picornaviruses/Cockenzie_virus/02_21/Cockenzie_virus_5BLASTp.aln.fas.treefile")
Cockenzie_tree<-midpoint.root(Cockenzie_tree)

Cockenzie_tree_info<-fortify(Cockenzie_tree)

##tidying the labels on the tree, removing underscores
#removing accession numbers and underscores - this is ugly code - fix!
Cockenzie_tree$tip.label<-gsub("[A-Z]{3}[0-9]{5}","",Cockenzie_tree$tip.label)
Cockenzie_tree$tip.label<-gsub("[A-Z]{2}_[0-9]{6}","",Cockenzie_tree$tip.label)
Cockenzie_tree$tip.label<-gsub("\\.[1-9]{1}","",Cockenzie_tree$tip.label)
Cockenzie_tree$tip.label<-gsub("_"," ",Cockenzie_tree$tip.label)
Cockenzie_tree$tip.label<-gsub("^ ","",Cockenzie_tree$tip.label)
Cockenzie_tree$tip.label<-gsub("  DCV variant","",Cockenzie_tree$tip.label)

#generating a vector of colours to colour the new virus in red
Cockenzie_cols<-ifelse(str_detect(Cockenzie_tree$tip.label,"Cockenzie"),"#56B4E9", 
                       ifelse(str_detect(Cockenzie_tree$tip.label,"Drosophila"),"#009E73","black")) #to get label order to label with black and red colours
Cockenzie_nodes<-Cockenzie_tree$node.label[str_detect(Cockenzie_tree$node.label,"[0-9]{2}")] #make subset of labels containing the bootstrap values > 59
Cockenzie_nodes<-Cockenzie_tree$node.label %in% Cockenzie_nodes[as.numeric(Cockenzie_nodes)>59]

#plotting the tree
Cockenzie_tree_plot<-ggtree(Cockenzie_tree, ladderize = TRUE, size = 2.3) + 
                      geom_tiplab(size=12.75, color=Cockenzie_cols, hjust = -0.02)+
                      geom_nodepoint(aes(subset=Cockenzie_nodes), color="#0072B2", size=7.5)+
                      geom_treescale(color = "grey", width = 0.5, x=0, y=7.5, linesize = 4, fontsize = 0)+
                      xlim(0, 6)

pdf(file = "../../../../thesis/ch2_metagenomic_virus_discovery/Cockenzie_tree_plot.pdf",height = 6,width = 16)
Cockenzie_tree_plot
dev.off()

##Crammond + Negevirus tree
Nege_tree<-read.tree(file = "Negev_like_viruses/02_21/Negeviruses_5BLASTp.aln.fas.treefile")
Nege_tree<-midpoint.root(Nege_tree)

Nege_tree_info<-fortify(Nege_tree)

##tidying the labels on the tree, removing underscores
#removing accession numbers and underscores - this is ugly code - fix!
Nege_tree$tip.label<-gsub("[A-Z]{3}[0-9]{5}","",Nege_tree$tip.label)
Nege_tree$tip.label<-gsub("[A-Z]{2}_[0-9]{9}","",Nege_tree$tip.label)
Nege_tree$tip.label<-gsub("[A-Z]{2}_[0-9]{6}","",Nege_tree$tip.label)
Nege_tree$tip.label<-gsub("\\.[1-9]{1}","",Nege_tree$tip.label)
Nege_tree$tip.label<-gsub("_"," ",Nege_tree$tip.label)
Nege_tree$tip.label<-gsub("^ ","",Nege_tree$tip.label)

#generating a vector of colours to colour the new virus in red
Nege_cols<-ifelse(str_detect(Nege_tree$tip.label,"Crammond|Hillwood|Almond"),"#56B4E9", 
                  ifelse(str_detect(Nege_tree$tip.label,"Drosophila"),"#009E73","black")) #to get label order to label with black and red colours
Nege_nodes<-Nege_tree$node.label[str_detect(Nege_tree$node.label,"[0-9]{2}")] #make subset of labels containing the bootstrap values > 59
Nege_nodes<-Nege_tree$node.label %in% Nege_nodes[as.numeric(Nege_nodes)>59]

#plotting the tree
Nege_tree_plot<-ggtree(Nege_tree, ladderize = TRUE, size = 2.3) + 
  geom_tiplab(size=12.75, color=Nege_cols, hjust = -0.02)+
  geom_nodepoint(aes(subset=Nege_nodes), color="#0072B2", size=7.5)+
  geom_treescale(color = "grey", width = 0.5, x=0, y=18, linesize = 4, fontsize = 0)+
  xlim(0, 7)

pdf(file = "../../../../thesis/ch2_metagenomic_virus_discovery/Nege_tree_plot.pdf",height = 13,width = 14)
Nege_tree_plot
dev.off()

##Midmar tombusvirus
Midmar_tree<-read.tree(file = "Dansoman_like_tombus_virus/02_21/Midmar_virus_5BLASTp.aln.fas.treefile")
Midmar_tree<-midpoint.root(Midmar_tree)
Midmar_tree_info<-fortify(Midmar_tree)

##tidying the labels on the tree, removing underscores
#removing accession numbers and underscores - this is ugly code - fix!
Midmar_tree$tip.label<-gsub("[A-Z]{3}[0-9]{5}","",Midmar_tree$tip.label)
Midmar_tree$tip.label<-gsub("[A-Z]{2}_[0-9]{9}","",Midmar_tree$tip.label)
Midmar_tree$tip.label<-gsub("[A-Z]{2}_[0-9]{6}","",Midmar_tree$tip.label)
Midmar_tree$tip.label<-gsub("\\.[1-9]{1}","",Midmar_tree$tip.label)
Midmar_tree$tip.label<-gsub("_"," ",Midmar_tree$tip.label)
Midmar_tree$tip.label<-gsub("^ ","",Midmar_tree$tip.label)
Midmar_tree$tip.label<-gsub("Dansoman","Drosophila Dansoman",Midmar_tree$tip.label)

#generating a vector of colours to colour the new virus in red
Midmar_cols<-ifelse(str_detect(Midmar_tree$tip.label,"Midmar"),"#56B4E9", 
                    ifelse(str_detect(Midmar_tree$tip.label,"Drosophila"),"#009E73","black")) #to get label order to label with black and red colours
Midmar_nodes<-Midmar_tree$node.label[str_detect(Midmar_tree$node.label,"[0-9]{2}")] #make subset of labels containing the bootstrap values > 59
Midmar_nodes<-Midmar_tree$node.label %in% Midmar_nodes[as.numeric(Midmar_nodes)>59]

#plotting the tree
Midmar_tree_plot<-ggtree(Midmar_tree, ladderize = TRUE, size = 2.3) + 
  geom_tiplab(size=12.75, color=Midmar_cols, hjust = -0.02)+
  geom_nodepoint(aes(subset=Midmar_nodes), color="#0072B2", size=7.5)+
  geom_treescale(color = "grey", width = 0.5, x=0, y=7.5, linesize = 4, fontsize = 0)+
  xlim(0, 12)

pdf(file = "../../../../thesis/ch2_metagenomic_virus_discovery/Midmar_tree_plot.pdf",height = 6,width = 14)
Midmar_tree_plot
dev.off()

##Gosford narnavirus
Gosford_tree<-read.tree(file = "Narnavirus/02_21/Gosford_narnavirus_5BLASTp.aln.fas.treefile")
Gosford_tree<-midpoint.root(Gosford_tree)
Gosford_tree_info<-fortify(Gosford_tree)

##tidying the labels on the tree, removing underscores
#removing accession numbers and underscores - this is ugly code - fix!
Gosford_tree$tip.label<-gsub("[A-Z]{3}[0-9]{5}","",Gosford_tree$tip.label)
Gosford_tree$tip.label<-gsub("[A-Z]{2}_[0-9]{9}","",Gosford_tree$tip.label)
Gosford_tree$tip.label<-gsub("[A-Z]{2}_[0-9]{6}","",Gosford_tree$tip.label)
Gosford_tree$tip.label<-gsub("\\.[1-9]{1}","",Gosford_tree$tip.label)
Gosford_tree$tip.label<-gsub("_"," ",Gosford_tree$tip.label)
Gosford_tree$tip.label<-gsub("^ ","",Gosford_tree$tip.label)

#generating a vector of colours to colour the new virus in red
Gosford_cols<-ifelse(str_detect(Gosford_tree$tip.label,"Gosford"),"#56B4E9","black") #to get label order to label with black and red colours
Gosford_nodes<-Gosford_tree$node.label[str_detect(Gosford_tree$node.label,"[0-9]{2}")] #make subset of labels containing the bootstrap values > 59
Gosford_nodes<-Gosford_tree$node.label %in% Gosford_nodes[as.numeric(Gosford_nodes)>59]

#plotting the tree
Gosford_tree_plot<-ggtree(Gosford_tree, ladderize = TRUE, size = 2.3) + 
  geom_tiplab(size=12.75, color=Gosford_cols, hjust = -0.02)+
  geom_nodepoint(aes(subset=Gosford_nodes), color="#0072B2", size=7.5)+
  geom_treescale(color = "grey", width = 0.5, x=0, y=7.5, linesize = 4, fontsize = 0)+
  xlim(0, 20)

pdf(file = "../../../../thesis/ch2_metagenomic_virus_discovery/Gosford_tree_plot.pdf",height = 5,width = 23)
Gosford_tree_plot
dev.off()

###################
## dsRNA viruses ##
###################

##Totiviruses
Toti_tree<-read.tree(file = "Toti_viruses/02_21/Totivirus_5BLASTp.aln.fas.treefile")
Toti_tree<-midpoint.root(Toti_tree)

Toti_tree_info<-fortify(Toti_tree)

##tidying the labels on the tree, removing underscores
#removing accession numbers and underscores - this is ugly code - fix!
Toti_tree$tip.label<-gsub("[A-Z]{3}[0-9]{5}","",Toti_tree$tip.label)
Toti_tree$tip.label<-gsub("[A-Z]{2}_[0-9]{9}","",Toti_tree$tip.label)
Toti_tree$tip.label<-gsub("[A-Z]{2}_[0-9]{6}","",Toti_tree$tip.label)
Toti_tree$tip.label<-gsub("\\.[1-9]{1}","",Toti_tree$tip.label)
Toti_tree$tip.label<-gsub("_"," ",Toti_tree$tip.label)
Toti_tree$tip.label<-gsub("^ ","",Toti_tree$tip.label)
Toti_tree$tip.label<-gsub(" genome type A","",Toti_tree$tip.label)

#generating a vector of colours to colour the new virus in red
Toti_cols<-ifelse(str_detect(Toti_tree$tip.label,"Craighall|Inverleith"),"#56B4E9", 
                  ifelse(str_detect(Toti_tree$tip.label,"Drosophila"),"#009E73",
                         ifelse(str_detect(Toti_tree$tip.label,"Trichinella"),"#F0E442","black")
                  )
) #to get label order to label with black and red colours
Toti_nodes<-Toti_tree$node.label[str_detect(Toti_tree$node.label,"[0-9]{2}")] #make subset of labels containing the bootstrap values > 59
Toti_nodes<-Toti_tree$node.label %in% Toti_nodes[as.numeric(Toti_nodes)>59]

#plotting the tree
Toti_tree_plot<-ggtree(Toti_tree, ladderize = TRUE, size = 2.3) + 
  geom_tiplab(size=12.75, color=Toti_cols, hjust = -0.02)+
  geom_nodepoint(aes(subset=Toti_nodes), color="#0072B2", size=7.5)+
  geom_treescale(color = "grey", width = 0.5, x=0, y=7.5, linesize = 4, fontsize = 0)+
  xlim(0, 16)

pdf(file = "../../../../thesis/ch2_metagenomic_virus_discovery/Toti_tree_plot.pdf",height = 8.5,width = 14)
Toti_tree_plot
dev.off()

###Reoviruses
Reo_tree<-read.tree(file = "Reoviruses/02_21/Reoviruses_5BLASTp.aln.fas.treefile")
Reo_tree<-midpoint.root(Reo_tree)
Reo_tree_info<-fortify(Reo_tree)

##tidying the labels on the tree, removing underscores
#removing accession numbers and underscores - this is ugly code - fix!
Reo_tree$tip.label<-gsub("[A-Z]{3}[0-9]{5}","",Reo_tree$tip.label)
Reo_tree$tip.label<-gsub("[A-Z]{2}_[0-9]{9}","",Reo_tree$tip.label)
Reo_tree$tip.label<-gsub("[A-Z]{2}_[0-9]{6}","",Reo_tree$tip.label)
Reo_tree$tip.label<-gsub("\\.[1-9]{1}","",Reo_tree$tip.label)
Reo_tree$tip.label<-gsub("_"," ",Reo_tree$tip.label)
Reo_tree$tip.label<-gsub("^ ","",Reo_tree$tip.label)
Reo_tree$tip.label<-gsub(" strain JKT-6423","",Reo_tree$tip.label)

#generating a vector of colours to colour the new virus in red
Reo_cols<-ifelse(str_detect(Reo_tree$tip.label,"Vogrie|Glencorse"),"#56B4E9", 
                  ifelse(str_detect(Reo_tree$tip.label,"Drosophila"),"#009E73","black")
                  )
Reo_nodes<-Reo_tree$node.label[str_detect(Reo_tree$node.label,"[0-9]{2}")] #make subset of labels containing the bootstrap values > 59
Reo_nodes<-Reo_tree$node.label %in% Reo_nodes[as.numeric(Reo_nodes)>59]

#plotting the tree
Reo_tree_plot<-ggtree(Reo_tree, ladderize = TRUE, size = 2.3) + 
  geom_tiplab(size=12.75, color=Reo_cols, hjust = -0.02)+
  geom_nodepoint(aes(subset=Reo_nodes), color="#0072B2", size=7.5)+
  geom_treescale(color = "grey", width = 0.5, x=0, y=16.5, linesize = 4, fontsize = 0)+
  xlim(0, 16)

pdf(file = "../../../../thesis/ch2_metagenomic_virus_discovery/Reo_tree_plot.pdf",height = 12,width = 14)
Reo_tree_plot
dev.off()
