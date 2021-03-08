####Playing with gggenes for making genome arrow plots

#######################################################
##To plot genome annotation of viruses using gggenes ##
#######################################################

##Tools for drawing gene arrow maps in ggplot2 (R)
##good documentation online (updated june-2019)
library(gggenes)
library(plyr)
library(dplyr)
library(tidyr)
library(scales)
library(ggplot2)

setwd("C:/Users/s1667991/Dropbox/PhD - 1st Year/Sequencing/05_19_totalRNASeq/genome_presentation")

###Combining the negative and positive sense viruses into a single plot so they are all on the same axis 
pns_genomes_newviruses <- read.csv("gene_positions_new_viruses_2.csv", header = TRUE)
pns_subgenes_newviruses <- read.csv("subgene_positions_new_viruses_2.csv", header = TRUE)

glimpse(pns_genomes_newviruses)
glimpse(pns_subgenes_newviruses)

pns_genomes_newviruses$start <- as.numeric(as.character(pns_genomes_newviruses$start))
pns_genomes_newviruses$end <- as.numeric(as.character(pns_genomes_newviruses$end))
pns_subgenes_newviruses$from <- as.numeric(as.character(pns_subgenes_newviruses$from))
pns_subgenes_newviruses$to <- as.numeric(as.character(pns_subgenes_newviruses$to))
pns_subgenes_newviruses$start <- as.numeric(as.character(pns_subgenes_newviruses$start))
pns_subgenes_newviruses$end <- as.numeric(as.character(pns_subgenes_newviruses$end))

# The palette:
#Now using the prot_cat col in subgene csv to colour the ORFs by cat of protein produced if known/ predicted
cbPalette <- c("glycoprotein"="#E69F00","other_protein"="#56B4E9","replicase"="#009E73","capsid_protein"="#F0E442","structural_protein"="#0072B2","nucleoprotein"="#D55E00")

##plotting segments/viruses as arrows, and then genes as subgene regions to illustrate ORFs

#for a pdf - needs to be quite tall...
pdf(width = 8.3, height = 14.7, file = "MW_05_19_seq_newvirusgenomes_gggene_3.pdf")

ggplot(pns_genomes_newviruses, aes(xmin = start, xmax = end, y = molecule)) +
  facet_wrap(~ molecule, scales = "free_y", ncol = 1) +
  geom_gene_arrow(fill = "white",arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"),arrow_body_height = unit(9, "mm")) +
  geom_subgene_arrow(data = pns_subgenes_newviruses,
                     aes(xmin = start, xmax = end, y = molecule, fill = prot_cat,
                         xsubmin = from, xsubmax = to), arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"),arrow_body_height = unit(7, "mm"), color="black", alpha=.7)+
  scale_fill_manual(values=cbPalette)+
  theme(panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size=9,colour = "black", face = "bold"),
        axis.text.x = element_text(size=12,colour = "black", face = "bold"),
        axis.title.x = element_text("Genome position",size = 14, face = "bold"),
        axis.title.y = element_blank(),
        axis.line.x = element_line(colour = "black", size = 1),
        axis.ticks.x = element_line(colour = "black", size = 1),
        strip.text = element_blank(),
        strip.background = element_blank())+
  geom_subgene_label(data = pns_subgenes_newviruses,
                     aes(xsubmin = from, xsubmax = to, label = subgene),
                     padding.y = grid::unit(0, "mm"),
                     align = "left",
                     grow = FALSE,
                     min.size = 0)
dev.off()

##And now the extended known virus genomes
pns_genomes_knownviruses <- read.csv("gene_positions_known_viruses.csv", header = TRUE)
pns_subgenes_knownviruses <- read.csv("subgene_positions_known_viruses.csv", header = TRUE)

glimpse(pns_genomes_knownviruses)
glimpse(pns_subgenes_knownviruses)

pns_genomes_knownviruses$start <- as.numeric(as.character(pns_genomes_knownviruses$start))
pns_genomes_knownviruses$end <- as.numeric(as.character(pns_genomes_knownviruses$end))
pns_subgenes_knownviruses$from <- as.numeric(as.character(pns_subgenes_knownviruses$from))
pns_subgenes_knownviruses$to <- as.numeric(as.character(pns_subgenes_knownviruses$to))
pns_subgenes_knownviruses$start <- as.numeric(as.character(pns_subgenes_knownviruses$start))
pns_subgenes_knownviruses$end <- as.numeric(as.character(pns_subgenes_knownviruses$end))

# The palette:
#Now using the prot_cat col in subgene csv to colour the ORFs by cat of protein produced if known/ predicted
cbPalette_known <- c("glycoprotein"="#E69F00","other_protein"="#56B4E9","replicase"="#009E73","capsid_protein"="#F0E442","structural_protein"="#0072B2","nucleoprotein"="#D55E00","polyprotein"="#CC79A7","polymerase_associated_protein"="#004949")

#for a pdf - needs to be quite tall...
pdf(width = 8.3, height = 10.7, file = "MW_05_19_seq_knownvirusgenomes_gggene.pdf")

ggplot(pns_genomes_knownviruses, aes(xmin = start, xmax = end, y = molecule)) +
  facet_wrap(~ molecule, scales = "free_y", ncol = 1) +
  geom_gene_arrow(fill = "white",arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"),arrow_body_height = unit(9, "mm")) +
  geom_subgene_arrow(data = pns_subgenes_knownviruses,
                     aes(xmin = start, xmax = end, y = molecule, fill = prot_cat,
                         xsubmin = from, xsubmax = to), arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"),arrow_body_height = unit(7, "mm"), color="black", alpha=.7)+
  scale_fill_manual(values=cbPalette_known)+
  theme(panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size=9,colour = "black", face = "bold"),
        axis.text.x = element_text(size=12,colour = "black", face = "bold"),
        axis.title.x = element_text("Genome position",size = 14, face = "bold"),
        axis.title.y = element_blank(),
        axis.line.x = element_line(colour = "black", size = 1),
        axis.ticks.x = element_line(colour = "black", size = 1),
        strip.text = element_blank(),
        strip.background = element_blank())
#+
  #geom_subgene_label(data = pns_subgenes_knownviruses,
                     #aes(xsubmin = from, xsubmax = to, label = subgene),
                     #padding.y = grid::unit(0, "mm"),
                     #align = "left",
                     #grow = FALSE,
                     #min.size = 0)
dev.off()


#Things to update - labelling of viruses, get the segments onto the same line, and re-name the proteins, genes and subgenes so that no "_"s appear in the fig
#currently editing in inkscape but could work out some way to get the plot to be generated like this