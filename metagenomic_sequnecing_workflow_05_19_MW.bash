####De novo Assembly of virus genomes - Sequencing from Edinburgh (Sep 16-Oct 18) - Collections by Megan Wallace#############
#############################################################################################################################

##########################
### 
### Trinity - Assembly
###
##########################

#Creating link between raw read data and 0_reads files for the three new datasets (2017 & 2018 data)

#in 0_reads file of DrosMS_JulOct17 
cd /mnt/drive1-6tb/Megan/DrosMS_JulOct17/0_reads
		ln -s ../190805_A00291_0195_AHKLMLDMXX_2_11793RL0001L01_1.fastq.gz DrosMS_JulOct17.1.fastq.gz 
		ln -s ../190805_A00291_0195_AHKLMLDMXX_2_11793RL0001L01_2.fastq.gz DrosMS_JulOct17.2.fastq.gz 

#to go within the 0_reads file of DrosMS_AprJun18 		
cd /mnt/drive1-6tb/Megan/DrosMS_AprJun18/0_reads	
		ln -s ../190805_A00291_0195_AHKLMLDMXX_2_11793RL0002L01_1.fastq.gz DrosMS_AprJun18.1.fastq.gz 
		ln -s ../190805_A00291_0195_AHKLMLDMXX_2_11793RL0002L01_2.fastq.gz DrosMS_AprJun18.2.fastq.gz  

#to go within the 0_reads file of DrosMS_JulOct18		
cd /mnt/drive1-6tb/Megan/DrosMS_JulOct18/0_reads	
		ln -s ../190805_A00291_0195_AHKLMLDMXX_2_11793RL0003L01_1.fastq.gz DrosMS_JulOct18.1.fastq.gz 
		ln -s ../190805_A00291_0195_AHKLMLDMXX_2_11793RL0003L01_2.fastq.gz DrosMS_JulOct18.2.fastq.gz 

#Also created soft link between new reads file in the SeptOct2016 folder under the directory 10_19_reanalyses, which will contain all re-analysis of this data usiung the new cd hit filtering
cd /mnt/drive1-6tb/Megan/SeptOct2016/10_19_reanalyses/0_reads
		ln -s /mnt/drive1-6tb/Megan/SeptOct2016/Raw_Data/170203_K00166_0176_BHGKTTBBXX_1_TP-D7-008_TP-D5-004_1.fastq.gz SeptOct2016.1.fastq.gz
		ln -s /mnt/drive1-6tb/Megan/SeptOct2016/Raw_Data/170203_K00166_0176_BHGKTTBBXX_1_TP-D7-008_TP-D5-004_2.fastq.gz SeptOct2016.2.fastq.gz

#And created soft link between new reads file in the Dec2016June2017 folder under the directory 10_19_reanalyses, which will contain all re-analysis of this data usiung the new cd hit filtering
cd /mnt/drive1-6tb/Megan/Dec2016June2017/10_19_reanalyses/0_reads
		ln -s /mnt/drive1-6tb/Megan/Dec2016June2017/raw/171024_K00300_0057_BHM2W3BBXX_1_TP-D7-005_TP-D5-007_1.fastq.gz Dec2016June2017.A.1.fastq.gz
		ln -s /mnt/drive1-6tb/Megan/Dec2016June2017/raw/171024_K00300_0057_BHM2W3BBXX_1_TP-D7-005_TP-D5-007_2.fastq.gz Dec2016June2017.A.2.fastq.gz
		ln -s /mnt/drive1-6tb/Megan/Dec2016June2017/raw/171024_K00300_0057_BHM2W3BBXX_2_TP-D7-005_TP-D5-007_1.fastq.gz Dec2016June2017.B.1.fastq.gz
		ln -s /mnt/drive1-6tb/Megan/Dec2016June2017/raw/171024_K00300_0057_BHM2W3BBXX_2_TP-D7-005_TP-D5-007_2.fastq.gz Dec2016June2017.B.2.fastq.gz
		ln -s /mnt/drive1-6tb/Megan/Dec2016June2017/raw/171024_K00300_0057_BHM2W3BBXX_3_TP-D7-005_TP-D5-007_1.fastq.gz Dec2016June2017.C.1.fastq.gz
		ln -s /mnt/drive1-6tb/Megan/Dec2016June2017/raw/171024_K00300_0057_BHM2W3BBXX_3_TP-D7-005_TP-D5-007_2.fastq.gz Dec2016June2017.C.2.fastq.gz

# Running Trinity with default parameters.
	# Trinity uses a fixed k-mer value of 25. 

# Creating a test file from the fastq.qz files to test running of Trinity - can't compress symbolic links so working with original files 

cd /localdisk/home/s1667991/MeganMetaGenomics/DrosMS_JulOct17
#using head to create file with only the first 200 sequences - 800 lines
zcat 190805_A00291_0195_AHKLMLDMXX_2_11793RL0001L01_1.fastq.gz | head -800 | gzip > test_DrosMS_JulOct17.1.fastq.gz 
zcat 190805_A00291_0195_AHKLMLDMXX_2_11793RL0001L01_2.fastq.gz | head -800 | gzip > test_DrosMS_JulOct17.2.fastq.gz
mv test_DrosMS_JulOct17.1.fastq.gz 0_reads/
mv test_DrosMS_JulOct17.2.fastq.gz 0_reads/

#Running the test file through Trinity to assess any errors in script
cd /localdisk/home/s1667991/MeganMetaGenomics
Trinity --trimmomatic --seqType fq --max_memory 20G --left DrosMS_JulOct17/0_reads/test_DrosMS_JulOct17.1.fastq.gz --right DrosMS_JulOct17/0_reads/test_DrosMS_JulOct17.2.fastq.gz --CPU 20 --output DrosMS_JulOct17/1_trinity/test_DrosMS_JulOct17_Trinity_10_1
##Works fine apart from final assembly of clusters of reads...because there's not enough reads

cd /mnt/drive1-6tb/Megan/	
		#DrosMS_JulOct17
		Trinity --trimmomatic --seqType fq --max_memory 20G --left DrosMS_JulOct17/0_reads/DrosMS_JulOct17.1.fastq.gz --right DrosMS_JulOct17/0_reads/DrosMS_JulOct17.2.fastq.gz --CPU 20 --output DrosMS_JulOct17/1_trinity/DrosMS_JulOct17_Trinity_10_19
		#DrosMS_AprJun18
		Trinity --trimmomatic --seqType fq --max_memory 20G --left DrosMS_AprJun18/0_reads/DrosMS_AprJun18.1.fastq.gz --right DrosMS_AprJun18/0_reads/DrosMS_AprJun18.2.fastq.gz --CPU 20 --output DrosMS_AprJun18/1_trinity/DrosMS_AprJun18_Trinity_10_19
		#DrosMS_JulOct18
		Trinity --trimmomatic --seqType fq --max_memory 20G --left DrosMS_JulOct18/0_reads/DrosMS_JulOct18.1.fastq.gz --right DrosMS_JulOct18/0_reads/DrosMS_JulOct18.2.fastq.gz --CPU 20 --output DrosMS_JulOct18/1_trinity/DrosMS_JulOct18_Trinity_10_19
		#SeptOct16 - permission denied to this file - changed permissions for user to rwX using chmod (chmod -v -R u+rwX dir)
		Trinity --trimmomatic --seqType fq --max_memory 20G --left SeptOct2016/10_19_reanalyses/0_reads/SeptOct2016.1.fastq.gz --right SeptOct2016/10_19_reanalyses/0_reads/SeptOct2016.2.fastq.gz --CPU 20 --output SeptOct2016/10_19_reanalyses/1_trinity/SeptOct2016_Trinity_10_19
		#Dec2016June2017 - three sets of fastq.gz for this single sequencing run....have combined output
		Trinity --trimmomatic --seqType fq --max_memory 20G --left Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.A.1.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.B.1.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.C.1.fastq.gz --right Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.A.2.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.B.2.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.C.2.fastq.gz --CPU 20 --output Dec2016June2017/10_19_reanalyses/1_trinity/Dec2016June2017_Trinity_10_19
		
#############################################################################
### 								
### 						Filtering with CD-HIT						
###																		
#############################################################################

#http://www.bioinformatics.org/cd-hit/cd-hit-user-guide.pdf

# Options
# -o output filename, required.
# -c sequence identity threshold, default 0.9, this is the default cd-hit's "global sequence identity" calculated as: number of identical amino acids or bases in alignment divided by the full length of the shorter sequence.
# -p 1 or 0, default 0, if set to 1, print alignment overlap in .clstr file.
# -d length of description in .clstr file, default 20, if set to 0, it takes the fasta defline and stops at first space.
# -b band_width of alignment, default 20.
# -T number of threads, default 1; with 0, all CPUs will be used.

##Example from tymeks code
	#cd-hit-est -o cdhit -c 0.85 -i Trinity.fasta -p 1 -d 0 -b 3 -T 10 
##Example from the cd-hit user guide 
	#cd-hit-est -i est_human -o est_human95 -c 0.95 -n 8 

cd /mnt/drive1-6tb/Megan
	cd-hit-est -o DrosMS_JulOct17/2_cd_hit/DrosMS_JulOct17_cd_hit -c 0.85 -i DrosMS_JulOct17/1_trinity/DrosMS_JulOct17_Trinity_10_19/Trinity.fasta -p 1 -d 10 -b 3 -T 10 
	
	cd-hit-est -o DrosMS_AprJun18/2_cd_hit/DrosMS_AprJun18_cd_hit -c 0.85 -i DrosMS_AprJun18/1_trinity/DrosMS_AprJun18_Trinity_10_19/Trinity.fasta -p 1 -d 10 -b 3 -T 10
	
	cd-hit-est -o DrosMS_JulOct18/2_cd_hit/DrosMS_JulOct18_cd_hit -c 0.85 -i DrosMS_JulOct18/1_trinity/DrosMS_JulOct18_Trinity_10_19/Trinity.fasta -p 1 -d 10 -b 3 -T 10

	cd-hit-est -o SeptOct2016/10_19_reanalyses/2_cd_hit/SeptOct2016_cd_hit -c 0.85 -i SeptOct2016/10_19_reanalyses/1_trinity/SeptOct2016_Trinity_10_19/Trinity.fasta -p 1 -d 10 -b 3 -T 10
	
	cd-hit-est -o Dec2016June2017/10_19_reanalyses/2_cd_hit/Dec2016June2017_cd_hit -c 0.85 -i Dec2016June2017/10_19_reanalyses/1_trinity/Dec2016June2017_Trinity_10_19/Trinity.fasta -p 1 -d 10 -b 3 -T 10

#-c 0.85 = chosen to keep the sequence identity threshold at 85%, as even the most similar viruses in my dataset - eg. prestney burn and grom, or the nora viruses...still have only ~71% identity at the nucleotide level
#-p 1 = printing alignment in the .clstr file
#-d 10 (default 20 = most detailed) = length of description in the clstr file, changed from tymeks 0...as wanted detial on the clusters created - maybe not neccesary?
#-b 3 = band_width of alignment, default = 20 - kept on 3 but maybe could increase? 
#-T 10 = kept to 10, to use 10 CPUs

######################################################
######################################################
######################################################
#
# Re-format post assembly and call translations
#
######################################################
######################################################
######################################################

#formatting into unwrapped fasta, replacing seq header with seq run name and removing spaces

##################DrosMS_JulOct17########################
#Copy the output files and re-format them. Format into unwrapped fasta; remove sequence header after path; swap spaces for _'s; swap TRINITY for Sequencing run name
cat /mnt/drive1-6tb/Megan/DrosMS_JulOct17/2_cd_hit/DrosMS_JulOct17_cd_hit | fasta_formatter | sed 's/ path.\+//g' | sed s'/ /_/g' | sed 's/TRINITY/DrosMS_JulOct17/g' > /mnt/drive1-6tb/Megan/DrosMS_JulOct17/3_fasta/DrosMS_JulOct17_cd_hit.fas
##################DrosMS_AprJun18########################
cat /mnt/drive1-6tb/Megan/DrosMS_AprJun18/2_cd_hit/DrosMS_AprJun18_cd_hit | fasta_formatter | sed 's/ path.\+//g' | sed s'/ /_/g' | sed 's/TRINITY/DrosMS_AprJun18/g' > /mnt/drive1-6tb/Megan/DrosMS_AprJun18/3_fasta/DrosMS_AprJun18_cd_hit.fas
##################DrosMS_JulOct18########################
cat /mnt/drive1-6tb/Megan/DrosMS_JulOct18/2_cd_hit/DrosMS_JulOct18_cd_hit | fasta_formatter | sed 's/ path.\+//g' | sed s'/ /_/g' | sed 's/TRINITY/DrosMS_JulOct18/g' > /mnt/drive1-6tb/Megan/DrosMS_JulOct18/3_fasta/DrosMS_JulOct18_cd_hit.fas
##################Dec2016June2017########################
cat /mnt/drive1-6tb/Megan/Dec2016June2017/10_19_reanalyses/2_cd_hit/Dec2016June2017_cd_hit | fasta_formatter | sed 's/ path.\+//g' | sed s'/ /_/g' | sed 's/TRINITY/Dec2016June2017/g' > /mnt/drive1-6tb/Megan/Dec2016June2017/10_19_reanalyses/3_fasta/Dec2016June2017_cd_hit.fas
##################SeptOct2016########################
cat /mnt/drive1-6tb/Megan/SeptOct2016/10_19_reanalyses/2_cd_hit/SeptOct2016_cd_hit | fasta_formatter | sed 's/ path.\+//g' | sed s'/ /_/g' | sed 's/TRINITY/SeptOct2016/g' > /mnt/drive1-6tb/Megan/SeptOct2016/10_19_reanalyses/3_fasta/SeptOct2016_cd_hit.fas


####and then creating concatenated protein sequences for BLAST#####

##This script creates three fasta files, 1) The DNA sequences of the longest ORF from each sequence, 2) The Protein from the longest ORF, and 3) a concatenated chimera of all protein sequences >=250 AAs (or all within 20% of the length of the longest, if none are >=250

##################DrosMS_JulOct17########################
#processing the output to creat concatenated putatative protein sequences for blasting
Rscript /localdisk/science/ProcessSPADES.Rscript /mnt/drive1-6tb/Megan/DrosMS_JulOct17/3_fasta/DrosMS_JulOct17_cd_hit.fas
##################DrosMS_AprJun18########################
Rscript /localdisk/science/ProcessSPADES.Rscript /mnt/drive1-6tb/Megan/DrosMS_AprJun18/3_fasta/DrosMS_AprJun18_cd_hit.fas
##################DrosMS_JulOct18########################
Rscript /localdisk/science/ProcessSPADES.Rscript /mnt/drive1-6tb/Megan/DrosMS_JulOct18/3_fasta/DrosMS_JulOct18_cd_hit.fas
##################Dec2016June2017########################
Rscript /localdisk/science/ProcessSPADES.Rscript /mnt/drive1-6tb/Megan/Dec2016June2017/10_19_reanalyses/3_fasta/Dec2016June2017_cd_hit.fas
####################SeptOct2016##########################
Rscript /localdisk/science/ProcessSPADES.Rscript /mnt/drive1-6tb/Megan/SeptOct2016/10_19_reanalyses/3_fasta/SeptOct2016_cd_hit.fas

##Have moved the _PROTBLASTCHIMERA.fasta and _PROT.fasta to a 4_proteins directory, newly made for each folder, and removed from 3_fasta dir

######################################################
######################################################
######################################################
#													
# 	Re-format post assembly and call translations	
#													
#####################################################
######################################################
######################################################

cd /mnt/drive1-6tb/Megan/DrosMS_JulOct17/4_proteins
#### identifying (using grep) the sequence names, and replacing the >>s with >s in the longest ORFs concatenated protein fasta
cat *.LongestORFS_PROTBLASTCHIMERA.fasta | grep --no-group-separator -B 1 '[A-WY]\{150,\}' | sed 's/>>/>/g' > DrosMS_JulOct17_LongORFs.fas
cd /mnt/drive1-6tb/Megan/DrosMS_JulOct17/3_fasta
#### and replacing the >>s with >s in the longest ORFs DNA fasta
cat *.LongestORFS_DNA.fasta | sed 's/>>/>/g' > DrosMS_JulOct17_Long_DNAs.fas

cd /mnt/drive1-6tb/Megan/DrosMS_AprJun18/4_proteins
cat *.LongestORFS_PROTBLASTCHIMERA.fasta | grep --no-group-separator -B 1 '[A-WY]\{150,\}' | sed 's/>>/>/g' > DrosMS_AprJun18_LongORFs.fas
cd /mnt/drive1-6tb/Megan/DrosMS_AprJun18/3_fasta
cat *.LongestORFS_DNA.fasta | sed 's/>>/>/g' > DrosMS_AprJun18_Long_DNAs.fas

cd /mnt/drive1-6tb/Megan/DrosMS_JulOct18/4_proteins
cat *.LongestORFS_PROTBLASTCHIMERA.fasta | grep --no-group-separator -B 1 '[A-WY]\{150,\}' | sed 's/>>/>/g' > DrosMS_JulOct18_LongORFs.fas
cd /mnt/drive1-6tb/Megan/DrosMS_JulOct18/3_fasta
cat *.LongestORFS_DNA.fasta | sed 's/>>/>/g' > DrosMS_JulOct18_Long_DNAs.fas

cd /mnt/drive1-6tb/Megan/Dec2016June2017/10_19_reanalyses/4_proteins
cat *.LongestORFS_PROTBLASTCHIMERA.fasta | grep --no-group-separator -B 1 '[A-WY]\{150,\}' | sed 's/>>/>/g' > Dec2016June2017_LongORFs.fas
cd /mnt/drive1-6tb/Megan/Dec2016June2017/10_19_reanalyses/3_fasta
cat *.LongestORFS_DNA.fasta | sed 's/>>/>/g' > Dec2016June2017_Long_DNAs.fas

cd /mnt/drive1-6tb/Megan/SeptOct2016/10_19_reanalyses/4_proteins
cat *.LongestORFS_PROTBLASTCHIMERA.fasta | grep --no-group-separator -B 1 '[A-WY]\{150,\}' | sed 's/>>/>/g' > SeptOct2016_LongORFs.fas
cd /mnt/drive1-6tb/Megan/SeptOct2016/10_19_reanalyses/3_fasta
cat *.LongestORFS_DNA.fasta | sed 's/>>/>/g' > SeptOct2016_Long_DNAs.fas
	
###################################################################
### 															
### Building a BLAST database, searching against this database 	###
###																
###################################################################

##Workflow based on script from Nathan, using diamond blastp against FlyMetagenomic_diamond database, looking for the ACCs in the virus database nr_viruses, and then joining the virus AACs to the nucleotide sequences to create a fasta file of virus like sequences
##Then creating a blast db from this file, and matching to a db of all of the known viruses to find known viruses in the metagenomic sequences

cd /mnt/drive1-6tb/Megan
#make a blast database to find those already known
#re-made this with the bits I know about, but haven't sent to darren ...Lasswade + Bunya L4
makeblastdb -in  Drosophila_Viruses_Nov2019.fas -out KnownDrosVirus -dbtype 'nucl'

######DrosMS_JulOct17#######
cd /mnt/drive1-6tb/Megan/DrosMS_JulOct17
#Using diamond to do a blastp search on the concatenated proteins found in the sequencing run (using e hit score of 0.01, -k parameter (no. of target seqs. to report alignments for) at 1, and use 20 CPUs) 
#now using new up to date database
diamond blastp  -e 0.01 -k 1 -p 20 -d /mnt/drive3-6tb/Darren/BLAST/NCBI_databases/FlyMetagenomic_diamond.dmnd -q 4_proteins/DrosMS_JulOct17_LongORFs.fas -f 6 > 5_BLAST/allSets_FlyMetagenomic.tsv 

##accessing the accession numbers from the 2nd field & uniquely sorting them
cut -f 2 5_BLAST/allSets_FlyMetagenomic.tsv | sort | uniq  > 5_BLAST/accs

#using blastdbcmd ask if those ACCs numbers are in the virus database (quick/dirty hack to see if they are viruses) - and if they are, output the results to a table
join -t $'\t' -1 1 -2 2  <(/mnt/drive3-6tb/Darren/BLAST/ncbi-blast-2.9.0+/bin/blastdbcmd -db /mnt/drive3-6tb/Darren/BLAST/NCBI_databases/nr_viruses_v5 -entry_batch 5_BLAST/accs -outfmt "%a#####%S" 2>/dev/null | sed 's/#####/\t/g' | sort -k 1,1 | uniq ) <(sort -k 2,2 5_BLAST/allSets_FlyMetagenomic.tsv) > 5_BLAST/DrosMS_JulOct17_NewVirusLikeHits.tsv

#experiment (and testing) suggests that Mimiviruses,Herpesviruses,Iridoviruses, and pox viruses are almost invariably unconvincing hits. Remember the objective is to find some viruses, not to fid all the viruses
cd /mnt/drive1-6tb/Megan/DrosMS_JulOct17
cat 5_BLAST/DrosMS_JulOct17_NewVirusLikeHits.tsv | grep -v '[Mm]imi\|[Pp]hage\|[Hh]erpes\|ourmia\|Acidianus\|phycodnavirus\|endogenous\|Marseillevirus\|Mollivirus\|phycodnavirus\|Ostreococcus\|moumouvirus\|Sulfolobus\|Pandoravirus\|Chlorella\|Phaeocystis\|megavirus\|Ectocarpus_siliculosus\|Feldmannia\|Klosneuvirus\|Lausannevirus\|Bathycoccus\|marseillevirus\|Cafeteria\|Archaea\|Carcinoma\|Acanthocystis\|[Mm]ito\|leukemia\|mamavirus\|Acanthamoeba\|Moumouvirus\|Cannes\|huxleyi\|Megavirus\|ichnovirus\|Port-Miou\|Chlorella_virus' | cut -f 3  | sort | uniq | grep -F -f - -A 1 --no-group-separator --no-filename 3_fasta/DrosMS_JulOct17_Long_DNAs.fas | paste  - - | sed 's/>//g' >  5_BLAST/DrosMS_JulOct17_VirusLikeSequences.tsv
##formatting as fasta
join -t $'\t' -1 3 -2 1  <(sort -t $'\t' -k 3,3 5_BLAST/DrosMS_JulOct17_NewVirusLikeHits.tsv) <(sort -t $'\t' -k 1,1 5_BLAST/DrosMS_JulOct17_VirusLikeSequences.tsv) | awk -F$'\t' '{print ">"$3"_"$4"%_p="$12"_"$1"\n"$14}' | sed 's/ /_/g' > 5_BLAST/DrosMS_JulOct17_NewVirusLikeHits.fas

## running a BLASTn against the new virus_like sequence file against the most up to date known drosophila virus databse tto find those viruses in this DB which are already known -- Known virus database replaced with more up to date one, with two extra viruses added from prev sequencing runs .... 
cat 5_BLAST/DrosMS_JulOct17_NewVirusLikeHits.fas | blastn -db ../KnownDrosVirus -num_threads 10 -outfmt 6 -perc_identity 90 -evalue 1e-50 | cut -f 1,2 | sort | uniq | more > 5_BLAST/DrosMS_JulOct17_already_known.tsv

##useful way of counting the number of sequences in a fasta file...eg. for the new virus like hits file
grep -c "^>" DrosMS_JulOct17_NewVirusLikeHits.fas #115

#####Trying to find the L gene of a new phasmavirus...
#closest BLAST HIT is Anopheles_triannulatus_orthophasmavirus for the M and S segments, so making a BLAST database from the L segment from this virus on BLAST, and then searching through the DNAs

makeblastdb -in MH822966.1_Anopheles_triannulatus_orthophasmavirus_isolate_AMA_segment_RNA1.fas -out Phasma_L_n -dbtype 'nucl'

cat ../3_fasta/DrosMS_JulOct17_Long_DNAs.fas | blastn -db Phasma_L_n -num_threads 10 -outfmt 6 -perc_identity 20 -evalue 1e-50 | cut -f 1,2 | sort | uniq | more > Phasma_L_n_matches.tsv
#no matches...

#trying to match the proteins - made a blastdp from the ORFs from metagenomic sequencing + queried the target p against it
makeblastdb -in ../4_proteins/DrosMS_JulOct17_LongORFs.fas -out JulOct17_LongORF -dbtype 'prot'

blastp -evalue 0.01 -num_threads 20 -max_target_seqs 15 -outfmt 7 -db JulOct17_LongORF -query QBK47217.1_RdRp_Anopheles_triannulatus_orthophasmavirus_p.fas -out Phasma_L_p_matches.txt

# BLASTP 2.7.1+
# Query: QBK47217.1 RNA-dependent RNA polymerase [Anopheles triannulatus orthophasmavirus]
# Database: JulOct17_LongORF
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 1 hits found
#QBK47217.1	DrosMS_JulOct17_DN19267_c0_g3_i1_len=492	45.732	164	88	1	363	525	1	164	1.65e-36	135
# BLAST processed 1 queries

##Only one match found...which is too small (492 bp), but does encode part of the RdRp

############################################################################################################

######DrosMS_AprJun18#######
cd /mnt/drive1-6tb/Megan/DrosMS_AprJun18
#Using diamond to do a blastp search on the concatenated proteins found in the sequencing run (using e hit score of 0.01, -k parameter (no. of target seqs. to report alignments for) at 1, and use 20 CPUs) 
#now using new up to date database
diamond blastp  -e 0.01 -k 1 -p 20 -d /mnt/drive3-6tb/Darren/BLAST/NCBI_databases/FlyMetagenomic_diamond.dmnd -q 4_proteins/DrosMS_AprJun18_LongORFs.fas -f 6 > 5_BLAST/allSets_FlyMetagenomic.tsv 

cut -f 2 5_BLAST/allSets_FlyMetagenomic.tsv | sort | uniq  > 5_BLAST/accs

#using blastdbcmd ask if those ACCs numbers are in the virus database (quick/dirty hack to see if they are viruses) - and if they are, output the results to a table
join -t $'\t' -1 1 -2 2  <(/mnt/drive3-6tb/Darren/BLAST/ncbi-blast-2.9.0+/bin/blastdbcmd -db /mnt/drive3-6tb/Darren/BLAST/NCBI_databases/nr_viruses_v5 -entry_batch 5_BLAST/accs -outfmt "%a#####%S" 2>/dev/null | sed 's/#####/\t/g' | sort -k 1,1 | uniq ) <(sort -k 2,2 5_BLAST/allSets_FlyMetagenomic.tsv) > 5_BLAST/DrosMS_AprJun18_NewVirusLikeHits.tsv

#experiment (and testing) suggests that Mimiviruses,Herpesviruses,Iridoviruses, and pox viruses are almost invariably unconvincing hits. Remember the objective is to find some viruses, not to fid all the viruses
cd /mnt/drive1-6tb/Megan/DrosMS_AprJun18
cat 5_BLAST/DrosMS_AprJun18_NewVirusLikeHits.tsv | grep -v '[Mm]imi\|[Pp]hage\|[Hh]erpes\|ourmia\|Acidianus\|phycodnavirus\|endogenous\|Marseillevirus\|Mollivirus\|phycodnavirus\|Ostreococcus\|moumouvirus\|Sulfolobus\|Pandoravirus\|Chlorella\|Phaeocystis\|megavirus\|Ectocarpus_siliculosus\|Feldmannia\|Klosneuvirus\|Lausannevirus\|Bathycoccus\|marseillevirus\|Cafeteria\|Archaea\|Carcinoma\|Acanthocystis\|[Mm]ito\|leukemia\|mamavirus\|Acanthamoeba\|Moumouvirus\|Cannes\|huxleyi\|Megavirus\|ichnovirus\|Port-Miou\|Chlorella_virus' | cut -f 3  | sort | uniq | grep -F -f - -A 1 --no-group-separator --no-filename 3_fasta/DrosMS_AprJun18_Long_DNAs.fas | paste  - - | sed 's/>//g' >  5_BLAST/DrosMS_AprJun18_VirusLikeSequences.tsv
##formatting as fasta
join -t $'\t' -1 3 -2 1  <(sort -t $'\t' -k 3,3 5_BLAST/DrosMS_AprJun18_NewVirusLikeHits.tsv) <(sort -t $'\t' -k 1,1 5_BLAST/DrosMS_AprJun18_VirusLikeSequences.tsv) | awk -F$'\t' '{print ">"$3"_"$4"%_p="$12"_"$1"\n"$14}' | sed 's/ /_/g' > 5_BLAST/DrosMS_AprJun18_NewVirusLikeHits.fas
##useful way of counting the number of sequences in a fasta file...eg. for the new virus like hits file
grep -c "^>" DrosMS_AprJun18_NewVirusLikeHits.fas #81

## running a BLASTn against the new virus_like sequence file against the most up to date known drosophila virus databse tto find those viruses in this DB which are already known -- Known virus database replaced with more up to date one, with two extra viruses added from prev sequencing runs .... 
cat 5_BLAST/DrosMS_AprJun18_NewVirusLikeHits.fas | blastn -db ../KnownDrosVirus -num_threads 10 -outfmt 6 -perc_identity 90 -evalue 1e-50 | cut -f 1,2 | sort | uniq | more > 5_BLAST/DrosMS_AprJun18_already_known.tsv

############################################################################################################

######DrosMS_JulOct18#######
cd /mnt/drive1-6tb/Megan/DrosMS_JulOct18
#Using diamond to do a blastp search on the concatenated proteins found in the sequencing run (using e hit score of 0.01, -k parameter (no. of target seqs. to report alignments for) at 1, and use 20 CPUs) 
#now using new up to date database
diamond blastp  -e 0.01 -k 1 -p 20 -d /mnt/drive3-6tb/Darren/BLAST/NCBI_databases/FlyMetagenomic_diamond.dmnd -q 4_proteins/DrosMS_JulOct18_LongORFs.fas -f 6 > 5_BLAST/allSets_FlyMetagenomic.tsv 

cut -f 2 5_BLAST/allSets_FlyMetagenomic.tsv | sort | uniq  > 5_BLAST/accs

#using blastdbcmd ask if those ACCs numbers are in the virus database (quick/dirty hack to see if they are viruses) - and if they are, output the results to a table
join -t $'\t' -1 1 -2 2  <(/mnt/drive3-6tb/Darren/BLAST/ncbi-blast-2.9.0+/bin/blastdbcmd -db /mnt/drive3-6tb/Darren/BLAST/NCBI_databases/nr_viruses_v5 -entry_batch 5_BLAST/accs -outfmt "%a#####%S" 2>/dev/null | sed 's/#####/\t/g' | sort -k 1,1 | uniq ) <(sort -k 2,2 5_BLAST/allSets_FlyMetagenomic.tsv) > 5_BLAST/DrosMS_JulOct18_NewVirusLikeHits.tsv

#experiment (and testing) suggests that Mimiviruses,Herpesviruses,Iridoviruses, and pox viruses are almost invariably unconvincing hits. Remember the objective is to find some viruses, not to fid all the viruses
cd /mnt/drive1-6tb/Megan/DrosMS_JulOct18
cat 5_BLAST/DrosMS_JulOct18_NewVirusLikeHits.tsv | grep -v '[Mm]imi\|[Pp]hage\|[Hh]erpes\|ourmia\|Acidianus\|phycodnavirus\|endogenous\|Marseillevirus\|Mollivirus\|phycodnavirus\|Ostreococcus\|moumouvirus\|Sulfolobus\|Pandoravirus\|Chlorella\|Phaeocystis\|megavirus\|Ectocarpus_siliculosus\|Feldmannia\|Klosneuvirus\|Lausannevirus\|Bathycoccus\|marseillevirus\|Cafeteria\|Archaea\|Carcinoma\|Acanthocystis\|[Mm]ito\|leukemia\|mamavirus\|Acanthamoeba\|Moumouvirus\|Cannes\|huxleyi\|Megavirus\|ichnovirus\|Port-Miou\|Chlorella_virus' | cut -f 3  | sort | uniq | grep -F -f - -A 1 --no-group-separator --no-filename 3_fasta/DrosMS_JulOct18_Long_DNAs.fas | paste  - - | sed 's/>//g' >  5_BLAST/DrosMS_JulOct18_VirusLikeSequences.tsv
##formatting as fasta
join -t $'\t' -1 3 -2 1  <(sort -t $'\t' -k 3,3 5_BLAST/DrosMS_JulOct18_NewVirusLikeHits.tsv) <(sort -t $'\t' -k 1,1 5_BLAST/DrosMS_JulOct18_VirusLikeSequences.tsv) | awk -F$'\t' '{print ">"$3"_"$4"%_p="$12"_"$1"\n"$14}' | sed 's/ /_/g' > 5_BLAST/DrosMS_JulOct18_NewVirusLikeHits.fas
##useful way of counting the number of sequences in a fasta file...eg. for the new virus like hits file
grep -c "^>" DrosMS_JulOct18_NewVirusLikeHits.fas #135

## running a BLASTn against the new virus_like sequence file against the most up to date known drosophila virus databse tto find those viruses in this DB which are already known -- Known virus database replaced with more up to date one, with two extra viruses added from prev sequencing runs .... 
cat 5_BLAST/DrosMS_JulOct18_NewVirusLikeHits.fas | blastn -db ../KnownDrosVirus -num_threads 10 -outfmt 6 -perc_identity 90 -evalue 1e-50 | cut -f 1,2 | sort | uniq | more > 5_BLAST/DrosMS_JulOct18_already_known.tsv

############################################################################################################

######SeptOct2016#######
cd /mnt/drive1-6tb/Megan/SeptOct2016/10_19_reanalyses
#Using diamond to do a blastp search on the concatenated proteins found in the sequencing run (using e hit score of 0.01, -k parameter (no. of target seqs. to report alignments for) at 1, and use 20 CPUs) 
#now using new up to date database
diamond blastp  -e 0.01 -k 1 -p 20 -d /mnt/drive3-6tb/Darren/BLAST/NCBI_databases/FlyMetagenomic_diamond.dmnd -q 4_proteins/SeptOct2016_LongORFs.fas -f 6 > 5_BLAST/allSets_FlyMetagenomic.tsv 

cut -f 2 5_BLAST/allSets_FlyMetagenomic.tsv | sort | uniq  > 5_BLAST/accs

#using blastdbcmd ask if those ACCs numbers are in the virus database (quick/dirty hack to see if they are viruses) - and if they are, output the results to a table
join -t $'\t' -1 1 -2 2  <(/mnt/drive3-6tb/Darren/BLAST/ncbi-blast-2.9.0+/bin/blastdbcmd -db /mnt/drive3-6tb/Darren/BLAST/NCBI_databases/nr_viruses_v5 -entry_batch 5_BLAST/accs -outfmt "%a#####%S" 2>/dev/null | sed 's/#####/\t/g' | sort -k 1,1 | uniq ) <(sort -k 2,2 5_BLAST/allSets_FlyMetagenomic.tsv) > 5_BLAST/SeptOct2016_NewVirusLikeHits.tsv

#experiment (and testing) suggests that Mimiviruses,Herpesviruses,Iridoviruses, and pox viruses are almost invariably unconvincing hits. Remember the objective is to find some viruses, not to fid all the viruses
cd /mnt/drive1-6tb/Megan/SeptOct2016/10_19_reanalyses
cat 5_BLAST/SeptOct2016_NewVirusLikeHits.tsv | grep -v '[Mm]imi\|[Pp]hage\|[Hh]erpes\|ourmia\|Acidianus\|phycodnavirus\|endogenous\|Marseillevirus\|Mollivirus\|phycodnavirus\|Ostreococcus\|moumouvirus\|Sulfolobus\|Pandoravirus\|Chlorella\|Phaeocystis\|megavirus\|Ectocarpus_siliculosus\|Feldmannia\|Klosneuvirus\|Lausannevirus\|Bathycoccus\|marseillevirus\|Cafeteria\|Archaea\|Carcinoma\|Acanthocystis\|[Mm]ito\|leukemia\|mamavirus\|Acanthamoeba\|Moumouvirus\|Cannes\|huxleyi\|Megavirus\|ichnovirus\|Port-Miou\|Chlorella_virus' | cut -f 3  | sort | uniq | grep -F -f - -A 1 --no-group-separator --no-filename 3_fasta/SeptOct2016_Long_DNAs.fas | paste  - - | sed 's/>//g' >  5_BLAST/SeptOct2016_VirusLikeSequences.tsv
##formatting as fasta
join -t $'\t' -1 3 -2 1  <(sort -t $'\t' -k 3,3 5_BLAST/SeptOct2016_NewVirusLikeHits.tsv) <(sort -t $'\t' -k 1,1 5_BLAST/SeptOct2016_VirusLikeSequences.tsv) | awk -F$'\t' '{print ">"$3"_"$4"%_p="$12"_"$1"\n"$14}' | sed 's/ /_/g' > 5_BLAST/SeptOct2016_NewVirusLikeHits.fas
##useful way of counting the number of sequences in a fasta file...eg. for the new virus like hits file
grep -c "^>" SeptOct2016_NewVirusLikeHits.fas #178

## running a BLASTn against the new virus_like sequence file against the most up to date known drosophila virus databse tto find those viruses in this DB which are already known -- Known virus database replaced with more up to date one, with two extra viruses added from prev sequencing runs .... 
cat 5_BLAST/SeptOct2016_NewVirusLikeHits.fas | blastn -db ../../KnownDrosVirus -num_threads 10 -outfmt 6 -perc_identity 90 -evalue 1e-50 | cut -f 1,2 | sort | uniq | more > 5_BLAST/SeptOct2016_already_known.tsv

############################################################################################################

######Dec2016June2017#######
cd /mnt/drive1-6tb/Megan/Dec2016June2017/10_19_reanalyses
#Using diamond to do a blastp search on the concatenated proteins found in the sequencing run (using e hit score of 0.01, -k parameter (no. of target seqs. to report alignments for) at 1, and use 20 CPUs) 
#now using new up to date database
diamond blastp  -e 0.01 -k 1 -p 20 -d /mnt/drive3-6tb/Darren/BLAST/NCBI_databases/FlyMetagenomic_diamond.dmnd -q 4_proteins/Dec2016June2017_LongORFs.fas -f 6 > 5_BLAST/allSets_FlyMetagenomic.tsv 

cut -f 2 5_BLAST/allSets_FlyMetagenomic.tsv | sort | uniq  > 5_BLAST/accs

#using blastdbcmd ask if those ACCs numbers are in the virus database (quick/dirty hack to see if they are viruses) - and if they are, output the results to a table
join -t $'\t' -1 1 -2 2  <(/mnt/drive3-6tb/Darren/BLAST/ncbi-blast-2.9.0+/bin/blastdbcmd -db /mnt/drive3-6tb/Darren/BLAST/NCBI_databases/nr_viruses_v5 -entry_batch 5_BLAST/accs -outfmt "%a#####%S" 2>/dev/null | sed 's/#####/\t/g' | sort -k 1,1 | uniq ) <(sort -k 2,2 5_BLAST/allSets_FlyMetagenomic.tsv) > 5_BLAST/Dec2016June2017_NewVirusLikeHits.tsv

#experiment (and testing) suggests that Mimiviruses,Herpesviruses,Iridoviruses, and pox viruses are almost invariably unconvincing hits. Remember the objective is to find some viruses, not to fid all the viruses
cd /mnt/drive1-6tb/Megan/Dec2016June2017/10_19_reanalyses
cat 5_BLAST/Dec2016June2017_NewVirusLikeHits.tsv | grep -v '[Mm]imi\|[Pp]hage\|[Hh]erpes\|ourmia\|Acidianus\|phycodnavirus\|endogenous\|Marseillevirus\|Mollivirus\|phycodnavirus\|Ostreococcus\|moumouvirus\|Sulfolobus\|Pandoravirus\|Chlorella\|Phaeocystis\|megavirus\|Ectocarpus_siliculosus\|Feldmannia\|Klosneuvirus\|Lausannevirus\|Bathycoccus\|marseillevirus\|Cafeteria\|Archaea\|Carcinoma\|Acanthocystis\|[Mm]ito\|leukemia\|mamavirus\|Acanthamoeba\|Moumouvirus\|Cannes\|huxleyi\|Megavirus\|ichnovirus\|Port-Miou\|Chlorella_virus' | cut -f 3  | sort | uniq | grep -F -f - -A 1 --no-group-separator --no-filename 3_fasta/Dec2016June2017_Long_DNAs.fas | paste  - - | sed 's/>//g' >  5_BLAST/Dec2016June2017_VirusLikeSequences.tsv
##formatting as fasta
join -t $'\t' -1 3 -2 1  <(sort -t $'\t' -k 3,3 5_BLAST/Dec2016June2017_NewVirusLikeHits.tsv) <(sort -t $'\t' -k 1,1 5_BLAST/Dec2016June2017_VirusLikeSequences.tsv) | awk -F$'\t' '{print ">"$3"_"$4"%_p="$12"_"$1"\n"$14}' | sed 's/ /_/g' > 5_BLAST/Dec2016June2017_NewVirusLikeHits.fas
##useful way of counting the number of sequences in a fasta file...eg. for the new virus like hits file
grep -c "^>" Dec2016June2017_NewVirusLikeHits.fas #173

## running a BLASTn against the new virus_like sequence file against the most up to date known drosophila virus databse tto find those viruses in this DB which are already known -- Known virus database replaced with more up to date one, with two extra viruses added from prev sequencing runs .... 
cat 5_BLAST/Dec2016June2017_NewVirusLikeHits.fas | blastn -db ../../KnownDrosVirus -num_threads 10 -outfmt 6 -perc_identity 90 -evalue 1e-50 | cut -f 1,2 | sort | uniq | more > 5_BLAST/Dec2016June2017_already_known.tsv

####################################################################
###Searching for new viruses across all viral sequencing datasets###
####################################################################

##Trying to identify reads of the new viruses (8) in the rest of the viral sequencing datasets using BLASTn

##made a fasta file containing the new virus sequences (whole genomes or segments) to make into a database (/mnt/drive1-6tb/Megan/11_19_put_new_viruses.fas)

###starting by running a BLASTn of these as queries against the already known viruses and saving the results as a text file...incase any of them are viruses we already know about
cd /mnt/drive1-6tb/Megan
cat 11_19_put_new_viruses.fas | blastn -db KnownDrosVirus -num_threads 10 -outfmt 6 -perc_identity 90 -evalue 1e-50 | cut -f 1,2 | sort | uniq | more > 11_19_put_new_viruses_already_known.tsv

##Results - only the potential new picornavirus matched to some bits of an unnamed discistrovirus in the known virus database
# Picronavirus_DrosMS_JulOct18_Drosophila_Nasuta_cf_goose_dicistrovirus_81%_72.5%_cf_Emperyat_virus_76.4% DrosophilaNasuta_len=1102_cf_Empeyrat_virus_at_76.4%
# Picronavirus_DrosMS_JulOct18_Drosophila_Nasuta_cf_goose_dicistrovirus_81%_72.5%_cf_Emperyat_virus_76.4% DrosophilaNasuta_len=2822_cf_Goose_dicistrovirus_at_72.5%
# Picronavirus_DrosMS_JulOct18_Drosophila_Nasuta_cf_goose_dicistrovirus_81%_72.5%_cf_Emperyat_virus_76.4% DrosophilaNasuta_len=3191_cf_Goose_dicistrovirus_at_81%

##now checking for the new viruses across all virus databases 
##made db from the newviruslikehit.fas files from all 5 of the sequencing datasets
cd /mnt/drive1-6tb/Megan/
makeblastdb -in SeptOct2016/10_19_reanalyses/5_BLAST/SeptOct2016_NewVirusLikeHits.fas -out SeptOct2016_hits_db -dbtype 'nucl'
makeblastdb -in Dec2016June2017/10_19_reanalyses/5_BLAST/Dec2016June2017_NewVirusLikeHits.fas -out Dec2016June2017_hits_db -dbtype 'nucl'
makeblastdb -in Dec2016June2017/10_19_reanalyses/5_BLAST/Dec2016June2017_NewVirusLikeHits.fas -out Dec2016June2017_hits_db -dbtype 'nucl'
makeblastdb -in DrosMS_JulOct17/5_BLAST/DrosMS_JulOct17_NewVirusLikeHits.fas -out DrosMS_JulOct17_hits_db -dbtype 'nucl'
makeblastdb -in DrosMS_AprJun18/5_BLAST/DrosMS_AprJun18_NewVirusLikeHits.fas -out DrosMS_AprJun18_hits_db -dbtype 'nucl'
makeblastdb -in DrosMS_JulOct18/5_BLAST/DrosMS_JulOct18_NewVirusLikeHits.fas -out DrosMS_JulOct18_hits_db -dbtype 'nucl'

##BLASTn - using the .fas file w the putative new viruses in it as a querie file agains each of the databases
cd /mnt/drive1-6tb/Megan
cat 11_19_put_new_viruses.fas | blastn -db SeptOct2016_hits_db -num_threads 10 -outfmt 6 -perc_identity 90 -evalue 1e-50 | cut -f 1,2 | sort | uniq | more > SeptOct2016/10_19_reanalyses/5_BLAST/11_19_put_new_viruses_SeptOct2016_matches.tsv
cat 11_19_put_new_viruses.fas | blastn -db Dec2016June2017_hits_db -num_threads 10 -outfmt 6 -perc_identity 90 -evalue 1e-50 | cut -f 1,2 | sort | uniq | more > Dec2016June2017/10_19_reanalyses/5_BLAST/11_19_put_new_viruses_Dec2016June2017_matches.tsv
cat 11_19_put_new_viruses.fas | blastn -db DrosMS_JulOct17_hits_db -num_threads 10 -outfmt 6 -perc_identity 90 -evalue 1e-50 | cut -f 1,2 | sort | uniq | more > DrosMS_JulOct17/5_BLAST/11_19_put_new_viruses_DrosMS_JulOct17_matches.tsv
cat 11_19_put_new_viruses.fas | blastn -db DrosMS_AprJun18_hits_db -num_threads 10 -outfmt 6 -perc_identity 90 -evalue 1e-50 | cut -f 1,2 | sort | uniq | more > DrosMS_AprJun18/5_BLAST/11_19_put_new_viruses_DrosMS_AprJun18_matches.tsv
cat 11_19_put_new_viruses.fas | blastn -db DrosMS_JulOct18_hits_db -num_threads 10 -outfmt 6 -perc_identity 90 -evalue 1e-50 | cut -f 1,2 | sort | uniq | more > DrosMS_JulOct18/5_BLAST/11_19_put_new_viruses_DrosMS_JulOct18_matches.tsv

###Doing directed searches to try to increase the completeness of the new found virus genomes

##Making protein databases from the hits files of all the 5 databases
cd /mnt/drive1-6tb/Megan/
makeblastdb -in SeptOct2016/10_19_reanalyses/4_proteins/SeptOct2016_LongORFs.fas -out SeptOct2016_hits_prot_db/SeptOct2016_hits_prot_db -dbtype 'prot'
makeblastdb -in Dec2016June2017/10_19_reanalyses/4_proteins/Dec2016June2017_LongORFs.fas -out Dec2016June2017_hits_prot_db/Dec2016June2017_hits_prot_db -dbtype 'prot'
makeblastdb -in DrosMS_JulOct17/4_proteins/DrosMS_JulOct17_LongORFs.fas -out DrosMS_JulOct17_hits_prot_db/DrosMS_JulOct17_hits_prot_db -dbtype 'prot'
makeblastdb -in DrosMS_AprJun18/4_proteins/DrosMS_AprJun18_LongORFs.fas -out DrosMS_AprJun18_hits_prot_db/DrosMS_AprJun18_hits_prot_db -dbtype 'prot'
makeblastdb -in DrosMS_JulOct18/4_proteins/DrosMS_JulOct18_LongORFs.fas -out DrosMS_JulOct18_hits_prot_db/DrosMS_JulOct18_hits_prot_db -dbtype 'prot'

##Inveresk virus
##closest BLAST hit from tBLASTx, complete n sequence from genbank = MG550265.1 Soybean cyst nematode nyami-like virus
##BLASTn
cd /mnt/drive1-6tb/Megan
cat SeptOct2016/10_19_reanalyses/5_BLAST/MG550265.1_Soybean_cyst_nematode_nyami-like_virus.fas | blastn -db SeptOct2016_hits_db -num_threads 10 -outfmt 6 -perc_identity 20 -evalue 1e-50 | cut -f 1,2 | sort | uniq | more > SeptOct2016/10_19_reanalyses/5_BLAST/Inveresk_close_blastn_SeptOct2016_matches.tsv
##no matches

##BLASTx
cd /mnt/drive1-6tb/Megan
blastx -evalue 0.01 -num_threads 20 -max_target_seqs 15 -outfmt 7 -db SeptOct2016_hits_prot_db/SeptOct2016_hits_prot_db -query SeptOct2016/10_19_reanalyses/5_BLAST/MG550265.1_Soybean_cyst_nematode_nyami-like_virus.fas -out SeptOct2016/10_19_reanalyses/5_BLAST/Inveresk_close_tblastx_SeptOct2016_matches.tsv
##Only one useful hit, but only an extension of ~80 bp to the each end of the genome, and w a lot of runs of Ts so not sure I trust it, codes protein is exactly the same length


##New Goose discistrovirus + emperyat like picornavirus (only in JulOct18 I think)
##Made fasta of closest BLASTn hits, eg. the TSA bit Darren found + bundaberg bee virus 
cd /mnt/drive1-6tb/Megan
##BLASTn
cat DrosMS_JulOct18/5_BLAST/picorna_bits_for_searching.fasta | blastn -db DrosMS_JulOct18_hits_db -num_threads 10 -outfmt 6 -perc_identity 20 -evalue 1e-50 | cut -f 1,2 | sort | uniq | more > DrosMS_JulOct18/5_BLAST/Picorna_close_blastn_JulOct18_matches.tsv
##no matches
##tBLASTx
cd /mnt/drive1-6tb/Megan
blastx -evalue 0.01 -num_threads 20 -max_target_seqs 15 -outfmt 7 -db DrosMS_JulOct18_hits_prot_db/DrosMS_JulOct18_hits_prot_db -query DrosMS_JulOct18/5_BLAST/picorna_bits_for_searching.fasta -out DrosMS_JulOct18/5_BLAST/Picorna_close_tblastx_JulOct18_matches.tsv
##some matches but none of any use - all were other viruses which happen to be ~25% similar at protein level to these but I know what they actually are

#########################################################################
##To Find Read Depth + map coverage across viral genomes - new viruses ##
#########################################################################

#make a bowtie2 database which contains all the viruses you want to map coverage for that you've built/found, and the CO1 sequences for the fly species whihc may or may not be in your sequencing pool, then map the raw metagenomic F + R reads to this database, to find the coverage of all of these viruses in the reads
##maybe include COI's for the species in the fasta file of viruses?
##Have made a fasta file containing the CO1 sequences of the Drosphila sp. I think are in my datasets, and some of those not in it
bowtie2-build DistinctiveDrosophilidCOIRegion.fas Dros_CO1_bowtie_db/Dros_CO1_bowtie_db
##Made a second bowtie database containing the viruses from all 5 databases I've discovered
bowtie2-build 11_19_new_viruses_5_datasets.fasta 11_19_new_viruses_5_datasets_bowtie_db/11_19_new_viruses_5_datasets_bowtie_db

##SeptOct2016
#mapping F and R reads against the Drosophila CO1 genes, and keep only the mapped reads in a sorted bam file
#update 01/20 now using more sensistive mapping method as was concerned about cross mapping between the CO1 genes of the different Drosophila species
bowtie2 -x Dros_CO1_bowtie_db/Dros_CO1_bowtie_db --very-sensitive -p 20 -q -1 SeptOct2016/10_19_reanalyses/0_reads/SeptOct2016.1.fastq.gz -2 SeptOct2016/10_19_reanalyses/0_reads/SeptOct2016.2.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > SeptOct2016/10_19_reanalyses/6_read_mapping/Dros_CO1_SeptOct2016_read_mapping2.bam
##trying to filter out low quality mapped reads
samtools view -b -q 10 SeptOct2016/10_19_reanalyses/6_read_mapping/Dros_CO1_SeptOct2016_read_mapping2.bam > SeptOct2016/10_19_reanalyses/6_read_mapping/Dros_CO1_SeptOct2016_read_mapping2.filtered.bam
#get the read depth at each point
cd /mnt/drive1-6tb/Megan/SeptOct2016/10_19_reanalyses/6_read_mapping
genomeCoverageBed -d -ibam Dros_CO1_SeptOct2016_read_mapping2.bam > SeptOct2016_Dros_CO1_ReadDepth2.tsv
genomeCoverageBed -d -ibam Dros_CO1_SeptOct2016_read_mapping2.filtered.bam > SeptOct2016_Dros_CO1_ReadDepth2.filtered.tsv

##Dec2016June2017
cd /mnt/drive1-6tb/Megan
#mapping F and R reads against the Drosophila CO1 genes, and keep only the mapped reads in a sorted bam file
bowtie2 -x Dros_CO1_bowtie_db/Dros_CO1_bowtie_db --very-sensitive -p 20 -q -1 Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.A.1.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.B.1.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.C.1.fastq.gz -2 Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.A.2.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.B.2.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.C.2.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > Dec2016June2017/10_19_reanalyses/6_read_mapping/Dros_CO1_Dec2016June2017_read_mapping2.bam
##trying to filter out low quality mapped reads
samtools view -b -q 10 Dec2016June2017/10_19_reanalyses/6_read_mapping/Dros_CO1_Dec2016June2017_read_mapping2.bam > Dec2016June2017/10_19_reanalyses/6_read_mapping/Dros_CO1_Dec2016June2017_read_mapping2.filtered.bam
#get the read depth at each point
cd /mnt/drive1-6tb/Megan/Dec2016June2017/10_19_reanalyses/6_read_mapping
genomeCoverageBed -d -ibam Dros_CO1_Dec2016June2017_read_mapping2.bam > Dec2016June2017_Dros_CO1_ReadDepth2.tsv
genomeCoverageBed -d -ibam Dros_CO1_Dec2016June2017_read_mapping2.filtered.bam > Dec2016June2017_Dros_CO1_ReadDepth2.filtered.tsv

##DrosMS_JulOct17
cd /mnt/drive1-6tb/Megan
#mapping F and R reads against the Drosophila CO1 genes, and keep only the mapped reads in a sorted bam file
bowtie2 -x Dros_CO1_bowtie_db/Dros_CO1_bowtie_db --very-sensitive -p 20 -q -1 DrosMS_JulOct17/0_reads/DrosMS_JulOct17.1.fastq.gz -2 DrosMS_JulOct17/0_reads/DrosMS_JulOct17.2.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_JulOct17/6_read_mapping/Dros_CO1_DrosMS_JulOct17_read_mapping2.bam
##trying to filter out low quality mapped reads
samtools view -b -q 10 DrosMS_JulOct17/6_read_mapping/Dros_CO1_DrosMS_JulOct17_read_mapping2.bam > DrosMS_JulOct17/6_read_mapping/Dros_CO1_DrosMS_JulOct17_read_mapping2.filtered.bam
#get the read depth at each point
cd /mnt/drive1-6tb/Megan/DrosMS_JulOct17/6_read_mapping
genomeCoverageBed -d -ibam Dros_CO1_DrosMS_JulOct17_read_mapping2.bam > DrosMS_JulOct17_Dros_CO1_ReadDepth2.tsv
genomeCoverageBed -d -ibam Dros_CO1_DrosMS_JulOct17_read_mapping2.filtered.bam > DrosMS_JulOct17_Dros_CO1_ReadDepth2.filtered.tsv
cd /mnt/drive1-6tb/Megan

##DrosMS_AprJun18
cd /mnt/drive1-6tb/Megan
#mapping F and R reads against the Drosophila CO1 genes, and keep only the mapped reads in a sorted bam file
bowtie2 -x Dros_CO1_bowtie_db/Dros_CO1_bowtie_db --very-sensitive -p 20 -q -1 DrosMS_AprJun18/0_reads/DrosMS_AprJun18.1.fastq.gz -2 DrosMS_AprJun18/0_reads/DrosMS_AprJun18.2.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_AprJun18/6_read_mapping/Dros_CO1_DrosMS_AprJun18_read_mapping2.bam
##trying to filter out low quality mapped reads
samtools view -b -q 10 DrosMS_AprJun18/6_read_mapping/Dros_CO1_DrosMS_AprJun18_read_mapping2.bam > DrosMS_AprJun18/6_read_mapping/Dros_CO1_DrosMS_AprJun18_read_mapping2.filtered.bam
#get the read depth at each point
cd /mnt/drive1-6tb/Megan/DrosMS_AprJun18/6_read_mapping
genomeCoverageBed -d -ibam Dros_CO1_DrosMS_AprJun18_read_mapping2.bam > DrosMS_AprJun18_Dros_CO1_ReadDepth2.tsv
genomeCoverageBed -d -ibam Dros_CO1_DrosMS_AprJun18_read_mapping2.filtered.bam > DrosMS_AprJun18_Dros_CO1_ReadDepth2.filtered.tsv
cd /mnt/drive1-6tb/Megan

##DrosMS_JulOct18
cd /mnt/drive1-6tb/Megan
#mapping F and R reads against the Drosophila CO1 genes, and keep only the mapped reads in a sorted bam file
bowtie2 -x Dros_CO1_bowtie_db/Dros_CO1_bowtie_db --very-sensitive -p 20 -q -1 DrosMS_JulOct18/0_reads/DrosMS_JulOct18.1.fastq.gz -2 DrosMS_JulOct18/0_reads/DrosMS_JulOct18.2.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_JulOct18/6_read_mapping/Dros_CO1_DrosMS_JulOct18_read_mapping2.bam
#filtering out low quality mapped reads
samtools view -b -q 10 DrosMS_JulOct18/6_read_mapping/Dros_CO1_DrosMS_JulOct18_read_mapping2.bam > DrosMS_JulOct18/6_read_mapping/Dros_CO1_DrosMS_JulOct18_read_mapping2.filtered.bam
#get the read depth at each point
cd /mnt/drive1-6tb/Megan/DrosMS_JulOct18/6_read_mapping
genomeCoverageBed -d -ibam Dros_CO1_DrosMS_JulOct18_read_mapping2.bam > DrosMS_JulOct18_Dros_CO1_ReadDepth2.tsv
genomeCoverageBed -d -ibam Dros_CO1_DrosMS_JulOct18_read_mapping2.filtered.bam > DrosMS_JulOct18_Dros_CO1_ReadDepth2.filtered.tsv
cd /mnt/drive1-6tb/Megan

#SeptOct2016
#mapping F + R reads agains the newly described genomes from all 5 datasets to look at coverage across the genomes
cd /mnt/drive1-6tb/Megan
bowtie2 -x 11_19_new_viruses_5_datasets_bowtie_db/11_19_new_viruses_5_datasets_bowtie_db -k 1 -p 20 -q -1 SeptOct2016/10_19_reanalyses/0_reads/SeptOct2016.1.fastq.gz -2 SeptOct2016/10_19_reanalyses/0_reads/SeptOct2016.2.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > SeptOct2016/10_19_reanalyses/6_read_mapping/11_19_new_viruses_5_datasets_SeptOct2016_read_mapping.bam
#get the read depth at each point
cd /mnt/drive1-6tb/Megan/SeptOct2016/10_19_reanalyses/6_read_mapping
genomeCoverageBed -d -ibam 11_19_new_viruses_5_datasets_SeptOct2016_read_mapping.bam > SeptOct2016_11_19_new_viruses_5_datasets_ReadDepth.tsv
#Dec2016June2017
cd /mnt/drive1-6tb/Megan
#mapping F + R reads agains the newly described genomes from all 5 datasets to look at coverage across the genomes
bowtie2 -x 11_19_new_viruses_5_datasets_bowtie_db/11_19_new_viruses_5_datasets_bowtie_db -k 1 -p 20 -q -1 Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.A.1.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.B.1.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.C.1.fastq.gz -2 Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.A.2.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.B.2.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.C.2.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > Dec2016June2017/10_19_reanalyses/6_read_mapping/11_19_new_viruses_5_datasets_Dec2016June2017_read_mapping.bam
#get the read depth at each point
cd /mnt/drive1-6tb/Megan/Dec2016June2017/10_19_reanalyses/6_read_mapping
genomeCoverageBed -d -ibam 11_19_new_viruses_5_datasets_Dec2016June2017_read_mapping.bam > Dec2016June2017_11_19_new_viruses_5_datasets_ReadDepth.tsv
#DrosMS_JulOct17
#mapping F + R reads agains the newly described genomes from all 5 datasets to look at coverage across the genomes
cd /mnt/drive1-6tb/Megan
bowtie2 -x 11_19_new_viruses_5_datasets_bowtie_db/11_19_new_viruses_5_datasets_bowtie_db -k 1 -p 20 -q -1 DrosMS_JulOct17/0_reads/DrosMS_JulOct17.1.fastq.gz -2 DrosMS_JulOct17/0_reads/DrosMS_JulOct17.2.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_JulOct17/6_read_mapping/11_19_new_viruses_5_datasets_DrosMS_JulOct17_read_mapping.bam
#get the read depth at each point
cd /mnt/drive1-6tb/Megan/DrosMS_JulOct17/6_read_mapping
genomeCoverageBed -d -ibam 11_19_new_viruses_5_datasets_DrosMS_JulOct17_read_mapping.bam > JulOct17_11_19_new_viruses_5_datasets_ReadDepth.tsv
#DrosMS_AprJun18
#mapping F + R reads agains the newly described genomes from all 5 datasets to look at coverage across the genomes
cd /mnt/drive1-6tb/Megan
bowtie2 -x 11_19_new_viruses_5_datasets_bowtie_db/11_19_new_viruses_5_datasets_bowtie_db -k 1 -p 20 -q -1 DrosMS_AprJun18/0_reads/DrosMS_AprJun18.1.fastq.gz -2 DrosMS_AprJun18/0_reads/DrosMS_AprJun18.2.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_AprJun18/6_read_mapping/11_19_new_viruses_5_datasets_DrosMS_AprJun18_read_mapping.bam
#get the read depth at each point
cd /mnt/drive1-6tb/Megan/DrosMS_AprJun18/6_read_mapping
genomeCoverageBed -d -ibam 11_19_new_viruses_5_datasets_DrosMS_AprJun18_read_mapping.bam > AprJun18_11_19_new_viruses_5_datasets_ReadDepth.tsv
gan
#DrosMS_JulOct18
#mapping F + R reads agains the newly described genomes from all 5 datasets to look at coverage across the genomes
cd /mnt/drive1-6tb/Megan
bowtie2 -x 11_19_new_viruses_5_datasets_bowtie_db/11_19_new_viruses_5_datasets_bowtie_db -k 1 -p 20 -q -1 DrosMS_JulOct18/0_reads/DrosMS_JulOct18.1.fastq.gz -2 DrosMS_JulOct18/0_reads/DrosMS_JulOct18.2.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_JulOct18/6_read_mapping/11_19_new_viruses_5_datasets_DrosMS_JulOct18_read_mapping.bam
#get the read depth at each point
cd /mnt/drive1-6tb/Megan/DrosMS_JulOct18/6_read_mapping
genomeCoverageBed -d -ibam 11_19_new_viruses_5_datasets_DrosMS_JulOct18_read_mapping.bam > JulOct18_11_19_new_viruses_5_datasets_ReadDepth.tsv

##################################################################################################
##investigating effect of --very-sensitive read mapping on mapping of reads to new virus genomes##
##################################################################################################

#SeptOct2016
cd /mnt/drive1-6tb/Megan
bowtie2 -x 11_19_new_viruses_5_datasets_bowtie_db/11_19_new_viruses_5_datasets_bowtie_db --very-sensitive -p 20 -q -1 SeptOct2016/10_19_reanalyses/0_reads/SeptOct2016.1.fastq.gz -2 SeptOct2016/10_19_reanalyses/0_reads/SeptOct2016.2.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > SeptOct2016/10_19_reanalyses/6_read_mapping/11_19_new_viruses_5_datasets_SeptOct2016_read_mapping_sensitive.bam
#get the read depth at each point
cd /mnt/drive1-6tb/Megan/SeptOct2016/10_19_reanalyses/6_read_mapping
genomeCoverageBed -d -ibam 11_19_new_viruses_5_datasets_SeptOct2016_read_mapping_sensitive.bam > SeptOct2016_11_19_new_viruses_5_datasets_ReadDepth_sensitive.tsv
#Dec2016June2017
cd /mnt/drive1-6tb/Megan
#mapping F + R reads agains the newly described genomes from all 5 datasets to look at coverage across the genomes
bowtie2 -x 11_19_new_viruses_5_datasets_bowtie_db/11_19_new_viruses_5_datasets_bowtie_db --very-sensitive -p 20 -q -1 Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.A.1.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.B.1.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.C.1.fastq.gz -2 Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.A.2.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.B.2.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.C.2.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > Dec2016June2017/10_19_reanalyses/6_read_mapping/11_19_new_viruses_5_datasets_Dec2016June2017_read_mapping_sensitive.bam
#get the read depth at each point
cd /mnt/drive1-6tb/Megan/Dec2016June2017/10_19_reanalyses/6_read_mapping
genomeCoverageBed -d -ibam 11_19_new_viruses_5_datasets_Dec2016June2017_read_mapping_sensitive.bam > Dec2016June2017_11_19_new_viruses_5_datasets_ReadDepth_sensitive.tsv
#DrosMS_JulOct17
#mapping F + R reads agains the newly described genomes from all 5 datasets to look at coverage across the genomes
cd /mnt/drive1-6tb/Megan
bowtie2 -x 11_19_new_viruses_5_datasets_bowtie_db/11_19_new_viruses_5_datasets_bowtie_db --very-sensitive -p 20 -q -1 DrosMS_JulOct17/0_reads/DrosMS_JulOct17.1.fastq.gz -2 DrosMS_JulOct17/0_reads/DrosMS_JulOct17.2.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_JulOct17/6_read_mapping/11_19_new_viruses_5_datasets_DrosMS_JulOct17_read_mapping_sensitive.bam
#get the read depth at each point
cd /mnt/drive1-6tb/Megan/DrosMS_JulOct17/6_read_mapping
genomeCoverageBed -d -ibam 11_19_new_viruses_5_datasets_DrosMS_JulOct17_read_mapping_sensitive.bam > JulOct17_11_19_new_viruses_5_datasets_ReadDepth_sensitive.tsv
#DrosMS_AprJun18
#mapping F + R reads agains the newly described genomes from all 5 datasets to look at coverage across the genomes
cd /mnt/drive1-6tb/Megan
bowtie2 -x 11_19_new_viruses_5_datasets_bowtie_db/11_19_new_viruses_5_datasets_bowtie_db --very-sensitive -p 20 -q -1 DrosMS_AprJun18/0_reads/DrosMS_AprJun18.1.fastq.gz -2 DrosMS_AprJun18/0_reads/DrosMS_AprJun18.2.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_AprJun18/6_read_mapping/11_19_new_viruses_5_datasets_DrosMS_AprJun18_read_mapping_sensitive.bam
#get the read depth at each point
cd /mnt/drive1-6tb/Megan/DrosMS_AprJun18/6_read_mapping
genomeCoverageBed -d -ibam 11_19_new_viruses_5_datasets_DrosMS_AprJun18_read_mapping_sensitive.bam > AprJun18_11_19_new_viruses_5_datasets_ReadDepth_sensitive.tsv
#DrosMS_JulOct18
#mapping F + R reads agains the newly described genomes from all 5 datasets to look at coverage across the genomes
cd /mnt/drive1-6tb/Megan
bowtie2 -x 11_19_new_viruses_5_datasets_bowtie_db/11_19_new_viruses_5_datasets_bowtie_db --very-sensitive -p 20 -q -1 DrosMS_JulOct18/0_reads/DrosMS_JulOct18.1.fastq.gz -2 DrosMS_JulOct18/0_reads/DrosMS_JulOct18.2.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_JulOct18/6_read_mapping/11_19_new_viruses_5_datasets_DrosMS_JulOct18_read_mapping_sensitive.bam
#get the read depth at each point
cd /mnt/drive1-6tb/Megan/DrosMS_JulOct18/6_read_mapping
genomeCoverageBed -d -ibam 11_19_new_viruses_5_datasets_DrosMS_JulOct18_read_mapping_sensitive.bam > JulOct18_11_19_new_viruses_5_datasets_ReadDepth_sensitive.tsv

##########################################################
## Finding depth of coverage across known viral genomes ##
##########################################################

#making a bowtie2 database which contains all the viruses to look for and map coverage over in the datasets 
cd mnt/drive1-6tb/Megan
mkdir Drosophila_Viruses_Nov2019_bowtie_db
bowtie2-build Drosophila_Viruses_Nov2019.fas Drosophila_Viruses_Nov2019_bowtie_db/Drosophila_Viruses_Nov2019_bowtie_db

##mapping the raw metagenomic F + R reads to this bowtie2 database, to find the coverage of all of these viruses in the reads - and number of total reads in order to identify which virus is in which dataset

###############
#SeptOct2016###
###############

#mapping F + R reads agains the newly described genomes from all 5 datasets to look at coverage across the genomes - currently using -k -1 mapping, as expecting more divergence in the virus mapping than the Drosophila genes, and to stay consistent with what I did for the new viruses - but could change to --very-sensistive, to filter for quality with samtools if need be
cd /mnt/drive1-6tb/Megan
bowtie2 -x Drosophila_Viruses_Nov2019_bowtie_db/Drosophila_Viruses_Nov2019_bowtie_db -k 1 -p 20 -q -1 SeptOct2016/10_19_reanalyses/0_reads/SeptOct2016.1.fastq.gz -2 SeptOct2016/10_19_reanalyses/0_reads/SeptOct2016.2.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > SeptOct2016/10_19_reanalyses/6_read_mapping/Known_Dros_Viruses_Nov2019_SeptOct2016_read_mapping.bam
##this is taking an awful long time to map...maybe I should only map reads to the viruses I've already identified as being present in the datasets using diamond searches against a BLAST database. 

#examined the files SeptOct2016_already_known.tsv
## want to 1) extract second col (names of viruses with similarity to contigs in the dataset)
## 2) remove replicates so only unique sequence names from the known virus contigs remain
cd /mnt/drive1-6tb/Megan
mkdir SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_already_known_viruses
cut -f2 SeptOct2016/10_19_reanalyses/5_BLAST/SeptOct2016_already_known.tsv | sort | uniq > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_already_known_viruses/SeptOct2016_already_known_unique_viruses.tsv
## 3) match unique names to the sequences with the same titles in the Drosophila_viruses_2019.fas file - to create a fas file with the sequences I want to map reads to
##making a tsv file version of Drosophila_Viruses_Nov2019
cd /mnt/drive1-6tb/Megan
#first linearising the file so sequences aren't multi-line any more, then putting IDs on separate lines, and tab delimiting - then removing ^M carriage returns
awk 'NR==1{printf $0"\t";next}{printf /^>/ ? "\n"$0"\t" : $0}' Drosophila_Viruses_Nov2019/Drosophila_Viruses_Nov2019.fas |  awk -F"\t" '{gsub("^>","",$0);print $0}' | sed -e "s/\r//g" > Drosophila_Viruses_Nov2019/Drosophila_Viruses_Nov2019.tsv
##joining by the common fields in the first cols of both tsv files, then formatting as fasta
join -t $'\t' -1 1 -2 1 <(sort -t $'\t' -k 1,1 SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_already_known_viruses/SeptOct2016_already_known_unique_viruses.tsv) <(sort -t $'\t' -k 1,1 Drosophila_Viruses_Nov2019/Drosophila_Viruses_Nov2019.tsv) | awk '{print ">"$1"\n"$2}' > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_already_known_viruses/SeptOct2016_already_known_virus_seqs.fas
##checking the no. of sequences in this final fas file for the bowtie_db
grep -c "^>" SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_already_known_viruses/SeptOct2016_already_known_virus_seqs.fas #67 viruses/segments of viruses
## 4) create bowtie2 database from this named sequence fas file - which contains only the known viruses I already suspect are in this sequencing data set based on contig similarity
cd SeptOct2016/10_19_reanalyses/6_read_mapping
mkdir SeptOct2016_already_known_viruses/SeptOct2016_already_known_viruses_bowtie_db
bowtie2-build SeptOct2016_already_known_viruses/SeptOct2016_already_known_virus_seqs.fas SeptOct2016_already_known_viruses/SeptOct2016_already_known_viruses_bowtie_db/SeptOct2016_already_known_viruses_bowtie_db
## ?? - maybe I should then add the missing segments of any viruses present...incase there's not large contigs but small reads map
## 5) Map raw reads to this bowtie2 database, so that I can see coverage across the genomes for all of the known viruses present
cd /mnt/drive1-6tb/Megan
bowtie2 -x SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_already_known_viruses/SeptOct2016_already_known_viruses_bowtie_db/SeptOct2016_already_known_viruses_bowtie_db -k 1 -p 20 -q -1 SeptOct2016/10_19_reanalyses/0_reads/SeptOct2016.1.fastq.gz -2 SeptOct2016/10_19_reanalyses/0_reads/SeptOct2016.2.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_already_known_viruses/SeptOct2016_known_virus_read_mapping.bam
##optional...seeing if mapping reads with the --very-sensistive option makes a difference to results
bowtie2 -x SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_already_known_viruses/SeptOct2016_already_known_viruses_bowtie_db/SeptOct2016_already_known_viruses_bowtie_db --very-sensitive -p 20 -q -1 SeptOct2016/10_19_reanalyses/0_reads/SeptOct2016.1.fastq.gz -2 SeptOct2016/10_19_reanalyses/0_reads/SeptOct2016.2.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_already_known_viruses/SeptOct2016_known_virus_read_mapping_sensitive.bam
##optional...seeing if filtering reads by quality using q-10 makes a difference to results 
samtools view -b -q 10 SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_already_known_viruses/SeptOct2016_known_virus_read_mapping.bam > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_already_known_viruses/SeptOct2016_known_virus_read_mapping.filtered.bam
samtools view -b -q 10 SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_already_known_viruses/SeptOct2016_known_virus_read_mapping_sensitive.bam > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_already_known_viruses/SeptOct2016_known_virus_read_mapping_sensitive.filtered.bam

#get the read depth at each point
cd /mnt/drive1-6tb/Megan/SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_already_known_viruses/
genomeCoverageBed -d -ibam SeptOct2016_known_virus_read_mapping.bam > SeptOct2016_known_virus_ReadDepth.tsv
##for filtered, sensitive + filtered or sensitive
genomeCoverageBed -d -ibam SeptOct2016_known_virus_read_mapping.filtered.bam > SeptOct2016_known_virus_ReadDepth.filtered.tsv ##seems to be exactly the same as the non-filtered reads in terms of total reads per virus
genomeCoverageBed -d -ibam SeptOct2016_known_virus_read_mapping_sensitive.bam > SeptOct2016_known_virus_ReadDepth_sensitive.tsv ##less known viruses have a non-zero number of reads mapped to them, seems to increase reads mapping for some viruses - makes sense as not only first match considered
genomeCoverageBed -d -ibam SeptOct2016_known_virus_read_mapping_sensitive.filtered.bam > SeptOct2016_known_virus_ReadDepth_sensitive.filtered.tsv ##again reduces number of viruses with non-zero read numbers mapped to them. Significantly reduces number of reads mappped to each virus/segment in some cases - think this mapping may be the most imformative

##########################################################
##Dec2016June2017 - mapping reads to known virus genomes##
##########################################################

##now using --very-sensistive mapping 

cd /mnt/drive1-6tb/Megan
mkdir Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_already_known_viruses
chmod +x Dec2016June2017/10_19_reanalyses/5_BLAST/Dec2016June2017_already_known.tsv
cut -f2 Dec2016June2017/10_19_reanalyses/5_BLAST/Dec2016June2017_already_known.tsv | sort | uniq > Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_already_known_viruses/Dec2016June2017_already_known_unique_viruses.tsv
## 3) match unique names to the sequences with the same titles in the Drosophila_viruses_2019.fas file - already made tsv version of this for SeptOct2016
cd /mnt/drive1-6tb/Megan
join -t $'\t' -1 1 -2 1 <(sort -t $'\t' -k 1,1 Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_already_known_viruses/Dec2016June2017_already_known_unique_viruses.tsv) <(sort -t $'\t' -k 1,1 Drosophila_Viruses_Nov2019/Drosophila_Viruses_Nov2019.tsv) | awk '{print ">"$1"\n"$2}' > Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_already_known_viruses/Dec2016June2017_already_known_virus_seqs.fas
##checking the no. of sequences in this final fas file for the bowtie_db
grep -c "^>" Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_already_known_viruses/Dec2016June2017_already_known_virus_seqs.fas # 43 viruses/segments of viruses
## 4) create bowtie2 database from this named sequence fas file 
cd Dec2016June2017/10_19_reanalyses/6_read_mapping
mkdir Dec2016June2017_already_known_viruses/Dec2016June2017_already_known_viruses_bowtie_db
bowtie2-build Dec2016June2017_already_known_viruses/Dec2016June2017_already_known_virus_seqs.fas Dec2016June2017_already_known_viruses/Dec2016June2017_already_known_viruses_bowtie_db/Dec2016June2017_already_known_viruses_bowtie_db
## 5) Map raw reads to this bowtie2 database, so that I can see coverage across the genomes for all of the known viruses present
cd /mnt/drive1-6tb/Megan+-
bowtie2 -x Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_already_known_viruses/Dec2016June2017_already_known_viruses_bowtie_db/Dec2016June2017_already_known_viruses_bowtie_db --very-sensitive -p 20 -q -1 Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.A.1.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.B.1.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.C.1.fastq.gz -2 Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.A.2.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.B.2.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.C.2.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_already_known_viruses/Dec2016June2017_known_virus_read_mapping_sensitive.bam
##this takes quite a while to run...
#get the read depth at each point
cd /mnt/drive1-6tb/Megan/Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_already_known_viruses/
genomeCoverageBed -d -ibam Dec2016June2017_known_virus_read_mapping_sensitive.bam > Dec2016June2017_known_virus_ReadDepth_sensitive.tsv

###################
##DrosMS_JulOct17##
###################
cd /mnt/drive1-6tb/Megan
mkdir DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_already_known_viruses
chmod +x DrosMS_JulOct17/5_BLAST/DrosMS_JulOct17_already_known.tsv
cut -f2 DrosMS_JulOct17/5_BLAST/DrosMS_JulOct17_already_known.tsv | sort | uniq > DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_already_known_viruses/DrosMS_JulOct17_already_known_unique_viruses.tsv
## 3) match unique names to the sequences with the same titles in the Drosophila_viruses_2019.fas file - already made tsv version of this for SeptOct2016
cd /mnt/drive1-6tb/Megan
join -t $'\t' -1 1 -2 1 <(sort -t $'\t' -k 1,1 DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_already_known_viruses/DrosMS_JulOct17_already_known_unique_viruses.tsv) <(sort -t $'\t' -k 1,1 Drosophila_Viruses_Nov2019/Drosophila_Viruses_Nov2019.tsv) | awk '{print ">"$1"\n"$2}' > DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_already_known_viruses/DrosMS_JulOct17_already_known_virus_seqs.fas
##checking the no. of sequences in this final fas file for the bowtie_db
grep -c "^>" DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_already_known_viruses/DrosMS_JulOct17_already_known_virus_seqs.fas # 45 viruses/segments of viruses
## 4) create bowtie2 database from this named sequence fas file 
cd DrosMS_JulOct17/6_read_mapping
mkdir DrosMS_JulOct17_already_known_viruses/DrosMS_JulOct17_already_known_viruses_bowtie_db
bowtie2-build DrosMS_JulOct17_already_known_viruses/DrosMS_JulOct17_already_known_virus_seqs.fas DrosMS_JulOct17_already_known_viruses/DrosMS_JulOct17_already_known_viruses_bowtie_db/DrosMS_JulOct17_already_known_viruses_bowtie_db
## 5) Map raw reads to this bowtie2 database, so that I can see coverage across the genomes for all of the known viruses present
cd /mnt/drive1-6tb/Megan
bowtie2 -x DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_already_known_viruses/DrosMS_JulOct17_already_known_viruses_bowtie_db/DrosMS_JulOct17_already_known_viruses_bowtie_db --very-sensitive -p 20 -q -1 DrosMS_JulOct17/0_reads/DrosMS_JulOct17.1.fastq.gz -2 DrosMS_JulOct17/0_reads/DrosMS_JulOct17.2.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_already_known_viruses/DrosMS_JulOct17_known_virus_read_mapping_sensitive.bam
##this takes quite a while to run...
#get the read depth at each point
cd /mnt/drive1-6tb/Megan/DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_already_known_viruses/
genomeCoverageBed -d -ibam DrosMS_JulOct17_known_virus_read_mapping_sensitive.bam > DrosMS_JulOct17_known_virus_ReadDepth_sensitive.tsv

###################
##DrosMS_AprJun18##
###################
cd /mnt/drive1-6tb/Megan
mkdir DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_already_known_viruses
chmod +x DrosMS_AprJun18/5_BLAST/DrosMS_AprJun18_already_known.tsv
cut -f2 DrosMS_AprJun18/5_BLAST/DrosMS_AprJun18_already_known.tsv | sort | uniq > DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_already_known_viruses/DrosMS_AprJun18_already_known_unique_viruses.tsv
## 3) match unique names to the sequences with the same titles in the Drosophila_viruses_2019.fas file - already made tsv version of this for SeptOct2016
cd /mnt/drive1-6tb/Megan
join -t $'\t' -1 1 -2 1 <(sort -t $'\t' -k 1,1 DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_already_known_viruses/DrosMS_AprJun18_already_known_unique_viruses.tsv) <(sort -t $'\t' -k 1,1 Drosophila_Viruses_Nov2019/Drosophila_Viruses_Nov2019.tsv) | awk '{print ">"$1"\n"$2}' > DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_already_known_viruses/DrosMS_AprJun18_already_known_virus_seqs.fas
##checking the no. of sequences in this final fas file for the bowtie_db
grep -c "^>" DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_already_known_viruses/DrosMS_AprJun18_already_known_virus_seqs.fas # 34 viruses/segments of viruses
## 4) create bowtie2 database from this named sequence fas file 
cd DrosMS_AprJun18/6_read_mapping
mkdir DrosMS_AprJun18_already_known_viruses/DrosMS_AprJun18_already_known_viruses_bowtie_db
bowtie2-build DrosMS_AprJun18_already_known_viruses/DrosMS_AprJun18_already_known_virus_seqs.fas DrosMS_AprJun18_already_known_viruses/DrosMS_AprJun18_already_known_viruses_bowtie_db/DrosMS_AprJun18_already_known_viruses_bowtie_db
## 5) Map raw reads to this bowtie2 database, so that I can see coverage across the genomes for all of the known viruses present
cd /mnt/drive1-6tb/Megan
bowtie2 -x DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_already_known_viruses/DrosMS_AprJun18_already_known_viruses_bowtie_db/DrosMS_AprJun18_already_known_viruses_bowtie_db --very-sensitive -p 20 -q -1 DrosMS_AprJun18/0_reads/DrosMS_AprJun18.1.fastq.gz -2 DrosMS_AprJun18/0_reads/DrosMS_AprJun18.2.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_already_known_viruses/DrosMS_AprJun18_known_virus_read_mapping_sensitive.bam
##this takes quite a while to run...
#get the read depth at each point
cd /mnt/drive1-6tb/Megan/DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_already_known_viruses/
genomeCoverageBed -d -ibam DrosMS_AprJun18_known_virus_read_mapping_sensitive.bam > DrosMS_AprJun18_known_virus_ReadDepth_sensitive.tsv

###################
##DrosMS_JulOct18##
###################
cd /mnt/drive1-6tb/Megan
mkdir DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_already_known_viruses
chmod +x DrosMS_JulOct18/5_BLAST/DrosMS_JulOct18_already_known.tsv
cut -f2 DrosMS_JulOct18/5_BLAST/DrosMS_JulOct18_already_known.tsv | sort | uniq > DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_already_known_viruses/DrosMS_JulOct18_already_known_unique_viruses.tsv
## 3) match unique names to the sequences with the same titles in the Drosophila_viruses_2019.fas file - already made tsv version of this for SeptOct2016
cd /mnt/drive1-6tb/Megan
join -t $'\t' -1 1 -2 1 <(sort -t $'\t' -k 1,1 DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_already_known_viruses/DrosMS_JulOct18_already_known_unique_viruses.tsv) <(sort -t $'\t' -k 1,1 Drosophila_Viruses_Nov2019/Drosophila_Viruses_Nov2019.tsv) | awk '{print ">"$1"\n"$2}' > DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_already_known_viruses/DrosMS_JulOct18_already_known_virus_seqs.fas
##checking the no. of sequences in this final fas file for the bowtie_db
grep -c "^>" DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_already_known_viruses/DrosMS_JulOct18_already_known_virus_seqs.fas # 82 viruses/segments of viruses
## 4) create bowtie2 database from this named sequence fas file 
cd DrosMS_JulOct18/6_read_mapping
mkdir DrosMS_JulOct18_already_known_viruses/DrosMS_JulOct18_already_known_viruses_bowtie_db
bowtie2-build DrosMS_JulOct18_already_known_viruses/DrosMS_JulOct18_already_known_virus_seqs.fas DrosMS_JulOct18_already_known_viruses/DrosMS_JulOct18_already_known_viruses_bowtie_db/DrosMS_JulOct18_already_known_viruses_bowtie_db
## 5) Map raw reads to this bowtie2 database, so that I can see coverage across the genomes for all of the known viruses present
cd /mnt/drive1-6tb/Megan
bowtie2 -x DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_already_known_viruses/DrosMS_JulOct18_already_known_viruses_bowtie_db/DrosMS_JulOct18_already_known_viruses_bowtie_db --very-sensitive -p 20 -q -1 DrosMS_JulOct18/0_reads/DrosMS_JulOct18.1.fastq.gz -2 DrosMS_JulOct18/0_reads/DrosMS_JulOct18.2.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_already_known_viruses/DrosMS_JulOct18_known_virus_read_mapping_sensitive.bam
##this takes quite a while to run...
#get the read depth at each point
cd /mnt/drive1-6tb/Megan/DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_already_known_viruses/
genomeCoverageBed -d -ibam DrosMS_JulOct18_known_virus_read_mapping_sensitive.bam > DrosMS_JulOct18_known_virus_ReadDepth_sensitive.tsv

############################################################################################
##To count reads from each strand of the new viruses to confirm replication + strandedness##
############################################################################################

#Doing for all 5 datasets, and all new viruses
#ugly non-standard approach to count reads from each strand

#SeptOct2016
cd /mnt/drive1-6tb/Megan
#map forward only
bowtie2 -x 11_19_new_viruses_5_datasets_bowtie_db/11_19_new_viruses_5_datasets_bowtie_db -k 1 -p 20 -q -U SeptOct2016/10_19_reanalyses/0_reads/SeptOct2016.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_ForwardReadMappings.bam

#Dec2016June2017
cd /mnt/drive1-6tb/Megan
#map forward only
bowtie2 -x 11_19_new_viruses_5_datasets_bowtie_db/11_19_new_viruses_5_datasets_bowtie_db -k 1 -p 20 -q -U Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.A.1.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.B.1.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.C.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_ForwardReadMappings.bam

#DrosMS_JulOct17
cd /mnt/drive1-6tb/Megan
#map forward only
bowtie2 -x 11_19_new_viruses_5_datasets_bowtie_db/11_19_new_viruses_5_datasets_bowtie_db -k 1 -p 20 -q -U DrosMS_JulOct17/0_reads/DrosMS_JulOct17.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_ForwardReadMappings.bam

#DrosMS_AprJun18
cd /mnt/drive1-6tb/Megan
#map forward only
bowtie2 -x 11_19_new_viruses_5_datasets_bowtie_db/11_19_new_viruses_5_datasets_bowtie_db -k 1 -p 20 -q -U DrosMS_AprJun18/0_reads/DrosMS_AprJun18.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_ForwardReadMappings.bam

#DrosMS_JulOct18
cd /mnt/drive1-6tb/Megan
#map forward only
bowtie2 -x 11_19_new_viruses_5_datasets_bowtie_db/11_19_new_viruses_5_datasets_bowtie_db -k 1 -p 20 -q -U DrosMS_JulOct18/0_reads/DrosMS_JulOct18.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_ForwardReadMappings.bam

##now to count the number of forward and reverse reads for each virus 

##SeptOct2016
samtools view -F 16 SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_Forward_reads_mapped_per_virus.tsv
samtools view -f 16 SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_Reverse_reads_mapped_per_virus.tsv

##Dec2016June2017
samtools view -F 16 Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_Forward_reads_mapped_per_virus.tsv
samtools view -f 16 Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_Reverse_reads_mapped_per_virus.tsv

##DrosMS_JulOct17
samtools view -F 16 DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_Forward_reads_mapped_per_virus.tsv
samtools view -f 16 DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_Reverse_reads_mapped_per_virus.tsv

##DrosMS_AprJun18
samtools view -F 16 DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_Forward_reads_mapped_per_virus.tsv
samtools view -f 16 DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_Reverse_reads_mapped_per_virus.tsv

##DrosMS_JulOct18
samtools view -F 16 DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_Forward_reads_mapped_per_virus.tsv
samtools view -f 16 DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_Reverse_reads_mapped_per_virus.tsv

##Doing a further investigation of Inveresk virus - have saved reverse complement nucleotide as I think this may be the true sequence...and the mRNA rc may be what I currently have in the 5 datasets database
##made a bowtie db from the reverse complement sequence
cd /mnt/drive1-6tb/Megan
bowtie2-build Inveresk_Virus_reverse_complement.fas Inveresk_Virus_reverse_complement_db/Inveresk_Virus_reverse_complement_db

#mapping forward reads only from SeptOct 2016 as a sample
bowtie2 -x Inveresk_Virus_reverse_complement_db/Inveresk_Virus_reverse_complement_db -k 1 -p 20 -q -U SeptOct2016/10_19_reanalyses/0_reads/SeptOct2016.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_Inveresk_virus_RC_ForwardReadMappings.bam

samtools view -F 16 SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_Inveresk_virus_RC_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_Inveresk_virus_RC_Forward_reads_mapped_per_virus.tsv
samtools view -f 16 SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_Inveresk_virus_RC_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_Inveresk_virus_RC_Reverse_reads_mapped_per_virus.tsv

##made a fas file w both orientations of the virus, and Dimm CO1 
##made a bowtie db 
cd /mnt/drive1-6tb/Megan
bowtie2-build Inveresk_virus_strandedness_check.fas Inveresk_Virus_strandedness_check_db/Inveresk_Virus_strandedness_check_db

#mapping forward reads only from Sept Oct 2016
bowtie2 -x Inveresk_Virus_strandedness_check_db/Inveresk_Virus_strandedness_check_db -k 1 -p 20 -q -U SeptOct2016/10_19_reanalyses/0_reads/SeptOct2016.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_Inveresk_virus_strandedness_check_ForwardReadMappings.bam

samtools view -F 16 SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_Inveresk_virus_strandedness_check_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_Inveresk_virus_strandedness_check_Forward_reads_mapped.tsv
samtools view -f 16 SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_Inveresk_virus_strandedness_check_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_Inveresk_virus_strandedness_check_Reverse_reads_mapped.tsv

###Now doing the same test with the other viruses I think will probably be negative sense...and across all 5 datasets
##made fasta file with Dimm CO1 and the negative sense viruses in it
##made a bowtie db 
cd /mnt/drive1-6tb/Megan
mkdir negative_ssRNA_viruses_strandedness_check_db
bowtie2-build negative_ssRNA_viruses_strandedness_check.fasta negative_ssRNA_viruses_strandedness_check_db/negative_ssRNA_viruses_strandedness_check_db
#mapping forward reads only from all databases
bowtie2 -x negative_ssRNA_viruses_strandedness_check_db/negative_ssRNA_viruses_strandedness_check_db -k 1 -p 20 -q -U SeptOct2016/10_19_reanalyses/0_reads/SeptOct2016.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings.bam

##other datasets
bowtie2 -x negative_ssRNA_viruses_strandedness_check_db/negative_ssRNA_viruses_strandedness_check_db -k 1 -p 20 -q -U Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.A.1.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.B.1.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.C.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings.bam

bowtie2 -x negative_ssRNA_viruses_strandedness_check_db/negative_ssRNA_viruses_strandedness_check_db -k 1 -p 20 -q -U DrosMS_JulOct17/0_reads/DrosMS_JulOct17.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings.bam

bowtie2 -x negative_ssRNA_viruses_strandedness_check_db/negative_ssRNA_viruses_strandedness_check_db -k 1 -p 20 -q -U DrosMS_AprJun18/0_reads/DrosMS_AprJun18.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings.bam

bowtie2 -x negative_ssRNA_viruses_strandedness_check_db/negative_ssRNA_viruses_strandedness_check_db -k 1 -p 20 -q -U DrosMS_JulOct18/0_reads/DrosMS_JulOct18.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings.bam

#map forward and reverse reads
samtools view -F 16 SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_negative_ssRNA_viruses_strandedness_check_Forward_reads_mapped.tsv
samtools view -f 16 SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_negative_ssRNA_viruses_strandedness_check_Reverse_reads_mapped.tsv

#Dec2016June2017
samtools view -F 16 Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_negative_ssRNA_viruses_strandedness_check_Forward_reads_mapped.tsv
samtools view -f 16 Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_negative_ssRNA_viruses_strandedness_check_Reverse_reads_mapped.tsv

##DrosMS_JulOct17
samtools view -F 16 DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_negative_ssRNA_viruses_strandedness_check_Forward_reads_mapped.tsv
samtools view -f 16 DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_negative_ssRNA_viruses_strandedness_check_Reverse_reads_mapped.tsv

##DrosMS_AprJun18
samtools view -F 16 DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_negative_ssRNA_viruses_strandedness_check_Forward_reads_mapped.tsv
samtools view -f 16 DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_negative_ssRNA_viruses_strandedness_check_Reverse_reads_mapped.tsv

##DrosMS_JulOct18
samtools view -F 16 DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_negative_ssRNA_viruses_strandedness_check_Forward_reads_mapped.tsv
samtools view -f 16 DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_negative_ssRNA_viruses_strandedness_check_Reverse_reads_mapped.tsv

#optional - investigating difference if using --very-sensitive mapping to map forward and reverse reads to the new virus genomes for strandedness
#SeptOct2016
cd /mnt/drive1-6tb/Megan
bowtie2 -x negative_ssRNA_viruses_strandedness_check_db/negative_ssRNA_viruses_strandedness_check_db --very-sensitive -p 20 -q -U SeptOct2016/10_19_reanalyses/0_reads/SeptOct2016.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings_sensitive.bam
samtools view -F 16 SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings_sensitive.bam | cut -f 3 | sort | uniq -c > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_negative_ssRNA_viruses_strandedness_check_Forward_reads_mapped_sensitive.tsv
samtools view -f 16 SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings_sensitive.bam | cut -f 3 | sort | uniq -c > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_negative_ssRNA_viruses_strandedness_check_Reverse_reads_mapped_sensitive.tsv

#Dec2016June2017
cd /mnt/drive1-6tb/Megan
bowtie2 -x negative_ssRNA_viruses_strandedness_check_db/negative_ssRNA_viruses_strandedness_check_db --very-sensitive -p 20 -q -U Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.A.1.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.B.1.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.C.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings_sensitive.bam

samtools view -F 16 Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings_sensitive.bam | cut -f 3 | sort | uniq -c > Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_negative_ssRNA_viruses_strandedness_check_Forward_reads_mapped_sensitive.tsv
samtools view -f 16 Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings_sensitive.bam | cut -f 3 | sort | uniq -c > Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_negative_ssRNA_viruses_strandedness_check_Reverse_reads_mapped_sensitive.tsv

##DrosMS_JulOct17
cd /mnt/drive1-6tb/Megan
bowtie2 -x negative_ssRNA_viruses_strandedness_check_db/negative_ssRNA_viruses_strandedness_check_db --very-sensitive -p 20 -q -U DrosMS_JulOct17/0_reads/DrosMS_JulOct17.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings_sensitive.bam

samtools view -F 16 DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings_sensitive.bam | cut -f 3 | sort | uniq -c > DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_negative_ssRNA_viruses_strandedness_check_Forward_reads_mapped_sensitive.tsv
samtools view -f 16 DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings_sensitive.bam | cut -f 3 | sort | uniq -c > DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_negative_ssRNA_viruses_strandedness_check_Reverse_reads_mapped_sensitive.tsv

##DrosMS_AprJun18
cd /mnt/drive1-6tb/Megan
bowtie2 -x negative_ssRNA_viruses_strandedness_check_db/negative_ssRNA_viruses_strandedness_check_db --very-sensitive -p 20 -q -U DrosMS_AprJun18/0_reads/DrosMS_AprJun18.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings_sensitive.bam

samtools view -F 16 DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings_sensitive.bam | cut -f 3 | sort | uniq -c > DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_negative_ssRNA_viruses_strandedness_check_Forward_reads_mapped_sensitive.tsv
samtools view -f 16 DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings_sensitive.bam | cut -f 3 | sort | uniq -c > DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_negative_ssRNA_viruses_strandedness_check_Reverse_reads_mapped_sensitive.tsv

##DrosMS_JulOct18
cd /mnt/drive1-6tb/Megan
bowtie2 -x negative_ssRNA_viruses_strandedness_check_db/negative_ssRNA_viruses_strandedness_check_db --very-sensitive -p 20 -q -U DrosMS_JulOct18/0_reads/DrosMS_JulOct18.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings_sensitive.bam

samtools view -F 16 DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings_sensitive.bam | cut -f 3 | sort | uniq -c > DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_negative_ssRNA_viruses_strandedness_check_Forward_reads_mapped_sensitive.tsv
samtools view -f 16 DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_negative_ssRNA_viruses_strandedness_check_ForwardReadMappings_sensitive.bam | cut -f 3 | sort | uniq -c > DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_negative_ssRNA_viruses_strandedness_check_Reverse_reads_mapped_sensitive.tsv

##Now doing the same test with the other viruses I think will probably be positive sense...and across all 5 datasets
##made fasta file with Dimm CO1 and the positive sense viruses in it
##made a bowtie db 
cd /mnt/drive1-6tb/Megan
mkdir positive_viruses_strandedness_check_db
bowtie2-build positive_viruses_strandedness_check.fasta positive_viruses_strandedness_check_db/positive_viruses_strandedness_check_db
#mapping forward reads only from all databases
bowtie2 -x positive_viruses_strandedness_check_db/positive_viruses_strandedness_check_db -k 1 -p 20 -q -U SeptOct2016/10_19_reanalyses/0_reads/SeptOct2016.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_positive_viruses_strandedness_check_ForwardReadMappings.bam

bowtie2 -x positive_viruses_strandedness_check_db/positive_viruses_strandedness_check_db -k 1 -p 20 -q -U Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.A.1.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.B.1.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.C.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_positive_viruses_strandedness_check_ForwardReadMappings.bam

bowtie2 -x positive_viruses_strandedness_check_db/positive_viruses_strandedness_check_db -k 1 -p 20 -q -U DrosMS_JulOct17/0_reads/DrosMS_JulOct17.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_positive_viruses_strandedness_check_ForwardReadMappings.bam

bowtie2 -x positive_viruses_strandedness_check_db/positive_viruses_strandedness_check_db -k 1 -p 20 -q -U DrosMS_AprJun18/0_reads/DrosMS_AprJun18.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_positive_viruses_strandedness_check_ForwardReadMappings.bam

bowtie2 -x positive_viruses_strandedness_check_db/positive_viruses_strandedness_check_db -k 1 -p 20 -q -U DrosMS_JulOct18/0_reads/DrosMS_JulOct18.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_positive_viruses_strandedness_check_ForwardReadMappings.bam

#map forward and reverse reads
samtools view -F 16 SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_positive_viruses_strandedness_check_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_positive_viruses_strandedness_check_Forward_reads_mapped.tsv
samtools view -f 16 SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_positive_viruses_strandedness_check_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_positive_viruses_strandedness_check_Reverse_reads_mapped.tsv

#Dec2016June2017
samtools view -F 16 Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_positive_viruses_strandedness_check_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_positive_viruses_strandedness_check_Forward_reads_mapped.tsv
samtools view -f 16 Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_positive_viruses_strandedness_check_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_positive_viruses_strandedness_check_Reverse_reads_mapped.tsv

##DrosMS_JulOct17
samtools view -F 16 DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_positive_viruses_strandedness_check_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_positive_viruses_strandedness_check_Forward_reads_mapped.tsv
samtools view -f 16 DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_positive_viruses_strandedness_check_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_positive_viruses_strandedness_check_Reverse_reads_mapped.tsv

##DrosMS_AprJun18
samtools view -F 16 DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_positive_viruses_strandedness_check_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_positive_viruses_strandedness_check_Forward_reads_mapped.tsv
samtools view -f 16 DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_positive_viruses_strandedness_check_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_positive_viruses_strandedness_check_Reverse_reads_mapped.tsv

##DrosMS_JulOct18
samtools view -F 16 DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_positive_viruses_strandedness_check_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_positive_viruses_strandedness_check_Forward_reads_mapped.tsv
samtools view -f 16 DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_positive_viruses_strandedness_check_ForwardReadMappings.bam | cut -f 3 | sort | uniq -c > DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_positive_viruses_strandedness_check_Reverse_reads_mapped.tsv


#Looking at impact of --very-sensitive read mapping on mapping of forward reads only to positive sense viruses 

#SeptOct2016
cd /mnt/drive1-6tb/Megan
bowtie2 -x positive_viruses_strandedness_check_db/positive_viruses_strandedness_check_db --very-sensitive -p 20 -q -U SeptOct2016/10_19_reanalyses/0_reads/SeptOct2016.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_positive_viruses_strandedness_check_ForwardReadMappings_sensitive.bam

samtools view -F 16 SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_positive_viruses_strandedness_check_ForwardReadMappings_sensitive.bam | cut -f 3 | sort | uniq -c > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_positive_viruses_strandedness_check_Forward_reads_mapped_sensitive.tsv
samtools view -f 16 SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_positive_viruses_strandedness_check_ForwardReadMappings_sensitive.bam | cut -f 3 | sort | uniq -c > SeptOct2016/10_19_reanalyses/6_read_mapping/SeptOct2016_positive_viruses_strandedness_check_Reverse_reads_mapped_sensitive.tsv

#Dec2016June2017
bowtie2 -x positive_viruses_strandedness_check_db/positive_viruses_strandedness_check_db --very-sensitive -p 20 -q -U Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.A.1.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.B.1.fastq.gz,Dec2016June2017/10_19_reanalyses/0_reads/Dec2016June2017.C.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_positive_viruses_strandedness_check_ForwardReadMappings_sensitive.bam

samtools view -F 16 Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_positive_viruses_strandedness_check_ForwardReadMappings_sensitive.bam | cut -f 3 | sort | uniq -c > Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_positive_viruses_strandedness_check_Forward_reads_mapped_sensitive.tsv
samtools view -f 16 Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_positive_viruses_strandedness_check_ForwardReadMappings_sensitive.bam | cut -f 3 | sort | uniq -c > Dec2016June2017/10_19_reanalyses/6_read_mapping/Dec2016June2017_positive_viruses_strandedness_check_Reverse_reads_mapped_sensitive.tsv

#DrosMS_JulOct17
bowtie2 -x positive_viruses_strandedness_check_db/positive_viruses_strandedness_check_db --very-sensitive -p 20 -q -U DrosMS_JulOct17/0_reads/DrosMS_JulOct17.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_positive_viruses_strandedness_check_ForwardReadMappings_sensitive.bam

samtools view -F 16 DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_positive_viruses_strandedness_check_ForwardReadMappings_sensitive.bam | cut -f 3 | sort | uniq -c > DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_positive_viruses_strandedness_check_Forward_reads_mapped_sensitive.tsv
samtools view -f 16 DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_positive_viruses_strandedness_check_ForwardReadMappings_sensitive.bam | cut -f 3 | sort | uniq -c > DrosMS_JulOct17/6_read_mapping/DrosMS_JulOct17_positive_viruses_strandedness_check_Reverse_reads_mapped_sensitive.tsv

#DrosMS_AprJun18
bowtie2 -x positive_viruses_strandedness_check_db/positive_viruses_strandedness_check_db --very-sensitive -p 20 -q -U DrosMS_AprJun18/0_reads/DrosMS_AprJun18.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_positive_viruses_strandedness_check_ForwardReadMappings_sensitive.bam

samtools view -F 16 DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_positive_viruses_strandedness_check_ForwardReadMappings_sensitive.bam | cut -f 3 | sort | uniq -c > DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_positive_viruses_strandedness_check_Forward_reads_mapped_sensitive.tsv
samtools view -f 16 DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_positive_viruses_strandedness_check_ForwardReadMappings_sensitive.bam | cut -f 3 | sort | uniq -c > DrosMS_AprJun18/6_read_mapping/DrosMS_AprJun18_positive_viruses_strandedness_check_Reverse_reads_mapped_sensitive.tsv

#DrosMS_JulOct18
bowtie2 -x positive_viruses_strandedness_check_db/positive_viruses_strandedness_check_db --very-sensitive -p 20 -q -U DrosMS_JulOct18/0_reads/DrosMS_JulOct18.1.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_positive_viruses_strandedness_check_ForwardReadMappings_sensitive.bam

samtools view -F 16 DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_positive_viruses_strandedness_check_ForwardReadMappings_sensitive.bam | cut -f 3 | sort | uniq -c > DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_positive_viruses_strandedness_check_Forward_reads_mapped_sensitive.tsv
samtools view -f 16 DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_positive_viruses_strandedness_check_ForwardReadMappings_sensitive.bam | cut -f 3 | sort | uniq -c > DrosMS_JulOct18/6_read_mapping/DrosMS_JulOct18_positive_viruses_strandedness_check_Reverse_reads_mapped_sensitive.tsv

##########################################################################
## Reads mapped to new and known viruses - heatmap					    ##
##########################################################################

##Mapping reads from each of the datasets to a single bowtie2 database, created from the new and known virus genomes...known virus genomes selected from the pilot mapping of each dataset to the whole Drosophila virus list.
###Fasta file of virus genomes to map to is in .fas

##Creating a bowtie2 database, also contains the CO1 genes of all species found in the dataset to normalise reads with 
cd 02_21_read_mapping
mkdir 02_21_read_mapping/new_and_known_viruses_bowtie_db
bowtie2-build /mnt/drive1-6tb/Megan/02_21_read_mapping/new_and_known_viruses.fas new_and_known_viruses_bowtie_db/new_and_known_viruses_bowtie_db
## 5) Map raw reads to this bowtie2 database, so that I can see coverage across the genomes for all of the known viruses present
cd /mnt/drive1-6tb/Megan
#trialing on one of the datasets, but use the read_mapping_across_datasets.sh script to map across all datasets
bowtie2 -x 02_21_read_mapping/new_and_known_viruses_bowtie_db/new_and_known_viruses_bowtie_db --very-sensitive -p 8 -q -1 DrosMS_JulOct18/0_reads/DrosMS_JulOct18.1.fastq.gz -2 DrosMS_JulOct18/0_reads/DrosMS_JulOct18.2.fastq.gz | samtools view -u -F 4 - | samtools sort -@ 6 -m 6G - > 02_21_read_mapping/DrosMS_JulOct18_new_and_known_virus_read_mapping.bam

##################################################################
##To create phylogenies of newly described viruses using iqtree ##
##################################################################

##Burdiehouse burn virus 
##alignment
cd /mnt/drive1-6tb/Megan/Phylogenies/Burdiehouse_burn_virus/02_21
t_coffee -seq Burdiehouse_burn_virus_5BLASTp.fasta -mode accurate -n_core 8 -email=megan.wallace@ed.ac.uk ##for some reason t-coffee isn't working...apparently something thats fixed in the latest commit?? Maybe get Darren to update?
##so now trying clustal omega 
clustalo --infile=Burdiehouse_burn_virus_5BLASTp.fasta --seqtype=Protein --infmt=fa --full --outfile=Burdiehouse_burn_virus_5BLASTp.aln.fas --outfmt=fa ##this works
##tree
iqtree -s Burdiehouse_burn_virus_5BLASTp.aln.fas -nt 8 -bb 1000	

##Cockenzie virus 
##alignment 
cd /mnt/drive1-6tb/Megan/Phylogenies/Cockenzie_virus/02_21
clustalo --infile=Cockenzie_virus_5BLASTp.fasta --seqtype=Protein --infmt=fa --full --outfile=Cockenzie_virus_5BLASTp.aln.fas --outfmt=fa 
##tree
iqtree -s Cockenzie_virus_5BLASTp.aln.fas -nt 8 -bb 1000	

##Craighall and Inverleith Toti viruses
##alignment
cd /mnt/drive1-6tb/Megan/Phylogenies/Toti_viruses/02_21
clustalo --infile=Totivirus_5BLASTp.fasta --seqtype=Protein --infmt=fa --full --outfile=Totivirus_5BLASTp.aln.fas --outfmt=fa 
##tree
iqtree -s Totivirus_5BLASTp.aln.fas -nt 8 -bb 1000

##Crammond virga-like virus + Nege-like viruses 
cd /mnt/drive1-6tb/Megan/Phylogenies/Virga_nege_viruses/02_21
clustalo --infile=Negeviruses_5BLASTp.fasta --seqtype=Protein --infmt=fa --full --outfile=Negeviruses_5BLASTp.aln.fas --outfmt=fa 
##tree
iqtree -s Negeviruses_5BLASTp.aln.fas -nt 8 -bb 1000

##Dalkeith virus 
cd /mnt/drive1-6tb/Megan/Phylogenies/Dalkeith_virus/02_21
clustalo --infile=Dalkeith_virus_5BLASTp.fasta --seqtype=Protein --infmt=fa --full --outfile=Dalkeith_virus_5BLASTp.aln.fas --outfmt=fa 
##tree
iqtree -s Dalkeith_virus_5BLASTp.aln.fas -nt 8 -bb 1000

##Combined reovirus phylogeny - Vogrie virus and Glencorse burn virus seg 1 / RdRp 
cd /mnt/drive1-6tb/Megan/Phylogenies/Reoviridae/02_21
clustalo --infile=Reoviruses_5BLASTp.fasta --seqtype=Protein --infmt=fa --full --outfile=Reoviruses_5BLASTp.aln.fas --outfmt=fa 
##tree
iqtree -s Reoviruses_5BLASTp.aln.fas -nt 8 -bb 1000

##Gosford narnavirus
cd /mnt/drive1-6tb/Megan/Phylogenies/Gosford_narnavirus/02_21
clustalo --infile=Gosford_narnavirus_5BLASTp.fasta --seqtype=Protein --infmt=fa --full --outfile=Gosford_narnavirus_5BLASTp.aln.fas --outfmt=fa 
##tree
iqtree -s Gosford_narnavirus_5BLASTp.aln.fas -nt 8 -bb 1000

##Inveresk virus
cd /mnt/drive1-6tb/Megan/Phylogenies/Inveresk_virus/02_21
clustalo --infile=Inveresk_Nyamivirus_5BLASTp.fasta --seqtype=Protein --infmt=fa --full --outfile=Inveresk_Nyamivirus_5BLASTp.aln.fas --outfmt=fa 
##tree
iqtree -s Inveresk_Nyamivirus_5BLASTp.aln.fas -nt 8 -bb 1000

##Lasswade virus
cd /mnt/drive1-6tb/Megan/Phylogenies/Lasswade_virus/02_21
clustalo --infile=Lasswade_virus_5BLASTp.fasta --seqtype=Protein --infmt=fa --full --outfile=Lasswade_virus_5BLASTp.aln.fas --outfmt=fa 
##tree
iqtree -s Lasswade_virus_5BLASTp.aln.fas -nt 8 -bb 1000

##Midmar virus
cd /mnt/drive1-6tb/Megan/Phylogenies/Midmar_virus/02_21
clustalo --infile=Midmar_virus_5BLASTp.fasta --seqtype=Protein --infmt=fa --full --outfile=Midmar_virus_5BLASTp.aln.fas --outfmt=fa
##tree
iqtree -s Midmar_virus_5BLASTp.aln.fas -nt 8 -bb 1000

##North esk virus
cd /mnt/drive1-6tb/Megan/Phylogenies/North_esk_virus/02_21
clustalo --infile=North_esk_virus_5BLASTp.fasta --seqtype=Protein --infmt=fa --full --outfile=North_esk_virus_5BLASTp.aln.fas --outfmt=fa
##tree
iqtree -s North_esk_virus_5BLASTp.aln.fas -nt 8 -bb 1000

##Phlebo_viruses (Phlebo-like viruses - Sighthill and Tranent)
cd /mnt/drive1-6tb/Megan/Phylogenies/Phleboviruses/02_21
clustalo --infile=Phleboviruses_5BLASTp.fasta --seqtype=Protein --infmt=fa --full --outfile=Phleboviruses_5BLASTp.aln.fas --outfmt=fa
##tree
iqtree -s Phleboviruses_5BLASTp.aln.fas -nt 8 -bb 1000

##Sunshine
cd /mnt/drive1-6tb/Megan/Phylogenies/Sunshine/02_21
clustalo --infile=Sunshine_virus_5BLASTp.fasta --seqtype=Protein --infmt=fa --full --outfile=Sunshine_virus_5BLASTp.aln.fas --outfmt=fa
##tree
iqtree -s Sunshine_virus_5BLASTp.aln.fas -nt 8 -bb 1000
