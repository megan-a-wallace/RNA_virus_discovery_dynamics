######################################################
## New and known virus read mapping across datasets ##
######################################################

#first make a txt table with the viruses and their lengths 
for v in $(cat new.and.known.virus.list.txt)
do
	#assigning the variable that will be used to filter the bam files by virus
	virusid=$(echo $v)
	
	#fetching the virus sequence of interest with header in single line form, then retreiving the sequence, and measuring its length 
	refseqlength=$(grep -A1 "$virusid" new.and.known.virus.linear.fasta | awk '/^>/{getline; getline; print}' | awk '{print length($1)}')
	
	echo adding length of $virusid to table, which is $refseqlength
	
	echo -e "$virusid \t $refseqlength" >> virus.lengths.txt
	
	echo created table of virus lengths
	
	##And now getting the number of reads mapped to each id across each dataset
	readcount_SO16=$(samtools view -b -@ 6 -c SeptOct16_new_and_known_virus_read_mapping.bam $virusid)
	readcount_DJ1617=$(samtools view -b -@ 6 -c Dec2016Jun2017_new_and_known_virus_read_mapping.bam $virusid)
	readcount_JO17=$(samtools view -b -@ 6 -c DrosMS_JulOct17_new_and_known_virus_read_mapping.bam $virusid)
	readcount_AJ18=$(samtools view -b -@ 6 -c DrosMS_AprJun18_new_and_known_virus_read_mapping.bam $virusid)
	readcount_JO18=$(samtools view -b -@ 6 -c DrosMS_JulOct18_new_and_known_virus_read_mapping.bam $virusid)
	
	echo adding reads mapped to $virusid to table, which are $readcount_SO16, $readcount_DJ1617, $readcount_JO17, $readcount_AJ18 and $readcount_JO18
	
	echo -e "$virusid \t $readcount_SO16 \t $readcount_DJ1617 \t $readcount_JO17 \t $readcount_AJ18 \t $readcount_JO18" >> new.and.known.reads.mapped.per.dataset.txt
	
done



