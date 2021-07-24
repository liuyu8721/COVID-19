#!/bin/bash
export PATH=/export/home/liuy683/anaconda3/opt/mummer-3.23:$PATH

cd /Folder_with_Sequence_Fasta_Files/
ref=NC_045512.2.fa # The reference SARS-CoV-2 Wuhan Genome

for input in $(ls *.fasta); do

	prefix=$(echo $input | cut -d'.' -f1)
	echo $input

	dos2unix ${input}
	sed -i '1s/ /_/g' ${input}
	nucmer --forward -p ${prefix} ${ref} ${input}
	show-coords -r -c -l ${prefix}.delta > ${prefix}.coords
	show-snps ${prefix}.delta -T -l > ${prefix}.snps

done

mkdir fasta; mv *.fasta ./fasta/
mkdir snps; mv *.snps ./snps/
mkdir others; mv EPI* ./others/
