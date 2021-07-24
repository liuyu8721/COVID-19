#!/bin/bash

cd /Folder_with_Merged_Sequence_in_Fasta_Format_from_GISAID/

infile=MergedSequence #Assume sequence file: MergedSequence.fasta
Rscript header.R ${infile} #Get the viral sequence header and the row range of nucleotide sequence

mkdir ./${infile}/
for i in $(cat ${infile}.header); do
	f=$(echo $i | cut -d"," -f1)
	s=$(echo $i | cut -d"," -f2)
	e=$(echo $i | cut -d"," -f3)
	echo $f $s $e
	sed -n "${s},${e}p" ${infile}.fasta > ./${infile}/${f}.fasta
done
