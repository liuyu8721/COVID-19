# Rapid SARS-CoV-2 Variant Tracking through Weighted Network Analysis of Frequency Trajectories of Mutations
Early identification of Severe Acute Respiratory Syndrome Coronavirus 2 (SARS-CoV-2) variants has been able to timely track the emergence of clinically important strains and inform the public health response. Here, we proposed a weighted network framework to model the Frequency Trajectories of Mutations (FTMs) for SARS-CoV-2 variant tracing.

## Sequence curation
Sequences are retrieved according to submission date but leveraged according to collection date. 
FastaSplit.sh is used to separate the merged sequence file from GISAID into single genome fasta files.

## Bioinformatic analysis for mutation calling
SNPs and INDELs are determined by a bioinformatic framework proposed by Massacci, et al. (2020). Here, MutationCalling.sh is used to extract the SNPs and INDELs. And a R script, Mercatelli_and_Giorgi.R, adapted from Mercatelli and Giorgi (2020) is used to summarize the mutation information and translate them into proteins.

