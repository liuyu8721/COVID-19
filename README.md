# Rapid SARS-CoV-2 Variant Tracking through Weighted Network Analysis of Frequency Trajectories of Mutations
Early identification of Severe Acute Respiratory Syndrome Coronavirus 2 (SARS-CoV-2) variants has been able to timely track the emergence of clinically important strains and inform the public health response. Here, we proposed a weighted network framework to model the Frequency Trajectories of Mutations (FTMs) for SARS-CoV-2 variant tracing.

## Sequence curation
Sequences are retrieved according to submission date but leveraged according to collection date. 
FastaSplit.sh is used to separate the merged sequence file from GISAID into single genome fasta files.

## Bioinformatic analysis for mutation calling
SNPs and INDELs are determined by a bioinformatic framework proposed by Massacci, et al. (2020). Here, MutationCalling.sh is used to extract the SNPs and INDELs. And a R script, Mercatelli_and_Giorgi.R, adapted from Mercatelli and Giorgi (2020) is used to summarize the mutation information and translate them into proteins.

## Generation and filtration of FTMs
A mutation frequency at a sampling week on a specific site was calculated as the fraction of genomes with the mutation of all genomes sampled at that week. Then the frequency trajectory of a mutation s(1≤s≤S) can be denoted as: y(s)={y(s, t): 1≤t≤T}, where t denotes the week number and t=1 represents the first complete calendar week of 2020 (from January 5 to January 11, 2020). When aggregating the mutation events for each mutation site, all possible mutation directions (e.g. C→T and C→G) were considered to allow the distinction of different variant branches.




