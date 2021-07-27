# Rapid SARS-CoV-2 Variant Tracking through Weighted Network Analysis of Frequency Trajectories of Mutations
Early identification of Severe Acute Respiratory Syndrome Coronavirus 2 (SARS-CoV-2) variants has been able to timely track the emergence of clinically important strains and inform the public health response. Here, we proposed a weighted network framework to model the Frequency Trajectories of Mutations (FTMs) for SARS-CoV-2 variant tracing.

## Sequence curation
Sequences are retrieved according to submission date but leveraged according to collection date. FastaSplit.sh is used to separate the merged sequence file from GISAID into single genome fasta files.

## Bioinformatic analysis for mutation calling
SNPs and INDELs are determined by a bioinformatic framework proposed by Massacci, et al. (2020). Here, MutationCalling.sh is used to extract the SNPs and INDELs. And a R script, Mercatelli_and_Giorgi.R, adapted from Mercatelli and Giorgi (2020) is used to summarize the mutation information and translate them into proteins.

## Generation and filtration of FTMs
A mutation frequency at a sampling week on a specific site was calculated as the fraction of genomes with the mutation of all genomes sampled at that week. Then the frequency trajectory of a mutation s(1≤s≤S) can be denoted as: y(s)={y(s, t): 1≤t≤T}, where t denotes the week number and t=1 represents the first complete calendar week of 2020 (from January 5 to January 11, 2020). When aggregating the mutation events for each mutation site, all possible mutation directions (e.g. C→T and C→G) were considered to allow the distinction of different variant branches. FTMCalculation.R shows an example to aggregate India mutation events to generate FTMs.  
Following, we assume that most of the mutation events are randomly introduced and could not accumulate in the population over time. A hierarchical clustering analysis using Ward's method was applied to group and exclude them before investigating the temporal clustering patterns. To accelerate the calculation, all mutation frequencies less than a threshold (e.g. 0.1%) were first excluded. FTMFiltration.R shows an example to FTM filtration.

## Weighted network construction
Mutations represented by FTMs are referred to as nodes and pairwise synchronous relationships are quantified using Pearson correlation. WeightedNetworkAnalysis.R shows how to create a weighted network model for FTMs.

## Module identification
The topological overlap measure based clustering is utilized to identify modules which may help us to detect and define variants.

## Core mutations for variant determination
Mutations with connection strengths to any other intra-modular nodes larger than a cutoff are chosen as core mutations to determine a variant.

## Phylogenetic branch assessment
Phylogenetic analysis might be necessary for further validation of genome sequences defined by module core mutations. 

## References
Massacci A, Sperandio E, D'Ambrosio L, Maffei M, Palombo F, Aurisicchio L, Ciliberto G, Pallocca M. 2020. Design of a companion bioinformatic tool to detect the emergence and geographical distribution of SARS-CoV-2 Spike protein genetic variants. J Transl Med. 18(1):494.  
Mercatelli D, Giorgi FM. 2020. Geographic and genomic distribution of SARS-CoV-2 mutations. Front Microbiol. 11:1800.

