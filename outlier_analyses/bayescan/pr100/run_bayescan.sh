#! /usr/bin/env bash

BayeScan2 $1 \
	-o bft_bscan \
	-n 15000 \
	-pilot 15000 \
	-burn 100000 \
	-pr_odds 100 \
	-threads $2 \


#
# alleles.txt  Name of the genotypes data input file 
# -d discarded Optional input file containing list of loci to discard
# -snp         Use SNP genotypes matrix
# --------------------------- 
# | Output                  | 
# --------------------------- 
# -od .        Output file directory, default is the same as program file
## -o alleles   Output file prefix, default is input file without the extension
# -fstat       Only estimate F-stats (no selection)
# --------------------------- 
# | Parameters of the chain | 
# --------------------------- 
# -n 5000      Number of outputted iterations, default is 5000 
# -thin 10     Thinning interval size, default is 10 
# -nbp 20      Number of pilot runs, default is 20 
# -pilot 5000  Length of pilot runs, default is 5000 
# -burn 50000  Burn-in length, default is 50000 
# --------------------------- 
# | Parameters of the model | 
# --------------------------- 
# -pr_odds 10  Prior odds for the neutral model, default is 10 
# -lb_fis 0    Lower bound for uniform prior on Fis (dominant data), default is 0
# -hb_fis 1    Higher bound for uniform prior on Fis (dominant data), default is 1
# -beta_fis    Optional beta prior for Fis (dominant data, m_fis and sd_fis need to be set)
# -m_fis 0.05  Optional mean for beta prior on Fis (dominant data with -beta_fis)
# -sd_fis 0.01 Optional std. deviation for beta prior on Fis (dominant data with -beta_fis)
# -aflp_pc 0.1 Threshold for the recessive genotype as a fraction of maximum band intensity, default is 0.1
# --------------------------- 
# | Output files            | 
# --------------------------- 
# -out_pilot   Optional output file for pilot runs
# -out_freq    Optional output file for allele frequencies
