#!/bin/bash

# filter variant calls to snps
# on MQ
# DP (median +- *1.33)
# where parent 1(2) is homozygous ref(alt)
# >10bp from an indel

module load bcftools/intel/1.9

VCF=$1

mdp=$(grep -v '#' $VCF | cut -f8 | grep "^DP" | cut -f1 -d';' | cut -f2 -d'=' |\
      python -c 'from sys import stdin; from numpy import median; l = [int(i) for i in stdin]; print(int(round(median(l))))')
filter='TYPE="snp" && DP >= '
filter+="$((mdp*2/3))"
filter+=' && DP <= '
filter+="$((mdp*4/3))"
filter+=' && QUAL > 50 && GT[0]="RR" && GT[1]="AA"'

bcftools filter -g 10 -i "$filter" $VCF

