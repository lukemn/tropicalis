#!/bin/bash

REF=$1
VCF=$2
DIR=`dirname $VCF`
# thinned, filtered SNPs
SITES=$DIR/snp.sites.bed

# variable sequences
cut -f1 $SITES | uniq > $DIR/vscaf
# ... and their lengths (>1kb)
scripts/lookup.py -k $REF.fai -rc 2 -s $DIR/vscaf | sort -k2,2gr | awk '($2>1000)' > $DIR/uvscaf
ln -s $DIR/uvscaf msg.chrLengths

NSEQ=$(cat $REF.fai | wc -l)
NVARSEQ=$(cat $DIR/vscaf | wc -l)
# filtered
NSNP=$(grep -v "#" $VCF | wc -l)
# thinned
NFSNP=$(cat $SITES | wc -l)
VLEN=$(awk '{S+=$2} END {print S/1e6}' $DIR/uvscaf)

printf "%s reference sequences\n%s variants pass filtering\n%s to be used after thinning\n%s variable seq spanning %.2f Mb\n" $NSEQ $NSNP $NFSNP $NVARSEQ $VLEN

mkdir -p hmm_fit hmm_data/refs/par1 hmm_data/refs/par2
# ref is first sample
i=1; awk '{OFS="\t"} {print $1,$3,$4}' $SITES > hmm_data/refs/par${i}/par${i}.alleles
i=2; awk '{OFS="\t"} {print $1,$3,$5}' $SITES > hmm_data/refs/par${i}/par${i}.alleles



# grep -v "#" $VCF | cut -f1,2,4,5 | awk '{OFS="\t"} {print $1,$2-1,$2,$3,$4}'> $SITES