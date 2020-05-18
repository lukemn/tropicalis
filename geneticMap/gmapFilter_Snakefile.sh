#!/bin/bash

# process the inferred ancestries at pseudomarkers to make the final genetic map.
# 1. run redundancy, filtering, and large-scale (linkage group) consistency checks
# 2. take the most likely map from a bunch of iterations of greedy marker ordering
# 3. maximise genetic and physical concordance within chromosomes
# 4. save the output, some stats, flag any problematic regions for followup (un-orientable without more information, suspect contigs)

from os import path
configfile: "config.yaml" 
gmaps = config['assembly']
print(gmaps)

rdas = ["{}/rqtl/{}_q90_bestOrders_LL.rda".format(i, i) for i in gmaps.split()]
print(rdas)

rule all:
    input:      
        rdas

rule makeInitialMap:
    input:
        "{assembly}/rqtl/pheno.csv",
        "{assembly}/rqtl/geno.riself.bin50.csv",   
        "{assembly}/rqtl/pmap.riself.bin50.csv",
        "{assembly}/rqtl/riself.rqtl.dat.bin50.rda",
        "{assembly}/rqtl/reorderedSeq.csv",
        "{assembly}/rqtl/bin50.markerpos.txt"
    output:
        "{assembly}/rqtl/RF_LOD.png",
        "{assembly}/rqtl/{assembly}_rimap.rda",
        "{assembly}/rqtl/{assembly}_q90_bestOrders_LL.rda"
    threads:
        20
    log:
        "{assembly}/crossFilter.log"
    shell:
        "scripts/crossFilter.R {wildcards.assembly}/rqtl {threads} > {log}"
       
