#!/usr/bin/env Rscript

## Required: 
##   individual ID(s) (-i, comma-delimited, default is all) 
##   working directory (-d, > hmm_data/i and hmm_fit/i for line i)
## Optional: 
##   contig(s) (-c, comma-delimited, default is all)
##   number of threads for parallel (-np, default is 1)
##   many others (see below), n.b. genotype priors, recrate & rfac.

## 1. filter and reformat pileup data by individual
##  from full sorted pileup : sprintf("%s/pileup/%s.p*", dir, indiv)
##  removes indels, converts read string to base counts
##  intersects with both parents using a single file for all contigs : sprintf('%s/hmm_data/refs/%s/%s.alleles', dir, sp, sp)
##  writes filtered data.frame : sprintf("%s/hmm_data/%s.hmmdata", dir, indiv)
## 2. fit hmm by individual
##  further site filtering. default is to sample a single variant per read
##  run forward/backward algorithm given priors, reads, base quals
##  dump probabilities and plots (for -y)
##  writes data.frame : sprintf("%s/hmm_prob/%s.hmmprob.txt", dir, indiv)
## sequence lengths are assumed to be in "msg.chrLengths" (no header; chr length in first two col)

# requires MSG and additional R code (below) and two C routines (countalleles, hmmprobs)
MSGD="scripts"
source(sprintf("%s/ded.R", MSGD))
source(sprintf("%s/hmmlib_lmn.R", MSGD))
library(R.methodsS3, quietly = T)
library(R.oo, quietly = T)
library(data.table, quietly = T)
library(parallel, quietly = T)
library(dplyr, quietly = T)
library(ggplot2)
Sys.setenv(OMP_NUM_THREADS = 1)
# options(error=quote(q("yes")))

opts <- getopts()
# opts = list()
indivs <- opts$i; stopifnot(!is.null(indivs)); indivs <- argsplit(indivs)
dir <- opts$d; if(is.null(dir)) dir = "."                         # main working directory
v = opts$v; if(is.null(v)) v=0 else v=as.integer(v)               # verbosity (0|1)
contigLengths <- getContigLengths(f="msg.chrLengths", v=v)
contigs=opts$c; if(is.null(contigs)) contigs <- contigLengths$chr else contigs <- argsplit(contigs)
np <- opts$t; if(is.null(np)) np = 1                              # N MPI threads
clean <- opts$o; if(is.null(clean)) clean=F else clean=T          # overwrite
job <- opts$j; if(is.null(job)) job=0                             # job: 1 = write hmm_data, 2 = hmm_fit, 3 = rqtl, plot indivs, 0=all
print(opts)
##########################################################################################
## optional args, fixed parameters
indir = sprintf("%s/hmm_data", dir); outdir = sprintf("%s/hmm_fit", dir)
## scaling factor for recombination likelihood. Sensitive (larger = more recombination).
rfac <- as.numeric(opts$r); if(!len(rfac)) opts$r=rfac=1e-12
## sample a single variant from each read, default = TRUE
one.site.per.read <- opts$u; if(is.null(one.site.per.read)) opts$u=one.site.per.read=T else one.site.per.read=F
## prior genotype frequencies (par1 hom, het, par2 hom)
priors <- opts$z; if(is.null(opts$z)) priors="0.5,0.0,0.5"; opts$z=priors=as.numeric(argsplit(priors))
## expected recombination events per ~genome (sum of pileup sequence lengths). Not very sensitive. Default ~1/chrom
recrate <- as.numeric(opts$a); if(!len(recrate)) recrate=6/sum(contigLengths$length)
# male/X ignored at present
sex <- opts$s; if(is.null(sex)) sex = 'female'
sex.chrom <- opts$x; if(!is.null(sex.chrom)) sex.chrom <- argsplit(sex.chrom)
## prob genotyping error / neither parental allele
deltapar1 <- as.numeric(opts$p); if(!len(deltapar1)) deltapar1=1e-4
deltapar2 <- as.numeric(opts$q); if(!len(deltapar2)) deltapar2=1e-4
ploidy = 2
minCoverage <- 0 # ignore positions with < minCoverage read depth
maxCoverage <- 30 # subsample bases above maxCoverage
theta <- 1 # scaling factor for read independence (0-1, 1 = fully independent)
species <- c("par1", "par2") # strain labels
pupsp <- "par1" # reference
toplot <- opts$y; if(!is.null(toplot)) toplot <- argsplit(toplot) # contigs to plot
alleles <- c("A","C","G","T","N")
##########################################################################################

write_hmm_data <- function(indiv, dir, indir, contigs, NP=1, v=0){
  
  for(sp in species) {
    reffile <- sprintf('%s/refs/%s/%s.alleles', indir, sp, sp)
    stopifnot(file.exists(reffile))
  }
  pupfile <- Sys.glob(sprintf("%s/pileup/%s.p*", dir, indiv))
  if(len(pupfile)!=1) {cat("MISSING or >1 ", pupfile,"\n"); return(1)}
  ## samtools-0.1.9 pileup with consensus
  # pups <- fread(sprintf("cut -f1-4,9,10 %s", pupfile), data.table=F, showProgress = F, verbose = F)
  ## replace with v1.6 mpileup (chrom, pos, ref, nreads, reads, quals). 
  ## consensus is only used to check ! NA
  pups <- fread(pupfile, data.table=F, showProgress = F, verbose = F, sep='\t')
  names(pups) <- c('contig', 'pos', 'ref', 'cons', 'reads', 'quals')
  todo = NULL
  for(i in contigs){
    if(!i %in% pups$contig) {
      if(v) cat("Contig", i, "not in pileup -- skipping\n")
    } else {
      todo = c(todo, i)
    }
  }
  stopifnot(len(todo)>0)
  cat(sprintf('PARSING PILEUP DATA FOR %s: %s SCAFFOLDS IN PILEUP > %s\n', indiv, len(todo), indir))
  contigs = todo
  pups = pups[pups$contig %in% todo,]
  hmm_data <- sprintf("%s/%s.hmmdata", indir, indiv)
  if(file.exists(hmm_data) & !clean) {cat("File exists -- skipping:", hmm_data, "\n"); return(1)}
  
  pup <- do.call(rbind, mclapply(split(pups, pups$contig), mc.cores=NP, function(pup) {
    
    contig = pup$contig[1]
    pup = pup[,-1]
    stopifnot(ncol(pup) == 5)
    
    if(v) cat("\nWriting HMM input data file for", indiv, contig, "\n")
    if(v) cat("\tStarting with a total of", nrow(pup), "positions.\n")
    ok=pup$ref %in% alleles; pup <- pup[ok,]
    if(v) cat("\tRemoving", sum(!ok), "positions at which ref is not [ACGT].\n")
    pup <- pup[grep("^\\*+$",as.vector(pup$reads),invert=T),]
    if(nrow(pup)>0){
      if(v) cat("\tDecoding read bases...")
      tab <- decode.pileup.bases.count(pup$reads, pup$ref, MSGD, alleles, v=v)
      if(nrow(tab) != nrow(pup)) {print(str(tab)); print(str(pup)); stop("nrows differ after decoding")}
      pup <- cbind(pup, tab); pup$bad=""
      for(sp in species) {
        reffile <- sprintf('%s/refs/%s/%s.alleles', indir, sp, sp)
        ref <- fread(cmd=sprintf('grep -w "%s" %s | cut -f2,3', contig, reffile), data.table=F, header=F, showProgress = F, verbose = F, sep='\t')
        names(ref)=c("pos","allele")
        ## remove positions that are indels in either par1 or par2 (based on ref)
        non_indels_pos <- intersect(ref$pos,pup$pos)
        pup <- pup[pup$pos %in% non_indels_pos,]
        ref <- ref[ref$pos %in% non_indels_pos,]
        if(v) cat("\tKeeping", len(non_indels_pos), "intersection...\n")      
        if(nrow(ref) != nrow(pup)) ref <- ref[match(pup$pos, ref$pos),]
        stopifnot(nrow(ref) == nrow(pup), ref$pos == pup$pos)
        spref <- sprintf("%sref", sp)
        pup[[spref]] <- ref$allele      
        ok <- !is.na(pup[[spref]])
        if(v) cat("\tMasking", sum(!ok), "positions at which", sp, "ref was NA or N.\n")
        pup$bad[!ok] <- paste(sp, "NA/N")      
        ok <- pup[[spref]] %in% c("A","C","G","T")
        if(v) cat("\tMasking", sum(!ok), "positions at which", sp, "ref is not [ACGT].\n")
        pup$bad[!ok] <- paste(sp, "not ACGT")
      }
      ok <- !is.na(pup$ref)
      if(v) cat("\tMasking", sum(!ok), "positions at which par1 ref (pileup) is NA\n")
      pup$bad[!ok] <- "par1 ref NA"
      ok <- pup$ref %in% c("A","C","G","T")
      if(v) cat("\tMasking", sum(!ok), "positions at which par1 ref (pileup) is not [ACGT].\n")
      pup$bad[!ok] <- "par1 ref not ACGT"
      ok <- pup$ref == pup$par1ref
      if(v) cat("\tMasking", sum(!ok), "positions at which par1 ref (pileup) disagrees with par1 ref.\n")
      pup$bad[!ok] <- "par1 ref disagree"    
      if(v) cat("\tAfter filtering, keeping", sum(pup$bad=="", na.rm=T), "positions.\n")
      cbind(contig, pup) 
    }
  }))
  write.table(pup, file=hmm_data, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
}

fit_hmm <- function(indiv, indir, outdir, contigs, NP=1, v=0){
  
  cat(sprintf('FITTING HMM FOR %s: %s SCAFFOLDS > %s\n', indiv, len(contigs), outdir))
  
  hmm_data=sprintf("%s/%s.hmmdata", indir, indiv)
  stopifnot(file.exists(hmm_data))
  all_data <- fread(hmm_data, data.table=F, na.strings="", showProgress = F, verbose = F, colClasses = list(character='bad'), sep='\t')
  hmm_prob <- sprintf("%s/%s.hmmprob.txt", outdir, indiv)
  if(file.exists(hmm_prob) & !clean) {cat("HMM fit for indiv", indiv, " already exists ... skipping\n"); return(1)}
  contigs <- contigs[contigs %in% all_data$contig]
  stopifnot(len(contigs)>0)
  
  probs <- do.call(rbind, mclapply(contigs, mc.cores=NP, mc.preschedule=F, function(contig) {
    phi <- priors
    data <- all_data[all_data$contig==contig,]
    if (nrow(data)==0) return()
    data$read <- factor.contiguous(data$pos) # within 5bp
    total_sites <- len(unique(data$read))
    if(v) cat("\tRound 2: Total number of contiguous regions", total_sites, "\n")      
    ok <- !is.na(data$bad) | !is.na(data$par1ref) & !is.na(data$par2ref) & !is.na(data$cons)
    if(v) cat("\tRound 2: Removing", sum(!ok), "sites at which par1/par2/cons allele unknown\n")
    data$bad[!ok] <- "par1/par2/cons unknown"  
    ok <- data$A + data$C + data$G + data$T > 0
    if(v) cat("\tRound 2: Removing", sum(!ok), "sites at which cons allele is known but reads are unknown\n")
    data$bad[!ok] <- "reads unknown"      
    ok <- !is.na(data$bad) | data$par1ref %in% alleles & data$par2ref %in% alleles
    data$bad[!ok] <- "par1/par2 not in ACGT"
    if(v) cat("\tRound 2: Removing", sum(!ok), "sites at which par1/par2 ref not %in% {", paste(alleles, collapse=", "), "}\n")      
    data <- data[is.na(data$bad),]      
    ok <- data$par1ref != data$par2ref
    if(v) cat("\tRemoving", sum(!ok), "sites at which par1 == par2\n")
    data <- data[ok,]
    data$count <- data$A + data$C + data$G + data$T #+ data$N
    ok <- data$count >= minCoverage
    if(v) cat("\tRemoving", sum(!ok), "sites at where coverage is < ",minCoverage,"\n")
    data <- data[ok,]
    if(one.site.per.read) {
      data$read <- factor(data$read)
      ok <- !duplicated(data$read)
      if(v) cat("\tRemoving", sum(!ok), "sites from same reads\n")
      data <- data[ok,,drop=F]
      if(v) cat("\tNumber of informative markers:", nrow(data), "\n")
    }
    if(v) cat("\tFinal total of", nrow(data), "sites at which par1 != par2\n")
    if (nrow(data)==0) return(1)
    
    ancestries <- c("par1/par1","par1/par2","par2/par2")
    L <- nrow(data)
    K <- len(ancestries)
    if(v) cat("\tTransition probabilities  ...\n")
    r <- recrate; if(v) cat("\trec rate set to ", r, "\n")      
    d <- c(NA, diff(data$pos))
    p <- 1 - exp(-r*d*rfac)
    Pi <- array(dim=c(L,K,K), dimnames=list(NULL, ancestries, ancestries))
    if(ploidy == 2) {
      Pi[,"par1/par1","par1/par1"] <- Pi[,"par1/par2","par1/par2"] <- Pi[,"par2/par2","par2/par2"] <- 1-p
      Pi[,"par1/par1","par1/par2"] <- Pi[,"par1/par2","par1/par1"] <- Pi[,"par1/par2","par2/par2"] <- Pi[,"par2/par2","par1/par2"] <- p
      Pi[,"par1/par1","par2/par2"] <- Pi[,"par2/par2","par1/par1"] <- 0
    } else {
      Pi[,"par1","par1"] <- Pi[,"par2","par2"] <- 1-p
      Pi[,"par1","par2"] <- Pi[,"par2","par1"] <- p
    }
    Pi[1,,] <- NA      
    
    ## Allele frequencies in parental backgrounds
    alleles = alleles[1:4]
    ppar1 <- ppar2 <- matrix(NA, nrow=4, ncol=4, dimnames=list(alleles, alleles))
    ppar1[] <- deltapar1/3
    diag(ppar1) <- 1-deltapar1
    ppar2[] <- deltapar2/3
    diag(ppar2) <- 1-deltapar2
    p1 <- ppar1[data$par1ref,,drop=F]
    p2 <- ppar2[data$par2ref,,drop=F]
    p12 <- array(c(p1,p2), dim=c(dim(p1),2))
    dimnames(p12) <- list(NULL, alleles, NULL)
    ds=sum((data$A+data$C+data$G+data$T+data$N)>maxCoverage)
    if(v) cat("\ttake <=", maxCoverage, "reads (downsample", ds, "sites), translate bases, quals ...\n")
    N <- min(max(data$A+data$C+data$G+data$T+data$N),maxCoverage) ## Total number of reads
    eps <- paste('eps',seq(1,N,by=1),sep='')
    read <- paste('read',seq(1,N,by=1),sep='')
    
    y <- data[,c(alleles,"reads","quals","par1ref"),drop=F]
    y$selected_allele <- NA
    y[,eps] <- rep(0.0,N)
    y[,read] <- rep(5,N)
    fun <- function (x) {if (x=='A') return (0); if (x=='C') return (1); if (x=='G')	return (2); if (x=='T')	return (3); if (x=='N') return (5)}
    for(i in 1:nrow(y)) {
      total.reads <- unlist(strsplit(cleanupReadPileup(y[i,"reads"],y[i,"par1ref"]),''))
      y[i,"selected_allele"] <- total.reads[sample(len(total.reads),1)] ## Sample one read for plotting
      qual <- qual_corrected <- NULL
      for (s in 1:min(len(total.reads),N)) {y[i,read[s]] <-lapply(total.reads[s],fun); qual<-c(qual,(charToInt(unlist(strsplit(y[i,"quals"],''))[s])-33))}
      for (g in 1:len(qual)) {qual_corrected[g]<-qual[g]*(theta^(rank(-qual)[g]-1)); y[i,eps[g]] <- 10^(-(qual_corrected[g])/10)}
    }
    data$read_allele <- as.vector(y[,"selected_allele"])
    
    if(v) cat("\tEmission probabilities ...\n")
    prob = Pr.y.given.z(y=y[,read,drop=F], p=p12, n=N, eps=y[,eps,drop=F], ploidy=ploidy, C=TRUE, dir=MSGD, chrom=contig, id=indiv)
    colnames(prob) <- paste("Pr(y|", ancestries, ")")
    if(v) cat("\tPosterior probability ...\n")
    hmm <- forwardback.ded(Pi=Pi, delta=phi, prob=prob)
    prob[prob < 1e-30] <- 1e-30
    data <- cbind(data, prob)
    data$est <- apply(prob, 1, which.max)
    Pr.z.given.y <- exp(hmm$logalpha + hmm$logbeta - hmm$LL)
    Pr.z.given.y[Pr.z.given.y < 1e-30] <- 1e-30
    colnames(Pr.z.given.y) <- paste("Pr(", ancestries, "|y)")
    cbind(data, Pr.z.given.y)
  }))
  write.table(probs, file=hmm_prob, quote=F, row.names = F, sep='\t')
}

plotAncestry <- function(indivs, outdir, contigs, np=1, v=0){
  
  span = sum(contigLengths$length[contigLengths$chr %in% contigs])/1e6
  cat(sprintf('PLOTTING GENOTYPES: %s LINES, %s SCAFFOLDS (%.2f Mb) > %s\n', len(indivs), len(contigs), span, outdir))
  
  o = mclapply(indivs, mc.cores = np, function(indiv) {
    data <- readProbs(indiv, outdir, contigs)
    p <- ggplot(data) + geom_line(aes(pos/1e6, y=hid, col = haplotype, size = prob)) + geom_rug(aes(pos/1e6), sides='b', size=1e-2, alpha=0.5) + 
      coord_cartesian(ylim = c(-1.5, 1.4)) + theme_classic() + theme(legend.position = 'none') + scale_y_continuous(breaks = c(-1, 0, 1), labels = unique(data$haplotype)) + 
      scale_x_continuous(expand=c(0,0), breaks = scales::pretty_breaks(n=2)) + labs(x = 'Mb', y='') + facet_grid(.~contig, space='free', scales = 'free') + 
      scale_size(range=c(1e-3, 12))
    if(len(contigs)==1){
      plot.file <- sprintf("%s/%s-%s.hmmprob.pdf", outdir, indiv, data$contig[1])
    } else {
      plot.file <- sprintf("%s/%s-%s-%s.hmmprob.pdf", outdir, indiv, data$contig[1], rev(data$contig)[1])
    }
    ggsave(p, filename = plot.file, h = 2, w = max(8, span/2))
  })
}

dumpRQTL <- function(indivs, hmmdir, outdir='rqtl', bin=50, minProb=0.8, maxNA=0.3, np=1, v=0, plot=F, reorderSeq=T){
  
  cat(sprintf('COLLAPSING GENOTYPES: %s LINES %s > %s\n', len(indivs), hmmdir, outdir))
  cat(sprintf('\tFiltering to > %s posterior probability\n
              \tBinning by %s markers\n
              \tFull interpolation between markers\n
              \tMaximum marker per cent NA %s\n\n', minProb, bin, maxNA))
  
  ## Merge sampled markers across lines and write R/qtl files.
  ## Generate (plot) some summary stats.
  ## Ancestries are thresholded and missing data is imputed from flanking markers.
  ## Pseudomarkers are created using binned calls (bin markers).
  ## Minimal filtering is done (QC, redundancy): TBD in Rqtl
  ## Pseudomarkers with > maxNA missing data are removed. 
  ## Files are dumped for both f2 and riself cross-types (the latter excluding all het calls)
  ## if reorderSeq, sequences are renamed in ascending order by physical size
  
  ## threshold PP and filter on physical size and number of markers per haplotype block (minBlock, minBlockMarkers)
  ## return list of per-individual lists containing thresholded data, number of breakpoints/contig, 
  ## and % markers in the grey zone (0.05 >= PP < minProb)
  
  # hmmdir="hmm_fit"; outdir='rqtl'; indivs = system("ls pileup/ | sed 's/.pup//'", intern=T); bin=50; minProb=0.8; maxNA=0.3; reorderSeq=T; species <- c("par1", "par2"); contigLengths <- getContigLengths(f="msg.chrLengths")
  
  setwd(outdir)
  if(!file.exists(hmmdir)) hmmdir=sprintf("../%s", hmmdir)
  print(indivs)
  
  fdata = mclapply(indivs, mc.cores = np, function(i) {
    probs = readProbs(i, hmmdir)
    pr = table(round(probs$prob, 1))
    gz = sum(pr[as.numeric(names(pr)) %% 1 > 0])/sum(pr)
    data = filterBreaks(probs, minpp=minProb, minBlock=1e4, minBlockMarkers=3)
    btab = aggregate(data=data, nb~contig, min)
    names(btab)[2] <- names(data)[12] <- i
    list(calls = data[,c(2,3,12)], breaks = btab, grey = gz)
  })
  
  ## full marker X line haplotype matrix. very sparse due to coverage/read sampling
  gtout = lapply(fdata, '[[', 1) %>% Reduce(function(d1,d2) full_join(d1,d2, by=c('contig', 'pos')), .)
  gtout <- gtout[order(gtout$contig, gtout$pos),]
  
  ## % missing data by contig, marker, line
  nac <- unlist(lapply(split(gtout[,-(1:2)], gtout$contig), function(x) sum(is.na(x))/prod(dim(x))))
  susc <- sus(nac)
  gtout = gtout[!gtout[,1] %in% names(susc),]
  nam <- apply(gtout[,-(1:2)], 1, function(x) sum(is.na(x))/len(x))
  nal <- apply(gtout[,-(1:2)], 2, function(x) sum(is.na(x)))/nrow(gtout)
  susl <- sus(nal, nsd=2, ix=T)
  
  if(reorderSeq){
    # rename Ctg%%% by size; assuming < 1000 seq
    contigLengths <- contigLengths[order(contigLengths$length, decreasing=T),]
    contigLengths$contig = contigLengths$chr
    contigLengths$seq = sprintf("Ctg%.3i", 1:nrow(contigLengths))
    write.table(contigLengths[,-1], file = 'reorderedSeq.csv', row.names = F, quote=F, sep = ',')
    gtout <- merge(contigLengths[,c('contig', 'seq')], gtout, sort=F)
    gtout <- gtout[order(gtout$seq, gtout$pos),-1]
    names(gtout)[1] = 'contig'
  }
  
  if(plot){
    ## % hets
    hetc <- unlist(lapply(split(gtout[,-(1:2)], gtout$contig), function(x) sum(x=='par1par2', na.rm=T)/sum(!is.na(x))))
    p <- qplot(hetc, xlab = '% het by sequence'); ggsave('hetsBySeq.pdf', h=4, w=4)
    sushc <- sus(hetc)
    hetm <- apply(gtout[,-(1:2)], 1, function(x) sum(x=='par1par2', na.rm=T)/sum(!is.na(x)))
    p <- qplot(hetm, xlab = '% het by marker'); ggsave('hetsByMarker.pdf', h=4, w=4)
    sushm <- sus(hetm, nsd=5, ix=T)
    hetl <- apply(gtout[,-(1:2)], 2, function(x) sum(x=='par1par2', na.rm=T)/sum(!is.na(x)))
    # drop lines that are majority het
    # disable for backcross
    # if(plot) p <- qplot(hetl, xlab = '% het by line'); ggsave('hetsByLine.pdf', h=4, w=4)
    # sushl = sus(hetl, nsd=2, ix=T)
    # gtout <- gtout[,!sushl]
    
    ## n breakpoints by contig, line
    lbycontig <- lapply(fdata, '[[', 2) %>% Reduce(function(d1,d2) full_join(d1,d2, by=c('contig')), .)
    bycontig = apply(lbycontig[,-1], 1, function(x) sum(x, na.rm=T))
    p <- qplot(bycontig, xlab = 'breakpoints per sequence'); ggsave('bpBySeq.pdf', h=4, w=4)
    byline = apply(lbycontig[,-1], 2, function(x) sum(x, na.rm=T))
    susbyline <- sus(byline, nsd=2)
    p <- qplot(byline, xlab = 'breakpoints per line'); ggsave('bpByLine.pdf', h=4, w=4)
    bycontigAny <- apply(lbycontig[,-1], 1, function(x) sum(x>0, na.rm=T))/(ncol(lbycontig)-1)
    p <- qplot(bycontigAny, xlab = '% >=1 breakpoint/line/sequence'); ggsave('bpBySeqAny.pdf', h=4, w=4)
    
    ## breakpoints by physical distance
    contigL <- getContigLengths(f="../msg.chrLengths")
    contigL$contig <- as.character(sapply(contigL$chr, function(x) rev(strsplit(x, '_')[[1]])[1]))
    bycontig <- merge(data.frame(contig = lbycontig$contig, nb = bycontig), contigL[,-1])
    bycontig$gd <- bycontig$nb/(bycontig$length/1e6)
    p <- qplot(bycontig$gd, xlab = 'breakpoint per Mb'); ggsave('bpByDist.pdf', h=4, w=4)
    p <- qplot(bycontig$nb, bycontig$length/1e6, xlab = 'breakpoints', ylab = 'Mb'); ggsave('bpByMb.pdf', h=4, w=4)
    
    ## parental proportions
    pprop <- do.call(rbind, apply(gtout[,-(1:2)], 2, function(x) as.data.frame(table(x)/len(x))))
    pprop$line <- sapply(rownames(pprop), function(x) strsplit(x, '.', fixed=T)[[1]][1])
    p <- ggplot(pprop, aes(Freq, fill=x)) + geom_histogram(position = 'dodge') + xlab('genotype frequencies by line'); ggsave('genoFreq.pdf', h=4, w=4)
    anc = c('par1par1', 'par1par2', 'par2par2')
    ppropm <- do.call(rbind, mclapply(split(gtout[,-(1:2)], gtout$contig), mc.cores = np, function(x) {
      do.call(rbind, lapply(1:nrow(x), function(i) sapply(anc, function(z) sum(x[i,]==z, na.rm=T))))
    }))
    ppropm <- cbind(gtout[,1:2], ppropm)
    ppropm <- melt(as.data.table(ppropm), id = c('contig', 'pos'))
    p <- ggplot(ppropm, aes(value/(ncol(gtout)-2), fill=variable)) + geom_histogram(position = 'dodge') + xlab('genotype frequencies by marker'); ggsave('genoFreqMarker.pdf', h=4, w=4)
    p <- ggplot(ppropm, aes(pos/1e6, value/(ncol(gtout)-2), col=variable)) + geom_point(alpha=0.3, size=0.3) + xlab('Mb') + ylab('%') + facet_grid(.~contig, space='free', scales = 'free'); ggsave('genoFreqSeq.pdf', h=4, w=20)
  }
  
  ## interpolate missing positions where flanking markers agree (any distance, ignore edges)
  ## NA 14797632 > 862329 (bin 30)
  ## 14583869 > 91850 (bin 50)
  gtinterpol <- as.data.frame(do.call(rbind, mclapply(split(gtout[,-(1:2)], gtout$contig), mc.cores = np, function(x) infill(x))))
  # gtinterpol <- mclapply(split(gtout[,-(1:2)], gtout$contig), mc.cores = np, function(x) as.data.frame(infill(x)))
  # clr = unlist(lapply(gtinterpol, nrow))
  # clc = unlist(lapply(gtinterpol, ncol))
  cat(sprintf('before/after interpolation NA = %s/%s\n', sum(is.na(gtout[,-(1:2)])),sum(is.na(gtinterpol))))
  nam <- apply(gtinterpol, 1, function(x) sum(is.na(x))/len(x))
  gtinterpol <- gtinterpol[!nam>maxNA,]
  gtout <- gtout[!nam>maxNA,]
  
  ## approx breakpoints within contigs, as evidence for misjoins
  ## ignore terminal windows, which commonly show spurious het calls
  binnedbp <- do.call(rbind, mclapply(split(cbind(gtout[,1:2], gtinterpol), gtout$contig), mc.cores = np, function(x) {
    if(nrow(x)>3){
      y = x[,1:2]
      x = as.matrix(x[,-(1:2)])
      a = data.frame(switch=unlist(lapply(2:nrow(x), function(i) sum(x[i,]!=x[i-1,], na.rm=T))))
      a$b = rep(seq(1, nrow(a)+bin, bin), each=bin)[1:nrow(a)]
      a$pos = y$pos[-1]
      if(nrow(a)>0) {
        sw = aggregate(data=a, switch~b, sum)
        sw$start = aggregate(data=a, pos~b, min)$pos
        sw$end = aggregate(data=a, pos~b, max)$pos
        sw$contig = y$contig[1]
        sw$switch[sw$b==1] <- NA
        sw$switch[sw$b==max(sw$b)] <- NA
        sw
      }
    }
  }))
  if(plot) {
    p <- qplot(binnedbp$switch, xlab = sprintf('sum breakpoints (%s marker bins)', bin)); ggsave('bpBin.pdf', h=4, w=4)
    p <- ggplot(subset(binnedbp, switch>0), aes(start/1e6, switch)) + geom_point() + labs(x='Mb', y = sprintf('sum breakpoints (%s marker bins >0)', bin)) + facet_grid(.~contig, space='free', scale='free'); ggsave('bpBinContig.pdf', h=3, w=30)
  }
  write.table(binnedbp, file = sprintf('bin%s.markerpos.txt', bin), row.names=F, quote=F, sep='\t')
  
  ## make rqtl cross files using binned data as pseudomarkers (or majority assignment if < bin)
  ## set bins with mixed parentage to NA
  bingt <- mclapply(split(gtinterpol, gtout$contig), mc.cores = np, function(x) {
    if(nrow(x)==1){
      cbind(V1=1, x)
    } else {
      spl = rep(seq(1, nrow(x)+bin, bin), each=bin)[1:nrow(x)]
      if(len(spl)>1){
        o = do.call(rbind, lapply(split(x, spl), function(y) apply(y, 2, function(z) {if(len(table(z))==1) names(table(z)) else NA})))
      } else {
        o = apply(x, 2, function(y) names(which.max(table(y))))
      }
      as.data.frame(cbind(unique(spl), o))
    }
  })
  # unlist(lapply(bingt, nrow))
  # unlist(lapply(bingt, ncol))
  bingto <- do.call(rbind, bingt)
  lines = names(gtout)[-(1:2)]
  
  ## F2 code.
  bmarkers = rownames(bingto)
  bmarkers[grep('.', bmarkers, fixed=T, invert=T)] <- paste0(bmarkers[grep('.', bmarkers, fixed=T, invert=T)], '.1')
  bcontig = unlist(lapply(strsplit(bmarkers, '.', fixed=T), '[[', 1))
  bpos = bingto[,1]
  gto = as.matrix(bingto[,-1])
  gto[gto=='par1par1'] <- 'AA' # ref
  gto[gto=='par1par2'] <- 'AB'
  gto[gto=='par2par2'] <- 'BB'
  # final NA filter
  nas <- apply(gto, 1, function(x) sum(is.na(x)))/len(lines)
  gto <- gto[nas < maxNA,]
  bcontig <- bcontig[nas < maxNA]
  bmarkers <- bmarkers[nas < maxNA]
  bpos <- bpos[nas < maxNA]
  gto <- data.frame(cbind(bmarkers, gto))
  names(gto) <- c('marker', lines)
  pmap <- data.frame(marker = bmarkers, chr = bcontig, position = bpos)
  pheno = as.data.frame(matrix(rep(0, len(lines)), 1))
  names(pheno) = lines
  pheno <- cbind('dummy', pheno)
  
  # require(qtl2)
  # write_control_file(output_file = 'cross.F2.yaml', crosstype = 'f2', geno_file = 'geno.F2.csv', pheno_file = 'pheno.csv', pmap_file = 'pmap.F2.csv', alleles = c('A','B'), geno_codes = c(AA=1L, AB=2L, BB=3L), na.strings = "NA", overwrite = T)
  write.csv(gto, sprintf('geno.F2.bin%s.csv', bin), row.names=F, quote=F)
  write.csv(pmap, sprintf('pmap.F2.bin%s.csv', bin), row.names=F, quote=F)
  write.csv(pheno, 'pheno.csv', row.names=F, quote=F)
  
  # riself code, assuming two-class priors
  # this will fail if the first 10 lines do not collectively show all haplotypes
  diplos = sort(levels(factor(unlist(bingto[,2:11]))))
  if(length(diplos)==3) diplos = diplos[c(1,3)]
  
  bingto <- do.call(rbind, bingt)
  bmarkers = rownames(bingto)
  bmarkers[grep('.', bmarkers, fixed=T, invert=T)] <- paste0(bmarkers[grep('.', bmarkers, fixed=T, invert=T)], '.1')
  bcontig = unlist(lapply(strsplit(bmarkers, '.', fixed=T), '[[', 1))
  bpos = bingto[,1]
  lines = names(gtout)[-(1:2)]
  gto = as.matrix(bingto[,-1])
  # dump proportions
  gtop <- data.frame(AA=apply(gto, 1, function(x) sum(x=='par1par1', na.rm=T)), AB=apply(gto, 1, function(x) sum(x=='par1par2', na.rm=T)), BB=apply(gto, 1, function(x) sum(x=='par2par2', na.rm=T)), nas=apply(gto, 1, function(x) sum(is.na(x))))
  gtop$marker = rownames(gtop)
  write.csv(gtop, sprintf('geno.props.bin%s.csv', bin), row.names=F, quote=F)
  gto[gto==diplos[1]] <- 'A'
  gto[gto==diplos[2]] <- 'B'
  gto[!gto %in% c('A','B')] = NA
  # final NA filter
  nas <- apply(gto, 1, function(x) sum(is.na(x)))/len(lines)
  gto <- gto[nas < maxNA,]
  bcontig <- bcontig[nas < maxNA]
  bmarkers <- bmarkers[nas < maxNA]
  bpos <- bpos[nas < maxNA]
  gto <- data.frame(cbind(bmarkers, gto))
  names(gto) <- c('marker', lines)
  pmap <- data.frame(marker = bmarkers, chr = bcontig, position = bpos)
  # write_control_file(output_file = 'cross.riself.yaml', crosstype = 'riself', geno_file = 'geno.riself.csv', pheno_file = 'pheno.csv', pmap_file = 'pmap.riself.csv', alleles = c('A', 'B'), geno_codes = c(A=1L, B=2L), na.strings = "NA", overwrite = T)
  write.csv(gto, sprintf('geno.riself.bin%s.csv', bin), row.names=F, quote=F)
  write.csv(pmap, sprintf('pmap.riself.bin%s.csv', bin), row.names=F, quote=F)
  
  save(gtout, gtinterpol, binnedbp, bingto, bcontig, bmarkers, bpos, contigLengths, file = sprintf('riself.rqtl.dat.bin%s.rda', bin))
}

for(i in indivs) {
  if(job %in% c(0,1)) write_hmm_data(i, dir, indir, contigs, NP=np, v=v)
  if(job %in% c(0,2)) fit_hmm(i, indir, outdir, contigs, NP=np, v=v)
}

if(job %in% c(0,3)){
  if(len(toplot)>0) plotAncestry(indivs, outdir, toplot, np=np, v=v)
  dumpRQTL(indivs, outdir, bin=50, np=np, v=v, plot=T)
}

