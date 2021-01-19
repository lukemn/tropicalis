#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = T)
WD   = args[1]
np   = as.integer(args[2])
pref = basename(dirname(WD))

require(qtl, quietly = T)
require(parallel, quietly = T)
require(ggplot2, quietly = T)
require(data.table, quietly = T)

# functions to:
# (1) filter genetic markers (binned pseudomarkers from RIL GBS)
# (2) make a genetic map
# (3) assess (and attempt to maximise) physical/genetic concordance
# (4) output some stats to stdout
#     a map (r/qtl)
#     genetically ordered and oriented (where possible) assembly contigs/scaffolds.
# The main functions are commented below.

# read tidy bug where single marker/scaf
source('scripts/my_read_cross.R')
nameSplitChar="." # markers are Ctg%%%.# (where # is bin number)
maxsim = 0.95     # maximum line similarity
maxLineNA = 0.6   # maximum missing data per line
maxMarkerNA = 0.3 # maximium missing data per marker
nsamples = 50     # wihin linkage group ordering iterations

### utils ###
updateErr <- function(cross, chrom){
  estErr <- function(i, cross, chrom, np=1){
    tempmap <- est.map(cross, chr=chrom, error.prob=i, map.function = 'morgan', maxit = 1e6, n.cluster=np)
    sum(sapply(tempmap, attr, "loglik"))
  }
  if(missing(chrom)) chrom=1:nchr(cross)
  optimise(f=estErr, interval=c(5e-2, 1e-4), cross=cross, chrom=chrom, np=1, maximum=T)$max
}

getTermini = function(s, y, ori){
  # add scaffold termini to a block from the fasta index
  # marked 'end' in dupe column.
  slen = subset(fai, scaf==s)$len
  if(ori){
    n = y[1,]; n$start = n$end+1; n$end = slen; if('ix' %in% names(n)) n$ix=NA
    n$dupe='end'
    y <- rbind(n, y)
    n = y[nrow(y),]; n$end = n$start-1; n$start = 0; if('ix' %in% names(n)) n$ix=NA
    n$dupe='end'
    y <- rbind(y, n)
  } else {
    n = y[1,]; n$end = n$start-1; n$start = 0; if('ix' %in% names(n)) n$ix=NA
    n$dupe='end'
    y <- rbind(n, y)
    n = y[nrow(y),]; n$start = n$end+1; n$end = slen; if('ix' %in% names(n)) n$ix=NA
    n$dupe='end'
    y <- rbind(y, n)
  }
  y
}

getOri = function(y, gte=F, tol=0, dominant=F){
  # vector orientation
  # ascending=F
  y <- y[!is.na(y)]
  if(length(y)<2) return(NA)
  if(gte) mdiff = diff(y) >= -tol else mdiff = diff(y) > -tol
  if(all(mdiff)) return(F)
  if(all(!mdiff)) return(T)
  if(dominant){
    ot = table(mdiff)
    return(as.logical(ifelse(ot[1]>ot[2], T, F)))
  }
  return(NA)
}

startingMap <- function(WD, pref, XO_MAD_quantiles=c(0.9, 0.95, 0.99, 1)){
  
  # Filter the data and form an initial map, iterating marker ordering within 
  # linkage groups and taking the most likely order across runs.
  # Input: 
  #   pseudomarkers and minimal rqtl files
  #   marker positions and fasta index
  # Output, saved to input directory:
  #   cross object with genetically unique markers ordered within 6 linkage groups.
  #   duplicate list
  #   marker and contig info
  #   ... over a range of stringencies (removing lines based on median abs deviation quantile in XO number)
  #   Orphan markers are saved in the initial dump after basic filtering and LG formation, 
  #   but are discarded in the final ordered maps.
  
  cat("##\nReading in cross, assembly and marker info\n##\n")
  print(date())
  # to rerun filtering/LG formation again after dropping lines, 
  # list individuals in `todrop`
  dropf = sprintf('%s/todrop', WD)
  if(file.exists(dropf)){
    todrop = scan(dropf, what='character', quiet = T)
    cat(sprintf('%s individuals to drop\n', length(todrop)))
  }
  cross <- my.read.cross(dir = WD, format='tidy', genfile='geno.riself.bin50.csv', mapfile = 'pmap.riself.bin50.csv', phefile = 'pheno.csv', na.strings = 'NA', genotypes = c('A', 'B'), alleles = c('A', 'B'), crosstype = 'riself')
  if(file.exists(dropf)){
    todrop = todrop[todrop %in% cross$pheno$id]
    if(length(todrop)>0) cross = subset(cross, ind=paste0('-', todrop))
  }
  
  fai <- fread(sprintf('%s/reorderedSeq.csv', WD))
  names(fai) <- c('len', 'fa', 'scaf')
  mpos <- fread(sprintf('%s/bin50.markerpos.txt', WD))
  mpos$marker <- mpos$dupe <- paste(mpos$contig, mpos$b, sep=nameSplitChar)
  
  # terminal markers are dropped in mpos
  mar = colnames(pull.geno(cross))
  mterm = mar[!mar %in% mpos$marker]
  cross <- drop.markers(cross, markers = mterm)
  cat(sprintf('dropping %s terminal markers\n', length(mterm)))
  
  #####
  cat("##\nDedupe and filtering\n##\n")
  dupes <- findDupMarkers(cross, exact.only = F)
  tagscaf <- unique(unlist(lapply(strsplit(names(dupes), nameSplitChar, fixed = T), '[[', 1)))
  taglen <- sum(fai$len[fai$scaf %in% tagscaf])
  dupscaf <- unique(unlist(lapply(strsplit(unlist(dupes), nameSplitChar, fixed = T), '[[', 1)))
  duplen <- sum(fai$len[fai$scaf %in% dupscaf[!dupscaf %in% tagscaf]])
  cscaf = unique(unlist(lapply(pull.map(cross), function(x) tstrsplit(names(x), nameSplitChar, fixed = T)[[1]])))
  allscaf = unique(c(tagscaf, dupscaf, cscaf))
  fullen = sum(fai$len[fai$scaf %in% allscaf])
  cat(sprintf('collapse %s dupes to %s tags (spanning %.2f Mb), %.2f covered by all %s marker scaffolds\n', length(unlist(dupes)), length(dupes), taglen/1e6, fullen/1e6, length(allscaf)))
  cross <- drop.markers(cross, unlist(dupes))
  summary(cross)
  
  ## remove similar lines
  cg <- comparegeno(cross)
  # hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes"); rug(cg[lower.tri(cg)])
  wh <- which(cg > maxsim, arr=TRUE)
  wh <- wh[wh[,1] < wh[,2],]
  g <- pull.geno(cross)
  # lapply(1:nrow(wh), function(x) table(g[wh[x,1],], g[wh[x,2],]))
  if(nrow(wh)>0){
    cat(sprintf('dropping %s lines due to high similarity (>%.3f)\n', nrow(wh), maxsim))
    wh = data.frame(wh)
    wh$rowNA = sapply(wh$row, function(i) sum(is.na(g[i,])))
    wh$colNA = sapply(wh$col, function(i) sum(is.na(g[i,])))
    wh$todrop = apply(wh, 1, function(x) x[1:2][which.max(x[3:4])])
    cross <- subset(cross, ind=-wh$todrop)
    print(summary(cross))
  }
  
  ## remove lines based on missing data
  tokeep = ntyped(cross)>(sum(nmar(cross))*(1-maxLineNA))
  if(length(tokeep) < nind(cross)){
    cat(sprintf('dropping %s lines due to missing data (>%.2f)', sum(!(tokeep)), maxLineNA))
    cross <- subset(cross, ind=tokeep)
    print(summary(cross))
  }
  
  ## marker segregation distortion
  gt <- geno.table(cross); gt <- gt[order(gt$P),]
  mingt = apply(gt[,3:4], 1, min)
  dropgt = gt[mingt <= 3 | gt$P.value < 1e-20,]
  if(nrow(dropgt)>0){
    cat(sprintf('dropping %s markers (%s bp) with strong segregation distortion (min class <3 | -log10 p > 20)\n', nrow(dropgt), sum(fai$len[fai$scaf %in% dropgt$chr])))
    cross <- drop.markers(cross, rownames(dropgt))
    print(head(dropgt, 10))
  }
  
  ## marker missing data
  mna <- ntyped(cross, 'mar')
  todrop = names(mna)[mna < (nind(cross)*(1-maxMarkerNA))]
  keepchr = sapply(strsplit(names(mna)[!names(mna) %in% todrop], nameSplitChar, fixed = T), '[[', 1)
  dropchr <- unique(sapply(strsplit(todrop, nameSplitChar, fixed = T), '[[', 1))
  if(length(dropchr)>0) {
    mrep = minrep = sum(!dropchr %in% unique(keepchr))
    minrep = min(sapply(dropchr, function(x) sum(keepchr==x)))
    lostd = sum(subset(fai, scaf%in%dropchr[!dropchr %in% keepchr])$len)
    sprintf('dropping %s markers due to missing data (>%.2f), %s sole scaffold markers (%.2f Kb)', length(todrop), maxMarkerNA, mrep, lostd/1000)
    cross <- drop.markers(cross, todrop)
    summary(cross)
  }
  
  #####
  cat("##\nForm linkage groups\n##\n")
  ## recombination fraction v LOD
  cross <- est.rf(cross)
  chk = checkAlleles(cross, threshold=2)
  if(length(chk)>0) {print(chk); cross <- drop.markers(cross, chk$marker)}
  
  # Plot RF ~ LOD
  i=names(cross$geno)
  # i=c(1:6)
  rf <- pull.rf(cross, chr = i)
  lod <- pull.rf(cross, what="lod", chr = i)
  f = as.numeric(lod) > 2
  rfd = data.frame(x = as.numeric(rf)[f], y = as.numeric(lod)[f])
  p <- ggplot(rfd, aes(x, y)) + geom_point(alpha=0.5, stroke=0) + theme_classic() + labs(x="Recombination fraction", y="LOD score") + geom_rug()
  ggsave(sprintf('%s/RF_LOD.png', WD))
  
  # check LOD ~ LGs
  lodr = range(round(quantile(rfd$y[rfd$x<0.25], c(0.5, 0.85), na.rm=T)))
  lodr = n = sixes = minlg = maxdross = lodr[1]:lodr[2]
  for(i in seq_along(lodr)){
    m=lodr[i]
    lg <- formLinkageGroups(cross, min.lod=m)
    lgt = table(lg[,2])
    n[i] = length(lgt)
    sixes[i] = sum(lgt[1:6])/sum(lgt)
    minlg[i] = min(lgt[1:6])
    maxdross[i] = ifelse(length(lgt)>6, max(lgt[-(1:6)]), NA)
  }
  
  # choose the lowest value that gives 6 LGs,
  # else the value that maximises the number of markers in the first 6 LGs
  lodlg = data.frame(lod=lodr, n, sixes, minlg, maxdross)
  print(lodlg)
  stopifnot(max(lodlg$n)>6)
  if(sum(lodlg$n==6)>0){
    minLOD = min(lodlg$lod[lodlg$n==6])  
  } else {
    lodlg = subset(lodlg, n>6 & minlg>20)
    minLOD = lodlg$lod[which.max(lodlg$sixes)]
  }
  cat(sprintf('\nusing minimum LOD %s\n', minLOD))
  
  # contig ~ LGs
  lg <- formLinkageGroups(cross, min.lod=minLOD)
  ctglg = subset(data.frame(table(lg), stringsAsFactors = F), Freq>0)
  lgt = table(ctglg$LG)
  
  cat(sprintf('%s markers in %s non-singleton LGs\n', sum(ctglg$Freq), length(lgt[lgt>1])))
  ctglg = table(ctglg$origchr, ctglg$LG)
  splitctg = apply(ctglg[,1:6], 1, sum)>1
  badctgs = sum(splitctg)
  if(badctgs>0){
    cat(sprintf('%s contigs span >1 primary LGs\nbad contigs:\n', badctgs))
    print(table(lg)[splitctg,1:10,drop=F])
  }
  
  # dump starting map -------------------------------------------------------
  cat("##\nSaving filtered cross data\n##\n")
  ## initial order markers and estimate map distances
  cross <- formLinkageGroups(cross, min.lod=minLOD, reorgMarkers=TRUE)
  summary(cross)
  
  print("ML error rate ... ")
  gterr = updateErr(cross, chr=1:6)
  print(gterr)
  
  cross = orderMarkers(cross, map.function = 'morgan', error.prob = gterr, window = 3)
  xos = countXO(cross)
  cat(sprintf('XOs (median %s)\n', median(xos)))
  print(summary(xos))
  mads = abs(xos-median(xos))
  cat('XOs > 99% MAD\n')
  print(xos[mads > quantile(mads, 0.99)])
  print(summaryMap(cross))
  
  save(cross, fai, mpos, dupes, gterr, file = sprintf('%s/%s_rimap_1.rda', WD, pref))
  
  
  # iterate marker ordering-------------------------------------------------------
  cat("##\nIterating within LG marker ordering\n##\n")
  # dropping suspect lines at 90, 95, 99% percentile in XO MAD
  
  iterateMarkerOrdering <- function(cross, gterr=5e-3, win=5){
    
    runiter <- function(cross){
      so = mclapply(1:nsamples, mc.cores = np, function(i) orderMarkers(cross, chr=1:6, error.prob=gterr, map.function='morgan', window=win))
      maps = lls = xos = lapply(1:6, function(x) list())
      linexos = data.frame()
      for(chrom in 1:6){
        lls[[chrom]] = unlist(lapply(so, function(i) attr(i$geno[[chrom]]$map, 'loglik')))
        xos[[chrom]] = unlist(lapply(so, function(i) sum(countXO(i, chr=chrom))))
        maps[[chrom]] = lapply(so, function(i) i$geno[[chrom]])
        linexos = rbind(linexos, do.call(rbind, lapply(so, function(i) {
          xo=countXO(i, chr=chrom)
          data.frame(line=names(xo), xo, chrom, row.names = NULL)
        })))
      }
      
      ix = unlist(lapply(lls, which.max))
      bestlls = unlist(lapply(1:6, function(i) unlist(lls[[i]])[ix[i]]))
      bestmaps = lapply(1:6, function(i) maps[[i]][[ix[i]]])
      lxo = aggregate(data = linexos, xo~line, mean)
      print(summary(lxo$xo))
      return(list(bestlls, bestmaps, lxo))
    }
    
    for(q in sort(XO_MAD_quantiles, decreasing = T)){
      
      o = runiter(cross)
      bestlls = o[[1]]
      bestmaps = o[[2]]
      lxo = o[[3]]
      
      if(q < 1) {
        mads = abs(lxo$xo - median(lxo$xo))
        dropl = as.character(lxo$line[mads > quantile(mads, q)])
        writeLines(dropl, con=sprintf('%s/%s_%s.drop', WD, pref, q*100))
        cross = subset(cross, ind=paste0('-', dropl))
        o = runiter(cross)
        bestlls = o[[1]]
        bestmaps = o[[2]]
        lxo = o[[3]]
      }
      
      rcross=cross
      cliks = sapply(1:6, function(i) attr(cross$geno[[i]]$map, 'loglik'))
      for(i in 1:6) if(bestlls[i] > cliks[i]) rcross$geno[[i]] <- bestmaps[[i]]
      bliks = sapply(1:6, function(i) attr(rcross$geno[[i]]$map, 'loglik'))
      cat(sprintf('improvement in LL after %s iterations (%s XO quantile, %s lines):\n', nsamples, q, nind(cross)))
      print(cliks-bliks)
      print(summaryMap(cross)[1:6,])
      print(summaryMap(rcross)[1:6,])
      save(rcross, dupes, fai, mpos, file = sprintf('%s/%s_q%s_bestOrders_LL.rda', WD, pref, q*100))
      
      cross = rcross
    }
  }
  
  iterateMarkerOrdering(cross)
  
}

geneticPhysicalConcordance <- function(cross, dupeList, fai, mpos, pref, gterr=5e-3, maxPoG_LOD= -3){
  
  # Check agreement between physical and genetic data, and
  # reorder the map where changes can be made without strongly reducing the likelihood (>= maxPoG_LOD)
  # Add in duplicate markers and attempt to order them.
  # Input: 
  #   cross object, duplicate list, assembly fasta index
  # Output:
  #   as above, with harmonised (where possible) data.frame of genetic map, 
  #   and genetic map with duplicate markers inserted.
  #   iterated as in startingMap for line XO filtering (XO_MAD_quantiles)
  
  # dupe list > df
  duped <- data.frame(dupe = unlist(dupeList), stringsAsFactors = F)
  duped$marker <- rep(names(dupeList), sapply(dupeList, len))
  duped$scaf<- unlist(lapply(strsplit(duped$dupe, nameSplitChar, fixed = T), '[[', 1))
  duped <- merge(duped, mpos[,c('dupe', 'start', 'end')], sort=F)
  duped <- duped[order(duped$marker, duped$start),]
  
  getMarkerPos <- function(cross, mpos, lgs=1:6){
    
    gmap <- pull.map(cross, chr = lgs)
    mapn = names(unlist(gmap))
    mapdf <- data.frame(genetic=as.numeric(unlist(gmap)), 
                        ix = unlist(sapply(gmap, function(i) 1:length(i))),
                        marker = gsub('^[1-9].', '', mapn),
                        lg = tstrsplit(mapn, '.', fixed = T)[[1]],
                        stringsAsFactors = F
    )
    mapdf <- merge(mapdf, mpos[,c('marker', 'start', 'end')], sort=F)
    mapdf$scaf <- unlist(lapply(strsplit(mapdf$marker, nameSplitChar, fixed = T), '[[', 1))
    mapdf
  }
  
  # order by physical distance within each scaffold, where not conflicting with genetic order
  # adds column 'ori', which is 1 if ordered unambiguously, 0 if unoriented, or if conflicted 
  # a proportion of markers supporting the dominant orientation..
  orderByPhysWithinScaf <- function(d){
    
    od = do.call(rbind, lapply(split(d, d$lg), function(lg){
      
      lgg = aggregate(data = lg, genetic ~ scaf, median)
      lgg = lgg[order(lgg$genetic),]
      lgg$scafo = 1:nrow(lgg)
      lg = plyr::join(lg, lgg[,c('scaf', 'scafo')])
      
      oo = do.call(rbind, lapply(split(lg, lg$scaf), function(x) {
        
        # input is sorted in ascending order by genetic distance, then marker start
        # x = subset(d, scaf=='Ctg019')
        
        scafx = subset(x, dupe==F)
        dupex = subset(x, dupe!=F)
        sx = x$scaf[1]
        nmark = nrow(scafx)
        ori = getOri(scafx$start)
        # print(sx)
        
        if(nmark>1){
          # can orient
          if(!is.na(ori)){
            if(ori){
              # rev, monotonic
              xo = do.call(rbind, lapply(split(x, x$genetic), function(y) y[order(y$start, decreasing = T),]))
            } else {
              # fwd, monotonic
              xo = x
            }
            ox = cbind(getTermini(sx, xo, ori), ori=1)
          } else {
            # non mono. try rev, else take most common
            xo = do.call(rbind, lapply(split(x, x$genetic), function(y) y[order(y$start, decreasing = T),]))
            ori = getOri(subset(xo, dupe==F)$start, gte=T)
            if(is.na(ori)){
              # still suspect
              ori1 = table(diff(scafx$start)>0)
              ori2 = table(diff(subset(xo, dupe==F)$start)>0)
              oris = ori1; if(max(ori2)>max(oris)) oris = ori2
              stopifnot(length(oris)==2)
              if(oris[1]>oris[2]){
                # rev
                xo = xo[order(xo$start, decreasing = T),]
                ox = cbind(getTermini(sx, xo, T), ori=oris[1]/sum(oris))
              } else {
                if(oris[2]>oris[1]){
                  # fwd
                  xo = x[order(x$start),]
                  ox = cbind(getTermini(sx, xo, F), ori=oris[2]/sum(oris))
                } else {
                  if(oris[2]==oris[1]){
                    # ??
                    ox = cbind(getTermini(sx, x, F), ori=0)
                  }
                }  
              }
            } else {
              ox = cbind(getTermini(sx, xo, ori), ori=1)
            }
          }
        } else {
          # <=1 marker in the map, cannot orient from genetic data
          ox = cbind(getTermini(sx, x, F), ori=0)
        }
        # discard termini for now
        ox = subset(ox, dupe!='end')
        ox$o = 1:nrow(ox)
        ox
      }))
      oo = oo[order(oo$scafo, oo$o),]
      oo$o = 1:nrow(oo)
      oo
    }))
    rownames(od) = NULL
    od
    
  }
  
  cat(sprintf('%s Test full reordering by physical distance\n', pref))
  for(i in 1:2){
    # genetic map data.frame with physical info and dupe status.
    mapdf <- getMarkerPos(cross, mpos)
    
    # Rounding of genetic positions gives preference to physical information.
    # This should be fine unless there are very large numbers of XOs and clean data
    d = mapdf; d$dupe=F; d$marker_scaf = d$scaf; d$genetic = round(d$genetic, 3)
    d = d[order(d$lg, d$genetic, d$start),]
    # initial order ('ix')
    od = orderByPhysWithinScaf(d)
    od = od[order(od$lg, od$ix),]
    sus = subset(od, ori<1)
    suscaf = unique(sus[,c('scaf', 'ori')])
    cat(sprintf('%s: Of %s contigs in the map, %s are unorientable, %s are discordant (spanning %s and %s bp)\n', 
                pref, length(unique(od$scaf)), sum(suscaf$ori==0), sum(suscaf$ori>0), sum(subset(fai, scaf %in% suscaf$scaf[suscaf$ori==0])$len), sum(subset(fai, scaf %in% suscaf$scaf[suscaf$ori>0])$len)))
    
    # check likelihoods for all changes
    # the oriented order is in column 'o'
    # ggplot(od, aes(ix, o)) + geom_point() + facet_grid(.~lg, scales='free')
    cat(sprintf('%s breaks: initial order vs. oriented order\n', pref))
    print(lapply(1:6, function(i) rle(diff(od$o[od$lg==i]))))
    suggested = sapply(split(od, od$lg), function(x) sum(x$ix!=x$o))
    cat(sprintf('%s likelihoods for new order\n', pref))
    newliks = mclapply(1:6, mc.cores = np, function(i) cbind(chr=i, compareorder(cross, chr=i, order = order(od$o[od$lg==i]), error.prob = gterr, map.function = 'morgan')) )
    ok = sapply(newliks, function(x) {
      if(suggested[as.numeric(x$chr[1])]>0){
        if(x$LOD[2] > maxPoG_LOD) {
          cat(sprintf("LG %s: LOD %.3f, map length %.3f, order accepted\n", x$chr[1], x$LOD[2], x$length[2]-x$length[1]))
          T
        } else {
          cat(sprintf("\tLG %s: LOD %.3f, map length %.3f, order rejected\n", x$chr[1], x$LOD[2], x$length[2]-x$length[1]))
          F
        }
      } else {
        cat(sprintf("LG %s: no change\n", x$chr[1]))
        F
      }
    })
    
    # apply any changes
    for(i in 1:6) {
      if(ok[[i]]){
        oi = switch.order(cross, chr=i, order=order(od$o[od$lg==i]), error.prob=gterr, map.function = 'morgan')
        cross$geno[[i]] <- oi$geno[[i]]
      }
    }
  }
  
  mapdf <- getMarkerPos(cross, mpos); mapdf$genetic = round(mapdf$genetic, 3)
  mapdf$o = mapdf$ix = 1:nrow(mapdf)
  
  cat(sprintf('%s Test case-by-case reordering by physical distance\n', pref))
  orderByScafWithinLG <- function(d){
    
    od = do.call(rbind, lapply(split(d, d$lg), function(x){
      
      # x = split(mapdf, mapdf$lg)[[6]]
      rownames(x) = NULL
      x$ix = 1:nrow(x)
      lgx = x$lg[1]
      likto = -0.5 # move markers even if likelihood is slightly negative (LOD)
      tol = 1e-2   # else, move even if against genetic distance (tol cM)
      rl = rle(x$scaf)
      rlt = table(rl$values)
      disco = rlt[rlt>1]
      jumps = sapply(names(disco), function(i) diff(range(which(rl$values==i))))
      dscaf = names(disco)[order(jumps)]
      k=1; it=1
      while(length(dscaf)>0){
        cix = which(rl$values==dscaf[k])
        if(length(cix)>2) cix=cix[1:2]
        inner = rl$value[(min(cix)+1):(max(cix)-1)]
        ixes = which(x$scaf==dscaf[k])
        ixes = ixes[ixes <= sum(rl$lengths[1:max(cix)])]
        ixg = diff(ixes)>1
        ixl = ixes[1:which(ixg)]; ixr = ixes[(which(ixg)+1):length(ixes)]
        l = rev(ixl)[1]; ll = ixl[1]
        r = ixr[1]; rr = rev(ixr)[1]
        inix = (l+1):(r-1)
        lg = x$genetic[l]; rg = x$genetic[r]; cg = x$genetic[l+1]
        igd = unique(x$genetic[l:r])
        # genetic distance to flanking markers
        lgen = cg-lg; rgen = rg-cg
        # genetic distance to end of flanking contigs
        lgenc = cg-min(x$genetic[ixl])
        rgenc = max(x$genetic[ixr])-cg
        # x[max(c(1,ll-2)):min(c(nrow(x),rr+2)),]
        # x[(l-1):(r+1),]
        
        # if an outer is terminal, then see if we can clump it, else deal with the inners
        if(length(inner)>1 & (min(ixl)==1 | max(ixr)==nrow(x))){
          cat(sprintf('Terminal case\n'))
          cat(sprintf('LG %s genetic conflict %s-%s-%s: %.3f-%.3f<->%.3f-%.3f cM\n', lgx, dscaf[k], inner, dscaf[k], lgenc, lgen, rgen, rgenc))
          # small > large
          nl = length(ixl)
          nr = length(ixr)
          if(nr>=nl){
            lo = c(inix, ixl, r:nrow(x)); while(min(lo)>1) lo = c(min(lo)-1, lo)
            lik = compareorder(cross, chr=lgx, order = lo, map.function = 'morgan')[2,1]
            if(lik > likto){
              cat(sprintf('\tjoin R lik %.3f\n', lik))
              x = x[lo,]
            } else {
              cat(sprintf('\tfailed: will not join R lik %.3f\n', lik))
            }
          } else {
            lo = c(1:l, ixr, inix); while(max(lo) < nrow(x)) lo = c(lo, max(lo)+1)
            lik = compareorder(cross, chr=lgx, order = lo, map.function = 'morgan')[2,1]
            if(lik > likto){
              cat(sprintf('\tjoin L lik %.3f\n', lik))
              x = x[lo,]
            } else {
              cat(sprintf('\tfailed: will not join L lik %.3f\n', lik))
            }
          }
        } else {
          if(length(igd)>1){
            # conflict with genetic order
            cat(sprintf('LG %s genetic conflict %s-%s-%s: %.3f-%.3f<->%.3f-%.3f cM\n', lgx, dscaf[k], inner, dscaf[k], lgenc, lgen, rgen, rgenc))
            
            if(sum(x$scaf %in% inner)>1){
              # not orphan, test moving toward other markers
              inix = (l+1):(r-1)
              oix = x$scaf %in% inner; oix[inix]=F
              go = mean(x$genetic[oix])
              ing = mean(x$genetic[inix])
              if(any(is.na(c(go, ing)))) {cat('\terr. bailing\n'); break}
              lo = 1:nrow(x)
              if(ing > go){
                # lo[ll:(r-1)] = c(inix, ll:l)
                lo[inix:rr] = rr:inix
                lik = compareorder(cross, chr=lgx, order = lo, map.function = 'morgan')[2,1]
                if(lik > likto){
                  cat(sprintf('\tmove R lik %.3f\n', lik))
                  x = x[lo,]
                } else {
                  cat(sprintf('\tfailed: will not move R lik %.3f\n', lik))
                }
              } else {
                lo[(l+1):rr] = c(r:rr, inix)
                lik = compareorder(cross, chr=lgx, order = lo, map.function = 'morgan')[2,1]
                if(lik > likto){
                  cat(sprintf('\tmove L lik %.3f\n', lik))
                  x = x[lo,]
                } else {
                  cat(sprintf('\tfailed: will not move L lik %.3f\n', lik))
                }
              }
            } else {
              cat(sprintf('\torphan %s %s-%s-%s\n', k, dscaf[k], inner, dscaf[k]))
              cat(sprintf('\t%s-%s : %s-%s\n', ll, l, r, rr))
              move=F
              if(lgenc<tol & rgenc>tol) {cat('\tmove L on G\n'); exc = rbind(x[(l+1):(r-1),], x[ll:l,], x[r:rr,]); move=T}
              if(rgenc<tol & lgenc>tol) {cat('\tmove R on G\n'); exc = rbind(x[ll:l,], x[r:rr,], x[(l+1):(r-1),]); move=T}
              if(move) {
                x[ll:rr,] <- exc
              } else {
                cat(sprintf('\tfailed: nearest end %s > override %s cM\n', min(c(lgenc, rgenc)), tol))
              }
            }
          }
        }
        rl <- rle(x$scaf)
        rlt = table(rl$values)
        disco = rlt[rlt>1]
        dscafi = dscaf
        dscaf = NULL
        if(length(disco)>0){
          jumps = sapply(names(disco), function(i) min(diff(which(rl$values==i))))
          dscaf = names(disco)[order(jumps)]
        }
        if(sum(dscafi %in% dscaf)==length(dscaf)) break
        it = it+1
        k = min(length(dscaf), k+1)
      }
      x$o = 1:nrow(x)
      x
    }))
    rownames(od) = NULL
    od
  }
  
  od = orderByScafWithinLG(mapdf)
  od = od[order(od$lg, od$ix),]
  
  # reorder the cross, check full lik
  do.call(rbind, mclapply(1:6, mc.cores = 6, function(i) cbind(chr=i, compareorder(cross, chr=i, order = order(od$o[od$lg==i]), error.prob = gterr, map.function = 'morgan'))))
  
  newo <- lapply(c(1:6), function(i) {
    oi = switch.order(cross, chr=i, order=order(od$o[od$lg==i]), error.prob=gterr, map.function = 'morgan')
    oi$geno[[i]]
  })
  for(i in 1:6) cross$geno[[i]] <- newo[[i]]
  
  mapdf = od; mapdf = mapdf[order(mapdf$lg, mapdf$o),]
  # cross <- est.rf(cross)
  # plotRF(cross)
  # plot.map(cross)
  print(summary.map(cross))
  # (gterr = updateErr(cross))
  handle = sprintf('%s/%s_rimap_2.rda', WD, pref)
  cat(sprintf('Saving cross and map df to %s\n', handle))
  save(cross, dupes, duped, gterr, fai, mpos, mapdf, file = handle)
  
  expandMap <- function(cross, dupeList){
    
    addDupes <- function(cross, dupeList, mpos){
      # merge map markers with duplicates
      
      # get marker : dupe physical coordinates from dupe list
      duped <- data.frame(dupe = unlist(dupeList), stringsAsFactors = F)
      duped$marker <- rep(names(dupeList), sapply(dupeList, len))
      duped$scaf<- unlist(lapply(strsplit(duped$dupe, nameSplitChar, fixed = T), '[[', 1))
      duped <- merge(duped, mpos[,c('dupe', 'start', 'end')], sort=F)
      duped <- duped[order(duped$marker, duped$start),]
      
      # get marker physical and genetic coordinates from cross
      mapdf <- getMarkerPos(cross, mpos)
      mapdf$marker_scaf <- mapdf$scaf
      
      # merge and sort
      mapd = merge(mapdf[,c('marker', 'genetic', 'lg', 'marker_scaf')], 
                   duped[,c('dupe', 'marker', 'start', 'end', 'scaf')], 
                   sort=F)
      mapd <- rbind(mapd, cbind(mapdf[,c('marker', 'genetic', 'lg', 'start', 'end', 'scaf', 'marker_scaf')], dupe=F))
      mapd$genetic = round(mapd$genetic, 3)
      mapd <- mapd[order(mapd$lg, mapd$genetic, mapd$scaf, mapd$start),]
      
      mapd
    }
    
    orderDupes <- function(d, addTermini=T){
      
      # order markers within scaffolds, where possible
      # label ori = 
      #   1 (inferred from >1 marker, physical/genetic orientation of MARKERS consistent)
      #   0 (<1 marker, unknown)
      #   0>&<1 (inferred from >1 marker, physical/genetic orientation not consistent, 
      #          float indicates the number of markers supporting the dominant order)
      
      o = do.call(rbind, lapply(split(d, d$lg), function(lg){
        
        # preserve genetic order
        lgg = aggregate(data = lg, genetic ~ scaf, median)
        lgg = lgg[order(lgg$genetic),]
        lgg$scafo = 1:nrow(lgg)
        lg = plyr::join(lg, lgg[,c('scaf', 'scafo')])
        
        oo = do.call(rbind, lapply(split(lg, lg$scaf), function(x) {
          
          # input is sorted in ascending order by genetic distance, then marker start
          # x = subset(d, scaf=='Ctg022')
          
          x$scaf = as.character(x$scaf)
          scafx = subset(x, dupe==F)
          sx = x$scaf[1]
          # may be able to orient from duplicate genetic data
          # scafx = x[!duplicated(x$genetic),]
          if(nrow(scafx)==1 & length(unique(x$genetic)>1)) {
            udupe = subset(x, marker_scaf!=sx)
            scafx <- rbind(scafx, udupe[!duplicated(udupe$genetic),])
          }
          dupex = subset(x, dupe!=F)
          nmark = nrow(scafx)
          ori = getOri(scafx$start)
          # print(sx)
          
          if(nmark>1){
            # can orient
            if(!is.na(ori)){
              if(ori){
                # rev, monotonic
                xo = do.call(rbind, lapply(split(x, x$genetic), function(y) y[order(y$start, decreasing = T),]))
              } else {
                # fwd, monotonic
                xo = x
              }
              ox = cbind(getTermini(sx, xo, ori), ori=1, strand=ifelse(ori, '-', '+'))
            } else {
              # non mono. try rev, else take most common
              xo = do.call(rbind, lapply(split(x, x$genetic), function(y) y[order(y$start, decreasing = T),]))
              ori = getOri(subset(xo, dupe==F)$start, gte=T)
              if(is.na(ori)){
                # still suspect
                ori1 = table(diff(scafx$start)>0)
                ori2 = table(diff(subset(xo, dupe==F)$start)>0)
                oris = ori1; if(max(ori2)>max(oris)) oris = ori2
                stopifnot(length(oris)==2)
                if(oris[1]>oris[2]){
                  # rev
                  ox = cbind(getTermini(sx, xo, T), ori=oris[1]/sum(oris), strand='-')
                } else {
                  if(oris[2]>oris[1]){
                    # fwd
                    ox = cbind(getTermini(sx, x, T), ori=oris[2]/sum(oris), strand='+')
                  } else {
                    if(oris[2]==oris[1]){
                      # ??
                      ox = cbind(getTermini(sx, x, F), ori=0, strand=NA)
                    }
                  }  
                }
              } else {
                ox = cbind(getTermini(sx, xo, ori), ori=1, strand=ifelse(ori, '-', '+'))
              }
            }
          } else {
            # <=1 marker in the map
            ox = cbind(getTermini(sx, x, F), ori=0, strand=NA)
          }
          # discard termini, for now
          if(!addTermini) ox = subset(ox, dupe!='end')
          cbind(ox, o = 1:nrow(ox))
        }))
        oo = oo[order(oo$genetic, oo$scafo, oo$o),]
        oo$o = 1:nrow(oo)
        oo
      }))
      rownames(o) = NULL
      
      # check ordering of marker scaffolds (genetic blocks) within LGs
      # move dupes that interrupt map scaffolds to nearest end based on genetic distance
      # or physical if tied.
      
      oo = do.call(rbind, lapply(split(o, o$lg), function(x){
        
        # x = split(o, o$lg)[[2]]
        x$diff = c(1, abs(diff(x$start)))
        x$diff[-1][x$scaf[2:nrow(x)] != x$scaf[1:(nrow(x)-1)]] = 0
        x$cpos = cumsum(x$diff)
        tol = 1e-2
        # rl = rle(x$marker_scaf)
        rl = rle(x$scaf)
        rlt = table(rl$values)
        disco = rlt[rlt>1]
        jumps = sapply(names(disco), function(i) diff(range(which(rl$values==i))))
        dscaf = names(disco)[order(jumps)]
        k=1; it=1
        while(length(dscaf)>0 & it<3){
          cix = which(rl$values==dscaf[k])
          if(length(cix)>2) cix=cix[1:2]
          inner = rl$value[(min(cix)+1):(max(cix)-1)]
          # compare genetic/physical distances of inner to container
          # ixes = which(x$marker_scaf==dscaf[k])
          ixes = which(x$scaf==dscaf[k])
          ixes = ixes[ixes <= sum(rl$lengths[1:max(cix)])]
          ixg = diff(ixes)>1
          ixl = ixes[1:which(ixg)]; ixr = ixes[(which(ixg)+1):length(ixes)]
          l = rev(ixl)[1]; ll = ixl[1]
          r = ixr[1]; rr = rev(ixr)[1]
          lg = x$genetic[l]; rg = x$genetic[r]
          igd = unique(x$genetic[(l+1):(r-1)])
          if(length(igd)==1){
            lgen = igd-lg; rgen = rg-igd
            lphy = sum(x$diff[ll:l]); rphy = sum(x$diff[r:rr])
            cat(sprintf('%s LG %s case %s %s-%s-%s: %.3f<cM>%.3f %.3f<kb>%.3f\n', 
                        pref, x$lg[1], k, dscaf[k], inner, dscaf[k], lgen, rgen, lphy/1e3, rphy/1e3))
            cat(sprintf('\t%s-%s : %s-%s\n', ll, l, r, rr))
            # x[(ll-1):(rr+1),]
            # x[(l-1):(r+1),]
            # x[(l-10):(r+10),]
            xx = x[(ll):(rr),]
            inrl = length(rle(xx$scaf)$values)
            ax = xx[order(xx$scaf),]; ax = ax[order(ax$genetic),]; arl = length(rle(ax$scaf)$values)
            dx = xx[order(xx$scaf, decreasing = T),]; dx = dx[order(dx$genetic),]; drl = length(rle(dx$scaf)$values)
            
            if(arl < inrl) {
              x[ll:rr,] <- ax
            } else {
              if(drl < inrl){
                x[ll:rr,] <- dx
              } else {
                cat(sprintf('\t%s failed: no simple resolution\n', pref))
              }
            } 
          } else {
            cat(sprintf('\t%s failed: no simple resolution\n', pref))
          }
          
          rl <- rle(x$scaf)
          rlt = table(rl$values)
          disco = rlt[rlt>1]
          jumps = sapply(names(disco), function(i) diff(range(which(rl$values==i))))
          dscaf = names(disco)[order(jumps)]        
          it = it+1
          k = min(length(dscaf), k+1)
        }
        x
      }))
      rownames(oo) = NULL
      
      # final check to make sure contig termini are correct, update cumsum
      ooo = do.call(rbind, lapply(split(oo, oo$lg), function(x){
        
        # x = subset(oo, lg==2); y = subset(oo, scaf=='Ctg029')
        xx = do.call(rbind, lapply(split(x, x$scaf), function(y){
          term = subset(y, dupe=='end')
          bod = subset(y, dupe!='end')
          if(nrow(term)!=2) {print(y$scaf[1]); break}
          tgl = subset(fai, scaf==y$scaf[1])$len+1
          if(min(term$start)!=0) y$end[which.min(y$start)] = 0
          if(max(term$end)!=tgl) y$end[which.max(y$end)] = tgl
          if(max(term$start)!=(max(bod$end)+1)) y$start[y$start==max(term$start)] = max(bod$end)+1
          if(min(term$end)!=(min(bod$start)-1)) y$end[y$end==min(term$end)] = min(bod$start)-1
          if(y$start[1]==0){
            y$diff = c(0, abs(diff(y$end)))
          } else {
            y$diff = abs(c(diff(y$end), 0))
          }
          
          # check marker orientation consistent with dupes
          # set any dupe mismatch diffs to NA
          mori = getOri(y$start, dominant = T)
          dori = getOri(subset(y, dupe!=F)$start)
          # plot(y$genetic[-1], diff(y$start))
          # plot(y$start[-c(1, nrow(y))], y$genetic[-c(1, nrow(y))], col = factor(diff(y$start)>0), pch=20)
          if(is.na(dori)) {
            cat(sprintf('LG %s: %s dupe split over %s junction\n', y$lg[1], y$scaf[1], ifelse(mori, sum(diff(y$start)>0), sum(diff(y$start)<0))))
            if(mori){
              bix = which(diff(y$start)>0)
            } else {
              bix = which(diff(y$start)<0)
            }
            for(i in bix){
              print(y[(i-1):(i+1),-c(9:12)])
              cat('\n')
              y[i, c('diff', 'cpos')] = NA
            }
            # check for remaining sus large diffs
            lix = which(y$diff>1e6)
            if(length(lix)>0){
              cat(sprintf('LG %s: %s setting %s sus diffs > 1Mb to NA\n', y$lg[1], y$scaf[1], length(lix)))
              for(i in lix){
                print(y[(i-1):(i+1),-c(9:12)])
                cat('\n')
                y[i, c('diff', 'cpos')] = NA
              }
            }
          }
          y
        }))
        
        # update cumsum
        xx$diff[-1][xx$scaf[2:nrow(xx)] != xx$scaf[1:(nrow(xx)-1)]] = 0
        ix = xx$dupe=='end'; xx$diff[ix] = xx$end[ix]-xx$start[ix]
        xx$cpos[!is.na(xx$diff)] = cumsum(xx$diff[!is.na(xx$diff)])
      
        # yy = subset(xx, scaf=='Ctg005'); max(yy$end)-fai$len[fai$scaf==yy$scaf[1]]; sum(yy$diff)-fai$len[fai$scaf==yy$scaf[1]]
        # yy = subset(xx, scaf=='Ctg013'); max(yy$end)-fai$len[fai$scaf==yy$scaf[1]]; sum(yy$diff)-fai$len[fai$scaf==yy$scaf[1]]
        
        # yy = subset(xx, scaf=='Ctg010'); max(yy$end)-fai$len[fai$scaf==yy$scaf[1]]; sum(yy$diff, na.rm=T)-fai$len[fai$scaf==yy$scaf[1]]
        # yy = subset(xx, scaf=='Ctg013'); max(yy$end)-fai$len[fai$scaf==yy$scaf[1]]; sum(yy$diff)-fai$len[fai$scaf==yy$scaf[1]]
        
        
        ctgspan = sum(fai$len[fai$scaf %in% xx$scaf])
        cat(sprintf("LG %s, span vs true span %s\n\n", xx$lg[1], ctgspan-max(xx$cpos, na.rm=T)))
        # stopifnot(ctgspan==diff(range(xx$cpos))+xx$cpos[1])
        xx
      }))
      rownames(ooo) = NULL
      
      # sapply(split(d, d$lg), function(x) {rlt = table(rle(x$scaf)$values); length(rlt[rlt>1])})
      # sapply(split(o, o$lg), function(x) {rlt = table(rle(x$scaf)$values); length(rlt[rlt>1])})
      # sapply(split(oo, oo$lg), function(x) {rlt = table(rle(x$scaf)$values); length(rlt[rlt>1])})
      ooo      
    }
    
    fixBadDiffs_NIC58_flye <- function(omap){
      
      # flye manual fixes
      # check full (marker+dupe) orientations
      cori = sapply(split(omap, omap$scaf), function(x) getOri(x$start, gte=T))
      # 6 cases with some trouble (marker ordering correct, but pattern of duplicates not)
      # 1 resolved, check alignments for the other 5 (all < 0.5cM changes)
      table(cori, useNA = 'al')
      sus = names(which(is.na(cori))); print(sus)
      # i=6; sus[i]; subi=subset(omap, scaf==sus[i]); rownames(subi)=NULL; plot(subi$cpos[-1], diff(subi$start)); ix = which(omap$scaf==sus[i])
      
      # Ctg003 :  - strand, dupes Ctg003.12701-Ctg003.12401 (5764539-5643694) out of place at 0cM (sibs at 0.172 cM)
      
      # Ctg007 : - strand, 2 isolated dupe bubbles. 
      #           6.223cM Ctg007.9651 (3552561 3605794) > 6.641cM Ctg007.9751 (3647579 3687064)
      #           13.952cM Ctg007.8651 > Ctg007.8451,Ctg007.8351 (2697111-2641698) > 14.377 cM Ctg007.8601
      
      # Ctg010: - strand, 23.257cM dupes Ctg010.17351,17301,17251 (terminal at 2810169-3289292) < 22.834cM Ctg010.1
      
      # Ctg011: - strand, 37.535cM dupe Ctg011.6401 < 37.010cM marker Ctg011.6351
      
      # Ctg015: - strand, 30.244cM dupe Ctg015.1201 < 29.825 Ctg015.1101
      
      # Ctg029: + strand, single marker, single dupe, resolvable
      i=6; sus[i]; subi=subset(omap, scaf==sus[i]); rownames(subi)=NULL; plot(subi$cpos[-1], diff(subi$start)); ix = which(omap$scaf==sus[i])
      exc = omap[ix,]; exc = exc[order(exc$start, decreasing = T),]; exc$genetic[c(1, nrow(exc))] = sort(exc$genetic[c(1, nrow(exc))])
      exc$o = exc$o[order(exc$o)]; exc$diff = c(abs(diff(c(exc$start[1],exc$end[1]))), abs(diff(exc$start)))
      omap[ix,] = exc
      
      # update cumsum
      omap = do.call(rbind, lapply(split(omap, omap$lg), function(x) {x$cpos[!is.na(x$diff)] = cumsum(x$diff[!is.na(x$diff)]); x}))
      omap
    }
    
    fixBadDiffs_NIC58_flye_stitch1 <- function(omap){
      
      # flye manual fixes (joins from other pacbio ass)
      # check full (marker+dupe) orientations
      cori = sapply(split(omap, omap$scaf), function(x) getOri(x$start, gte=T))
      table(cori, useNA = 'al')
      sus = names(which(is.na(cori))); print(sus)
      
      # stitched flye manual fixes
      # 3 cases, all < 0.5cM
      i=2; sus[i]; subi=subset(omap, scaf==sus[i]); rownames(subi)=NULL; plot(subi$cpos[-1], diff(subi$start)); ix = which(omap$scaf==sus[i])
      # Ctg001: 2 dupes markers < 0.5cM out of order, leave
      # Ctg012: 2 terminal dupes switched, but < 0.5cM. fix end pos, but otherwise leave
      exc = omap[ix,]; exc = rbind(exc[1:162,], exc[166,], exc[163:165,])
      exc$o = exc$o[order(exc$o)]; exc$diff = c(abs(diff(c(exc$start[1],exc$end[1]))), abs(diff(exc$start)))
      exc$diff[164] = diff(c(exc$start[164], exc$end[164]))
      omap[ix,] = exc
      # Ctg019 large block of dupe markers OOO, again < 0.5cM
      
      # update cumsum
      omap = do.call(rbind, lapply(split(omap, omap$lg), function(x) {x$cpos[!is.na(x$diff)] = cumsum(x$diff[!is.na(x$diff)]); x}))
      omap
    }
    
    fixBadDiffs_NIC58_flye_stitch2 <- function(omap){
      
      # flye manual fixes (pacbio ass + 0cM gap spanning reads)
      # check full (marker+dupe) orientations
      cori = sapply(split(omap, omap$scaf), function(x) getOri(x$start, gte=T))
      table(cori, useNA = 'al')
      sus = names(which(is.na(cori))); print(sus)
      
      # stitched flye manual fixes
      # 3 cases, all < 0.5cM
      # same as before, presumably
      # i=1; sus[i]; subi=subset(omap, scaf==sus[i]); rownames(subi)=NULL; plot(subi$cpos[-1], diff(subi$start)); ix = which(omap$scaf==sus[i])
      # Ctg003: 1 marker (Ctg003.4051, 3360717 3425185) placed after 2 dupes (3426298 - 3526455), 0.42cM 
      
      # Ctg012: terminal dupes switched, 0.42cm. fix end pos, but otherwise leave
      # dupes up to 540959 come after final marker up to 3292160
      i=2; sus[i]; ix = which(omap$scaf==sus[i]); subi=omap[ix,]; rownames(subi)=NULL
      # plot(subi$cpos[-1], diff(subi$start))
      exc = subi; exc = rbind(exc[1:161,], exc[163,], exc[162,], exc[164:166,])
      exc$o = exc$o[order(exc$o)]; exc$diff = c(abs(diff(c(exc$start[1],exc$end[1]))), abs(diff(exc$end)))
      exc$diff[162:166] = exc$cpos[162:166] = NA
      # plot(exc$diff)
      # plot(cumsum(exc$diff)[-1], diff(exc$start))
      # plot(exc$cpos, exc$genetic)
      # set dominant ori and add scaforder
      omap$strand = as.character(omap$strand)
      exc$strand = as.character(factor(getOri(exc$start, dominant = T), levels=c(T,F), labels=c('+', '-')))
      exc$ori = 1
      omap[ix,] = exc
      
      # Ctg015 large block of dupe markers OOO, 0.42cM
      # force into order
      i=3; ix = which(omap$scaf==sus[i]); sus[i]; subi=omap[ix,]; rownames(subi)=NULL
      # plot(subi$cpos[-1], diff(subi$start))
      i=which(omap$marker=='Ctg015.201')
      omap$diff[i] = omap$end[i]-omap$start[i]
      omap$diff[i+1] = omap$end[i+1]-omap$start[i+1]
      # update cumsum
      rownames(omap)=NULL
      omap = do.call(rbind, lapply(split(omap, omap$lg), function(x) {
        lgg = aggregate(data = x, genetic ~ scaf, median)
        lgg = lgg[order(lgg$genetic),]
        lgg$scafo = 1:nrow(lgg)
        x = merge(x[,-which(names(x)=='scafo')], lgg[,c('scaf', 'scafo')], sort=F)
        x = x[order(x$scafo, x$o),]
        x$cpos[is.na(x$diff)] = NA
        # x$diff[!is.na(x$diff)] = c(x$end[1], abs(diff(x$start[!is.na(x$diff)])))
        x$cpos[!is.na(x$diff)] = cumsum(x$diff[!is.na(x$diff)])
        x
      }))
      omap
    }
    
    dumpJunctions_NIC58_flye <- function(omap){
      # dump unoriented contigs for 
      # 1. checking against other assemblies
      # 2. local reassembly if required
      # all are on tips flanked by larger sequences (2 for Ctg008 and Ctg011, 1 for Ctg023)
      buffcM = 1
      unori = merge(merge(data.frame(scaf=sort(unique(subset(omap, ori!=1)$scaf))), fai), unique(omap[,c('lg', 'scaf', 'genetic', 'marker_scaf')]))
      names(unori)[c(1,2,3,6)] = c('uno', 'uno_len', 'uno_in', 'scaf')
      unori = merge(unori, fai)
      names(unori)[c(1,7,8)] = c('ori', 'ori_len', 'ori_in')
      unori = unori[,c(1,2,8,4,3,7,5,6)]
      unori = unori[order(unori$lg, unori$genetic),]
      lapply(unori$uno, function(i) {
        con = unique(subset(omap, scaf==i)[,c('lg', 'genetic', 'marker_scaf')])
        unique(subset(omap, lg==con$lg[1] & abs(genetic-con$genetic[1])<buffcM)[,c('genetic', 'marker_scaf')])
      })
      fwrite(unori, file = sprintf('%s/%s_unorderedContigs.csv', WD, pref))
      
      # and all oriented contigs
      umap = merge(unique(omap[,c(3,8,11)]), fai)
      umap <- merge(umap, aggregate(data = omap, genetic~scaf, min)); names(umap)[names(umap)=='genetic'] = 'genetic_start'
      umap <- merge(umap, aggregate(data = omap, genetic~scaf, max)); names(umap)[names(umap)=='genetic'] = 'genetic_end'
      umap <- umap[order(umap$lg, umap$genetic_start, umap$genetic_end),]
      umap$gap_cM = unlist(lapply(split(umap, umap$lg), function(x) c(x$genetic_start[-1]-x$genetic_end[-nrow(x)], NA)))
      fwrite(umap, file = sprintf('%s/%s_orderedContigs.csv', WD, pref))
    }
    
    # add duplicates, sort ascending by LG, physical marker position in contig
    cat(sprintf('%s Add duplicate markers\n', pref))
    d = mapd = addDupes(cross, dupeList, mpos)
    
    # order duplicates first within contigs, then LGs
    # add terminal contig positions and cumulative physical positions
    # along the linkage groups
    # (this can be > assembly span if there are map errors)
    cat(sprintf('%s Order duplicate markers\n', pref))
    omap <- orderDupes(mapd)
    
    ## NIC58 round 1: 5/6 cases irreconcilable phys/genetic discordance, but all < 0.5cM
    # omap <- fixBadDiffs_NIC58_flye(omap)
    ## NIC58 round 3: 3 cases irreconcilable phys/genetic discordance, but all < 0.5cM
    # omap <- fixBadDiffs_NIC58_flye_stitch1(omap)
    # omap <- fixBadDiffs_NIC58_flye_stitch2(omap)
    
    osum = aggregate(data = omap, diff~ori, sum)
    if(!0 %in% osum$ori) osum = rbind(osum, data.frame(ori=0, diff=0))
    if(sum(osum$ori %% 1)==0) osum = rbind(osum, data.frame(ori=0.5, diff=0))
    orit = table(factor(c(0,1,0.5))); orit[]=0
    ot = table(unique(omap[,c('scaf', 'ori')])$ori)
    if('0' %in% names(ot)) orit[1] = ot[names(ot)=='0']
    if(length(ot[as.numeric(names(ot)) %% 1])>0) orit[2] = ot[as.numeric(names(ot)) %% 1 > 0]
    if('1' %in% names(ot)) orit[3] = ot[names(ot)=='1']
    
    cat(sprintf('%s After duplicate and termini addition:\n%.2f Mb ordered (%s contig)\n%.2f Mb genetically unorderable (%s contig)\n%.2f Mb ambiguous (%s contig)\nvs. assembly span %.2f Mb (%.2f off)\n%s contigs', 
                pref,
                osum$diff[osum$ori==1]/1e6,
                orit[3],
                osum$diff[osum$ori==0]/1e6, 
                orit[1],
                osum$diff[osum$ori%%1>0]/1e6, 
                orit[2],
                sum(fai$len[fai$scaf %in% omap$scaf])/1e6,
                (sum(osum$diff)-sum(fai$len[fai$scaf %in% omap$scaf]))/1e6, 
                length(unique(omap$scaf))
    ))
    
    # write out any unoriented junctions and estimated genetic gaps
    # dump a stranded bed file for converting to pseudochromosomes
    dumpJunctions_NIC58_flye(omap)
    
    omap
  }
  
  mm <- expandMap(cross, dupeList)
  
  save(cross, dupes, duped, gterr, fai, mpos, mapdf, mm, file = sprintf('%s/%s_rimap_3.rda', WD, pref))
}

startingMap(WD, pref)

load(sprintf('%s/%s_rimap_1.rda', WD, pref))
for(q in c(1, 0.99, 0.95, 0.9)){
  opref = sprintf('%s_q%s', pref, q*100)
  load(sprintf('%s/%s_bestOrders_LL.rda', WD, opref))
  geneticPhysicalConcordance(rcross, dupes, fai, mpos, opref)
}

gatherStats <- function(WD, pref, q=1){
  
  # compile some useful stats across assemblies
  load(sprintf('%s/%s/rqtl/%s_q%s_rimap_3.rda', WD, pref, pref, q*100), verbose = T)
  
  PD=sprintf('%s/%s', WD, pref)
  ofai = fread(sprintf('%s/ref/genome.fa.fai', PD))[,1:2]
  names(ofai) = c('fa', 'len')
  
  refBamStats = fread(cmd = sprintf('grep "^SN" %s/bam/parents/NIC58.stats | cut -f2-3', PD))
  altBamStats = fread(cmd = sprintf('grep "^SN" %s/bam/parents/JU1373.stats | cut -f2-3', PD))
  
  # n seq, n diallelic SNPs, thinned SNPs, n variable seq
  dSNPS = fread(cmd=sprintf('cut -f1 -d" " %s/variants.log', PD), header = F)$V1
  
  osum = aggregate(data = mm, diff~ori, sum)
  if(!0 %in% osum$ori) osum = rbind(osum, data.frame(ori=0, diff=0))
  if(sum(osum$ori %% 1)==0) osum = rbind(osum, data.frame(ori=0.5, diff=0))
  
  mapBubbles = sapply(split(mm, mm$lg), function(x) {
    x = subset(x, dupe==F)
    rl = rle(x$scaf)
    rlt = table(rl$values)
    sum(rlt>1)
  })
  
  splitCtg = table(mm$scaf, mm$lg)
  splitCtg = splitCtg[apply(splitCtg, 1, function(x) sum(x>0)>1),]
  
  csum = summaryMap(cross)
  
  list(map = sprintf('%s/%s', WD, pref),
       xoQuantile = q,
       n_markers = csum$n.mar[7],
       map_span_cM = csum$length[7],
       map_span_Mb = sum(fai$len[fai$scaf %in% mm$scaf])/1e6,
       ordered_span_Mb = osum$diff[osum$ori==1]/1e6,
       unordered_span_Mb = osum$diff[osum$ori==0]/1e6,
       ambiguous_span_Mb = osum$diff[osum$ori%%1>0]/1e6,
       bubbles = sum(mapBubbles),
       bubbled_LG = sum(mapBubbles>0),
       LG_splits = nrow(splitCtg),
       n_map_contigs = length(unique(mm$scaf)),
       map_fai = list(fai),
       original_fai = list(ofai),
       nSNP = dSNPS[2],
       NIC58_Illumina_bam_stats = list(refBamStats),
       JU1373_Illumina_bam_stats = list(altBamStats), 
       rqtl_cross = cross, 
       map_df = mm
  )
}

home=F
if(home){
  WD='/Volumes/scratch/ctrop/gmap/NIC58_ref/cfass/'
  ass = sort(c('ra', 'flye', 'canu', 'wtdbg2'))
  mapStats = lapply(ass, function(i) gatherStats(WD, i))
  save(mapStats, file = '~/Documents/ctrop/papers/1_geneticMap/tropicalis_geneticMaps.rda')
  
  load('~/Documents/ctrop/papers/1_geneticMap/tropicalis_geneticMaps.rda')
  
  # preen flye map
  cross = mapStats[[2]]$rqtl_cross
  mm = mapStats[[2]]$map_df
  fai = mapStats[[2]]$map_fai[[1]]
  
  # no apparent bubbles
  summary(cross)
  summaryMap(cross)
  cross <- est.rf(cross)
  plotMap(cross)
  plotRF(cross)
  
 

}
