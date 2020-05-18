######################################################################
#
# read.cross.tidy.R
#
# copyright (c) 2005-2015, Karl W Broman
# last modified Aug, 2015
# first written Aug, 2014
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
#
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
#
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Part of the R/qtl package
# Contains: read.cross.tidy
#           [See read.cross.R for the main read.cross function.]
#
######################################################################

######################################################################
#
# read.cross.tidy
#
# read data in comma-delimited format, with separate files for phenotype,
# genotype, and map data
#
######################################################################

# dir=WD; genfile='geno.riself.bin50.csv'; mapfile = 'pmap.riself.bin50.csv'; phefile = 'pheno.csv'; na.strings = 'NA';genotypes = c('A', 'B')

my.read.cross <- function (format = c("csv", "csvr", "csvs", "csvsr", "mm", "qtx", 
                                      "qtlcart", "gary", "karl", "mapqtl", "tidy"), dir = "", file, 
                           genfile, mapfile, phefile, chridfile, mnamesfile, pnamesfile, 
                           na.strings = c("-", "NA"), genotypes = c("A", "H", "B", "D", 
                                                                    "C"), alleles = c("A", "B"), estimate.map = TRUE, convertXdata = TRUE, 
                           error.prob = 1e-04, map.function = c("haldane", "kosambi", 
                                                                "c-f", "morgan"), BC.gen = 0, F.gen = 0, crosstype, ...) 
{
  if (format == "csvrs") {
    format <- "csvsr"
    warning("Assuming you mean 'csvsr' rather than 'csvrs'.\n")
  }
  format <- match.arg(format)
  if (format == "csv" || format == "csvr") {
    cross <- read.cross.csv(dir, file, na.strings, genotypes, 
                            estimate.map, rotate = (format == "csvr"), ...)
  }
  else if (format == "csvs" || format == "csvsr") {
    if (missing(phefile) && !missing(file) && !missing(genfile)) {
      phefile <- genfile
      genfile <- file
    }
    else if (missing(genfile) && !missing(file) && !missing(phefile)) {
      genfile <- file
    }
    cross <- read.cross.csvs(dir, genfile, phefile, na.strings, 
                             genotypes, estimate.map, rotate = (format == "csvsr"), 
                             ...)
  }
  else if (format == "qtx") {
    cross <- read.cross.qtx(dir, file, estimate.map)
  }
  else if (format == "qtlcart") {
    if (missing(mapfile) && !missing(genfile)) 
      mapfile <- genfile
    cross <- read.cross.qtlcart(dir, file, mapfile)
  }
  else if (format == "karl") {
    if (missing(genfile)) 
      genfile <- "gen.txt"
    if (missing(mapfile)) 
      mapfile <- "map.txt"
    if (missing(phefile)) 
      phefile <- "phe.txt"
    cross <- read.cross.karl(dir, genfile, mapfile, phefile)
  }
  else if (format == "mm") {
    if (missing(mapfile) && !missing(genfile)) 
      mapfile <- genfile
    cross <- read.cross.mm(dir, file, mapfile, estimate.map)
  }
  else if (format == "gary") {
    if (missing(genfile)) 
      genfile <- "geno.dat"
    if (missing(mnamesfile)) 
      mnamesfile <- "mnames.txt"
    if (missing(chridfile)) 
      chridfile <- "chrid.dat"
    if (missing(phefile)) 
      phefile <- "pheno.dat"
    if (missing(pnamesfile)) 
      pnamesfile <- "pnames.txt"
    if (missing(mapfile)) 
      mapfile <- "markerpos.txt"
    cross <- read.cross.gary(dir, genfile, mnamesfile, chridfile, 
                             phefile, pnamesfile, mapfile, estimate.map, na.strings)
  }
  else if (format == "mapqtl") {
    cross <- read.cross.mq(dir = dir, locfile = genfile, 
                           mapfile = mapfile, quafile = phefile)
  }
  else if (format == "tidy") {
    if (!missing(file) && !missing(genfile) && !missing(mapfile) && 
        missing(phefile)) {
      phefile <- mapfile
      mapfile <- genfile
      genfile <- file
    }
    if (missing(genfile)) 
      genfile <- "gen.csv"
    if (missing(phefile)) 
      phefile <- "phe.csv"
    if (missing(mapfile)) 
      mapfile <- "map.csv"
    cross <- read.tidy(dir = dir, genfile = genfile, 
                             phefile = phefile, mapfile = mapfile, na.strings = na.strings, 
                             genotypes = genotypes)
  }
  estimate.map <- cross[[2]]
  cross <- cross[[1]]
  chrnam <- names(cross$geno)
  if (all(regexpr("^[Cc][Hh][Rr]", chrnam) > 0)) {
    chrnam <- substr(chrnam, 4, nchar(chrnam))
    if (all(regexpr("^[Oo][Mm][Oo][Ss][Oo][Mm][Ee]", chrnam) > 
            0)) 
      chrnam <- substr(chrnam, 8, nchar(chrnam))
  }
  if (sum(chrnam == "x") > 0) 
    chrnam[chrnam == "x"] <- "X"
  names(cross$geno) <- chrnam
  for (i in 1:length(cross$geno)) if (names(cross$geno)[i] == 
                                      "X") 
    class(cross$geno[[i]]) <- "X"
  chrtype <- sapply(cross$geno, class)
  if (any(chrtype == "X") && convertXdata) {
    if (class(cross)[1] == "bc") 
      cross <- fixXgeno.bc(cross)
    if (class(cross)[1] == "f2") {
      if (missing(alleles)) 
        alleles <- c("A", "B")
      cross <- fixXgeno.f2(cross, alleles)
    }
  }
  cross <- read.cross.bcsft(cross = cross, BC.gen = BC.gen, F.gen = F.gen, ...)
  if (estimate.map) {
    cat(" --Estimating genetic map\n")
    map.function <- match.arg(map.function)
    newmap <- est.map(cross, error.prob = error.prob, map.function = map.function)
    cross <- replace.map(cross, newmap)
  }
  for (i in 1:nchr(cross)) storage.mode(cross$geno[[i]]$data) <- "integer"
  if (class(cross)[1] != "4way") {
    if (length(alleles) > 2) {
      warning("length of arg alleles should be 2")
      alleles <- alleles[1:2]
    }
    if (length(alleles) < 2) 
      stop("length of arg alleles should be 2")
  }
  else {
    if (missing(alleles)) 
      alleles <- c("A", "B", "C", "D")
    if (length(alleles) > 4) {
      warning("length of arg alleles should be 4 for a 4-way cross")
      alleles <- alleles[1:4]
    }
    if (length(alleles) < 4) 
      stop("length of arg alleles should be 4 for a 4-way cross")
  }
  if (any(nchar(alleles)) != 1) {
    warning("Each item in arg alleles should be a single character")
    alleles <- substr(alleles, 1, 1)
  }
  attr(cross, "alleles") <- alleles
  type <- class(cross)[1]
  if (!missing(crosstype)) {
    if (crosstype == "risib") 
      cross <- convert2risib(cross)
    else if (crosstype == "riself") 
      cross <- convert2riself(cross)
    else class(cross)[1] <- crosstype
  }
  summary(cross)
  cat(" --Cross type:", class(cross)[1], "\n")
  cross
}

read.tidy <- function(dir, genfile, phefile, mapfile, na.strings=c("-","NA"),
                      genotypes=c("A","H","B","D","C"), ...)
{
  # create file names
  if(missing(genfile)) genfile <- "gen.csv"
  if(missing(phefile)) phefile <- "phe.csv"
  if(missing(mapfile)) mapfile <- "map.csv"
  
  if(!missing(dir) && dir != "") {
    genfile <- file.path(dir, genfile)
    phefile <- file.path(dir, phefile)
    mapfile <- file.path(dir, mapfile)
  }
  
  args <- list(...)
  
  if("" %in% na.strings) {
    na.strings <- na.strings[na.strings != ""]
    warning("Including \"\" in na.strings will cause problems; omitted.")
  }
  
  # if user wants to use comma for decimal point, we need
  dec <- ifelse("dec" %in% names(args), args[["dec"]], ".")
  
  # "sep" not in the "..." argument and so take sep=","
  sep <- ifelse("sep" %in% names(args), args[["sep"]], ",")
  
  # read the data files
  args <- c(args, list(sep = sep, na.strings = na.strings, row.names = 1, header = TRUE, stringsAsFactors = FALSE))
  
  gen   <- do.call("read.table", c(args, list(file = genfile)))
  pheno <- do.call("read.table", c(args, list(file = phefile)))
  map   <- do.call("read.table", c(args, list(file = mapfile)))
  
  # Check individual IDs
  mp <- setdiff(colnames(gen), colnames(pheno))
  if (length(mp) > 0) {
    warning(length(mp), " individuals with genotypes but no phenotypes\n    ", paste(mp, collapse="|"), "\n")
    pheno[mp] <- NA
  }
  
  mg <- setdiff(colnames(pheno), colnames(gen))
  if (length(mg) > 0) {
    warning(length(mg), " individuals with phenotypes but no genotypes\n    ", paste(mg, collapse="|"), "\n")
    gen[mg] <- NA
  }
  
  # ensure individual order is consistent
  ids   <- colnames(pheno)
  pheno <- pheno[, ids]
  gen   <- gen[,   ids]
  
  # Check markers
  genm   <- rownames(gen)
  mapm   <- rownames(map)
  mnames   <- intersect(genm, mapm)
  
  mg <- setdiff(mapm, genm)
  if (length(mg) > 0) warning("Removing ", length(mg)," genotyped markers with missing positions\n")
  
  mm <- setdiff(genm, mapm)
  if (length(mm) > 0) warning("Removing", length(mm), " mapped markers with no genotypes\n")
  
  map[[2]] <- asnumericwithdec(map[[2]], dec)
  
  # replace empty cells with NA
  add_na <- function(x) replace(x, !is.na(x) & x == "", NA)
  pheno <- add_na(pheno)
  gen   <- add_na(gen)
  
  
  # check chromosomes
  chr <- map[[1]]
  
  if(any(is.na(chr))) {stop("There are missing chromosome IDs. Check row(s) ", paste(which(is.na(chr)), collapse=","), sep = "")}
  
  # look for strange entries in the genotype data
  if(length(genotypes) > 0) {
    temp <- Filter(Negate(is.na), unique(unlist(gen)))
    wh <- !(temp %in% genotypes)
    if(any(wh)) {
      warn <- "The following unexpected genotype codes were treated as missing.\n    "
      ge <- paste("|", paste(temp[wh],collapse="|"),"|",sep="")
      warn <- paste(warn,ge,"\n",sep="")
      warning(warn)
    }
    
    # convert genotype data
    gen <- as.matrix(gen)
    allgeno <- matrix(match(gen, genotypes), ncol = ncol(gen), dimnames = dimnames(gen))
  } else {
    genotypes <- Filter(Negate(is.na), unique(unlist(gen)))
    gen <- as.data.frame(lapply(gen, factor, levels = genotypes), rownames(gen))
    allgeno <- data.matrix(gen)
  }
  
  # convert phenotype data
  # pheno must be rotated to allow for numeric and factor variables
  pheno <- data.frame(t(pheno))
  pheno <- data.frame(lapply(pheno, sw2numeric, dec = dec), row.names = NULL, stringsAsFactors = TRUE)
  
  # add id column if informative identifiers are provided
  default.ids <- make.names(seq_len(nrow(pheno)))
  if (!all(ids %in% default.ids)) pheno$id <- ids
  
  # re-order the markers by chr and position
  # try to figure out the chr labels
  if(all(chr %in% c(1:999,"X","x"))) { # 1...19 + X
    tempchr <- chr
    tempchr[chr=="X" | chr=="x"] <- 1000
    tempchr <- as.numeric(tempchr)
    neworder <- order(tempchr, map[[2]])
  } else {
    # prevent reordering of chromosomes
    tempchr <- factor(chr, levels=unique(chr))
    neworder <- order(tempchr, map[[2]])
  }
  chr <- chr[neworder]
  map <- map[neworder, ]
  allgeno <- allgeno[rownames(map), , drop = FALSE]
  mnames <- mnames[neworder]
  
  # fix up map information
  # number of chromosomes
  uchr <- unique(chr)
  n.chr <- length(uchr)
  geno <- vector("list", n.chr)
  names(geno) <- uchr
  min.mar <- 1
  allautogeno <- NULL
  
  for(i in 1:n.chr) { # loop over chromosomes
    print(i)
    # create map
    temp.map <- map[[2]][chr==uchr[i]]
    names(temp.map) <- mnames[chr==uchr[i]]
    
    # pull out appropriate portion of genotype data
    data <- t(allgeno[names(temp.map), ])
    if(dim(data)[1]==1) {data=t(data); colnames(data)=names(temp.map)}
    
    # drop genotype rownames
    rownames(data) <- NULL
    
    geno[[i]] <- list(data=data, map=temp.map)
    if(uchr[i] == "X" || uchr[i] == "x")
      class(geno[[i]]) <- "X"
    else {
      class(geno[[i]]) <- "A"
      if(is.null(allautogeno)) allautogeno <- data
      else allautogeno <- cbind(allautogeno,data)
    }
  }
  
  if(is.null(allautogeno)) allautogeno <- allgeno
  
  # check that data dimensions match
  n.mar1 <- sapply(geno,function(a) ncol(a$data))
  n.mar2 <- sapply(geno,function(a) length(a$map))
  n.phe <- ncol(pheno)
  n.ind1 <- nrow(pheno)
  n.ind2 <- sapply(geno,function(a) nrow(a$data))
  if(any(n.ind1 != n.ind2)) {
    cat(n.ind1,n.ind2,"\n")
    stop("Number of individuals in genotypes and phenotypes do not match.");
  }
  if(any(n.mar1 != n.mar2)) {
    cat(n.mar1,n.mar2,"\n")
    stop("Numbers of markers in genotypes and marker names files do not match.");
  }
  
  # print some information about the amount of data read
  cat(" --Read the following data:\n");
  cat("\t", n.ind1, " individuals\n");
  cat("\t", sum(n.mar1), " markers\n");
  cat("\t", n.phe, " phenotypes\n");
  
  if(all(is.na(allgeno)))
    warning("There is no genotype data!\n")
  
  # determine map type: f2 or bc or 4way?
  if(all(is.na(allgeno))) warning("There is no genotype data!\n")
  if(all(is.na(allautogeno)) || max(allautogeno,na.rm=TRUE)<=2) type <- "bc"
  else if(max(allautogeno,na.rm=TRUE)<=5) type <- "f2"
  else type <- "4way"
  cross <- list(geno=geno,pheno=pheno)
  class(cross) <- c(type,"cross")
  
  # check that nothing is strange in the genotype data
  cross.type <- class(cross)[1]
  if(cross.type=="f2") max.gen <- 5
  else if(cross.type=="bc") max.gen <- 2
  else max.gen <- 14
  
  # check that markers are in proper order
  #     if not, fix up the order
  for(i in 1:n.chr) {
    if(any(diff(cross$geno[[i]]$map)<0)) {
      o <- order(cross$geno[[i]]$map)
      cross$geno[[i]]$map <- cross$geno[[i]]$map[o]
      cross$geno[[i]]$data <- cross$geno[[i]]$data[,o,drop=FALSE]
    }
  }
  
  # if 4-way cross, make the maps matrices
  if(type=="4way") {
    for(i in 1:n.chr)
      cross$geno[[i]]$map <- rbind(cross$geno[[i]]$map, cross$geno[[i]]$map)
  }
  
  # return cross + indicator of whether to run est.map
  list(cross,FALSE)
}

# end of read.cross.tidy.R

sw2numeric <- function(x, dec) {
  x[x == ""] <- NA
  wh1 <- is.na(x)
  n <- sum(!is.na(x))
  y <- suppressWarnings(asnumericwithdec(as.character(x), dec))
  wh2 <- is.na(y)
  m <- sum(!is.na(y))
  if(n==m || (n-m) < 2 || (n-m) < n*0.05) {
    if(sum(!wh1 & wh2) > 0) {
      u <- unique(as.character(x[!wh1 & wh2]))
      if(length(u) > 1) {
        themessage <- paste("The phenotype values", paste("\"", u, "\"", sep="", collapse=" "))
        themessage <- paste(themessage, " were", sep="")
      }
      else {
        themessage <- paste("The phenotype value \"", u, "\" ", sep="")
        themessage <- paste(themessage, " was", sep="")
      }
      themessage <- paste(themessage, "interpreted as missing.")
      warning(themessage)
      
    }
    return(y)
  }
  else return(x)
}

asnumericwithdec <- function(x, dec=".") {
  if(dec!=".") x <- gsub(paste0("\\", dec), ".", x)
  as.numeric(x)
}

read.cross.bcsft <- function(..., BC.gen = 0, F.gen = 0, cross = NULL, force.bcsft = FALSE, estimate.map=TRUE){
  ## Must specify s = BC.gen and t = F.gen.
  ## Later: Could import in clever way from qtlcart? See qtlcart_io.R and their software.
  
  ## Make sure we only estimate map once!
  if(is.null(cross)) # Estimate map at end of this routine (called read.cross.bcsft directly).
    cross <- read.cross(..., estimate.map = FALSE)
  else # Estimate map in parent read.cross() call (read.cross.bcsft is pass-through from read.cross).
    estimate.map <- FALSE
  
  force.bcsft <- force.bcsft | (BC.gen > 0 | F.gen > 0)
  if((class(cross)[1] %in% c("bc","f2")) & force.bcsft) {
    
    # deal with ... args
    dots <- list(...)
    
    if("verbose" %in% names(dots)) verbose <- dots$verbose
    else verbose <- TRUE
    
    if("error.prob" %in% names(dots)) error.prob <- dots$error.prob
    else error.prob <- 0.0001
    
    if("map.function" %in% names(dots)) map.function <- dots$map.function
    else map.function <- "haldane"
    
    cross <- convert2bcsft(cross, BC.gen, F.gen, estimate.map = estimate.map,
                           error.prob=error.prob, map.function=map.function,
                           verbose=verbose)
  }
  
  cross
}
