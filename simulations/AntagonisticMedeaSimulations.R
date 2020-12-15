
# Simulations with perfectly linked antagonistic drivers 

# This file starts with a function to perform simulations, as described in the Methods section of the paper.
# The function starts by forming the starting population, with genotypes drawn from the equilibrium probability distribution but sexes assigned deterministically according to the equilibrium frequencies.
# Then it loops through ngen generations and returns data.frame(GenotypeFreqs, Sexes, PopSize, S), where each row is a generation.
# Selection operates in a Wright-Fisher framework if population regulation Reg = "fixed" and in an exponential growth regime if Reg = "exponential". 
# V2 is the viability of homozygotes for the allele with initial frequency p with heterozygous mothers, and V1 is the viability of the analogous homozygotes for the allele with initial frequency 1-p. 
# In other words, if p = 0.2, the resident Medea, with frequency 0.8, kills with penetrance 1-V1, and the invading Medea kills with penetrance 1-V2. 

###########################################################

patchdrive <- function(p = 0.05, S = 0.95, N = 100, V1= 0.05, V2 = 0.05, B = 50, Him = 0.005, ngen = 3,Reg = "exponential"){

	countGenos <- function(x,y){
		counts <- apply(rbind(x,y), 1, sum)
		cbind(sum(counts ==0), sum(counts == 1), sum(counts ==2))
	}

	#rounding function via JLutils
	round_preserve_sum <- function(x, digits = 0) {
	  up <- 10^digits
	  x <- x * up
	  y <- floor(x)
	  indices <- tail(order(x - y), round(sum(x)) - sum(y))
	  y[indices] <- y[indices] + 1
	  y / up
	}

	# Starting genotype frequencies from extended HWE
	Fhat <- S/(2-S) # equilibrium inbreeding coeffecient given S
	
	JJ <- 	((1-Fhat)*p^2 + p*Fhat)
	JN <- 	(1-Fhat)*2*p*(1-p)
	NN  <-	((1-Fhat)*(1-p)^2 + (1-p)*Fhat)
	
	hermprob <- (1+Fhat)/2
	maleprob <- (1-Fhat)/2
	sexcount <- round_preserve_sum(c(hermprob,maleprob)*N, 0)
	
	StartingGenotypes <- c(rmultinom(1, sexcount[1], c(JJ, JN,NN)), rmultinom(1, sexcount[2], c(JJ,JN,NN)))
	
	herms <- matrix(nrow = sum(StartingGenotypes[c(1:3)]), ncol = 2)
	herms[,1] <- c(rep(0,StartingGenotypes[1]), rep(0, StartingGenotypes[2]), rep(1, StartingGenotypes[3]))
	herms[,2] <- c(rep(0,StartingGenotypes[1]), rep(1, StartingGenotypes[2]), rep(1, StartingGenotypes[3]))
	
	males <- matrix(nrow = sum(StartingGenotypes[c(4:6)]), ncol = 2)
	males[,1] <- c(rep(0,StartingGenotypes[4]), rep(0, StartingGenotypes[5]), rep(1, StartingGenotypes[6]))
	males[,2] <- c(rep(0,StartingGenotypes[4]), rep(1, StartingGenotypes[5]), rep(1, StartingGenotypes[6]))
	
	PopSize <- N
	
	GenotypeFreqs <- matrix(nrow = ngen, ncol = 3)
	GenotypeFreqs[1,] <- countGenos(herms,males)
	
	Sexes <- matrix(nrow = ngen, ncol = 2)
	Sexes[1,] <- c(dim(herms)[1], dim(males)[1])	
	
	
	for( gen in 2: ngen){
		
	offspring <- NULL; 
	viability <- NULL
	sex <- NULL
	# In the case where males exist:
	if(!is.na(males[1])){
	for(i in 1:length(herms[,1])){
		kid <- matrix(ncol = 2, nrow = B)
		kidviability <- NULL
		kidsex <- NULL
		for(j in 1:(round(S*B))){
			kid[j,] <- sample(herms[i,], 2, replace = T)
				if(sum(kid[j,]) == 2 & sum(herms[i,]) == 1){
					kidviability[j] <- rbinom(1,1,V2)
				}else if(sum(kid[j,]) == 0 & sum(herms[i,]) == 1){
					kidviability[j] <- rbinom(1,1,V1)
				}else{
					kidviability[j] <- 1
				}
				kidsex[j] <- rbinom(1,1,Him)
			}
		dad <- sample(seq(1:length(males[,1])), 1)	
		if(round(S*B) < B){
		for(k in ((round(S*B))+1): B){
			kid[k,] <- c(sample(herms[i,], 1), sample(males[dad,], 1))		
				if(sum(kid[k,]) == 2 & sum(herms[i,]) == 1){
				kidviability[k] <- rbinom(1,1,V2)
				}else if(sum(kid[k,]) == 0 & sum(herms[i,]) == 1){
					kidviability[k] <- rbinom(1,1,V1)
				}else{
				kidviability[k] <- 1
				}
				kidsex[k] <- rbinom(1,1, 0.5)
			}
		}	
		offspring <- rbind(offspring, kid)
		viability <- c(viability, kidviability)
		sex <- c(sex,kidsex)
		}
	}
	
	
	# In the case where there are no males:
	if(is.na(males[1])){
	for(i in 1:length(herms[,1])){
		kid <- matrix(ncol = 2, nrow = B)
		kidviability <- NULL
		kidsex <- NULL
		for(j in 1:B){
			kid[j,] <- sample(herms[i,], 2, replace = T)
				if(sum(kid[j,]) == 2 & sum(herms[i,]) == 1){
					kidviability[j] <- rbinom(1,1,V2)
				}else if(sum(kid[j,]) == 0 & sum(herms[i,]) == 1){
					kidviability[j] <- rbinom(1,1,V1)
				}else{
					kidviability[j] <- 1
				}
				kidsex[j] <- rbinom(1,1,Him)
			}
	
		offspring <- rbind(offspring, kid)
		viability <- c(viability, kidviability)
		sex <- c(sex, kidsex)
		}
		}
	
	# Select offspring for the next generation
	if(Reg == "fixed"){	
	survivors <- offspring[viability == 1,]	
	survivorsex <- sex[viability == 1]	
	nextgen <- sample(seq(1:length(survivors[,1])), N)
	nextgengenos <- survivors[nextgen,]
	nextgensex <- survivorsex[nextgen]
	herms <- nextgengenos[nextgensex ==0,]
	if(length(nextgensex == 1) >0){
	males <- matrix(nextgengenos[nextgensex ==1,], ncol =2)
	}else{
	males <- matrix(nrow=0, ncol = 2)
		}
	GenotypeFreqs[gen,] <- countGenos(herms,males)
	Sexes[gen,] <- c(dim(herms)[1], dim(males)[1])	
		}
			
	else if(Reg == "exponential"){
	survivors <- offspring[viability == 1,]	
	survivorsex <- sex[viability == 1]
	herms <- survivors[survivorsex == 0,]
	if(length(survivorsex == 1) >0){
	males <- matrix(survivors[survivorsex ==1,], ncol =2)
	}else{
	males <- matrix(nrow=0, ncol = 2)
		}
	GenotypeFreqs[gen,] <- countGenos(herms,males)
	Sexes[gen,] <- c(dim(herms)[1], dim(males)[1])	
	PopSize[gen] <- dim(survivors)[1]
	}	

	} # Loop through generations. 

	data.frame(GenotypeFreqs, Sexes, PopSize, S)
}


###########################################################

# Now we can call this function to simulate data

# Lots of patches, exponential growth
# Figure 12C top
n.reps <- 250
n.patches <- 250
Selfing <- c(0,0.1, 0.25,0.5,0.75,0.9,1)
for(m in 1:7){
	GlobalFreqs <- matrix(nrow = (n.reps), ncol = 8)
	for(x in 1:n.reps){
	PatchFreqs <- matrix(nrow = n.patches, ncol = 9)
		for(z in 1:n.patches){
			drive.out <- as.matrix(patchdrive(N = 4, S = Selfing[m], B = 50, ngen = 3, Reg = "exponential", p = 0.2))
			lastgen <- dim(drive.out)[1]
			PatchFreqs[z,] <-c(drive.out[lastgen,], drive.out[1,c(4,5)])
			}	
		GlobalFreqs[x,] <- c(apply(PatchFreqs, 2, sum)[1:6], PatchFreqs[1,c(8,9)])			
	}
	write.csv(GlobalFreqs, paste0("S=",Selfing[m], ",N=4,250patches,B=50,ngen=3,Reg=exp,p=0.2,AntagonisticDrivers.csv"))
}
pS.exp <- matrix(nrow = 250, ncol = 7)
for(m in 1:7){infile <- read.csv(file = paste0("S=",Selfing[m], ",N=4,250patches,B=50,ngen=3,Reg=exp,p=0.2,AntagonisticDrivers.csv"), head = T);
	sim.freq <- (2*infile[,2] + infile[,3] )/(2*infile[,7])
	pS.exp[,m] <- sim.freq
	}

# Fixed population size, no patches
# Figure 12C bottom
Selfing <- c(0,0.1, 0.25,0.5,0.75,0.9,1)
for(m in 1:7){
	n.reps <- 250
	FWFreqs <- matrix(nrow = (n.reps), ncol = 6)
	for(x in 1:n.reps){
			drive.out <- as.matrix(patchdrive(N = 1000, S = Selfing[m], B = 50, ngen = 3, Reg = "fixed", p = 0.8))
			lastgen <- dim(drive.out)[1]		
			FWFreqs[x,] <- drive.out[lastgen,1:6]			
			}	
	write.csv(FWFreqs, paste0("S=",Selfing[m],"N=1000,B=50,ngen=3,Reg=fixed,p=0.2,AntagonisticDrivers.csv"))
}
pS.fixed <- matrix(nrow = 250, ncol = 7)
for(m in 1:7){infile <- read.csv(file = paste0("S=",Selfing[m], ",N=1000,B=50,ngen=3,Reg=fixed,p=0.2,AntagonisticDrivers.csv"), head = T);
	sim.freq <- (2*infile[,2] + infile[,3] )/(2*infile[,7])
	pS.fixed[,m] <- sim.freq
	}

# Small populations, fixed N, to assay drift
# Figure S10C

Selfing <- c(0,0.5,1)
for(x in 1:length(Selfing)){
	for(rep in 1:200){
drive.out <- patchdrive(S = Selfing[x], N = 1000, p = 0.2, B =50,ngen = 200, Reg = "fixed");
write.csv(drive.out, paste0("S=",Selfing[x], ",Antagonistic,rep", rep, ".csv"))
points(apply(drive.out, 1, function(x){(2*x[1]+ x[2])/(2*sum(x[1:3]))}), ylim = c(0,0.5), col = rainbow(200)[rep])
	}
}

# Large populations, fixed N
# Model unequal penetrances
# Case 1: resident allele has penetrance 0.6, invader has 0.95
# this for each S value:
Selfing <- c(0,0.05, 0.25,0.5,0.75,0.95,1)
plot(x = c(0,100), y = c(0,1), col = 0)
#antV <- seq(0,1, 0.05)
date();
for(i in 1:7){
drive.out <- as.matrix(patchdrive(N = 20000, S = Selfing[i], B = 50, ngen = 100, Reg = "fixed", p = 0.2, V2 = 0.05, V1 = 0.4))
write.csv(drive.out, paste0("V1.0.4_V2.0.05_N20000.S",Selfing[i],".csv"))
points((drive.out[,1]*2 + drive.out[,2])/40000, col = rainbow(7)[i])
print(date())
}


# Case 2: resident allele has penetrance 0.95, invader has 0.6
# this for each S value:
Selfing <- c(0,0.05, 0.25,0.5,0.75,0.95,1)
plot(x = c(0,100), y = c(0,1), col = 0)
#antV <- seq(0,1, 0.05)
date();
for(i in 1:7){
drive.out <- as.matrix(patchdrive(N = 20000, S = Selfing[i], B = 50, ngen = 100, Reg = "fixed", p = 0.2, V2 = 0.4, V1 = 0.05))
write.csv(drive.out, paste0("V1.0.05_V2.0.4_N20000.S",Selfing[i],".csv"))
points((drive.out[,1]*2 + drive.out[,2])/40000, col = rainbow(7)[i])
print(date())
}

# to plot allele frequencies (Fig S10A):
for(i in 1){
	infile <- read.csv(file = paste0("Antagonistic.N=20000,B=50,S=",Selfing[i],",ngen=100,fixed,p=0.2,V=0.05.csv"))
	plot(apply(infile, 1, function(x){(2*x[2]+ x[3])/(2*sum(x[2:4]))}), ylim = c(0,.3), col = 1, ty = "l", ylab = "Minor drive haplotype frequency", xlab = "Generation")
}
for(i in 2:7){
	infile <- read.csv(file = paste0("Antagonistic.N=20000,B=50,S=",Selfing[i],",ngen=100,fixed,p=0.2,V=0.05.csv"))
	points(apply(infile, 1, function(x){(2*x[2]+ x[3])/(2*sum(x[2:4]))}), ylim = c(0,.3), col = i, ty = "l")
}

# to plot heterozygote frequency (Fig S10B):
for(i in 1){
	infile <- read.csv(file = paste0("Antagonistic.N=20000,B=50,S=",Selfing[i],",ngen=100,fixed,p=0.2,V=0.05.csv"))
	plot(apply(infile, 1, function(x){(x[3])/(sum(x[2:4]))}), ylim = c(0,.5), col = 1, ty = "l", ylab = "Heterozygote frequency", xlab = "Generation")
}
for(i in 2:7){
	infile <- read.csv(file = paste0("Antagonistic.N=20000,B=50,S=",Selfing[i],",ngen=100,fixed,p=0.2,V=0.05.csv"))
	points(apply(infile, 1, function(x){(x[3])/(sum(x[2:4]))}), ylim = c(0,.5), col = i, ty = "l")
}


# to plot allele frequencies for unequal penetrances (Fig S10D):
for(i in 1){ #these infiles have B =50, p = 0.2, Reg = "fixed".
	infile <- read.csv(file = paste0("V1.0.4_V2.0.05_N20000.S",Selfing[i],".csv"))
	plot(apply(infile, 1, function(x){(2*x[2]+ x[3])/(2*sum(x[2:4]))}), ylim = c(0,1), col = 1, ty = "l", ylab = "Strong Medea invader frequency", xlab = "Generation")
}
for(i in 2:7){
	infile <- read.csv(file = paste0("V1.0.4_V2.0.05_N20000.S",Selfing[i],".csv"))
	points(apply(infile, 1, function(x){(2*x[2]+ x[3])/(2*sum(x[2:4]))}), ylim = c(0,.3), col = i, ty = "l")
}
abline(h = 0.2, lty =2)


# to plot heterozygote frequency for unequal penetrances (Fig S10E):
for(i in 1){
	infile <- read.csv(file = paste0("V1.0.4_V2.0.05_N20000.S",Selfing[i],".csv"))
	plot(apply(infile, 1, function(x){(x[3])/(sum(x[2:4]))}), ylim = c(0,1), col = 1, ty = "l", ylab = "Heterozygote frequency", xlab = "Generation")
}
for(i in 2:7){
	infile <- read.csv(file = paste0("V1.0.4_V2.0.05_N20000.S",Selfing[i],".csv"))
	points(apply(infile, 1, function(x){(x[3])/(sum(x[2:4]))}), ylim = c(0,1), col = i, ty = "l")
}





