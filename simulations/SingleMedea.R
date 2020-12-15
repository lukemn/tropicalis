
# Simulations with a single driver

# This file starts with a function to perform simulations, as described in the Methods section of the paper.
# The function starts by forming the starting population, with genotypes and sexes drawn from the equilibrium probability distribution.
# Then it loops through ngen generations and returns data.frame(GenotypeFreqs, Sexes), where each row is a generation.
# Selection operates in a Wright-Fisher framework. 

###########################################################

drive <- function(p = 0.05, S = 0.95, N = 100, V = 0.05, B = 50, Him = 0.005, ngen = 100, Reg = "fixed"){
	
	countGenos <- function(x,y){
		counts <- apply(rbind(x,y), 1, sum)
		cbind(sum(counts ==0), sum(counts == 1), sum(counts ==2))
	}
		
	# Starting genotype frequencies from extended HWE
	Fhat <- S/(2-S) # equilibrium inbreeding coeffecient given S
	
	JJh <- 	((1-Fhat)*p^2 + p*Fhat)*(1+Fhat)/2
	JNh <- 	(1-Fhat)*2*p*(1-p)*(1+Fhat)/2
	NNh  <-	((1-Fhat)*(1-p)^2 + (1-p)*Fhat)*(1+Fhat)/2
	JJm <- 	((1-Fhat)*p^2 + p*Fhat)*(1-Fhat)/2
	JNm <- 	(1-Fhat)*2*p*(1-p)*(1-Fhat)/2
	NNm <-	((1-Fhat)*(1-p)^2 + (1-p)*Fhat)*(1-Fhat)/2
	c(JJh,JNh,NNh,JJm,JNm,NNm)

	StartingGenotypes <- rmultinom(1,N, c(JJh,JNh,NNh,JJm,JNm,NNm))
	
	herms <- matrix(nrow = sum(StartingGenotypes[c(1:3)]), ncol = 2)
	herms[,1] <- c(rep(0,StartingGenotypes[1]), rep(0, StartingGenotypes[2]), rep(1, StartingGenotypes[3]))
	herms[,2] <- c(rep(0,StartingGenotypes[1]), rep(1, StartingGenotypes[2]), rep(1, StartingGenotypes[3]))
	
	males <- matrix(nrow = sum(StartingGenotypes[c(4:6)]), ncol = 2)
	males[,1] <- c(rep(0,StartingGenotypes[4]), rep(0, StartingGenotypes[5]), rep(1, StartingGenotypes[6]))
	males[,2] <- c(rep(0,StartingGenotypes[4]), rep(1, StartingGenotypes[5]), rep(1, StartingGenotypes[6]))
	
	
	GenotypeFreqs <- matrix(nrow = ngen, ncol = 3)
	GenotypeFreqs[1,] <- countGenos(herms,males)
	
	Sexes <- matrix(nrow = ngen, ncol = 2)
	Sexes[1,] <- c(dim(herms)[1], dim(males)[1])	
	
	for( gen in 2: ngen){
	  
	  if (gen %% (round(ngen/10))==0) cat(sprintf('gen %s / %s\n', gen, ngen))
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
	          kidviability[j] <- rbinom(1,1,V)
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
	            kidviability[k] <- rbinom(1,1,V)
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
	          kidviability[j] <- rbinom(1,1,V)
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
	  
	  
	} # Loop through generations. 
	
	data.frame(GenotypeFreqs, Sexes)
}

###########################################################

# Now we can call this function to simulate data

# Large populations (Figure 12A)

Selfing <- c(0,0.05, 0.25,0.5,0.75,0.95,1)
Selfing <- 0.95
for(x in 1:length(Selfing)){
  drive.out <- drive(p = 0.05, S = Selfing[x], N = 20000, ngen=1000);
  # drive.out_95 <- mclapply(1:10, mc.cores=5, function(i) drive(p = 0.05, S = 0.95, N = 1000, ngen=100))
  # drive.out_1 <- mclapply(1:10, mc.cores=5, function(i) drive(p = 0.05, S = 1, N = 1000, ngen=100))
  # dro = as.data.frame(do.call(rbind, lapply(1:5, function(i) cbind(i, g = 1:nrow(drive.out_95[[i]]), f = 1-apply(drive.out_95[[i]], 1, function(x){(2*x[3]+ x[2])/(2*sum(x[1:3]))})))))
  # ggplot(dro, aes(g, f, col = factor(i))) + geom_point()
  write.csv(drive.out, paste0("S=",Selfing[x], ".csv"))
  plot(1-apply(drive.out, 1, function(x){(2*x[3]+ x[2])/(2*sum(x[1:3]))}), ylim = c(0,1))
}

# Small populations (Figure 12B)

Selfing <- c(0,0.05, 0.25,0.5,0.75,0.95,1)
for(x in 1:length(Selfing)){
	for(rep in 1:250){
drive.out <- drive(S = Selfing[x], N = 1000, p = 0.05);
write.csv(drive.out, paste0("S=",Selfing[x], "rep", rep, ".csv"))
plot(apply(drive.out, 1, function(x){(2*x[3]+ x[2])/(2*sum(x[1:3]))}), ylim = c(0,0.5))
	}
}

S1 <- matrix(nrow = 100, ncol = 250)
for(i in 1:250){infile <- read.csv(file = paste0("S=1rep",i,".csv"), head = T);
	sim.freq <- (2*infile[,2] + infile[,3] )/2000
	S1[,i] <- sim.freq}
S0.95 <- matrix(nrow = 100, ncol = 250)
for(i in 1:250){infile <- read.csv(file = paste0("S=0.95rep",i,".csv"), head = T);
	sim.freq <- (2*infile[,2] + infile[,3] )/2000
	S0.95[,i] <- sim.freq}
S0.75 <- matrix(nrow = 100, ncol = 250)
for(i in 1:250){infile <- read.csv(file = paste0("S=0.75rep",i,".csv"), head = T);
	sim.freq <- (2*infile[,2] + infile[,3] )/2000
	S0.75[,i] <- sim.freq}
S0.5 <- matrix(nrow = 100, ncol = 250)
for(i in 1:250){infile <- read.csv(file = paste0("S=0.5rep",i,".csv"), head = T);
	sim.freq <- (2*infile[,2] + infile[,3] )/2000
	S0.5[,i] <- sim.freq}
S0.25 <- matrix(nrow = 100, ncol = 250)
for(i in 1:250){infile <- read.csv(file = paste0("S=0.25rep",i,".csv"), head = T);
	sim.freq <- (2*infile[,2] + infile[,3] )/2000
	S0.25[,i] <- sim.freq}
S0.05 <- matrix(nrow = 100, ncol = 250)
for(i in 1:250){infile <- read.csv(file = paste0("S=0.05rep",i,".csv"), head = T);
	sim.freq <- (2*infile[,2] + infile[,3] )/2000
	S0.05[,i] <- sim.freq}
S0 <- matrix(nrow = 100, ncol = 250)
for(i in 1:250){infile <- read.csv(file = paste0("S=0rep",i,".csv"), head = T);
	sim.freq <- (2*infile[,2] + infile[,3] )/2000
	S0[,i] <- sim.freq}


# Plot the results (Fig 12B)
classes <- list(S1, S0.95, S0.75, S0.5, S0.25, S0.05, S0)
par(mfrow= c(4,2))
for(q in 1:7){hist(classes[[q]][100,], br = seq(-0.01, 1.01,.01), ylim = c(0,250))}
