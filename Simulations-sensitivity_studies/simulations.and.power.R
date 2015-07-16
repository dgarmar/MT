
## Test adonis vs asymptotic approximation
## A, B, C(A)  model

#options(scipen=0)

library(vegan)

set.seed(123456789)

hellingerDist <- function (x1, x2) {
  # Function to calculate the hellinger distance between two points
  a <- (sqrt(x1) - sqrt(x2))
  b <- sqrt(sum(a*a))
  return(b)
}

DistH <- function(data) {
  # Function to calculate the interdistances matrix
  k <- ncol(data)
  interdist <- matrix(data=0,nrow=k,ncol=k)
  for(i in 2:k) {
    for(j in 1:(i-1)) {
      interdist[i,j] <- hellingerDist(data[,i],data[,j])
    }
  }
  return(interdist)
}

eigenG <- function (interdist,tol=10^-12) {
  # Function to calculate the eigenvalues
  A <- (- 0.5) * interdist^2
  A <-  A + t(A)
  
  n <- ncol(A)
  I <- diag (1,nrow=n)
  J <- matrix(1/n,ncol=n,nrow=n)
  Aux <- I - J
  
  G <- Aux %*% A %*% Aux
  
  e <- eigen(G,symmetric=T,only.values=T)$values
  
  index <- abs(e) > tol
  
  return(e[index])
  
}

generator <- function (dat,scenario,t=0.5,ha=1){
  # Function to simulate genes according to the proposed scenarios
  if (scenario == 0){ # Null hypothesis case
    
    return(dat)
    
  } else if (scenario == 1){ # 1 factor (A), 1 level 
    
    s <- (a-1)*b*n # 80 in the example case
    dat.alt <- apply(dat[,(s+1):ncol(dat)],2,function(x){HoHa(x,t=t,ha=ha)}) # x is p0
    dat[,(s+1):ncol(dat)] <- dat.alt[,sample(ncol(dat.alt))] # suffle
    return (dat)
    
  } else if (scenario == 2){ # 1 factor (A), 2 levels
    
    s1.1 <- (a-2)*b*n # 40 in the example case
    s1.2 <- (a-1)*b*n # 80 in the example case
    
    dat.alt1.1 <- apply(dat[,(s1.1+1):s1.2],2,function(x){HoHa(x,t=t,ha=ha)}) 
    dat.alt1.2 <- apply(dat[,(s1.2+1):ncol(dat)],2,function(x){HoHa(x,t=t*2,ha=ha)}) # # t changes between levels (*1.5, arbitrary)
    
    dat[,(s1.1+1):s1.2] <- dat.alt1.1[,sample(ncol(dat.alt1.1))] 
    dat[,(s1.2+1):ncol(dat)] <- dat.alt1.2[,sample(ncol(dat.alt1.2))] 
    return (dat)
    
  } else if (scenario == 3){ # 2 factors (A,B), 1 level in both
    
    s1.1 <- (a-1)*b*n # 80 in the example case
    s1.1b <- (a-2)*b*n # 40 in our example case
    s2.1 <- n*(a*b-1) # 110 in the example case
    s2.1b <- n*((a-1)*b-1) # 70 in the example case
    s2.1c <- n*((a-2)*b-1) # 30 in the example case
    if (ha != 6) {ha2 <- ha+1} else { ha2=1 } 
    
    dat.alt1.1 <- apply(dat[,(s1.1+1):s2.1],2,function(x){HoHa(x,t=t,ha=ha)}) 
    dat.alt2.1 <- apply(dat[,(s2.1+1):ncol(dat)],2,function(x){HoHa(x,t=t,ha=ha2)}) # ha changes between factors
    
    dat[,(s1.1+1):s2.1] <- dat.alt1.1[,sample(ncol(dat.alt1.1))] 
    dat[,(s2.1+1):ncol(dat)] <- dat.alt2.1[,sample(ncol(dat.alt2.1))]
    dat[,(s2.1b+1):s1.1] <- dat.alt2.1[,sample(ncol(dat.alt2.1))]
    dat[,(s2.1c+1):s1.1b] <- dat.alt2.1[,sample(ncol(dat.alt2.1))]
    return (dat)
  }
}

HoHa <- function (p0, ha, t) {
  # Function to obtain Ha from Ho
  points <- list(c(1,0,0),c(0,1,0),c(0,0,1),c(2/3,1/3,0),c(0,1/3,2/3),c(2/3,0,1/3))
  p <- points[[ha]] # To select one of the points to build the straigth line
  return ( (p-p0)*t + p0 ) # Straigth line equation. Given two points and a t value, returns a point over it.
  # The corresponding angles to the different t values and p's can be obtained from the showAngle function
  
}

showAngle <- function (p0=c(0.4,0.3,0.2), ha=1, t=0.5){
  # Function to compute the corresponding angle given different t and p values
  
  points <- list(c(1,0,0),c(0,1,0),c(0,0,1),c(2/3,1/3,0),c(0,1/3,2/3),c(2/3,0,1/3))
  p <- points[[ha]] # To select one of the points to build the straigth line
  p1 <- (p-p0)*t + p0
  
  c11 <- t(p0)%*%p0
  c22 <- t(p1)%*%p1
  c12 <- t(p0)%*%p1
  
  return(as.numeric(acos(c12/sqrt(c11*c22))*180/pi))
  
}

# Simulation scenario: truncated normal 

sim.tnorm <- function(nb.dim = 3,nb.perm = 10000, m = 0.4){ 
  # Function to perturb Ho
  sim = NULL
  sumJ = rep(0,nb.perm)
  sdSim <- 0.01*rnorm(n=1,mean=2,sd=0.25)
  for(i in 1:(nb.dim-1)){
    simI = NULL
    
    if (i!=1) m = (1-m)/2  
    
    for(j in 1:nb.perm){
      simIJ = rnorm(n=1,mean=m,sd=sdSim)
      if (simIJ < 0) simIJ <- 0
      if (simIJ > 1-sumJ[j]) simIJ <- 1-sumJ[j]
      sumJ[j] = sumJ[j] + simIJ
      simI = c(simI,simIJ)
    }
    sim = rbind(sim,simI)
  }
  sim = rbind(sim,simI)
  for(j in 1:nb.perm){
    sim[nb.dim,j] <- 1- sum(sim[(1:(nb.dim-1)),j])
  }
  return(sim)
}


nG   <- 1000    # number of genes
nS   <- 3           # number of splicing isoforms

a    <- 3     # phenotypes / conditions
b    <- 4     # tissues
n    <- 50    # replicates (individuals) nested to phenotypes


nb.mont        <- 10^4  # number of Montecarlo generations per gene

nb.perm = 5000


labelA <- gl(a, b*n)                    # df = a-1
labelB <- gl(b, n, length = a*b*n)      # df = b-1
labelI <- gl(n,1,length=a*b*n)          # df = a*(c-1)

#                                    AB   df = (a-1)*(b-1)
#                                   Error df = a*(b-1)*(n-1)


eigenStats <- matrix(ncol=3,nrow=nG+1)
pvals      <- matrix(ncol=3,nrow=nG)
scores     <- matrix(ncol=3,nrow=nG)


# Parameters to change 
m.values   <- c(0.4,0.9)       # m, most abundant isoform's relative abundance in Ho 
ha.values  <- c(1,6)               # ha, type of Ha: Select one sample case
t.values   <- seq(-1,1,0.01)   # t, degree of "alternativity". Related to the angle
sc.values  <- 1:3                  # sc, scenario: number of factors (and levels) that change:
# scenario = 0  -->  Null hypothesis case
# scenario = 1  -->  1 factor (A), 1 level 
# scenario = 2  -->  1 factor (A), 2 levels
# scenario = 3  -->  2 factors (A,B), 1 level in both

iterations <- length(sc.values)*length(m.values)*length(ha.values)*length(t.values) 
storage <- NULL
error <- 0

j = 1 # Counter

for (sc in sc.values){
  for (m in m.values){
    for (ha in ha.values){
      for (t in t.values){
        
        for (i in 1:(nG+1)) { 
          dat <- sim.tnorm(nS,(n*a*b),m)  # dat contains the simulated data under Ho
          if (i != nG+1) dat <- generator (dat, sc,t, ha) # Generate Ha from Ho 
          
          nb.mont <- 10^4
          if (i == nG+1) nb.mont <- 10^6 # Only for last iteration, from which the null distribution is obtained
          
          if (sum(dat<0)>0) {
            error <- 1 
            break
          }
          
          d = DistH(dat) # d contains the interdistances matrix
          e <- eigenG(d) # e contains the eigenvalues
          
          eigenStats [i,1] <- length (e) # Summary of eigenvalues
          eigenStats [i,2] <- sum(e>0)
          eigenStats [i,3] <- sum(e<0)
          
          if (eigenStats [i,3]>0)  e <- abs(e) # If eigenvalues <0 |eigenvalues|
          
          # Asymptotic Montecarlo distribution 
          # (approximation of the statistic distribution under Ho from chi2 variables)
          
          randomChisqNA  <- matrix(
            rchisq(nb.mont*eigenStats[i,1],df=a-1),
            nrow=eigenStats[i,1],ncol=nb.mont)
          
          randomChisqNB  <- matrix(
            rchisq(nb.mont*eigenStats[i,1],df=b-1),
            nrow=eigenStats[i,1],ncol=nb.mont)
          
          randomChisqNAB <- matrix(
            rchisq(nb.mont*eigenStats[i,1],df=((a-1)*(b-1))),
            nrow=eigenStats[i,1],ncol=nb.mont)
          
          randomChisqInd <- matrix(
            rchisq(nb.mont*eigenStats[i,1],df=(a*(n-1))),
            nrow=eigenStats[i,1],ncol=nb.mont)
          
          randomChisqError <- matrix(
            rchisq(nb.mont*eigenStats[i,1],df=(a*(b-1)*(n-1))),
            nrow=eigenStats[i,1],ncol=nb.mont)
          
          asymptMSA  <- e %*% randomChisqNA # Here asymptMSA refers to SSA
          asymptMSB  <- e %*% randomChisqNB
          asymptMSAB <- e %*% randomChisqNAB
          asymptMSE  <- e %*% randomChisqError
          asymptMSI  <- e %*% randomChisqInd
          
          asymptFA  <- asymptMSA / asymptMSI * (a*(n-1)) / (a-1) # Obtaining the F scores
          asymptFB  <- asymptMSB / asymptMSE * (a*(b-1)*(n-1)) / (b-1)
          asymptFAB <- asymptMSAB / asymptMSE * (a*(b-1)*(n-1))/(a-1)/(b-1)
          
          if (i == nG+1){break}
          
          ado <- adonis(as.dist(d) ~ labelA*labelB + labelI %in% labelA, permutations=1) 
          
          # always balanced groups
          pseudoFA  <- ado$aov.tab[1,4]     # pseudo F score
          pseudoFB  <- ado$aov.tab[2,4]     # pseudo F score
          pseudoFAB <- ado$aov.tab[3,4]     # pseudo F score
          
          #  we have to change the F numerator for factor A, individuals are random  
          pseudoFA  <- ado$aov.tab[1,4]/ado$aov.tab[4,4] # pseudo F score for A
          
          
          # F scores storage
          scores[i,1] <- pseudoFA           # equal groups
          scores[i,2] <- pseudoFB           # equal groups
          scores[i,3] <- pseudoFAB          # equal groups
          
          
          # pvalues computation and storage
          
          # pvalue asymptotic for A
          #           pvals[i,1] <- sum((asymptFA>pseudoFA)/nb.mont)
          #           #pvalue asymptotic for B
          #           pvals[i,2] <- sum((asymptFB>pseudoFB)/nb.mont)
          #           #pvalue asymptotic for AB interaction
          #           pvals[i,3] <- sum((asymptFAB>pseudoFAB)/nb.mont)
          
        }
        
        # Compute sensitivity
        if (error != 1){
          alpha <- 0.001 # significance level
          
          asymptFA <- sort(asymptFA[1,],decreasing=T) # Corresponding to the last iteration(Ho)
          asymptFB <- sort(asymptFAB[1,],decreasing=T)
          asymptFAB <- sort(asymptFAB[1,],decreasing=T)
          
          #  Critical values at the given significance level
          
          CaA <- asymptFA[round(alpha*nb.mont)]
          CaB <- asymptFB[round(alpha*nb.mont)]
          CaAB <- asymptFAB[round(alpha*nb.mont)] 
          
          # ErrII and potency
          
          potA <- (sum(scores[,1] >= CaA)/nG)
          potB <- (sum(scores[,2] >= CaB)/nG)
          potAB <- (sum(scores[,3] >= CaAB)/nG)
          
          # Compute angle for storing and store in "storage" variable
          ho <- c(m,(1-m)/2,(1-m)/2)
          angle <- showAngle(ho,ha,t)
          alt <- HoHa (ho,ha,t)
          storage <- rbind(storage,c(sc,m,ha,t,angle,potA,potB,potAB,alt[1],alt[2],alt[3]))
          cat(j,"\n")
          j = j + 1
          
        }else{
          error <- 0
          ho <- c(m,(1-m)/2,(1-m)/2)
          angle <- showAngle(ho,ha,t)
          alt <- HoHa (ho,ha,t)
          storage <- rbind(storage,c(sc,m,ha,t,angle,NA,NA,NA,NA,NA,NA))
        }
        
      }    
    }
  }
}


colnames(storage) <- c("scenario","ho","ha","t","angle","potA","potB","potAB","alt1","alt2","alt3")
storage <- data.frame(storage)

save.image("Ha.dgarrido.RData")
