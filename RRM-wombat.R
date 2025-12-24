###
rm(list=ls())
options(digits=3, width=65)

#### Import data from WOMBAT run
Data <- read.table("mrrtst.dat", header=FALSE)
head(Data)
colnames(Data) <- c("number","animal","sire","dam","cgroup","subject","weight","age")
dim(Data)
##### ==> from : https://github.com/Rostamabd/Random-Regression-Analysis/tree/master

#### Unique ages (control variable for RR) 
Age <- sort(unique(Data[,8]))

### Standardized time values (-1,1)
w <- -1 + 2 * ((Age - min(Age)) / (max(Age) - min(Age)))

### Compute Legendre polynomials (order 0-2, i.e. 3 coefficients for quadratic) 
source("legendre_Coeff.R")  # Ensure Legdre(p) returns at least 3 values (order 0,1,2)

#### Phi matrix (Legendre basis evaluated at standardized ages)
n_coeff <- 3  # Quadratic Legendre: orders 0-2
Phi <- matrix(0, nrow=length(w), ncol=n_coeff)
for(i in 1:length(w)){
  p <- w[i]
  Phi[i,] <- Legdre(p)[1:n_coeff]  # Take first 3 values
}
Phi <- round(Phi, 3)
head(Phi)
dim(Phi)

### Optional: save Legendre matrix
write.table(Phi, file = "plyleg.txt", quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)

### Read WOMBAT solutions for additive genetic RR coefficients (animal effect)
Solutions <- read.table("RnSoln_animal.dat", header=TRUE)
dim(Solutions)
head(Solutions)

### Unique animal IDs
animID <- sort(unique(Solutions$Original_ID))
length(animID)

### Matrix for predicted BV trajectory per age
BV_All_Day <- matrix(NA, nrow=length(animID), ncol=length(Age))

### Loop over animals to compute trajectory
for(i in 1:length(animID)) {
  cat("animal", i, "\n")
  
  # Extract coefficients for this animal, ordered by Tr (1=order 0, 2=order 1, 3=order 2)
  animal_rows <- Solutions[Solutions$Original_ID == animID[i], ]
  coeffs <- animal_rows$Solution[order(animal_rows$Tr)]
  
  # Safety check
  if(length(coeffs) != n_coeff) {
    stop(paste("Animal", animID[i], "has", length(coeffs), "coefficients, expected", n_coeff))
  }
  
  BV_All_Day[i,] <- Phi %*% coeffs
}

### Total BV as sum over all ages (matching original script's summation logic)
BV <- rowSums(BV_All_Day)

### Combine ID and total BV
ID_RR_BV <- cbind(animID, BV)
colnames(ID_RR_BV) <- c("ID", "BV")
head(ID_RR_BV)

### Optional: write results
write.table(ID_RR_BV, file="ID_RR_BV_wombat.txt", quote=FALSE, sep=" ", row.names=FALSE, col.names=TRUE)
