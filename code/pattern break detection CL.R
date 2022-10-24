nDevYears <- 9
paidMatrix <- matrix(NA, nrow = nDevYears, ncol = nDevYears)

paidMatrix[1:9, 1] <- rpois(9, 2000)
paidMatrix[1:8, 2] <- paidMatrix[1:8, 1] + rpois(8, 800)
paidMatrix[1:7, 3] <- paidMatrix[1:7, 2] + rpois(7, 400)
paidMatrix[1:6, 4] <- paidMatrix[1:6, 3] + rpois(6, 200)
paidMatrix[1:5, 5] <- paidMatrix[1:5, 4] + rpois(5, 100)
paidMatrix[1:4, 6] <- paidMatrix[1:4, 5] + rpois(4, 100)
paidMatrix[1:3, 7] <- paidMatrix[1:3, 6] + rpois(3, 80)
paidMatrix[1:2, 8] <- paidMatrix[1:2, 7] + rpois(2, 40)
paidMatrix[1, 9] <- paidMatrix[1:1, 8] + rpois(1, 10)

#adapted triangle

paidMatrix[9, 1] <- paidMatrix[9, 1]*1.1
paidMatrix[8, 2] <- paidMatrix[8, 2]*1.1
paidMatrix[7, 3] <- paidMatrix[7, 3]*1.1
paidMatrix[6, 4] <- paidMatrix[6, 4]*1.1
paidMatrix[5, 5] <- paidMatrix[5, 5]*1.1
paidMatrix[4, 6] <- paidMatrix[4, 6]*1.1
paidMatrix[3, 7] <- paidMatrix[3, 7]*1.1
paidMatrix[2, 8] <- paidMatrix[2, 8]*1.1
paidMatrix[1, 9] <- paidMatrix[1, 9]*1.1

devFac <- rep(NA, nDevYears - 1)
for(iRun in 1:(nDevYears - 1)){
  devFac[iRun] <- sum(paidMatrix[1:(nDevYears - iRun), iRun + 1])/sum(paidMatrix[1:(nDevYears - iRun), iRun])
  #plot(paidMatrix[1:(nDevYears - iRun), iRun], paidMatrix[1:(nDevYears - iRun), iRun + 1], pch = 18, cex = 1.5, col = 4)
  #lines(c(-20000, 20000), c(devFac[iRun]*-20000, devFac[iRun]*20000), lty = 3, lwd = 1.5) 
}

devFac

completeTriangle <- paidMatrix
for(iRun in 1:(nDevYears - 1)){
  completeTriangle[(nDevYears - iRun + 1):nDevYears, iRun + 1] <- 
    devFac[iRun]*completeTriangle[(nDevYears - iRun + 1):nDevYears, iRun]
}
completeTriangle

#Extraction of the ultimate values per accident year.

ultimates <- completeTriangle[, nDevYears] 
reserves <- rep(0, nDevYears)
for(iActual in 1:nDevYears){
  reserves[iActual] <- ultimates[iActual] - completeTriangle[iActual, (nDevYears - iActual + 1)] 
}
round(reserves, digits = 2)
sum(reserves)


#Calculation of the different types of development factors.

devFacMean <- rep(NA, nDevYears - 1)
devFacWMean <- rep(NA, nDevYears - 1)
devFacMin <- rep(NA, nDevYears - 1)
devFacMax <- rep(NA, nDevYears - 1)

weights <- sample(rnorm(nrow(paidMatrix) - 1, mean = 1, sd = 0.5))

for(iRun in 1:(nDevYears - 1)){
  
  devFacMean[iRun] <- mean(paidMatrix[1:(nDevYears - iRun), iRun + 1]/paidMatrix[1:(nDevYears - iRun), iRun], na.rm = T)
  devFacWMean[iRun] <- sum(weights[1:(nrow(paidMatrix) - iRun)]*(paidMatrix[1:(nDevYears - iRun), iRun + 1]/paidMatrix[1:(nDevYears - iRun), iRun]), na.rm = T)/sum(weights[1:(nrow(paidMatrix) - iRun)], na.rm = T)
  devFacMin[iRun] <- min(paidMatrix[1:(nDevYears - iRun), iRun + 1]/paidMatrix[1:(nDevYears - iRun), iRun], na.rm = T)
  devFacMax[iRun] <- max(paidMatrix[1:(nDevYears - iRun), iRun + 1]/paidMatrix[1:(nDevYears - iRun), iRun], na.rm = T)
}

plot(devFac, ylim = c(1, 1.4))
lines(devFacMean)
lines(devFacWMean)
lines(devFacMin)
lines(devFacMax)

lambda <- rev(cumprod(rev(devFac)))
#lambda <- rev(cumprod(rev(devFacMax)))
alpha <- 1/lambda

library(ChainLadder)
library(ggplot2)
library(plotly)
library(stringr)
library(plyr)

# library(gdata)
# paidVec <- unmatrix(paidMatrix, byrow=T)
# 
# #The names of the elements are a bit awkward, let's make them more recognizable
# 
# names(paidVec) <- laply(str_split(names(paidVec), ":"), function(xx) str_replace_all(xx[1], "r", ""))
# 
# #In order to make the ChainLadder functions work, the missing values need to be 
# #removed from the vector (i.e only the observed triangle needs to be stored).
# 
# paidCumVec <- paidVec[!is.na(paidVec)]
# paidCumTr <- triangle(paidCumVec)
# colnames(paidCumTr) <- as.character(1:ncol(paidCumTr))
# 
# #This enables us to plot the dev patterns quite easily, using the usual functions:
# 
# plot(paidCumTr, ylab = 'Observed Payments')
# plot(paidCumTr, lattice = TRUE)

resultMack <- MackChainLadder(paidMatrix)

f_vec <- resultMack$f
sigma_vec <- resultMack$sigma
F_list <- list()
C_list <- list()

#next some different definitions for the pearsonReds object can be found, please just select one only!

#Definition 1: take all residuals into account 

pearsonReds <- c()

for (iCol in 1:(nDevYears - 1)) {
    F_list[[iCol]] <- paidMatrix[1:(nDevYears - iCol), iCol + 1] / paidMatrix[1:(nDevYears - iCol), iCol]
    C_list[[iCol]] <- paidMatrix[1:(nDevYears - iCol), iCol]
    pearsonReds <- c(pearsonReds, (F_list[[iCol]] - f_vec[iCol]) * sqrt(C_list[[iCol]]) / sqrt(sigma_vec[iCol]))
}

pearsonReds <- sample(pearsonReds, 28)

#Definition 2: not most recent year

pearsonReds <- c()

for(iCol in 1:(nDevYears - 2)){
  F_list[[iCol]] <- paidMatrix[1:(nDevYears - iCol - 1), iCol + 1]/paidMatrix[1:(nDevYears - iCol - 1), iCol]
  C_list[[iCol]] <- paidMatrix[1:(nDevYears - iCol - 1), iCol]
  pearsonReds <- c(pearsonReds, (F_list[[iCol]] - f_vec[iCol])*sqrt(C_list[[iCol]])/sqrt(sigma_vec[iCol]))
}

#Definition 3: not 2 most recent years

pearsonReds <- c()

for(iCol in 1:(nDevYears - 3)){
  F_list[[iCol]] <- paidMatrix[1:(nDevYears - iCol - 2), iCol + 1]/paidMatrix[1:(nDevYears - iCol - 2), iCol]
  C_list[[iCol]] <- paidMatrix[1:(nDevYears - iCol - 2), iCol]
  pearsonReds <- c(pearsonReds, (F_list[[iCol]] - f_vec[iCol])*sqrt(C_list[[iCol]])/sqrt(sigma_vec[iCol]))
}

#Definition 4: not second most recent year

pearsonReds <- c()

for(iCol in 1:(nDevYears - 1)){
  temp <- paidMatrix[1:(nDevYears - iCol), iCol + 1]/paidMatrix[1:(nDevYears - iCol), iCol]
  if(length(temp) > 1){
    F_list[[iCol]] <- temp[-(length(temp) - 1)]
  } else {
    F_list[[iCol]] <- temp
  }
  temp <- paidMatrix[1:(nDevYears - iCol), iCol]
  if(length(temp) > 1){
    C_list[[iCol]] <- temp[-(length(temp) - 1)]
  } else {
    C_list[[iCol]] <- temp
  }
  pearsonReds <- c(pearsonReds, (F_list[[iCol]] - f_vec[iCol])*sqrt(C_list[[iCol]])/sqrt(sigma_vec[iCol]))
}

#Definition 5: not third most recent year

pearsonReds <- c()
bound <- 3

for(iCol in 1:(nDevYears - 1)){
  temp <- paidMatrix[1:(nDevYears - iCol), iCol + 1]/paidMatrix[1:(nDevYears - iCol), iCol]
  if(length(temp) > (bound - 1)){
    F_list[[iCol]] <- temp[-(length(temp) - bound + 1)]
  } else {
    F_list[[iCol]] <- temp
  }
  temp <- paidMatrix[1:(nDevYears - iCol), iCol]
  if(length(temp) > (bound - 1)){
    C_list[[iCol]] <- temp[-(length(temp) - bound + 1)]
  } else {
    C_list[[iCol]] <- temp
  }
  pearsonReds <- c(pearsonReds, (F_list[[iCol]] - f_vec[iCol])*sqrt(C_list[[iCol]])/sqrt(sigma_vec[iCol]))
}

#Once you have selected some definition for the pearsonReds object, let's start the simulation study

obsDiag <- rep(NA, nDevYears)
for(iRow in 1:(nDevYears)){
  iCol <- nDevYears - iRow + 1
  obsDiag[iRow] <- paidMatrix[iRow, iCol]
}

nBoot <- 1000
nObs <- length(pearsonReds)
reserveBoot <- rep(0, nBoot)

for(iBoot in 1:nBoot){
  bootReds <- sample(pearsonReds, nObs)
  bootReds_list <- list()
  counter <- 1
  for(iCol in 1:(nDevYears - 1)){
    bootReds_list[[iCol]] <- bootReds[counter:(counter + (nDevYears - iCol) - 1)]
    counter <- counter + (nDevYears - iCol)
  }
  
  F_boot_list <- list()
  f_boot_vec <- rep(NA, nDevYears - 1)
  sigma_boot_vec <- rep(NA, nDevYears - 1)
  for(iCol in 1:(nDevYears - 1)){
    F_boot_list[[iCol]] <- (bootReds[[iCol]]*sqrt(sigma_vec[iCol]))/sqrt(C_list[[iCol]]) + f_vec[iCol]
    f_boot_vec[iCol] <- sum(F_boot_list[[iCol]]*C_list[[iCol]])/sum(C_list[[iCol]])
    sigma_boot_vec[iCol] <- mean((F_boot_list[[iCol]] - f_vec[iCol])^2*C_list[[iCol]])
  }
  
  paidMatrix_boot <- paidMatrix
  diag <- 0
  for(iDiag in 1:(nDevYears - 1)){
    for(iRow in (1 + iDiag):(nDevYears)){
      selCol <- nDevYears - iRow + iDiag + 1
      paidMatrix_boot[iRow, selCol] <-
          rnorm(1,
              mean = paidMatrix_boot[iRow, selCol - 1] * f_boot_vec[selCol - 1],
              sd = sqrt(paidMatrix_boot[iRow, selCol - 1] * sigma_boot_vec[selCol - 1]))
      print(paidMatrix_boot[iRow, selCol])
    }
  }
  reserveBoot[iBoot] <- sum(paidMatrix_boot[, nDevYears] - obsDiag)
}

summary(reserveBoot)
quantile(reserveBoot, 0.025)
quantile(reserveBoot, 0.05)
quantile(reserveBoot, 0.95)
quantile(reserveBoot, 0.975)

> summary(reserveBoot)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
4648    9008   11118   11587   13521   24958 #Definition 1: all

> summary(reserveBoot)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
5346    7421    8121    8017    8677   10517 #Definition 2: not most recent year

> summary(reserveBoot)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
5535    7686    8322    8352    8956   10916  #Definition 3: not 2 most recent years

> summary(reserveBoot)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
5420   10133   12316   12706   14610   26003 #Definition 4: not 2 most recent year

> summary(reserveBoot)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
5476    9644   12044   12308   14488   25595 #Definition 5: not 3 most recent year
