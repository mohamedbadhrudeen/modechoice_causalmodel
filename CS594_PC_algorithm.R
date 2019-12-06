library(data.table)
library(caret)
library(dimRed)
library(pcalg)
#library(varimpact)

#Copula Estimate

## For details about the arguements and outputs, refer to function 'sbgcop.mcmc' in R package 'sbgcop',
## https://cran.r-project.org/web/packages/sbgcop/index.html.
copula.estimate <- function (Y, n0 = dim(Y)[2] + 1, S0 = diag(dim(Y)[2])/n0, nsamp = 100, 
                             odens = max(1, round(nsamp/1000)), impute = any(is.na(Y)), 
                             plugin.threshold = 100, plugin.marginal = (apply(Y, 2, function(x) {
                               length(unique(x))
                             }) > plugin.threshold), seed = NULL, verb = TRUE) 
{
  require(sbgcop)
  ok_S0 <- all(eigen(S0)$val > 0) & dim(S0)[1] == dim(Y)[2] & 
    dim(S0)[2] == dim(Y)[2]
  ok_n0 <- (n0 >= 0)
  if (!ok_S0) {
    stop("Error: S0 must be a positive definite p x p matrix \n")
  }
  if (!ok_n0) {
    stop("Error: n0 must be positive \n")
  }
  
  vnames <- colnames(Y)
  Y <- as.matrix(Y)
  colnames(Y) <- vnames
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  set.seed(seed)
  R <- NULL
  for (j in 1:p) {
    R <- cbind(R, match(Y[, j], sort(unique(Y[, j]))))
  }
  Rlevels <- apply(R, 2, max, na.rm = TRUE)
  Ranks <- apply(Y, 2, rank, ties.method = "max", na.last = "keep")
  N <- apply(!is.na(Ranks), 2, sum)
  U <- t(t(Ranks)/(N + 1))
  Z <- qnorm(U)
  Zfill <- matrix(rnorm(n * p), n, p)
  Z[is.na(Y)] <- Zfill[is.na(Y)]
  S <- cov(Z)
  Y.pmean <- Y
  if (impute) {
    Y.pmean <- matrix(0, nrow = n, ncol = p)
  }
  LPC <- NULL
  C.psamp <- array(dim = c(p, p, floor(nsamp/odens)))
  Y.imp <- NULL
  if (impute) {
    Y.imp <- array(dim = c(n, p, floor(nsamp/odens)))
  }
  dimnames(C.psamp) <- list(colnames(Y), colnames(Y), 
                            1:floor(nsamp/odens))
  for (ns in 1:nsamp) {
    for (j in sample(1:p)) {
      Sjc <- S[j, -j] %*% solve(S[-j, -j])
      sdj <- sqrt(S[j, j] - S[j, -j] %*% solve(S[-j, 
                                                 -j]) %*% S[-j, j])
      muj <- Z[, -j] %*% t(Sjc)
      if (!plugin.marginal[j]) {
        for (r in 1:Rlevels[j]) {
          ir <- (1:n)[R[, j] == r & !is.na(R[, j])]
          lb <- suppressWarnings(max(Z[R[, j] == r - 
                                         1, j], na.rm = TRUE))
          ub <- suppressWarnings(min(Z[R[, j] == r + 
                                         1, j], na.rm = TRUE))
          Z[ir, j] <- qnorm(runif(length(ir), pnorm(lb, 
                                                    muj[ir], sdj), pnorm(ub, muj[ir], sdj)), 
                            muj[ir], sdj)
        }
      }
      ir <- (1:n)[is.na(R[, j])]
      Z[ir, j] <- rnorm(length(ir), muj[ir], sdj)
    }
    # relocate the mean to zero
    # added by Ruifei Cui
    Z = t( (t(Z)-apply(Z,2,mean)))
    
    S <- solve(rwish(solve(S0 * n0 + t(Z) %*% Z), n0 + 
                       n))
    if (ns%%odens == 0) {
      C <- S/(sqrt(diag(S)) %*% t(sqrt(diag(S))))
      lpc <- ldmvnorm(Z %*% diag(1/sqrt(diag(S))), 
                      C)
      LPC <- c(LPC, lpc)
      C.psamp[, , ns/odens] <- C
      if (impute) {
        Y.imp.s <- Y
        for (j in 1:p) {
          Y.imp.s[is.na(Y[, j]), j] <- quantile(Y[, 
                                                  j], pnorm(Z[is.na(Y[, j]), j], 0, sqrt(S[j, 
                                                                                           j])), na.rm = TRUE, type = 1)
        }
        Y.imp[, , ns/odens] <- Y.imp.s
        Y.pmean <- ((ns/odens - 1)/(ns/odens)) * Y.pmean + 
          (1/(ns/odens)) * Y.imp.s
      }
    }
    if (verb == TRUE & (ns%%(odens * 10)) == 0) {
      cat(round(100 * ns/nsamp), "percent done ", 
          date(), "\n")
    }
  }
  G.ps <- list(C.psamp = C.psamp, Y.pmean = Y.pmean, Y.impute = Y.imp, 
               LPC = LPC)
  class(G.ps) <- "psgc"
  return(G.ps)
}



setwd("C:/Users/WineGlass/Desktop/CS594_Causal_Inference/Term Project/Data_NHTS")
csv_files <- list.files(pattern = "csv$")

#Reading the data
data_ <- fread('trippub.csv', header = T, sep = ',')
#myvars <- c("TRVLCMIN", "TRPMILES", "GASPRICE")
data_1 <- data_[, .(R_AGE, R_SEX, EDUC, HHSIZE,  URBRUR, HHFAMINC, TRIPPURP, GASPRICE, TRVLCMIN, TRPTRANS)]

s <- data_1[TRVLCMIN> 0 & GASPRICE > 0 & HHFAMINC > 0 & HHSIZE > 0 & TRIPPURP > 0 & R_AGE > 0&  R_SEX > 0 & EDUC > 0 & TRPTRANS > 0] 

s$AGEC[s$R_AGE <18] = 1
s$AGEC[s$R_AGE >=18 & s$R_AGE < 25] = 2
s$AGEC[s$R_AGE >=25 & s$R_AGE < 35] = 3
s$AGEC[s$R_AGE >=35 & s$R_AGE < 50] = 4
s$AGEC[s$R_AGE >=50 & s$R_AGE < 70] = 5
s$AGEC[s$R_AGE >= 70] = 6

dt = c(3, 4, 5, 6, 11, 16)
s1 <- s[TRPTRANS %in% dt]
s1 <- s1[TRPTRANS %in% c(3, 4,5,6), TRPTRANS := 0] 
s1 <- s1[TRPTRANS %in% c(11,16), TRPTRANS := 1] 

s1$TRIPPURP <- as.factor(s1$TRIPPURP)
s1$URBRUR <- as.factor(s1$URBRUR)
s1$R_SEX <- as.factor(s1$R_SEX)
s1$EDUC <- as.factor(s1$EDUC)
s1$HHFAMINC <- as.factor(s1$HHFAMINC)
s1$AGEC <- as.factor(s1$AGEC)
s1$HHSIZE <- as.numeric(s1$HHSIZE)
s1$TRVLCMIN <- as.numeric(s1$TRVLCMIN)
s1$TRPTRANS <- as.numeric(s1$TRPTRANS)

s1 <- s1[,!"R_AGE"]
s1 <- setcolorder(s1, c("AGEC", "R_SEX", "EDUC", "HHSIZE",  "URBRUR", "HHFAMINC", "TRIPPURP", "GASPRICE", "TRVLCMIN", "TRPTRANS"))
s1 <- na.omit(s1)

#Y <- s1[,10]
#X <- s1[,1:9]

#vim<- varimpact(Y = s1[,10], data = s1[,1:9])

#s$TRIPPURP <- as.numeric(as.character(s$TRIPPURP))  
x <- copula.estimate(s1)
j = 1:100
corr_mean <- 0

for (i in j)
{
  m <- x$C.psamp[seq(i, length(x$C.psamp), 100)]
  corr_mean[i] <- mean(m)
}

V <- colnames(s1)
q<-matrix(unlist(corr_mean), ncol = length(s1), byrow = TRUE)

rownames(q) <- V
colnames(q) <- V

constra <- diag(10) * 0
rownames(constra) <- colnames(s1)
colnames(constra) <- colnames(s1) 
constra[c(2,11)] = 1

#suffStat <- list(C = cor(s), n = nrow(s))
suffStat <- list(C = q, n = nrow(s1))


pc.s1 <- pc(suffStat, indepTest = gaussCItest,
              labels = V, alpha = 0.01, fixedGaps = constra, skel.method = "stable")

 

fci.s1 <- fci(suffStat, indepTest=gaussCItest, alpha = 0.01, labels = V, skel.method = "stable", maj.rule = TRUE)

if (require(Rgraphviz)) {
  par(mfrow = c(1,2))
  
  plot(fci.s1, main = "FIC")
  plot(pc.s1, main = "Copula_PC")
}


#data(gmD)
#V <- colnames(gmD$x)
#suffStat <- list(dm = gmD$x, nlev = c(3,2,3,4,2), adaptDF = FALSE)

#pc.D <- pc(suffStat, indepTest = disCItest, alpha = 0.01, labels = V, verbose = TRUE)



