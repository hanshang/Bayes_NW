source("semimetric.R")
source("CoDa_NFR.R")
source("CoDa_recon.R")
source("LQDT_FPCA.R")

############
# read data
############

# set up a working directory

dir = "~/Dropbox/Todos/Nonparametric functional regression for a time series of density/"
setwd(paste(dir, "/data", sep=""))

# read life-table death counts

female_USA_qx = t(matrix(read.table("USA_female_lt.txt", skip = 2, header = TRUE)$qx, 111, 91))
male_USA_qx   = t(matrix(read.table("USA_male_lt.txt", skip = 2, header = TRUE)$qx, 111, 91))

radix = 10^5

female_USA_pop = male_USA_pop = matrix(NA, nrow(female_USA_qx), ncol(female_USA_qx))
for(ij in 1:nrow(female_USA_qx))
{
    start_pop_female = start_pop_male = radix
    for(ik in 1:ncol(female_USA_qx))
    {
        female_USA_pop[ij,ik] = female_USA_qx[ij,ik] * start_pop_female
        start_pop_female = start_pop_female - female_USA_pop[ij,ik]
        
        male_USA_pop[ij,ik] = male_USA_qx[ij,ik] * start_pop_male
        start_pop_male = start_pop_male - male_USA_pop[ij,ik]
    }
}
rownames(female_USA_pop) = rownames(male_USA_pop) = 1933:2023
colnames(female_USA_pop) = colnames(male_USA_pop) = 0:110

###############
# save figures
###############

savepdf("Fig_12a", width = 12, height = 10, toplines = 0.8)
plot(fts(0:110, t(female_USA_pop), xname = "Age", yname = "Life-table death counts"), 
     main = "USA: female data (1933-2023)", ylim = c(0, 7500))
dev.off()

savepdf("Fig_12b", width = 12, height = 10, toplines = 0.8)
plot(fts(0:110, t(male_USA_pop), xname = "Age", yname = "Life-table death counts"), 
     main = "USA: male data (1933-2023)", ylim = c(0, 7500))
dev.off()

######################
# auxiliary functions
######################

ffunopare.knn.gcv <- function(RESPONSES, CURVES, PRED, ..., Knearest=NULL, kind.of.kernel = "quadratic", 
                              semimetric = "deriv")
{
    if(is.vector(RESPONSES)) RESPONSES <- as.matrix(RESPONSES)
    if(is.vector(PRED)) PRED <- as.matrix(t(PRED))
    testfordim <- sum(dim(CURVES)==dim(PRED))==2
    twodatasets <- T
    if(testfordim) twodatasets <- sum(CURVES==PRED)!=prod(dim(CURVES))
    sm <- get(paste("semimetric.", semimetric, sep = ""))
    if(semimetric == "mplsr") stop("semimetric option mlpsr not allowed!")
    SEMIMETRIC1 <- sm(CURVES, CURVES, ...)
    kernel <- get(kind.of.kernel)
    n1 <- ncol(SEMIMETRIC1)
    if(is.null(Knearest))
    {
      step <- ceiling(n1/100)
      if(step == 0) step <- 1
      Knearest <- seq(from = 5, to = n1 %/% 2, by = step)	
      # the vector Knearest contains the sequence of the 
      # k-nearest neighbours used for computing the optimal bandwidth
    }
    kmax <- max(Knearest)	
    p <- ncol(CURVES)
    RESPONSES.estimated <- matrix(0, nrow = n1, ncol = p)
    Bandwidth.opt <- 0
    HAT.RESP <- matrix(0, nrow = n1, ncol = length(Knearest))
    BANDWIDTH <- matrix(0, nrow = n1, ncol = kmax)
    lKnearest <- length(Knearest)
    HAT.RESP <- array(0,c(n1,lKnearest,p))
    for(i in 1:n1) 
    {
      Norm.diff <- SEMIMETRIC1[, i]	
      Norm.order <- order(Norm.diff)
      zz <- sort(Norm.diff)[2:(kmax + 2)]
      BANDWIDTH[i,  ] <- 0.5 * (zz[-1] + zz[ - (kmax + 1)])
      z <- zz[ - (kmax + 1)]
      ZMAT <- matrix(rep(z, kmax), nrow = kmax, byrow = T)
      UMAT <- ZMAT/BANDWIDTH[i,  ]
      KNUM <- kernel(UMAT)
      KNUM[col(KNUM) > row(KNUM)] <- 0
      Kdenom <- apply(KNUM[Knearest,  ], 1, sum)
      WEIGHTS <- KNUM[Knearest,  ]/Kdenom
      Ind.curves <- Norm.order[2:(kmax + 1)]
      HAT.RESP[i,,] <- WEIGHTS %*% RESPONSES[Ind.curves,]
    }
    CRITARR <- array(0,c(n1,p,lKnearest))
    for(i in 1:n1)
    {
      CRITARR[i,,] <- (t(HAT.RESP[i,,]) - RESPONSES[i,])^2
    }
    Criterium <- apply(CRITARR, 3, sum)
    index.opt <- order(Criterium)[1]
    RESPONSES.estimated <- HAT.RESP[, index.opt,]
    knearest.opt <- Knearest[index.opt]
    Bandwidth.opt <- BANDWIDTH[, knearest.opt]
    Cv.estimated <- sum((RESPONSES.estimated - RESPONSES)^2)/(n1*p)
    if(twodatasets) 
    {
      SEMIMETRIC2 <- sm(CURVES, PRED, ...)
      Bandwidth2 <- 0
      n2 <- ncol(SEMIMETRIC2)
      for(k in 1:n2) {
        Sm2k <- SEMIMETRIC2[, k]
        Bandwidth2[k] <- sum(sort(Sm2k)[knearest.opt:(knearest.opt+1)])*0.5
      }
      KERNEL <- kernel(t(t(SEMIMETRIC2)/Bandwidth2))
      KERNEL[KERNEL < 0] <- 0
      KERNEL[KERNEL > 1] <- 0
      Denom <- apply(KERNEL, 2, sum)
      NUM <- t(KERNEL) %*% RESPONSES
      RESPONSES.predicted <- NUM/Denom
      return(list(Estimated.values = RESPONSES.estimated, 
                  Predicted.values = RESPONSES.predicted, Bandwidths = 
                    Bandwidth.opt, knearest.opt = knearest.opt, Cv = 
                    Cv.estimated))
    }
    else 
    {
      return(list(Estimated.values = RESPONSES.estimated, Bandwidths
                  = Bandwidth.opt, knearest.opt = knearest.opt, Cv = 
                    Cv.estimated))
    }
}

ffunopare.knn <- function(RESPONSES, CURVES, PRED, neighbour,..., kind.of.kernel = "quadratic", 
                          semimetric = "deriv")
{
    if(is.vector(RESPONSES)) RESPONSES <- as.matrix(RESPONSES)
    if(is.vector(PRED)) PRED <- as.matrix(t(PRED))
    testfordim <- sum(dim(CURVES)==dim(PRED))==2
    twodatasets <- T
    if(testfordim) twodatasets <- sum(CURVES==PRED)!=prod(dim(CURVES))
    sm <- get(paste("semimetric.", semimetric, sep = ""))
    if(semimetric == "mplsr") stop("semimetric option mlpsr not allowed!")
    SEMIMETRIC1 <- sm(CURVES, CURVES, ...)
    kernel <- get(kind.of.kernel)
    p1 <- ncol(SEMIMETRIC1)
    n1 <- ncol(SEMIMETRIC1)
    if(neighbour >= n1)
      stop(paste("try a smaller number of neighbour \n than ", neighbour))
    bandwidth.knn1 <- 0
    for(j in 1:p1) {
      Sem <- SEMIMETRIC1[, j]
      knn.to.band <- Sem[order(Sem)[neighbour:(neighbour + 1)]]
      bandwidth.knn1[j] <- 0.5 * sum(knn.to.band)
    }
    KERNEL1 <- kernel(t(t(SEMIMETRIC1)/bandwidth.knn1))
    KERNEL1[KERNEL1 < 0] <- 0
    KERNEL1[KERNEL1 > 1] <- 0
    diag(KERNEL1) <- 0
    Denom1 <- apply(KERNEL1, 2, sum)
    NUM1 <- t(KERNEL1) %*% RESPONSES
    RESPONSES.estimated <- NUM1/Denom1
    Cv.estimated <- sum((RESPONSES.estimated - RESPONSES)^2)/(n1*p1)
    if(twodatasets)
    {
        SEMIMETRIC2 <- sm(CURVES, PRED, ...)
        Bandwidth2 <- 0
        p2 <- ncol(SEMIMETRIC2)
        bandwidth.knn2 <- 0
        for(j in 1:p2)
        {
            Sem <- SEMIMETRIC2[, j]
            knn.to.band <- Sem[order(Sem)[neighbour:(neighbour + 1)]]
            bandwidth.knn2[j] <- 0.5 * sum(knn.to.band)
        }
        KERNEL2 <- kernel(t(t(SEMIMETRIC2)/bandwidth.knn2))
        KERNEL2[KERNEL2 < 0] <- 0
        KERNEL2[KERNEL2 > 1] <- 0
        Denom2 <- apply(KERNEL2, 2, sum)
        NUM2 <- t(KERNEL2) %*% RESPONSES
        RESPONSES.predicted <- NUM2/Denom2
        return(list(Estimated.values = RESPONSES.estimated, Predicted.values = RESPONSES.predicted, Cv = Cv.estimated))
    }
    else
    {
        return(list(Estimated.values = RESPONSES.estimated, Cv = Cv.estimated))
    }
}

# Quadratic kernel function

quadratic <- function(u)
{
    return(1 - (u)^2)
}

##################
# females & males
##################

# training_data: 1933-2000
# testing_data: 2001-2023

# Bayes NW

n_test = 23
n_age = 111
FoF_regre_list_female_mat = FoF_regre_list_male_mat = matrix(NA, n_age, n_test)
for(iwk in 1:n_test)
{
    FoF_regre_list_female_mat[,iwk] = CoDa_NFR(dat = female_USA_pop[1:(67+iwk),]/radix)
    FoF_regre_list_male_mat[,iwk]   = CoDa_NFR(dat = male_USA_pop[1:(67+iwk),]/radix)
    print(iwk); rm(iwk)
}
year = 1933:2023
n_year = length(year)
colnames(FoF_regre_list_female_mat) = colnames(FoF_regre_list_male_mat) = year[69:91]

# CoDa

CoDa_female_val = CoDa_male_val = 
CoDa_female_val_normalization = CoDa_male_val_normalization = 
CoDa_female_val_arima = CoDa_male_val_arima = 
CoDa_female_val_arima_normalization = CoDa_male_val_arima_normalization = matrix(NA, n_age, n_test)
for(iwk in 1:n_test)
{
    # forecasting_method = "ETS"
  
    CoDa_female_val[,iwk] = CoDa_recon(dat = female_USA_pop[1:(67+iwk),], normalize = "TRUE", fore_method = "ETS",
                                       fh = 1, varprop = 0.95, constant = radix)$d_x_t_star_fore
    CoDa_male_val[,iwk]   = CoDa_recon(dat = male_USA_pop[1:(67+iwk),],   normalize = "TRUE", fore_method = "ETS",
                                       fh = 1, varprop = 0.95, constant = radix)$d_x_t_star_fore
    
    CoDa_female_val_normalization[,iwk] = CoDa_recon(dat = female_USA_pop[1:(67+iwk),], normalize = "FALSE", fore_method = "ETS",
                                                     fh = 1, varprop = 0.95, constant = radix)$d_x_t_star_fore
    CoDa_male_val_normalization[,iwk]   = CoDa_recon(dat = male_USA_pop[1:(67+iwk),],   normalize = "FALSE", fore_method = "ETS",
                                                     fh = 1, varprop = 0.95, constant = radix)$d_x_t_star_fore

    # forecasting_method = "ARIMA"
    
    CoDa_female_val_arima[,iwk] = CoDa_recon(dat = female_USA_pop[1:(67+iwk),], normalize = "TRUE", fore_method = "ARIMA",
                                             fh = 1, varprop = 0.95, constant = radix)$d_x_t_star_fore
    CoDa_male_val_arima[,iwk]   = CoDa_recon(dat = male_USA_pop[1:(67+iwk),],   normalize = "TRUE", fore_method = "ARIMA",
                                             fh = 1, varprop = 0.95, constant = radix)$d_x_t_star_fore
    
    CoDa_female_val_arima_normalization[,iwk] = CoDa_recon(dat = female_USA_pop[1:(67+iwk),], normalize = "FALSE", fore_method = "ARIMA",
                                                           fh = 1, varprop = 0.95, constant = radix)$d_x_t_star_fore
    CoDa_male_val_arima_normalization[,iwk]   = CoDa_recon(dat = male_USA_pop[1:(67+iwk),],   normalize = "FALSE", fore_method = "ARIMA",
                                                           fh = 1, varprop = 0.95, constant = radix)$d_x_t_star_fore
    print(iwk); rm(iwk)
}

## LQD

# some years produce NA values, as it depends on the variable adj_const

LQD_female_val = matrix(NA, n_age, n_test)
for(iwk in 1:n_test)
{
    LQD_female_val[,iwk] = lqd_fun_HS(data = female_USA_pop[1:(67+iwk),], fh = 1, var_prop = 0.95, fmethod = "ets", 
                                      adj_const = 60)$female_dens_forecast
    print(iwk); rm(iwk)
}

LQD_male_val = matrix(NA, n_age, n_test)
for(iwk in 1:n_test)
{
    LQD_male_val[,iwk]   = lqd_fun_HS(data = male_USA_pop[1:(67+iwk),],   fh = 1, var_prop = 0.95, fmethod = "ets", 
                                      adj_const = 60)$female_dens_forecast
    print(iwk); rm(iwk)
}

# Horta-Ziegelman

HZ_fore_female = HZ_fore_male = matrix(NA, n_age, n_test)
for(iwk in 1:n_test)
{
    ncomp = which(cumsum(ftsm(fts(0:(n_age-1), t(female_USA_pop[1:(67+iwk),])), order = 10)$varprop) >= 0.95)[1]
    HZ_fore = forecast(ftsm(fts(0:(n_age-1), t(female_USA_pop[1:(67+iwk),])), order = ncomp), h = 1)$mean$y
    HZ_fore_female[,iwk] = HZ_fore/sum(HZ_fore)*radix
    rm(ncomp); rm(HZ_fore)
    
    ncomp = which(cumsum(ftsm(fts(0:(n_age-1), t(male_USA_pop[1:(67+iwk),])), order = 10)$varprop) >= 0.95)[1]
    HZ_fore = forecast(ftsm(fts(0:(n_age-1), t(male_USA_pop[1:(67+iwk),])), order = ncomp), h = 1)$mean$y
    HZ_fore_male[,iwk] = HZ_fore/sum(HZ_fore)*radix
    rm(ncomp); rm(HZ_fore)
    print(iwk); rm(iwk)
}

###############
# save figures
###############

savepdf("Fig_13a", width = 12, height = 10, toplines = 0.8)
plot(fts(0:(n_age-1), t(female_USA_pop[69:n_year,])), xlab = "", ylab = "Life-table death count", 
     main = "USA: female data (2001-2023)", ylim = c(5, 4035))
dev.off()

savepdf("Fig_13b", width = 12, height = 10, toplines = 1.5)
plot(fts(0:(n_age-1), FoF_regre_list_female_mat*radix), xlab = "", ylab = "", 
     main = "Bayes NW estimator", ylim = c(5, 4035))
dev.off()

savepdf("Fig_13c", width = 12, height = 10, toplines = 1.5)
plot(fts(0:(n_age-1), CoDa_female_val), xlab = "", ylab = "", 
     main = "CoDa method", ylim = c(5, 4035))
dev.off()

savepdf("Fig_13d", width = 12, height = 10, toplines = 1.5)
plot(fts(0:(n_age-1), CoDa_female_val), type = "n", axes = FALSE, xlab = "", ylab = "")
dev.off()

savepdf("Fig_13e", width = 12, height = 10, toplines = 1.5)
plot(fts(0:(n_age-1), LQD_female_val), xlab = "", ylab = "", 
     main = "LQDT method", ylim = c(5, 4035))
dev.off()

savepdf("Fig_13f", width = 12, height = 10, toplines = 1.5)
plot(fts(0:(n_age-1), HZ_fore_female), xlab = "", ylab = "", 
     main = "HZ method", ylim = c(5, 4035))
dev.off()


savepdf("Fig_13g", width = 12, height = 10, toplines = 0.8)
plot(fts(0:(n_age-1), t(male_USA_pop[67:n_year,])), xlab = "Age", ylab = "Life-table death count", 
     main = "USA: male data (2001-2023)", ylim = c(0, 3475))
dev.off()

savepdf("Fig_13h", width = 12, height = 10, toplines = 1.5)
plot(fts(0:(n_age-1), FoF_regre_list_male_mat*radix), xlab = "", ylab = "", 
     main = "Bayes NW estimator", ylim = c(0, 3475))
dev.off()

savepdf("Fig_13i", width = 12, height = 10, toplines = 1.5)
plot(fts(0:(n_age-1), CoDa_male_val), xlab = "", ylab = "", 
     main = "CoDa method", ylim = c(0, 3475))
dev.off()

savepdf("Fig_13j", width = 12, height = 10, toplines = 1.5)
plot(fts(0:(n_age-1), LQD_male_val), xlab = "Age", ylab = "", 
     main = "LQDT method", ylim = c(0, 3475))
dev.off()

savepdf("Fig_13k", width = 12, height = 10, toplines = 1.5)
plot(fts(0:(n_age-1), HZ_fore_male), xlab = "Age", ylab = "", 
     main = "HZ method", ylim = c(0, 3475))
dev.off()


###############
# density norm
###############

# Bayes NW estimator

KLD_female_LT = KLD_male_LT = vector("numeric", n_test)
for(ik in 1:n_test)
{
    data_c = cbind(True = female_USA_pop[(68+ik),],
                   forecast = FoF_regre_list_female_mat[,ik])
    KLD_female_LT[ik] = mean(as.numeric(KLdiv(data_c, eps = 1e-16))[2:3])
    rm(data_c)
    
    data_c = cbind(True = male_USA_pop[(68+ik),], 
                   forecast = FoF_regre_list_male_mat[,ik])
    KLD_male_LT[ik] = mean(as.numeric(KLdiv(data_c, eps = 1e-16))[2:3])
    rm(data_c)
    print(ik); rm(ik)
}

# CoDa

KLD_CoDa_female_LT = KLD_CoDa_male_LT = KLD_CoDa_female_LT_normalization = KLD_CoDa_male_LT_normalization = 
KLD_CoDa_female_LT_arima = KLD_CoDa_male_LT_arima = KLD_CoDa_female_LT_arima_normalization = KLD_CoDa_male_LT_arima_normalization = vector("numeric", n_test)
for(ik in 1:n_test)
{
    # ETS
  
    data_c = cbind(True = female_USA_pop[(68+ik),],
                   forecast = CoDa_female_val[,ik])
    KLD_CoDa_female_LT[ik] = mean(as.numeric(KLdiv(data_c, eps = 1e-16))[2:3])
    rm(data_c)
     
    data_c = cbind(True = male_USA_pop[(68+ik),],
                   forecast = CoDa_male_val[,ik])
    KLD_CoDa_male_LT[ik] = mean(as.numeric(KLdiv(data_c, eps = 1e-16))[2:3])
    rm(data_c)
    
    # ETS (normalization)
    
    data_c = cbind(True = female_USA_pop[(68+ik),],
                   forecast = CoDa_female_val_normalization[,ik])
    KLD_CoDa_female_LT_normalization[ik] = mean(as.numeric(KLdiv(data_c, eps = 1e-16))[2:3])
    rm(data_c)
    
    data_c = cbind(True = male_USA_pop[(68+ik),],
                   forecast = CoDa_male_val_normalization[,ik])
    KLD_CoDa_male_LT_normalization[ik] = mean(as.numeric(KLdiv(data_c, eps = 1e-16))[2:3])
    rm(data_c)

    # ARIMA

    data_c = cbind(True = female_USA_pop[(68+ik),],
                   forecast = CoDa_female_val_arima[,ik])
    KLD_CoDa_female_LT_arima[ik] = mean(as.numeric(KLdiv(data_c, eps = 1e-16))[2:3])
    rm(data_c)
      
    data_c = cbind(True = male_USA_pop[(68+ik),],
                   forecast = CoDa_male_val_arima[,ik])
    KLD_CoDa_male_LT_arima[ik] = mean(as.numeric(KLdiv(data_c, eps = 1e-16))[2:3])
    rm(data_c)
    
    # ARIMA (normalization)
    
    data_c = cbind(True = female_USA_pop[(68+ik),],
                   forecast = CoDa_female_val_arima_normalization[,ik])
    KLD_CoDa_female_LT_arima_normalization[ik] = mean(as.numeric(KLdiv(data_c, eps = 1e-16))[2:3])
    rm(data_c)  
    
    data_c = cbind(True = male_USA_pop[(68+ik),],
                   forecast = CoDa_male_val_arima_normalization[,ik])
    KLD_CoDa_male_LT_arima_normalization[ik] = mean(as.numeric(KLdiv(data_c, eps = 1e-16))[2:3])
    rm(data_c); print(ik); rm(ik)    
}

# LQT

KLD_LQD_female_index = which(!apply(LQD_female_val == 0, 2, any))
KLD_LQD_female_LT = vector("numeric", length(KLD_LQD_female_index))
for(ik in KLD_LQD_female_index)
{
    data_c = cbind(True = female_USA_pop[(68+ik),],forecast = LQD_female_val[,ik])
    KLD_LQD_female_LT[ik] = mean(as.numeric(KLdiv(data_c, eps = 1e-16))[2:3])
    rm(data_c); print(ik); rm(ik)
}

KLD_LQD_male_index = which(!apply(LQD_male_val == 0, 2, any))
KLD_LQD_male_LT = vector("numeric", length(KLD_LQD_male_index))
for(ik in KLD_LQD_male_index)
{
    data_c = cbind(True = male_USA_pop[(68+ik),], forecast = LQD_male_val[,ik])
    KLD_LQD_male_LT[ik] = mean(as.numeric(KLdiv(data_c, eps = 1e-16))[2:3])
    rm(data_c); print(ik); rm(ik)
}
  
# Horta-Ziegelman

KLD_HZ_female_LT = KLD_HZ_male_LT = vector("numeric", n_test)
for(ik in 1:n_test)
{
    data_c = cbind(True = female_USA_pop[(68+ik),], forecast = HZ_fore_female[,ik])
    KLD_HZ_female_LT[ik] = mean(as.numeric(KLdiv(data_c, eps = 1e-16))[2:3])
    rm(data_c)
    
    data_c = cbind(True = male_USA_pop[(68+ik),], forecast = HZ_fore_male[,ik])
    KLD_HZ_male_LT[ik] = mean(as.numeric(KLdiv(data_c, eps = 1e-16))[2:3])
    rm(data_c); print(ik); rm(ik)    
}

# RW

KLD_RW_female_LT = KLD_RW_male_LT = vector("numeric", n_test)
for(ik in 1:n_test)
{
    data_c = cbind(True = female_USA_pop[(68+ik),], forecast = female_USA_pop[(67+ik),])
    KLD_RW_female_LT[ik] = mean(as.numeric(KLdiv(data_c, eps = 1e-16))[2:3])
    rm(data_c)
    
    data_c = cbind(True = male_USA_pop[(68+ik),], forecast = male_USA_pop[(67+ik),])
    KLD_RW_male_LT[ik] = mean(as.numeric(KLdiv(data_c, eps = 1e-16))[2:3])
    rm(data_c); print(ik); rm(ik)    
}

# RWD

RWD_female_LT = RWD_male_LT = matrix(NA, n_test, n_age)
for(ik in 1:n_test)
{
    for(ij in 1:n_age)
    {
        RWD_female_LT[ik,ij] = rwf(female_USA_pop[1:(67+ik),ij], drift = TRUE, h = 1)$mean
        RWD_male_LT[ik,ij]   = rwf(male_USA_pop[1:(67+ik),ij], drift = TRUE, h = 1)$mean
    }
}

KLD_RWD_female_LT = KLD_RWD_male_LT = vector("numeric", n_test)
for(ik in 1:n_test)
{
    data_c = cbind(True = female_USA_pop[(68+ik),], forecast = RWD_female_LT[ik,])
    KLD_RWD_female_LT[ik] = mean(as.numeric(KLdiv(data_c, eps = 1e-16))[2:3])
    rm(data_c)
    
    data_c = cbind(True = male_USA_pop[(68+ik),], forecast = RWD_male_LT[ik,])
    KLD_RWD_male_LT[ik] = mean(as.numeric(KLdiv(data_c, eps = 1e-16))[2:3])
    rm(data_c); print(ik); rm(ik)
}

############
# MSPE ratio
############

KLD_table_female_LT = list(Bayes_NW = KLD_female_LT,
                           CoDa = KLD_CoDa_female_LT,
                           LQD = KLD_LQD_female_LT[1:21],
                           HZ = KLD_HZ_female_LT)

KLD_table_male_LT = list(Bayes_NW = KLD_male_LT,
                          CoDa = KLD_CoDa_male_LT,
                          LQD = KLD_LQD_male_LT[1:16],
                          HZ = KLD_HZ_male_LT)

# save figures

savepdf("KLD_table_female_LT", width = 12, height = 10, toplines = 0.8)
boxplot(KLD_table_female_LT, notch = FALSE, ylim = c(0, 0.015), ylab = "KLD", main = "Female data")
dev.off()

savepdf("KLD_table_male_LT", width = 12, height = 10, toplines = 0.8)
boxplot(KLD_table_male_LT, notch = FALSE, ylim = c(0, 0.015), main = "Male data")
dev.off()

