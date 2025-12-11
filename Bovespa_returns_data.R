source("Bovespa_CoDa_NFR.RData")

#########################
# Descriptive statistics
#########################

# 5-minute returns

data.list = Bovespa_data.list
N = length(data.list)
ret = lapply(1:N, function(t){m = length(data.list[[t]]$Fech);
sapply(1:(m-1), function(j) log(data.list[[t]]$Fech[j+1]) - log(data.list[[t]]$Fech[j]))})
h.hat_5m = sapply(1:N, function(t) 2.34*sd(ret[[t]])*(length(ret[[t]])^(-(1/5))))
m = 5001 # number of grid points (this is arbitrary)
a = min(sapply(1:N, function(ik) min(ret[[ik]])))
b = max(sapply(1:N, function(ik) max(ret[[ik]])))
u = seq(from = a, to = b, length = m)

# Interval length
du = u[2] - u[1]

Y = sapply(1:N, function(t) density(ret[[t]], bw = h.hat_5m[t], kernel = 'epanechnikov', from = a, to = b, n = m)$y)
# correcting to ensure integral Y_t du = 1
for(t in 1:N)
{
    Y[,t] = Y[,t]/(sum(Y[,t]) * du)
}

Y_selected = Y[,c(which(date=="09/02/2010"),
                  which(date=="29/10/2009"),
                  which(date=="30/10/2009"))]

# plot

matplot(u, Y_selected, type = "l", ylab = "Density", xlab = "Grid point")
legend("topright", c("09/02/2010", "29/10/2009", "30/10/2009"), col = 1:3, lty = 1:3, cex = 0.8)

# Bayes NW method
        
Bovespa_CoDa_NFR <- function(data.list, band_choice = c("Silverman", "DPI"), kernel = c("gaussian", "epanechnikov"))
{
    kernel = match.arg(kernel)
    # Sample size
    N = length(data.list)
    
    # 5-minute returns
    ret = lapply(1:N, function(t){m = length(data.list[[t]]$Fech); 
    sapply(1:(m-1), function(j) log(data.list[[t]]$Fech[j+1]) - log(data.list[[t]]$Fech[j]))})
    
    # Creating the 'rule of thumb' bandwidths h.hat, to be used to estimate the daily densities. h.hat is a n-vector with h.hat[t] being Silverman's rule of thumb bandwidth for day t
    if(band_choice == "Silverman")
    {
        if(kernel == "gaussian")
        {
            h.hat_5m = sapply(1:N, function(t) 1.06*sd(ret[[t]])*(length(ret[[t]])^(-(1/5))))
        }
        if(kernel == "epanechnikov")
        {
            h.hat_5m = sapply(1:N, function(t) 2.34*sd(ret[[t]])*(length(ret[[t]])^(-(1/5)))) 
        }        
    }
    if(band_choice == "DPI")
    {
        if(kernel == "gaussian")
        {
            h.hat_5m = sapply(1:N, function(t) dpik(ret[[t]], kernel = "normal"))
        }
        if(kernel == "epanechnikov")
        {
            h.hat_5m = sapply(1:N, function(t) dpik(ret[[t]], kernel = "epanech"))
        }
    }

    # 2. Discretization
    # Evaluation points
    m = 5001 # number of grid points (this is arbitrary)
    a = min(sapply(1:N, function(ik) min(ret[[ik]])))
    b = max(sapply(1:N, function(ik) max(ret[[ik]])))
    u = seq(from = a, to = b, length = m)
    
    # Interval length
    du = u[2] - u[1]
    
    # Creating an (m x n) matrix which represents the observed densities. Y[j,t] is the density at date t evaluated at u[j]
    if(kernel == "gaussian")
    {
        Y = sapply(1:N, function(t) density(ret[[t]], bw = h.hat_5m[t], kernel = 'gaussian', from = a, to = b, n = m)$y)
    }
    if(kernel == "epanechnikov")
    {
        Y = sapply(1:N, function(t) density(ret[[t]], bw = h.hat_5m[t], kernel = 'epanechnikov', from = a, to = b, n = m)$y)
    }
    
    # correcting to ensure integral Y_t du = 1
    for(t in 1:N)
    {
        Y[,t] = Y[,t]/(sum(Y[,t])*du)
    }
    
    # adjustment
    
    radix = colSums(Y)[1]
    
    CoDa_NFR_den = CoDa_NFR(t(Y + 0.001))
    Yhat.fix_den = (CoDa_NFR_den$d_x_t_star_fore * radix) - 0.001
    
    Yhat.fix_den[Yhat.fix_den < 0] = 0
    Yhat.fix_den = Yhat.fix_den/(sum(Yhat.fix_den)*du)
    return(Yhat.fix_den)
}

Bovespa_ret_fun_CoDa_NFR <- function(ik, choice_band, choice_kernel)
{
    data.list_new = list()
    for(iw in 1:(202+ik))
    {
        data.list_new[[iw]] = Bovespa_data.list[[iw]]
    }
    dum = Bovespa_CoDa_NFR(data.list = data.list_new, band_choice = choice_band, kernel = choice_kernel)
    return(dum)
}

################################
# run parallel for testing data
################################

require(doMC)
registerDoMC(10)
result_Silverman_epan = foreach(iwk = 1:n_test) %dopar% Bovespa_ret_fun_CoDa_NFR(ik = iwk,
                                                        choice_band = "Silverman", choice_kernel = "epanechnikov")

Bovespa_CoDa_NFR_Silverman_epan = matrix(NA, 5001, n_test)
for(ik in 1:n_test)
{
    Bovespa_CoDa_NFR_Silverman_epan[,ik] = result_Silverman_epan[[ik]]
    print(ik); rm(ik)
}  

##############################
# Kullback-Leibler divergence
##############################

# kernel = "epan"

Bovespa_KLdiv_CoDa_NFR_Silverman_epan = matrix(NA, n_test, 2)
for(ik in 1:n_test)
{
    ## Silverman
    
    Bovespa_CoDa_NFR_compar_epan = cbind(Bovespa_test_density_epan[,ik], Bovespa_CoDa_NFR_Silverman_epan[,ik])
    colnames(Bovespa_CoDa_NFR_compar_epan) = c("True", "Estimate")
    Bovespa_KLdiv_CoDa_NFR_Silverman_epan[ik,] = as.numeric(KLdiv(Bovespa_CoDa_NFR_compar_epan, eps=1e-16))[2:3]
    rm(Bovespa_CoDa_NFR_compar_epan)
    print(ik); rm(ik)
}

Bovespa_KLdiv_CoDa_NFR_Silverman_epan_summary = round(sum(colMeans(Bovespa_KLdiv_CoDa_NFR_Silverman_epan)), 4)

# random walk

Bovespa_KLdiv_rw_epan = matrix(NA, n_test, 2)
for(ik in 2:n_test)
{
    ## Silverman epan
    
    Bovespa_rw_compar = cbind(Bovespa_test_density_epan[,ik], Bovespa_test_density_epan[,(ik-1)])
    colnames(Bovespa_rw_compar) = c("True", "Estimate")
    Bovespa_KLdiv_rw_epan[ik,] = as.numeric(KLdiv(Bovespa_rw_compar, eps=1e-16))[2:3]
    rm(Bovespa_rw_compar)
    print(ik); rm(ik)
}
                   
round(colMeans(Bovespa_KLdiv_rw[2:n_test,]), 4) # 0.1912 0.1977
Bovespa_KLdiv_rw_summary_epan = round(sum(colMeans(Bovespa_KLdiv_rw_epan[2:n_test,])), 4)

#####################
# comparison of KLDs
#####################

# Horta-Ziegelman

Bovespa_compar_summary <- cbind(rowMeans(Bovespa_KLdiv_Horta_Ziegelman_epan),

# LQT

rowMeans(Bovespa_KLdiv_lqd_arima_epan),

# CoDa

rowMeans(Bovespa_KLdiv_CoDa_epan),

# CoDa no normalization

rowMeans(Bovespa_KLdiv_CoDa_epan_norm),

# skew-t distribution

rowMeans(Bovespa_KLdiv_skew_t_epan),

# CoDa NFR

rowMeans(Bovespa_KLdiv_CoDa_NFR_Silverman_epan),

# RW

rowMeans(Bovespa_KLdiv_rw_epan))

rownames(Bovespa_compar_summary) = 1:102
colnames(Bovespa_compar_summary) = c("HZ", "LQDT", "CoDa normalised", "CoDa", "Skewed-t", "Bayes NW", "random walk")

boxplot(Bovespa_compar_summary, ylab = "KLD", xlab = "Method", notch = FALSE, outline = FALSE)

                   
