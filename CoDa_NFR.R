source("auxiliary.R")
source("density_norm.R")

##############################
# Compositional data analysis
##############################

# dat: n by p matrix

CoDa_NFR <- function(dat)
{
    n_year = nrow(dat)
    size_grid = n_age = ncol(dat)
    
    # standardize life table death to sum to 1
    
    dat_center = sweep(dat, 1, apply(dat, 1, sum), "/")
    
    # geometric mean
    
    alpha_x = vector("numeric", n_age)
    for(ik in 1:n_age)
    {
        alpha_x[ik] = geometric.mean(dat_center[,ik])
    }
    
    # standardization (closure operation)
    
    f_x_t = matrix(NA, n_year, n_age)
    for(ik in 1:n_year)
    {
            f_x_t[ik,] = dat[ik,]/alpha_x
    }
    
    # geometric mean and log-ratio transformation
    
    h_x_t = matrix(NA, n_year, n_age)
    for(ik in 1:n_year)
    {
        h_x_t[ik,] = log(f_x_t[ik,])
    }
  
    # NFR estimation and forecasting
    
    fore_val = as.numeric(ffunopare.knn.gcv(RESPONSES = h_x_t[2:n_year,], 
                                            CURVES = h_x_t[1:(n_year-1),],
                                            PRED = h_x_t[n_year,], semimetric = "L2")$Predicted.values)
    
    # Inverse clr transformation
    
    f_x_t_star_fore = exp(fore_val)/sum(exp(fore_val))
    d_x_t_star_fore = (f_x_t_star_fore * alpha_x)/sum((f_x_t_star_fore * alpha_x))
    return(d_x_t_star_fore)
}

##########################################
# compute errors of NFR density forecasts
##########################################

# dat: n by p data matrix
# den_fore: density forecast
# horizon: forecast horizon
# criterion: KLD or JSD

CoDa_NFR_err <- function(dat, den_fore, horizon, criterion)
{
    n = nrow(dat)
    true_dat = replace(dat, which(dat == 10^-5), NA)
  
    err = vector("numeric", (31 - horizon))
    for(ik in 1:(31 - horizon))
    {
        data_c = cbind(True = true_dat[(n - 31 + horizon + ik),], forecast = as.numeric(den_fore[,ik]))
        data_c = data_c[!is.na(data_c[,1]),]
        if(criterion == "density_norm")
        {
          err[ik] = density_norm(d1 = data_c[,2], d2 = data_c[,1], time_vector = 1:length(data_c[,2]))
        }
        else if(criterion == "KL_div")
        {
          err[ik] = mean(as.numeric(KLdiv(data_c, eps = 1e-16))[2:3])
        }
        else if(criterion == "JS_div_simple")
        {
          M = rowMeans(data_c)
          P_M = cbind(data_c[,1], M)
          E_M = cbind(as.numeric(data_c[,2]), M)
          colnames(E_M) = colnames(P_M) = c("True", "M")
          err[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
          rm(M); rm(P_M); rm(E_M)
        }
        else if(criterion == "JS_div_geo")
        {
          M = apply(data_c, 1, geometric.mean)
          P_M = cbind(data_c[,1], M)
          E_M = cbind(as.numeric(data_c[,2]), M)
          colnames(E_M) = colnames(P_M) = c("True", "M")
          err[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
          rm(M); rm(P_M); rm(E_M)
        }
        else
        {
          warning("Please specify a criterion from the list.")
        }
        rm(data_c)
    }
    rm(den_fore)
    return(mean(err, na.rm = TRUE))
}
