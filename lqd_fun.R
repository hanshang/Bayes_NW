# LQD to Density, Density to LQD

lqd2dens = ftsa:::lqd2dens
dens2lqd = ftsa:::dens2lqd

# data: n by p data matrix
# gridpoints: dense grid points
# h_scale: scaling parameter
# M: number of grid points
#

LQDT_FPCA <- function(data, gridpoints, h_scale = 1, M = 3001, m = 5001,
          lag_maximum = 4, no_boot = 1000, alpha_val = 0.1, p = 5,
          band_choice = c("Silverman", "DPI"), kernel = c("gaussian", "epanechnikov"),
          forecasting_method = c("uni", "multi"),
          varprop = 0.85, fmethod, VAR_type)
{
    N = nrow(data)
    if(getmode(trunc(diff(apply(data, 1, sum))) == 0))
    {
        Y = t(data)
        u = gridpoints
        du = u[2] - u[1]
    }
    else
    {
        kernel = match.arg(kernel)
        forecasting_method = match.arg(forecasting_method)
        if(!exists("h_scale"))
            h_scale = 1
        if(band_choice == "Silverman")
        {
            if(kernel == "gaussian")
            {
                h.hat_5m = sapply(1:N, function(t) 1.06 * sd(data[t,]) * (length(data[t, ])^(-(1/5))))
            }
            if(kernel == "epanechnikov")
            {
                h.hat_5m = sapply(1:N, function(t) 2.34 * sd(data[t,]) * (length(data[t, ])^(-(1/5))))
            }
            h.hat_5m = h_scale * h.hat_5m
        }
        if(band_choice == "DPI")
        {
            if(kernel == "gaussian")
            {
                h.hat_5m = sapply(1:N, function(t) dpik(data[t,], kernel = "normal"))
            }
            if(kernel == "epanechnikov")
            {
                h.hat_5m = sapply(1:N, function(t) dpik(data[t,], kernel = "epanech"))
            }
            h.hat_5m = h_scale * h.hat_5m
        }
        u = seq(from = min(data), to = max(data), length = m)
        du = u[2] - u[1]
        if(kernel == "gaussian")
        {
            Y = sapply(1:N, function(t) density(data[t, ], bw = h.hat_5m[t], kernel = "gaussian",
                        from = min(data), to = max(data), n = m)$y)
        }
        if(kernel == "epanechnikov")
        {
            Y = sapply(1:N, function(t) density(data[t, ], bw = h.hat_5m[t], kernel = "epanechnikov",
                        from = min(data), to = max(data), n = m)$y)
        }
    }
    n = ncol(Y)
    N = length(u)
    any(apply(Y, 2, function(z) any(z < 0)))
    dens = sapply(1:n, function(i) Y[, i]/trapzRcpp(X = u, Y = Y[, i]))
    lqd = matrix(0, nrow = M, ncol = n)
    c = vector("numeric", n)
    t = seq(0, 1, length.out = M)
    t0 = u[which.min(abs(u))]
    for(i in 1:n)
    {
        tmp = dens2lqd(dens = dens[, i], dSup = u, lqdSup = t, t0 = t0, verbose = FALSE)
        lqd[, i] = tmp$lqd
        c[i] = tmp$c
    }
    cut = res2 = list()
    for(i in 1:n)
    {
        res2[[i]] = lqd2dens(lqd = lqd[, i], lqdSup = t, t0 = t0, c = c[i], useSplines = FALSE, verbose = FALSE)
    }
    dens2 = sapply(res2, function(r) approx(x = r$dSup, y = r$dens, xout = u, yleft = 0, yright = 0)$y)
    totalMass = apply(dens2, 2, function(d) trapzRcpp(X = u, Y = d))
    L2Diff = sapply(1:n, function(i) sqrt(trapzRcpp(X = u, Y = (dens[, i] - dens2[, i])^2)))
    unifDiff = sapply(1:n, function(i){
                        interior = which(dens2[, i] > 0)
                        max(abs(dens[interior, i] - dens2[interior, i]))})
    c_fore = as.numeric(inv.logit(forecast(auto.arima(logit(c)), h = 1)$mean))
    t_new = t[2:(M - 1)]
    lqd_new = lqd[2:(M - 1), ]
    if(forecasting_method == "uni")
    {
        ftsm_fitting = ftsm(y = fts(x = t_new, y = lqd_new), order = 10)
        ftsm_ncomp = head(which(cumsum((ftsm_fitting$varprop^2)/sum((ftsm_fitting$varprop^2))) >= varprop), 1)
        den_fore = forecast(object = ftsm(y = fts(x = t_new, y = lqd_new), order = ftsm_ncomp), h = 1, method = fmethod)$mean$y
    }
    if(forecasting_method == "multi")
    {
        dt = diff(t)[1]
        foo_out = super_fun(Y = lqd_new, lag_max = lag_maximum, B = no_boot, alpha = alpha_val, du = dt, p = p,
                            m = (M - 2), u = t_new, select_ncomp = "TRUE")
        Ybar_est = foo_out$Ybar
        psihat_est = foo_out$psihat
        etahat_est = matrix(foo_out$etahat, ncol = n)
        selected_d0 = foo_out$d0
        score_object = t(etahat_est)
        colnames(score_object) = 1:selected_d0
        if(selected_d0 == 1)
        {
            etahat_pred_val = forecast(auto.arima(as.numeric(score_object)), h = 1)$mean
        }
        else
        {
            VAR_mod = VARselect(score_object)
            etahat_pred = predict(VAR(y = score_object, p = min(VAR_mod$selection[3], 3), type = VAR_type), n.ahead = 1)
            etahat_pred_val = as.matrix(sapply(1:selected_d0, function(t) (etahat_pred$fcst[[t]])[1]))
        }
        den_fore = Ybar_est + psihat_est %*% etahat_pred_val
    }
    den_fore = as.matrix(c(8, den_fore, 8))
    res_fore = lqd2dens(lqd = den_fore, lqdSup = t, t0 = t0, c = c_fore, useSplines = FALSE, verbose = FALSE)
    dens_fore = approx(x = res_fore$dSup, y = res_fore$dens, xout = u, yleft = 0, yright = 0)$y
    return(list(L2Diff = L2Diff, unifDiff = unifDiff, density_reconstruct = dens2,
              density_original = dens, dens_fore = dens_fore, totalMass = range(totalMass),
              u = u))
}

lqd_fun_HS <- function(data, fh, var_prop = 0.85, fmethod, adj_const)
{
    const = rowSums(data)[1]
    n = nrow(data)
    u = 0:110
    if(is.null(adj_const))
    {
      adj_const = mean(u)
    }
    u = u - adj_const
    du = diff(u)[1]
    N = length(u)
    
    # Try Forward transformation
    
    M = 5501
    lqd = matrix(0, nrow = M, ncol = n)
    c = rep(0, n)
    t = seq(0, 1, length.out = M)
    t0 = u[which.min(abs(u))] # closest value to 0
    
    female_dens = sapply(1:n, function(i){data[i,]/trapzRcpp(X = u, Y = data[i,])})
    
    for(i in 1:n)
    {
      tmp = dens2lqd(dens = female_dens[,i], dSup = u, lqdSup = t,
                     t0 = t0, verbose = FALSE)
      lqd[,i] = tmp$lqd
      c[i] = tmp$c
    }
    colnames(lqd) = 1:ncol(lqd)
    
    # Try cutting off boundary points with large LQD values (effectively setting the density to zero rather than trying to compute it numerically)
    res2 = list()
    for(i in 1:n)
    {
      res2[[i]] = lqd2dens(lqd = lqd[,i], lqdSup = t, t0 = t0, c = c[i],
                           useSplines = FALSE, cut = c(0, 0), verbose = FALSE)
    }
    
    
    # Make ad hoc correction to support of densities near the boundaries
    
    dens2 = sapply(res2, function(r){
      tmp = r$dSup + adj_const
      ind = which(tmp < 1)
      tmp[ind] = tmp[max(ind)]*(tmp[ind] - tmp[min(ind)])/diff(range(tmp[ind]))
      
      ind = which(tmp > 109)
      tmp[ind] = tmp[min(ind)] + (110 - tmp[min(ind)])*(tmp[ind] - tmp[min(ind)])/diff(range(tmp[ind]))
      
      approx(x = tmp, y = r$dens, xout = u + adj_const)$y
      
    })
    female_dens_recon = sapply(1:n, function(i){dens2[,i]/trapzRcpp(X = u, Y = dens2[,i])})
    female_dens_recon_norm = sapply(1:n, function(i){female_dens_recon[,i]/sum(female_dens_recon[,i])})
    
    # Plot LQD Functions
    
    #savepdf("Rplot_AP", width = 12, height = 10, toplines = 0.8)
    # pdf(file = "Rplot_AP.pdf", width = 12, height = 10)
    # par(mfrow = c(1,3))
    # plot(fts(0:110, t(data)), xlab = "Age", ylab = "Lifetable death count")
    # plot(fts(t, lqd))
    # plot(fts(0:110,female_dens_recon*const), xlab = "Age", ylab = "Re-constructed lifetable death count")
    # dev.off()
    
    ######################
    # Forecasting density
    ######################
    
    # point forecast of c and lqd
    
    if(fmethod == "arima")
    {
      c_fore = inv.logit(forecast(auto.arima(logit(c)), h = fh)$mean)
      
      ftsm_fitting = ftsm(y = fts(x = t, y = lqd), order = 10)
      ftsm_ncomp = head(which(cumsum((ftsm_fitting$varprop^2)/sum((ftsm_fitting$varprop^2))) >= var_prop),1)
      den_fore = forecast(object = ftsm(y = fts(x = t, y = lqd), order = ftsm_ncomp), h = fh, method = "arima")$mean$y
    }
    if(fmethod == "ets")
    {
      c_fore = inv.logit(forecast(ets(logit(c)), h = fh)$mean)
      
      ftsm_fitting = ftsm(y = fts(x = t, y = lqd), order = 10)
      ftsm_ncomp = head(which(cumsum((ftsm_fitting$varprop^2)/sum((ftsm_fitting$varprop^2))) >= var_prop),1)
      den_fore = forecast(object = ftsm(y = fts(x = t, y = lqd), order = ftsm_ncomp), h = fh, method = "ets")$mean$y
    }
    
    res_fore = list()
    for(i in 1:fh)
    {
      res_fore[[i]] = lqd2dens(lqd = as.matrix(den_fore)[,i], lqdSup = t, t0 = t0,
                               c = as.numeric(c_fore)[i], useSplines = FALSE, cut = c(0, 0), verbose = FALSE)
    }
    
    dens2_fore = sapply(res_fore, function(r){
      tmp = r$dSup + adj_const
      ind = which(tmp < 1)
      tmp[ind] = tmp[max(ind)]*(tmp[ind] - tmp[min(ind)])/diff(range(tmp[ind]))
      
      ind = which(tmp > 109)
      tmp[ind] = tmp[min(ind)] + (110 - tmp[min(ind)])*(tmp[ind] - tmp[min(ind)])/diff(range(tmp[ind]))
      
      approx(x = tmp, y = r$dens, xout = u + adj_const)$y
    })
    
    if(any(is.na(dens2_fore)))
    {
        dens2_fore_rm_zero <- replace(dens2_fore, which(is.na(dens2_fore)), 0)
    }
    else
    {
        dens2_fore_rm_zero = dens2_fore
    }
    
    female_dens_fore = sapply(1:fh, function(i){dens2_fore_rm_zero[,i]/trapzRcpp(X = u, Y = dens2_fore_rm_zero[,i])})
    female_dens_fore_norm = sapply(1:fh, function(i){female_dens_fore[,i]/sum(female_dens_fore[,i])})
    colnames(female_dens_fore) = 1:fh
    
    return(list(female_dens_recon = female_dens_recon_norm * const,
                female_dens_forecast = female_dens_fore_norm * const,
                c = c, c_fore = c_fore, lqd = lqd, lqd_fore = den_fore))
}

# compute distance between two densities

density_dist <- function(d1, d2, Grid)
{
    Positive = (d1 > 1e-6) & (d2 > 1e-6)
    Grid = Grid[Positive]
    d1 = d1[Positive]
    d2 = d2[Positive]
    grid.size=length(Grid)
    first_integral=c()
    for(t in 1:grid.size){
      first_integral[t]=trapz(Grid,(log(d1/d1[t])-log(d2/d2[t]))^2)
    }
    return(trapz(Grid,first_integral)/2)
}

# find mode

getmode <- function(v)
{
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
}
