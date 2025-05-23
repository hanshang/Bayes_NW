# dat: n by p data matrix
# normalize: normalization step
# fore_method: univariate time-series forecasting method
# drift_term: when RWF is used
# fh: forecast horizon
# varprop: explained variance
# constant: radix (in the case of life-table death counts, it's 10^5)

CoDa_recon <- function(dat, normalize = c("TRUE", "FALSE"), fore_method = c("RWF", "ETS", "ARIMA"),
                       drift_term, fh, varprop, constant)
{
  fore_method = match.arg(fore_method)
  n_year = nrow(dat)
  n_age = ncol(dat)
  dat_center = sweep(dat, 1, apply(dat, 1, sum), "/")
  alpha_x = vector("numeric", n_age)
  for (ik in 1:n_age) {
    alpha_x[ik] = geometric.mean(dat_center[, ik])
  }
  f_x_t = matrix(NA, n_year, n_age)
  if (normalize == "TRUE") {
    for (ik in 1:n_year) {
      f_x_t[ik, ] = (dat[ik, ]/alpha_x)/sum(dat[ik, ]/alpha_x)
    }
  }
  if (normalize == "FALSE") {
    for (ik in 1:n_year) {
      f_x_t[ik, ] = dat[ik, ]/alpha_x
    }
  }
  g_t = vector("numeric", n_year)
  h_x_t = matrix(NA, n_year, n_age)
  if (normalize == TRUE) {
    for (ik in 1:n_year) {
      g_t[ik] = geometric.mean(f_x_t[ik, ])
      h_x_t[ik, ] = log(f_x_t[ik, ]/g_t[ik])
    }
  }
  if (normalize == FALSE) {
    for (ik in 1:n_year) {
      h_x_t[ik, ] = log(f_x_t[ik, ])
    }
  }
  SVD_decomp = svd(h_x_t)
  ncomp = head(which(cumsum(SVD_decomp$d^2/sum(SVD_decomp$d^2)) >= 
                       varprop), 1)
  basis_fore = basis = SVD_decomp$v[, 1:ncomp]
  score_forecast = score = t(basis) %*% t(h_x_t)
  recon = basis %*% score
  f_x_t_star_recon = d_x_t_star_recon = matrix(NA, n_age, n_year)
  for (ik in 1:n_year) {
    f_x_t_star_recon[, ik] = exp(recon[, ik])/sum(exp(recon[, 
                                                            ik]))
    d_x_t_star_recon[, ik] = (f_x_t_star_recon[, ik] * alpha_x)/sum(f_x_t_star_recon[, 
                                                                                     ik] * alpha_x)
  }
  score_fore = matrix(NA, ncomp, fh)
  for (ik in 1:ncomp) {
    if (fore_method == "ARIMA") {
      score_fore[ik, ] = forecast(auto.arima(as.numeric(score_forecast[ik, 
      ])), h = fh)$mean
    }
    if (fore_method == "ETS") {
      score_fore[ik, ] = forecast(ets(as.numeric(score_forecast[ik, 
      ])), h = fh)$mean
    }
    if (fore_method == "RWF") {
      score_fore[ik, ] = rwf(as.numeric(score_forecast[ik, 
      ]), h = fh, drift = drift_term)$mean
    }
  }
  fore_val = basis_fore %*% score_fore
  f_x_t_star_fore = d_x_t_star_fore = matrix(NA, n_age, fh)
  for (ik in 1:fh) {
    f_x_t_star_fore[, ik] = exp(fore_val[, ik])/sum(exp(fore_val[, 
                                                                 ik]))
    d_x_t_star_fore[, ik] = (f_x_t_star_fore[, ik] * alpha_x)/sum((f_x_t_star_fore[, 
                                                                                   ik] * alpha_x))
  }
  return(list(d_x_t_star_recon = d_x_t_star_recon * constant, 
              d_x_t_star_fore = d_x_t_star_fore * constant))
}

# data: n by p data matrix
# normalization: TRUE or FALSE
# varprop: explained variance
# fmethod: forecasting method
# drift_term: random walk with or without drift

CoDa_FPCA <- function(data, normalization, varprop = 0.99, fmethod, drift_term)
{
    c = colSums(data)[1]
    dum = CoDa_recon(dat = data, normalize = normalization,
                     fore_method = fmethod, drift_term = drift_term, fh = 1,
                     varprop = varprop, constant = c)
    return(dum$d_x_t_star_fore)
}
