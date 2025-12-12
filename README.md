# Bayes_NW
Nonparametric density-on-density regression in Bayes Hilbert space

1. CoDa_NFR.R: centered log-ratio transformation + nonparametric functional regression for modeling and forecasting transformed data

2. CoDa_PCA.R: centered log-ratio transformation + functional principal component analysis

3. lqd_fun.R: log quantile density transformation of Petersen and Muller (2016)

4. Horta_Ziegelmann_FPCA.R: Functional principal component analysis ignoring density-valued objects, perform normalization in the end

5. load_packages.R: load R packages

6. save_function.R: internal functions for saving R figures

7. auxiliary.R: auxiliary functions for implementing nonparametric function-on-function regression

8. d-on-d_simulation.html: R markdown file for simulation DGP

9. Bovespa_Bayes_NW.RData: five existing method, proposed Bayes NW estimator, benchmark random walk

10. Bovespa_returns_data.R: Bovespa return density forecasting, same data set has previously been analysed in Kokoszka et al. (2019, IJF)

11. US_LT_density.R: Age-specific life-table death counts in the United States

# order of loading files

1. source("load_packages.R")
2. source("save_function.R")
3. source("CoDa_NFR.R")
4. source("CoDa_PCA.R")
5. source("lqd_fun.R")
6. load("Bovespa_Bayes_NW.RData")
7. file.edit("Bovespa_returns_data.R")
8. file.edit("US_LT_density.R")
