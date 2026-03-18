<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/hanshang/CLR_vs_CDF_transformation">
    <img src="MQ.png" alt="Logo" width="250" height="250">
    <img src="Toulouse.png" alt="Logo" width="200" height="200">    
  </a>

<h3 align="center">Density-valued time series: 
  <br />
  Nonparametric density-on-density regression</h3>


<!-- ABOUT THE PROJECT -->
## Abstract
This paper is concerned with forecasting probability density functions. Density functions are nonnegative and have a constrained integral; thus, they do not constitute a vector space. Implementing unconstrained functional time-series forecasting methods is problematic for such nonlinear and constrained data. A novel forecasting method is developed based on a nonparametric function-on-function regression, where both the response and the predictor are probability density functions. The asymptotic properties of our nonparametric regression estimator are established, as well as its finite-sample performance, through a series of Monte-Carlo simulation studies. Using Bovespa intraday 5-minute returns and age-specific period life tables from the United States, we assess and compare the finite-sample forecast accuracy of the proposed method with several existing methods.

### Main Results
The R script files in the `R Code` folder should be used in the following order:
1. CoDa_NFR.R: centered log-ratio transformation + nonparametric functional regression for modeling and forecasting transformed data
2. CoDa_PCA.R: centered log-ratio transformation + functional principal component analysis
3. lqd_fun.R: log quantile density transformation of Petersen and Muller (2016)
4. Horta_Ziegelmann_FPCA.R: Functional principal component analysis ignoring density-valued objects, perform normalization in the end
5. load_packages.R: load R packages
6. save_function.R: internal functions for saving R figures
7. auxiliary.R: auxiliary functions for implementing nonparametric function-on-function regression
8. density_norm.R: internal functions for computing norm of densities and obtaining mode of densities
9. d-on-d_simulation.html: R markdown file for simulation DGP
10. Bovespa_Bayes_NW.RData: five existing method, proposed Bayes NW estimator, benchmark random walk
11. Bovespa_returns_data.R: Bovespa return density forecasting, same data set has previously been analysed in Kokoszka et al. (2019, IJF)
12. US_LT_density.R: Age-specific life-table death counts in the United States

# order of loading files

1. source("load_packages.R")
2. source("save_function.R")
3. source("density_norm.R")
4. source("CoDa_NFR.R")
5. source("CoDa_PCA.R")
6. source("lqd_fun.R")
7. load("Bovespa_Bayes_NW.RData")
8. file.edit("Bovespa_returns_data.R")
9. file.edit("US_LT_density.R")

## Contact
arXiv link: [https://arxiv.org/abs/2507.04303](https://arxiv.org/abs/2507.04303)

Frederic Ferraty - ferraty@math.univ-toulouse.fr

Han Lin Shang - hanlin.shang@mq.edu.au

<br />
    <a href="https://github.com/hanshang/Bayes_NW"><strong>Explore R code »</strong></a>
<br />
