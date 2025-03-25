##################
# load R packages
##################

# List of required packages
packages <- c("MCS", "vars", "maps", "doMC", "ftsa", "MASS", "boot", "psych", "pdist", "pracma", "plotly",
              "fdapace", "fda.usc", "flexmix", "fdapace", "ggplot2", "truncnorm", "colorRamps", "matrixStats")

# Install missing packages
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

# Load all packages
lapply(packages, library, character.only = TRUE)

