################
# map of France
################

map_fr <- map_data("france")
p_map_fr =
  ggplot(data = map_fr, aes(x = long, y = lat, group = group, fill = region)) +
  geom_polygon() + coord_equal() + scale_fill_discrete(guide = "none") +
  ggtitle("France divided into 96 départements") +
  theme(
    plot.title = element_text(face = "bold", colour = "black", size = 18, hjust=0.5),
    axis.title.x = element_text(face = "bold", colour = "black", size = 18, hjust=0.5),
    axis.title.y = element_text(face = "bold", colour = "black", size = 18, hjust=0.5),
    axis.text.x = element_text(face = "bold", colour = "black", size = 14, hjust=0.5),
    axis.text.y = element_text(face = "bold", colour = "black", size = 14, hjust=0.5)
  ) +
  theme(aspect.ratio = 1.1)

# save the figure

ggsave("map_fr.png")
p_map_fr
dev.off()

####################
# build the dataset
####################

Popbydep = c(656955, 526050, 331315, 165197, 141756, 1079396, 326875, 265531, 152398, 309907, 372705, 278360, 2034469, 691453,
             142811, 348180, 647080, 296404, 240336, 532886, 596186, 116270, 408393, 539449, 520560, 600687, 429425,
             906554, 162421, 182258, 748468, 1400935, 190040, 1633440, 1176145, 1082073, 217139, 605380, 1264979, 257849, 411979, 327835, 764737, 226901,
             1437137, 682890, 173166, 330336, 76286, 815881, 490669, 563823, 169250, 305365, 730398, 181641, 755566, 1035866, 199596, 2588988,
             825077, 276903, 1452778, 660240, 683169, 226839, 479000, 1132607, 763204, 1876051, 233194, 547824, 560227, 432548, 828405, 2148271,
             1243788, 1423607, 1448625, 372627, 569769, 387898, 262618, 1073836, 560997, 683187, 437398, 370774, 359520, 332096, 140145, 1319401,
             1613762, 1670149, 1406041, 1248354, 376879, 358749, 290691, 859959, 279471)

# department names
Dep.names = c("Ain", "Aisne", "Allier", "Alpes-de-Haute-Provence", "Hautes-Alpes", "Alpes-Maritimes", "Ardèche", "Ardennes", "Ariège",
              "Aube", "Aude", "Aveyron", "Bouches-du-Rhône", "Calvados", "Cantal", "Charente", "Charente-Maritime", "Cher", "Corrèze", "Côte-d'Or", "Côtes-d'Armor", "Creuse", "Dordogne", "Doubs", "Drôme", "Eure","Eure-et-Loir", "Finistère", "Corse-du-Sud",
              "Haute-Corse", "Gard",
              "Haute-Garonne", "Gers", "Gironde", "Hérault", "Ille-et-Vilaine", "Indre", "Indre-et-Loire", "Isère", "Jura", "Landes", "Loir-et-Cher",
              "Loire", "Haute-Loire", "Loire-Atlantique", "Loiret", "Lot", "Lot-et-Garonne", "Lozère", "Maine-et-Loire", "Manche", "Marne", "Haute-Marne",
              "Mayenne", "Meurthe-et-Moselle", "Meuse", "Morbihan", "Moselle", "Nièvre", "Nord", "Oise", "Orne", "Pas-de-Calais", "Puy-de-Dôme",
              "Pyrénées-Atlantiques", "Hautes-Pyrénées", "Pyrénées-Orientales", "Bas-Rhin", "Haut-Rhin", "Rhône", "Haute-Saône", "Saône-et-Loire",
              "Sarthe", "Savoie", "Haute-Savoie", "Paris", "Seine-Maritime", "Seine-et-Marne", "Yvelines", "Deux-Sèvres", "Somme", "Tarn",
              "Tarn-et-Garonne", "Var", "Vaucluse", "Vendée", "Vienne", "Haute-Vienne", "Vosges", "Yonne", "Territoire de Belfort", "Essonne",
              "Hauts-de-Seine", "Seine-Saint-Denis", "Val-de-Marne", "Val-d'Oise", "Guadeloupe", "Martinique", "Guyane", "La Réunion", "Mayotte")
Dep.code = c(paste0("0", 1:9), as.character(10:19), "2A", "2B", as.character(21:95), as.character(971:974), "976")

#######################################
# upload covid19 dataset from data.gov
#######################################

options(timeout=100)
link = "https://protect-au.mimecast.com/s/7IPkCoV1Y2S4w5BXIVHCky?domain=data.gouv.fr"
invisible(file(link, open = "r"))
RAWDATA = read.csv(link, sep = ";", header = TRUE)
date.start = as.Date(RAWDATA[1, "jour"])
day.start = as.numeric(substr(date.start, start = 9, stop = 10))
month.start = as.numeric(substr(date.start, start = 6, stop = 7))
if(month.start < 10){month.start.code = paste0("0", month.start)}else{ month.start.code = as.character(month.start)}
year.start = as.numeric(substr(date.start, start = 1, stop = 4))
nbrow = nrow(RAWDATA)
date.last = as.Date(RAWDATA[nbrow, "jour"])
day.last = as.numeric(substr(date.last, start = 9, stop = 10))
month.last = as.numeric(substr(date.last, start = 6, stop = 7))
if(month.last < 10){month.last.code = paste0("0", month.last)}else{month.last.code = as.character(month.last)}
year.last = as.numeric(substr(date.last, start = 1, stop = 4))
Months = c("janvier", "février", "mars", "avril", "mai", "juin", "juillet", "août", "septembre", "octobre", "novembre", "décembre")
Dates = seq(as.Date(date.start), as.Date(date.last), by="days")
RAWDATA1 = RAWDATA[RAWDATA$sexe == 0, -2]
# reorganize the data (one department by row)
nbrow1 = nrow(RAWDATA1)
nbdep = length(Dep.names)
#nbdates = nbrow1 %/% nbdep
#nbdates = ceiling(nbrow1 / nbdep)

DATA = matrix(NA, nbdep, length(Dates))
varname = "hosp"
ylimit = c(0, 0.5)
for(ii in 1:nbdep){
  Indices = seq(ii, nbrow1, by = nbdep + 1)
  DATA[ii, ] = RAWDATA1[Indices, varname]
  print(ii); rm(ii)
}
# drop territories outside Europe
DATA = DATA[-((nbdep-4):nbdep), ]
# nb of hospitalizations per 100,000 inhabitants
DATA0 = 100000 * DATA / Popbydep[-((nbdep-4):nbdep)]
dimnames(DATA0)[[1]] = c(paste0("0", 1:9), as.character(10:19), "2A", "2B", as.character(21:29),  as.character(c(30:95)))
nbdates = ncol(DATA0)
df = data.frame(dates = as.Date(Dates), Nb_cas_hosp = round(as.vector(t(DATA0)), digits = 2), dept=rep(dimnames(DATA0)[[1]], each=nbdates))

ggp <- ggplot(df, aes(dates, Nb_cas_hosp, col = dept)) +
  geom_line() +
  ggtitle("nb of hospitalizations per 100,000 inhabitants") +
  theme(
    plot.title = element_text(face = "bold", colour = "black", size = 18, hjust=0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(face = "bold", colour = "black", size = 14, hjust=0.5),
    axis.text.y = element_text(face = "bold", colour = "black", size = 14, hjust=0.5)
  ) +
  theme(
    legend.title = element_text(face = "bold", colour = "black", size = 15, hjust=0.5)
  )
#theme(legend.position = "bottom")
ggplotly(ggp)

# (Dep.names[1:96])[which.max(rowSums(exp(Z)-0.1))] # "Territoire de Belfort"
# (Dep.names[1:96])[which.min(rowSums(exp(Z)-0.1))] # "Finistère"

# varname = "hosp"
# ylimit = c(0, 0.5)
# for(ii in 1:nbdep){
#  Indices = seq(ii, nbrow1, by = nbdep)
#  DATA[ii, ] = RAWDATA1[Indices, varname]
# }
# hospitalization ratios
# DATA0 = 100000 * DATA / Popbydep
# dimnames(DATA0)[[1]] = c(paste0("0", 1:9), as.character(10:19), "2A", "2B", as.character(21:29),  as.character(c(30:95, 971:974, 976)))
# nbdates = ncol(DATA0)

# Probability Density Functions of transformed data (PFDZ)
n.xgrid = 500
PDFZ = matrix(0, nbdates, n.xgrid)
#
Z = log(DATA0 + 0.1, base = exp(1))
for(ii in 1:nbdates){
    # kernel density estimate
    res.density = density(Z[, ii], from = min(Z)-1, to = max(Z)+1, n = n.xgrid)
    PDFZ[ii, ] = res.density$y
}
xgrid = res.density$x
rownames(PDFZ) = Dates
colnames(PDFZ) = xgrid

### create a N-colorpalette with "colorRampPalette" R function from "cold" to "warm where N = nbdates"
# colfunc = colorRampPalette(c("orange", "springgreen", "royalblue", "black"))
colfunc = colorRampPalette(rainbow(7))
Colors = colfunc(nbdates)
DENSITIESOFZ = t(PDFZ)
startpoint = min(Z)-1
endpoint = max(Z)+1
laps = (endpoint - startpoint)/(n.xgrid - 1)
DENSITIESOFZ.ts = ts(DENSITIESOFZ, start = startpoint, end = endpoint, deltat = laps)
ts.plot(DENSITIESOFZ.ts, gpars = list(col = Colors, xlab = "", bty = "l"))

# save the figure

savepdf("france_densities", width = 12, height = 10, toplines = 0.8)
plot(fts(xgrid, t(PDFZ)), xlab = "Grid point", ylab = "Estimated kernel densities")
dev.off()

############
# densities
############

PDFZ = replace(PDFZ, which(PDFZ==0), 10^-6)
n_train = ceiling(nrow(PDFZ)/3 * 2)
n_test = floor(nrow(PDFZ)/3)
size_grid = ncol(PDFZ)

source("CoDa_NFR.R")
source("CoDa_recon.R")
source("lqd_fun.R")

###########
# CoDa NFR 
###########

radix = mean(rowSums(PDFZ))

# expanding window

CoDa_NFR_den_fore = matrix(NA, n_test, size_grid)
for(iwk in 1:n_test)
{
    CoDa_NFR_den_fore[iwk,] = CoDa_NFR(dat = PDFZ[1:(n_train - 1 + iwk),]) 
    print(iwk); rm(iwk)
}
    
# graphical displays

par(mfrow = c(1,2))
plot(fts(xgrid, t(PDFZ[(n_train+1):nrow(PDFZ),])), xlab = "Grid point", ylab = "Holdout densities", 
     ylim = range(c(range(PDFZ[(n_train+1):nrow(PDFZ),]), range(CoDa_NFR_den_fore))))  

plot(fts(xgrid, t(CoDa_NFR_den_fore * radix)), xlab = "Grid point", ylab = "Estimated densities",
     ylim = range(c(range(PDFZ[(n_train+1):nrow(PDFZ),]), range(CoDa_NFR_den_fore))))


###################
# CoDa_FPCA
# expanding window
###################

## ETS & ARIMA & RWD & RW

# normalize == "FALSE"

CoDa_FPCA_den_fore = CoDa_FPCA_den_fore_arima = 
CoDa_FPCA_den_fore_rwf = CoDa_FPCA_den_fore_rw = matrix(NA, n_test, size_grid)
for(iwk in 1:n_test)
{
    # ETS
    CoDa_FPCA_den_fore[iwk,] = CoDa_recon(dat = PDFZ[1:(n_train - 1 + iwk),]/radix, normalize = "FALSE", 
                                          fore_method = "ETS", fh = 1, varprop = 0.9, const = radix)$d_x_t_star_fore

    # ARIMA
    CoDa_FPCA_den_fore_arima[iwk,] = CoDa_recon(dat = PDFZ[1:(n_train - 1 + iwk),]/radix, normalize = "FALSE", 
                                                fore_method = "ARIMA", fh = 1, varprop = 0.9, const = radix)$d_x_t_star_fore
    
    # RWD
    CoDa_FPCA_den_fore_rwf[iwk,] = CoDa_recon(dat = PDFZ[1:(n_train - 1 + iwk),]/radix, normalize = "FALSE", 
                                              fore_method = "RWF", drift_term = "TRUE", fh = 1, varprop = 0.9, const = radix)$d_x_t_star_fore
    
    # RW
    CoDa_FPCA_den_fore_rw[iwk,] = CoDa_recon(dat = PDFZ[1:(n_train - 1 + iwk),]/radix, normalize = "FALSE", 
                                             fore_method = "RWF", drift_term = "FALSE", fh = 1, varprop = 0.9, const = radix)$d_x_t_star_fore
    print(iwk); rm(iwk)    
}

# normalize == "TRUE"

CoDa_FPCA_den_normalize_fore = CoDa_FPCA_den_normalize_fore_arima = 
CoDa_FPCA_den_normalize_fore_rwf = CoDa_FPCA_den_normalize_fore_rw = matrix(NA, n_test, size_grid)
for(iwk in 1:n_test)
{
    # ETS
    CoDa_FPCA_den_normalize_fore[iwk,] = CoDa_recon(dat = PDFZ[1:(n_train - 1 + iwk),], normalize = "TRUE", 
                                          fore_method = "ETS", fh = 1, varprop = 0.9, const = radix)$d_x_t_star_fore
    
    # ARIMA
    CoDa_FPCA_den_normalize_fore_arima[iwk,] = CoDa_recon(dat = PDFZ[1:(n_train - 1 + iwk),], normalize = "TRUE", 
                                                fore_method = "ARIMA", fh = 1, varprop = 0.9, const = radix)$d_x_t_star_fore
    
    # RWD
    CoDa_FPCA_den_normalize_fore_rwf[iwk,] = CoDa_recon(dat = PDFZ[1:(n_train - 1 + iwk),], normalize = "TRUE", 
                                              fore_method = "RWF", drift_term = "TRUE", fh = 1, varprop = 0.9, const = radix)$d_x_t_star_fore
    
    # RW
    CoDa_FPCA_den_normalize_fore_rw[iwk,] = CoDa_recon(dat = PDFZ[1:(n_train - 1 + iwk),], normalize = "TRUE", 
                                             fore_method = "RWF", drift_term = "FALSE", fh = 1, varprop = 0.9, const = radix)$d_x_t_star_fore
    print(iwk); rm(iwk)    
}

##########
# HZ_FPCA 
##########

setwd("~/Dropbox/Todos/Nonparametric functional regression for a time series of density/code/Auxiliary_functions/")
source("Horta_Ziegelmann_FPCA.R")

## expanding window

# ncomp_select = "FALSE"

HZ_FPCA_den_fore = matrix(NA, n_test, size_grid)
for(iwk in 1:n_test)
{
    HZ_FPCA_den_fore[iwk,] = Horta_Ziegelmann_FPCA(data = PDFZ[1:(n_train - 1 + iwk),], 
                                          gridpoints = xgrid, m = 500, ncomp_select = "FALSE", D_val = 6)$Yhat.fix_den
    print(iwk); rm(iwk)
}

# ncomp_select = "TRUE"

HZ_FPCA_den_fore_ncomp_select = matrix(NA, n_test, size_grid)
for(iwk in 1:n_test)
{
    HZ_FPCA_den_fore_ncomp_select[iwk,] = Horta_Ziegelmann_FPCA(data = PDFZ[1:(n_train - 1 + iwk),], 
                                                   gridpoints = xgrid, m = 500, ncomp_select = "TRUE")$Yhat.fix_den
    print(iwk); rm(iwk)
}

 
############
# LQDT_FPCA 
############

# expanding window

LQDT_FPCA_den_fore = LQDT_FPCA_den_fore_arima = matrix(NA, n_test, size_grid)
for(iwk in 1:n_test)
{
    # ETS
  
    LQDT_FPCA_den_fore[iwk,] = LQDT_FPCA(data = PDFZ[1:(n_train - 1 + iwk),], 
                                         gridpoints = xgrid, m = size_grid,
                                         forecasting_method = "uni", fmethod = "ets")$dens_fore

    # ARIMA
    
    LQDT_FPCA_den_fore_arima[iwk,] = LQDT_FPCA(data = PDFZ[1:(n_train - 1 + iwk),], 
                                         gridpoints = xgrid, m = size_grid,
                                         forecasting_method = "uni", fmethod = "arima")$dens_fore
    print(iwk); rm(iwk)
}

# plot density forecasts

# savepdf("french_density_fore", width = 18, height = 10, toplines = 0.8)
savepdf("Fig_10", width = 18, height = 10, toplines = 0.8)
par(mfrow = c(2,3))
plot(fts(xgrid, t(PDFZ[(n_train+1):(n_train + n_test),])), xlab = "Grid", ylab = "Holdout density", ylim = c(0, 1.3))
plot(fts(xgrid, t(CoDa_NFR_den_fore * radix)), xlab = "", ylab = "Density forecasts", main = "Bayes NW", ylim = c(0, 1.3))
plot(fts(xgrid, t(CoDa_FPCA_den_fore_arima)), xlab = "", ylab = "Density forecasts", main = "CoDa", ylim = c(0, 1.3))

plot(fts(xgrid, t(PDFZ[(n_train+1):(n_train + n_test),])), xlab = "", type = "n", ylab = "", axes = F, ylim = c(0, 1.3))
plot(fts(xgrid, t(LQDT_FPCA_den_fore)), xlab = "Grid", ylab = "Density forecasts", main = "LQDT", ylim = c(0, 1.3))
plot(fts(xgrid, t(HZ_FPCA_den_fore)), xlab = "Grid", ylab = "Density forecasts", main = "HZ", ylim = c(0, 1.3))
#plot(fts(xgrid, t(PDFZ[(n_train):(n_train + n_test - 1),])), xlab = "Grid", ylab = "Density forecasts", main = "RW", ylim = c(0, 1.3))
dev.off()



##############################
# Kullback-Leibler divergence
##############################

KLdiv_NFR = KLdiv_CoDa = KLdiv_LQDT = KLdiv_HZ = KLdiv_RW = vector("numeric", n_test)
for(ik in 1:n_test)
{
    dat = cbind(true = PDFZ[(n_train + ik),], forecast = CoDa_NFR_den_fore[ik,]*radix)
    colnames(dat) = c("True", "Estimate")
    KLdiv_NFR[ik] = mean(as.numeric(KLdiv(dat, eps = 1e-16))[2:3])
    rm(dat)
    
    dat = cbind(true = PDFZ[(n_train + ik),], forecast = CoDa_FPCA_den_fore_arima[ik,])
    colnames(dat) = c("True", "Estimate")
    KLdiv_CoDa[ik] = mean(as.numeric(KLdiv(dat, eps = 1e-16))[2:3])
    rm(dat)
    
    dat = cbind(true = PDFZ[(n_train + ik),], forecast = LQDT_FPCA_den_fore[ik,])
    colnames(dat) = c("True", "Estimate")
    KLdiv_LQDT[ik] = mean(as.numeric(KLdiv(dat, eps = 1e-16))[2:3])
    rm(dat)
    
    dat = cbind(true = PDFZ[(n_train + ik),], forecast = HZ_FPCA_den_fore[ik,])
    colnames(dat) = c("True", "Estimate")
    KLdiv_HZ[ik] = mean(as.numeric(KLdiv(dat, eps = 1e-16))[2:3])
    rm(dat)
    
    dat = cbind(true = PDFZ[(n_train + ik),], forecast = PDFZ[(n_train + ik - 1),])
    colnames(dat) = c("True", "Estimate")
    KLdiv_RW[ik] = mean(as.numeric(KLdiv(dat, eps = 1e-16))[2:3])
    rm(dat); print(ik); rm(ik)
}

#####################
# summary statistics
#####################

KLD_output = cbind(KLdiv_NFR, KLdiv_CoDa, KLdiv_LQDT, KLdiv_HZ, KLdiv_RW)
colnames(KLD_output) = c("Bayes NW", "CoDa", "LQDT", "HZ", "RW")
rownames(KLD_output) = (n_train + 1):(n_train + n_test)

require(xtable)
xtable(apply(KLD_output, 2, summary), digits = 4)

