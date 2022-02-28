## Open libraries ----
library (tidyverse)
library(vegan)
#devtools::install_github ('zdealveindy/weimea')
library(weimea)
library(car)
library(FD)
library(dplyr)
library(MuMIn)
library(DHARMa)
library(rsq)

library(viridis)
library(VennDiagram)
library(cowplot)
#devtools::install_github("gavinsimpson/ggvegan")
library(ggvegan)
library(png)
library(viridis)

library(ape)
library(scales)
library(ecospat)
library(lmodel2)
#new functions from 'DecomposingFD package' ----

#theoretical papers: Scheiner, S.M., Kosman, E., Presley, S.J., Willig, M.R., 2017. Decomposing functional diversity. Methods Ecol. Evol. 8, 809-820. https://doi.org/10.1111/2041-210X.12696
#                    Scheiner, S.M., 2019. A compilation of and typology for abundance-, phylogenetic- And functional-based diversity metrics. bioRxiv 1-29. https://doi.org/10.1101/530782

#case-study paper that wrote the script/package: Malavasi, M., Bartak, V., Carranza, M.L., Simova, P., Acosta, A.T.R., 2018. Landscape pattern and plant biodiversity in Mediterranean coastal dune ecosystems: Do habitat loss and fragmentation really matter? J. Biogeogr. 45, 1367-1377. https://doi.org/10.1111/jbi.13215



## FTD: A function to calculate Scheiner et al. (2016)'s functional
## trait dispersion, written by S.A. Kothari
## with a few lines of code shamelessly stolen from J.J. Grossman

## contains two functions:
## 1. FTD uses a trait distance matrix to calculate functional diversity
## 2. FTD.comm applies FTD across the rows of a community x species matrix
## to calculate FTD for each community

## FTD requires:
## matrix or dist containing species functional distances (tdmat)
## tdmat should be rescaled to 0<=dij<=1
## Hill coefficient q
## weights are proportions of each species in each plot

FTD<-function(tdmat,weights=NULL,q=1){
  
  ## contingency for one-species communities
  if(length(tdmat)==1 && tdmat==0){
    tdmat<-as.matrix(tdmat)
  }
  
  ## is the input a (symmetric) matrix or dist? if not...
  if(!(class(tdmat) %in% c("matrix","dist"))){
    stop("distances must be class dist or class matrix")
  } else if(class(tdmat)=="matrix" && !isSymmetric(unname(tdmat))){
    warning("trait matrix not symmetric")
  } else if(class(tdmat)=="dist"){
    tdmat<-as.matrix(tdmat)
  }
  
  if(!isTRUE(all.equal(sum(diag(tdmat)),0))){
    warning("non-zero diagonal; species appear to have non-zero trait distance from themselves")
  }
  
  if(max(tdmat)>1 || min(tdmat)<0){
    tdmat<-(tdmat-min(tdmat))/(max(tdmat)-min(tdmat))
    warning("trait distances must be between 0 and 1; rescaling")
  }
  
  ## if no weights are provided, abundances are assumed equal
  if(is.null(weights)){
    nsp<-nrow(tdmat)
    weights<-rep(1/nsp,nsp)
  } else {
    nsp<-sum(weights>0)
  }
  
  if(!isTRUE(all.equal(sum(weights),1))){
    weights<-weights/sum(weights)
    warning("input proportional abundances do not sum to 1; summation to 1 forced")
  }
  
  tdmat.abund<-diag(weights) %*% tdmat %*% diag(weights)
  ## Here, because sum(weights)=1, the sum is here already adjusted by dividing by the 
  ## square of the number of species (if weights=NULL)
  ## or by multiplying by proportional abundances
  M<-sum(tdmat.abund)
  ## M equals Rao's Q in abundance-weighted calculations
  M.prime<-ifelse(nsp==1,0,M*nsp/(nsp-1))
  fij<-tdmat.abund/M
  
  ## calculating qHt
  ## fork -- if q=1, 1/(1-q) is undefined, so we use an analogue
  ## of the formula for Shannon-Weiner diversity
  ## if q!=1, we can calculate explicitly
  ## by definition, qHt=0 when all trait distances are zero
  if(isTRUE(all.equal(M,0))){
    qHt<-0
  } else if(q==1){
    fijlog<-fij*log(fij)
    fijlog[is.na(fijlog)]<-0
    qHt<-exp(-1*sum(fijlog))
  } else if(q==0){
    qHt<-sum(fij>0)
  } else {
    qHt<-sum(fij^q)^(1/(1-q))
  }
  
  ## getting qDT, qDTM, and qEt from qHt
  qDT<-(1+sqrt(1+4*qHt))/2
  qDTM<-1+qDT*M
  qEt<-qDT/nsp
  
  list(nsp=nsp,q=q,M=M,M.prime=M.prime,qHt=qHt,qEt=qEt,qDT=qDT,qDTM=qDTM)
}

## wrapper for the above function across multiple communities
## requires distance matrix w/ all species, scaled to 0-1
## community data matrix (communities x species)
## with species in same order as distance matrix
## if abund=T, values in community data matrix are treated as abundances
## otherwise, all are converted to presence-absence
## if match.names=T, the code will match species names across the
## trait distance matrix and comm data matrix and rearrange the latter

FTD.comm<-function(tdmat,spmat,q=1,abund=F,match.names=F){
  
  ## is the input a (symmetric) matrix or dist? if not...
  if(!(class(tdmat) %in% c("matrix","dist"))){
    stop("distances must be class dist or class matrix")
  } else if(class(tdmat)=="matrix" && !isSymmetric(unname(tdmat))){
    warning("trait matrix not symmetric")
  } else if(class(tdmat)=="dist"){
    tdmat<-as.matrix(tdmat)
  }
  
  if(!isTRUE(all.equal(sum(diag(tdmat)),0))){
    warning("non-zero diagonal; species appear to have non-zero trait distance from themselves")
  }
  
  if(max(tdmat)>1 || min(tdmat)<0){
    tdmat<-(tdmat-min(tdmat))/(max(tdmat)-min(tdmat))
    warning("trait distances must be between 0 and 1; rescaling")
  }
  
  if(abund==F){
    spmat[spmat>0]<- 1
    spmat<-spmat/rowSums(spmat)
  }
  
  n.comm<-nrow(spmat)
  if(match.names==T){
    sp.arr<-match(rownames(as.matrix(tdmat)),colnames(spmat))
    spmat<-spmat[,sp.arr]
  }
  
  ## apply FTD to each community in turn
  out<-apply(spmat,1,function(x) unlist(FTD(tdmat=tdmat,weights=x,q=q)))
  df.out<-data.frame(t(out))
  rownames(df.out)<-rownames(spmat)
  ## warning for zero-species communities
  if(sum(df.out$nsp==0)>0){
    warning("at least one community has no species")
  }
  
  nsp<-sum(colSums(spmat>0))
  ## this is Sw -- always an arithmetic mean, according to Evsey Kosman
  u.nsp<-mean(df.out$nsp)
  ## calculate mean richness, dispersion, evenness, FTD
  u.M<-sum(df.out$nsp*df.out$M)/sum(df.out$nsp)
  
  if(q==1){
    ## geometric mean -- limit of generalized mean as q->1
    u.qDT<-prod(df.out$qDT)^(1/n.comm)
  } else {
    ## generalized mean with m=1-q
    u.qDT<-(sum(df.out$qDT^(1-q))/n.comm)^(1/(1-q))
  }
  u.M.prime<-u.M*u.nsp/(u.nsp-1)
  
  ## calculate mean FTD and evenness
  u.qDTM<-1+u.qDT*u.M
  u.qEt<-u.qDT/u.nsp
  
  ## list more things
  list(com.FTD=df.out,nsp=nsp,u.nsp=u.nsp,u.M=u.M,u.M.prime=u.M.prime,u.qEt=u.qEt,u.qDT=u.qDT,u.qDTM=u.qDTM)
}

#for the 100m2 plots (L = plotxsp, e = plotxenv, t_w = spxwoody_traits, t_f = spxfern_traits)
#load dataframes ----
L <- read.delim (file = 'lala_100_spe.txt', row.names = 1)
e <- read.delim (file = 'lala_100_env.txt', row.names = 1)
t_w <- read.delim (file = 'lala_woody_traits.txt', row.names = 1)
rownames (t_w) <- t_w [,2]
t_f <- read.delim (file = 'lala_fern_traits.txt', row.names = 1)
rownames (t_f) <- t_f [,2]

#separate the species composition data in tree layer (L_w) for woody species and understory (L_u) for ferns. 
L_w <- L[, grepl("_T", names(L))]
L_u <- L[, grepl("_U", names(L))]

#remove Alsopodo (a tree fern) from the woody species composition dataset!
L_w <- subset(L_w, select=-Alsopodo_T)

#covert L_w to presence/absence (L_w01), to match data from understory
L_w01 <- L_w
L_w01[L_w01>0] <- 1

#separate the understory species composition in terrestrial (L_ft)
sp_list100 <- read.delim (file = 'lala_100_checklist.txt', row.names = 1, encoding = 'UTF-8')
fern_list <- sp_list100[sp_list100$fern == 'TRUE', ]
names(L_u) <- sub("_U", "", names(L_u))
L_f <- L_u[, names(L_u) %in% fern_list$species_abbrev ]

fern_list_ter <- fern_list[fern_list$life_form!='E',] #so ALSO includes 2 understory tree ferns!!!
L_ft <- L_f[, names(L_f) %in% fern_list_ter$species_abbrev ]

#calculate the percentage of ferns in the understory (average number of species per plot)
terr_list <- sp_list100[sp_list100$life_form!='E', ]
L_ut <- L_u[, names(L_u) %in% terr_list$species_abbrev ]
S_understory <- rowSums(L_ut[,])
S_ferns <- rowSums(L_ft[,]) 
Percentage_ferns <- S_ferns/S_understory*100 
mean(Percentage_ferns)  #25.6%...

terr_list2 <- terr_list[terr_list$life_form=='H', ]  #excluding tree & shrub & liana (seedlings), so only herbaceous species!
L_ut2 <- L_u[, names(L_u) %in% terr_list2$species_abbrev ]
S_understory2 <- rowSums(L_ut2[,])
Percentage_ferns2 <- S_ferns/S_understory2*100 
mean(Percentage_ferns2)  #66.1%...
rm(L, L_u, L_f, fern_list, L_ut, S_understory, S_ferns, Percentage_ferns, terr_list2, L_ut2, S_understory2, Percentage_ferns2)

#calculate CWM trait matrices (T) and FD  ----
#only use traits measured for both trees & ferns (Note: EWT, not succulence. leaf_N, not leaf_CN, log_LA, not LA)
#evaluate sp x trait values
hist(t_w$LDMC)
hist(t_w$Lth)
hist(t_w$Chl_SPAD)
hist(t_w$SLA)
hist(t_w$EWT)
hist(t_w$d13C)
hist((t_w$d15N))
hist(t_w$leaf_N)
hist(log(t_f$LA)) #log!

hist(t_f$LDMC)
hist(t_f$Lth)
hist(t_f$Chl_SPAD)
hist(t_f$SLA)
hist(t_f$EWT)
hist(t_f$d13C)
hist((t_f$d15N))
hist(t_f$leaf_N)
hist(log(t_w$LA)) #log!

t_f$log_LA <- log(t_f$LA)
t_w$log_LA <- log(t_w$LA)

traits <- c('LDMC', 'Lth', 'Chl_SPAD', 'log_LA', 'SLA', 'EWT', 'd13C', 'd15N', 'leaf_N')

#filter trait dataframes to include only those species present in the L dataframes. 
rownames(t_f) <- sub("_U", "", rownames(t_f))
t_ft <- t_f[rownames(t_f) %in% fern_list_ter$species_abbrev, ] #select only traits for (terrestrial) species that occur in the 100m2 plots

t_w <- t_w[rownames(t_w) %in% colnames(L_w), ] #select only traits for woody species that occur in the 100m2 plots
rownames(t_w) <- sub("_T", "", rownames(t_w))

#calculate woody species CWMs (presabs)
colnames(L_w) <- sub("_T", "", colnames(L_w))
L_w.t <- L_w[, colnames(L_w) %in% rownames(t_w)] #select only species for which traits are measured (needed for FD function to work)
L_w.t <- L_w.t[, order(colnames(L_w.t))] 
t_w <- t_w[order(rownames(t_w)), ] 

#NEW: fill in MV!! (3 chlorophyll & 6 deltaC/N/soil N)
t_w_noMV <- t_w[, traits]
for(i in 1:ncol(t_w_noMV)){
  t_w_noMV[is.na(t_w_noMV[,i]), i] <- mean(t_w_noMV[,i], na.rm = TRUE)}
t_w[, traits] <- t_w_noMV

#D x <- dbFD(t_w[, traits], L_w.t, w.abun = FALSE, stand.x = FALSE, calc.CWM = T, 
#D           calc.FDiv = F, calc.FGR = F, calc.FRic = F, corr='none')
#D T_w01 <- x$CWM
T_w01 <- cwm (com = L_w.t, traits = t_w[,traits])
#D rm(x)

#evaluate CWM values
hist(T_w01$LDMC)
hist(T_w01$Lth)
hist(T_w01$Chl_SPAD)
hist(T_w01$log_LA) 
hist(T_w01$SLA)
hist(T_w01$EWT)
hist(T_w01$d13C)
hist((T_w01$d15N+4)**2) #transform!!!
hist(T_w01$leaf_N)

# T_w01$d15N2 <- (T_w01$d15N + 4)**2

traits_w01 <- c('LDMC', 'Lth', 'Chl_SPAD', 'log_LA', 'SLA', 'EWT', 'd13C', 'd15N', 'leaf_N') #final trait list for analysis, not using transformed d15N data
#T_w01b <- T_w01 #for making PCA of merged trait datasets
T_w01graph <- T_w01 #retain unstandardized dataset for graphing
#D T_w01 <- as.data.frame(scale(T_w01)) #standardize all data! (mean=0, sd=1) #D scaling was moved close to the CWM-RDa analysis
# hist(T_w01$d15N2)

#woody species richness & functional diversity
#Also note that we standardize traits for FD calculation! (but NOT for CWM calculations)
colnames(L_w01) <- sub("_T", "", colnames(L_w01))
L_w01.t <- L_w01[, colnames(L_w01) %in% rownames(t_w)] #select only species for which traits are measured (needed for FD function to work)
L_w01.t <- L_w01.t[, order(colnames(L_w01.t))] 
t_w <- t_w[order(rownames(t_w)), ] 

testx2 <- as.data.frame(scale(t_w[, traits])) #standardize all data! (mean = 0, sd = 1)

Dt_w_gower <- vegdist(testx2, method="gower", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)
pcoa1 <- pcoa(Dt_w_gower, correction="none", rn=NULL)
biplot(pcoa1)
Dt_w_gower2 <- rescale(Dt_w_gower) #rescale distances to 0-1 range
pcoa2 <- pcoa(Dt_w_gower2, correction="none", rn=NULL) #check if rescaling did not alter distances
biplot(pcoa2)

FTW_w <- FTD.comm(Dt_w_gower2,L_w01.t,q=1,abund=F,match.names=F)
warnings()  
FD_w01 <- FTW_w$com.FTD  
FD_w01 <- rename(FD_w01, 'M.prime_g' = 'M.prime', 'qEt_g' = 'qEt', 'qDTM_g' = 'qDTM')
#the calculated species richness (S) only counts species with trait values --> replace by full S richness
FD_w01$S <- rowSums(L_w01[,1:122]) #species richness

#new FD measures to use: qET = functional trait evenness or Hill pairwise equatibility (qE(TP)), M.prime/M' = mean trait dispersion (M(TP))
rm( FTW_w, pcoa1, pcoa2, testx2, Dt_w_gower, Dt_w_gower2)
rm(x, x2, testx)

#ferns CWMs
L_ft.t <- L_ft[, colnames(L_ft) %in% rownames(t_ft)] #select only species for which traits are measured (needed for FD function to work)
L_ft.t <- L_ft.t[, order(colnames(L_ft.t))] 
t_ft <- t_ft[order(rownames(t_ft)), ] 

#NO missing values in trait matrix for ferns, so NO need to change the trait matrix

x <- dbFD(t_ft[, traits], L_ft.t, w.abun = FALSE, stand.x = FALSE, calc.CWM = T, 
          calc.FDiv = F, calc.FGR = F, calc.FRic = F, corr = 'none')
T_ft <- x$CWM
rm(x)

#evaluate CWM values
hist(T_ft$LDMC)
hist(T_ft$Lth)
hist(T_ft$Chl_SPAD)
hist(T_ft$log_LA)
hist(T_ft$SLA)
hist(T_ft$EWT)
hist(T_ft$d13C)
hist(log(T_ft$d15N+2.5))  #transform!!!
hist(T_ft$leaf_N)

T_ft$log_d15N <- log(T_ft$d15N + 2.5)

traits_ft <- c('LDMC', 'Lth', 'Chl_SPAD', 'log_LA', 'SLA', 'EWT', 'd13C', 'd15N', 'leaf_N') #final trait list for analysis, without log-transformed d15N variable
#T_ftb <- T_ft #for making PCA of merged trait datasets
T_ftgraph <- T_ft #keep unstanderdized traits for graphing
T_ft <- as.data.frame(scale(T_ft)) #standardize all data! (mean=0, sd=1)
#hist(T_ft$log_d15N)

#Ferns species richness & functional diversity

#Also note that we standardize traits before FD calculation! (but NOT for CWM calculations)
L_ft.t <- L_ft[, colnames(L_ft) %in% rownames(t_ft)] #select only species for which traits are measured (needed for FD function to work)
L_ft.t <- L_ft.t[, order(colnames(L_ft.t))] 
t_ft <- t_ft[order(rownames(t_ft)), ] 

testx3 <- t_ft[, traits]
testx3 <- as.data.frame(scale(testx3)) #standardize all data! (mean = 0, sd = 1)

Dt_ft_gower <- vegdist(testx3, method="gower", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)
pcoa1 <- pcoa(Dt_ft_gower, correction="none", rn=NULL)
biplot(pcoa1)
Dt_ft_gower2 <- rescale(Dt_ft_gower) #rescale distances to 0-1 range
pcoa2 <- pcoa(Dt_ft_gower2, correction="none", rn=NULL) #check if rescaling did not alter distances
biplot(pcoa2)

FTW_ft <- FTD.comm(Dt_ft_gower2,L_ft.t,q=1,abund=F,match.names=F)
warnings()  
FD_ft <- FTW_ft$com.FTD  
FD_ft <- rename(FD_ft, 'M.prime_g' = 'M.prime', 'qEt_g' = 'qEt', 'qDTM_g' = 'qDTM')

#the calculated species richness (S) only counts species with trait values --> replace by full S richness
FD_ft$S <- rowSums(L_ft[,1:76])

rm(FTW_ft, pcoa1, pcoa2, testx3, Dt_ft_gower, Dt_ft_gower2)

#explore data in e matrix ----

hist(e$elevation) #OK = logical since this follows the sampling scheme
#NOTE: heatload is calculated from slope & folded aspect, so collinearity likely an issue
hist(e$slope)
hist(e$folded_aspect_NE) #use only folded_aspect, not normal aspect
hist(e$heatload**2) #use square (**2) transformation!
e$heatload2 <- e$heatload**2
hist(e$fog_freq)

plot(e$elevation, e$fog_freq)
plot(e$folded_aspect_NE, e$heatload)
plot(e$slope, e$heatload)

#use soil depth as independent variable
hist(e$soildpt_mean)
hist(e$soildpt_median)
plot(e$soildpt_mean, e$soildpt_median)
e <- dplyr::rename(e,soil_depth = soildpt_mean) #use the mean soil depth measures as final 'soil depth' variable
hist(e$soil_depth)

hist(sqrt(e$rock_soil)) #many zero's...

hist(e$pH_H2O)
hist(log(e$Ca)) #log
hist(e$Cu)
hist(e$Fe)
hist(e$K)
hist(log(e$Mg)) #log
hist(log(e$Mn)) #log
hist(log(e$P))  #log
hist(log(e$Zn)) #log
hist(log(e$CN_ratio)) #log
hist(e$total_N)
#hist(log(e$total_C)) #log

e$transCa <- log(e$Ca)
e$transMg <- log(e$Mg) 
e$transMn <- log(e$Mn) 
e$transP <- log(e$P)
e$transZn <- log(e$Zn)
e$transCN_ratio <- log(e$CN_ratio)
#e$transtotal_C <- log(e$total_C)
e$trans_rock_soil <- sqrt(e$rock_soil)
e$trans_rock_soil <- as.numeric(e$trans_rock_soil)

soil <- c("pH_H2O", "transCa", "Cu", "Fe", "K", "transMg", "transMn", "transP", "transZn", "transCN_ratio", "total_N") #NOT including total_C 
e[, soil] <- scale(e[, soil]) #standardize all data! (mean = 0, sd = 1)

#PCA on transformed soil nutrient data for data reduction ----
pcasoil <- rda(e[, soil] ~ 1)
plotsoil = e[, soil]
plotsoil <- dplyr::rename(plotsoil, pH = pH_H2O, Ca = transCa, Mg = transMg, Mn = transMn, P = transP, Zn = transZn, 'C/N' = transCN_ratio, N = total_N)
pcasoil2 <- prcomp(plotsoil)

# biplot(pcasoil)
# summary(pcasoil2)
scores(pcasoil, choices = c(1:3), display="species", scaling = 0)

# Importance of components:
#                         PC1    PC2    PC3     PC4    PC5     PC6    PC7     PC8      PC9     PC10     PC11
# Eigenvalue            5.7329 2.0025 1.2089 0.68129 0.5390 0.32658 0.1836 0.14915 0.109520 0.041634 0.024866
# Proportion Explained  0.5212 0.1820 0.1099 0.06194 0.0490 0.02969 0.0167 0.01356 0.009956 0.003785 0.002261
# Cumulative Proportion 0.5212 0.7032 0.8131 0.87505 0.9241 0.95374 0.9704 0.98400 0.993955 0.997739 1.000000

# Species loadings
#                      PC1         PC2         PC3
# pH_H2O        -0.2576526 -0.44630955  0.07646534
# transCa        0.2842770 -0.26488032 -0.37671408
# Cu             0.1837649 -0.28203996  0.58535738
# Fe            -0.1179140  0.39806770 -0.51188853
# K              0.3610352 -0.02394300  0.02309165
# transMg        0.4045079  0.04673847  0.01743987
# transMn        0.1202223 -0.58747073 -0.29778098
# transP         0.3805703 -0.06191440 -0.14262217
# transZn        0.3872585  0.16800868  0.04747734
# transCN_ratio  0.3075829  0.32826572  0.29027486
# total_N        0.3275544 -0.07030982 -0.22627608

axes <- scores(pcasoil, 1:3, display = 'site')
axes <- as.data.frame(axes)
axes <- dplyr::rename(axes, PC1s = PC1, PC2s = PC2, PC3s = PC3)
axes$PC2s = -axes$PC2s #REVERT axis 2!!
e <- cbind(e, axes)

rm(axes)

#check for collinearity issues ----
model <- lm(cover_rock ~ elevation + heatload2 + trans_rock_soil + slope +
              folded_aspect_NE + soil_depth + PC1s + PC2s + PC3s + fog_freq, data = e)
vif(model) #EXCLUDE folded_aspect_NE (collinearity with heatload & slope)
model <- lm(cover_rock ~ elevation + heatload2 + trans_rock_soil +
              soil_depth + PC1s + PC2s + PC3s + fog_freq, data = e)
vif(model) #EXCLUDE folded_aspect_NE (collinearity with heatload & slope)
#also exclude slope, because no clear link with climate. 
rm(model)
#plot(e$folded_aspect_NE, e$heatload2)

#RDA (& varpart) on presence/absence L matrices ----
#1. woody species
rda.w01.0 <- rda(L_w01 ~ 1, e)                                         #intercept model
rda.w01.1 <- rda(L_w01 ~ 1 + elevation + heatload2 + fog_freq + 
                   trans_rock_soil + soil_depth + PC1s + PC2s + PC3s, e)      #full model


# #test the use of a DB-RDA to overcome horseshoe (use square root bray curtis distance)  --> did NOT improve model (so not used)
# #more info: https://archetypalecology.wordpress.com/2018/02/21/distance-based-redundancy-analysis-db-rda-in-r/
# #more info: https://www.davidzeleny.net/anadat-r/doku.php/en:rda_cca
# evars = c( 'elevation',"fog_freq", "heatload2", "PC1s", "PC3s", "soil_depth", "PC2s", "trans_rock_soil")
# L_w01_plot <- L_w01
# rankindex(e[,evars], L_w01_plot, indices = c("euc", "man", "gow","bra", "kul"), stepacross= FALSE, method = "spearman")
# testtt = capscale(L_w01_plot ~ 1 + elevation + fog_freq + heatload2 + PC1s + PC3s + soil_depth + PC2s + trans_rock_soil, e, dist="kul", sqrt.dist=TRUE)
# testtt = capscale(L_w01_plot ~ 1 + elevation + fog_freq + PC2s + PC3s + PC1s + heatload2, e, dist="kul", sqrt.dist=TRUE)
# 
# plot(testtt) # use base plot, might be done with ggplot2
# anova(testtt) 

#the COMBINED model
ordiR2step(rda.w01.0, scope = formula(rda.w01.1), direction = "forward", perm.max = 9999, R2scope = TRUE) #automatic stepwise model selection based on R2
rda.w01 <- rda(L_w01 ~ 1 + elevation + fog_freq + PC2s + PC3s + PC1s + heatload2, e)   #final model 
summary(rda.w01)
plot(rda.w01)

RsquareAdj(rda.w01) #model adj R2 0.2921899
anova.cca(rda.w01)
anova.cca(rda.w01, by = "axis", permutations = 9999)
anova.cca(rda.w01, by = "margin", permutations = 9999)

#variation partitioning on the FULL model (so NOT on the reduced model)
vpart.w01 <- varpart(L_w01, ~ elevation + heatload2 + fog_freq, ~ trans_rock_soil + soil_depth + PC1s + PC2s + PC3s, data = e)
vpart.w01

p2_label <- paste('total~R^2 == ', round(vpart.w01$part$fract$Adj.R.squared[3]*100, 1))
p2 <- draw.pairwise.venn(
  round((vpart.w01$part$indfract$Adj.R.squared[1] + vpart.w01$part$indfract$Adj.R.squared[2])*100/ vpart.w01$part$fract$Adj.R.squared[3], digits=1),
  round((vpart.w01$part$indfract$Adj.R.squared[3] + vpart.w01$part$indfract$Adj.R.squared[2])*100/ vpart.w01$part$fract$Adj.R.squared[3], digits=1), 
  round(vpart.w01$part$indfract$Adj.R.squared[2]*100/ vpart.w01$part$fract$Adj.R.squared[3], digits = 1),
  c("climate proxies", "soil"), fill = c(viridis(3)[2], viridis(3)[3]), col = c(viridis(3)[2], viridis(3)[3]), alpha = c(0.6, 0.6), ind = FALSE, cat.pos = c(330, 30), 
  scaled = T, cat.dist = -0.08,cat.cex = 1.5)
p2 <- ggdraw(p2) + 
  annotate(geom = 'text', x = 0.65, y = 0.02, label = p2_label, parse = TRUE, size = 5) +
  theme(plot.background = element_rect(fill = NA), plot.margin = margin(15, 15, 15, 15)) + theme(aspect.ratio=1/1)
p2


#ferns
rda.ft.0 <- rda(L_ft ~ 1, e)                                         #intercept model
rda.ft.1 <- rda(L_ft ~ 1 + elevation + heatload2 + fog_freq + 
                  trans_rock_soil + soil_depth + PC1s + PC2s + PC3s, e)      #full model

#the COMBINED model
ordiR2step(rda.ft.0, scope = formula(rda.ft.1), direction = "forward", perm.max = 9999, R2scope = TRUE) #automatic stepwise model selection based on R2
rda.ft <- rda(L_ft ~ 1 + elevation + fog_freq + PC2s + heatload2, e)   #final model
summary(rda.ft)
plot(rda.ft)

RsquareAdj(rda.ft) #model adj R2 0.2981856
anova.cca(rda.ft)
anova.cca(rda.ft, by = "axis", permutations = 9999)
anova.cca(rda.ft, by = "margin", permutations = 9999)

#variation partitioning on the FULL model (so NOT on the reduced model)
vpart.ft <- varpart(L_ft, ~ elevation + heatload2 + fog_freq, ~ trans_rock_soil + soil_depth + PC1s + PC2s + PC3s, data = e)
vpart.ft

p1_label <- paste('total~R^2 == ', round(vpart.ft$part$fract$Adj.R.squared[3]*100, 1))
p1 <- draw.pairwise.venn(
  round((vpart.ft$part$indfract$Adj.R.squared[1] + vpart.ft$part$indfract$Adj.R.squared[2])*100/ vpart.ft$part$fract$Adj.R.squared[3], digits = 1),
  round((vpart.ft$part$indfract$Adj.R.squared[3] + vpart.ft$part$indfract$Adj.R.squared[2])*100/ vpart.ft$part$fract$Adj.R.squared[3], digits = 1),
  round(vpart.ft$part$indfract$Adj.R.squared[2]*100/ vpart.ft$part$fract$Adj.R.squared[3], digits = 1),
  c("climate proxies", "soil"), fill = c(viridis(3)[2], viridis(3)[3]), col = c(viridis(3)[2], viridis(3)[3]), alpha = c(0.6, 0.6), ind = F, cat.pos = c(330, 30),
  scaled = T, cat.dist = -0.08,cat.cex = 1.5)
p1 <- ggdraw(p1) + 
  annotate(geom = 'text', x = 0.65, y = 0.02, label = p1_label, parse = TRUE, size = 5) +
  theme(plot.background = element_rect(fill = NA), plot.margin = margin(15, 15, 15, 15))+ theme(aspect.ratio=1/1)
p1

#explore if woody species composition has additional effect on fern species composition (residuals of RDA) ----
#environmental vars as 'conditional' = partialling effect out
summary(rda.w01.0)
sum(rda.w01.0$CA$eig)/length(rda.w01.0$CA$eig) #need to include all axes with eig > 0.1805055  --> PC1 to 15! (this is 72.755% of total variation)

rda.ft_pr.0 <- rda(L_ft ~ 1 + Condition(elevation + heatload2 + fog_freq + trans_rock_soil + soil_depth + PC1s + PC2s + PC3s), e)  
rda.ft_prt <- rda(L_ft ~ 1 + rda.w01.0$CA$u[,1] + rda.w01.0$CA$u[,2] + rda.w01.0$CA$u[,3] + rda.w01.0$CA$u[,4] + rda.w01.0$CA$u[,5] + rda.w01.0$CA$u[,6] + rda.w01.0$CA$u[,7] +
                    rda.w01.0$CA$u[,8] + rda.w01.0$CA$u[,9] + rda.w01.0$CA$u[,10] + rda.w01.0$CA$u[,11] + rda.w01.0$CA$u[,12] + rda.w01.0$CA$u[,13] + rda.w01.0$CA$u[,14] +
                    rda.w01.0$CA$u[,15] + Condition(elevation + heatload2 + fog_freq + trans_rock_soil + soil_depth + PC1s + PC2s + PC3s), e) #NOT including tree_dens!
anova.cca(rda.ft_prt, permutations = 9999) #significant! F=1.38, p=0.0002
RsquareAdj(rda.ft_prt) #model adj R2 0.07083862
anova.cca(rda.ft_prt, by = "axis", permutations = 9999)
anova.cca(rda.ft_prt, by = "margin", permutations = 9999)
plot(rda.ft_prt)

# ordiR2step(rda.ft_pr.0, scope = formula(rda.ft_prt), direction = "both", perm.max = 9999, R2scope = TRUE) #automatic stepwise model selection based on R2
# rda(formula = L_ft ~ Condition(elevation + heatload2 + fog_freq + trans_rock_soil +
#                                  soil_depth + PC1s + PC2s + PC3s) + rda.w01.0$CA$u[, 2] + rda.w01.0$CA$u[, 5] + rda.w01.0$CA$u[, 4]
#     + rda.w01.0$CA$u[, 10], data = e) #final model includes axes 2, 5, 4 & 10

#RDA (& varpart) on presence/absence T (trait) matrices ----
#no need for transformation

#woody species
#note that CWM's were Z-transformed AND that d15N (and LA) were transformed before Z-transformation!

T_w01s <- scale (T_w01)

rda.Tw01.0 <- rda(T_w01s[, traits_w01] ~ 1, e)                                         #intercept model
rda.Tw01.1 <- rda(T_w01s[, traits_w01] ~ 1 + elevation + heatload2 + fog_freq + 
                    trans_rock_soil + soil_depth + PC1s + PC2s + PC3s, e)      #full model

#the COMBINED model
ordiR2step(rda.Tw01.0, scope = formula(rda.Tw01.1), direction = "forward", perm.max = 9999, R2scope = F) #automatic stepwise model selection based on R2


rda.Tw01 <- rda(T_w01[, traits_w01] ~ 1 + elevation + PC2s + PC3s + PC1s + fog_freq, e)   #final model 1 
summary(rda.Tw01)
plot(rda.Tw01)

RsquareAdj(rda.Tw01.1) 
anova.cca(rda.Tw01)
anova.cca(rda.Tw01, by = "axis", permutations = 9999)
anova.cca(rda.Tw01, by = "margin", permutations = 9999)

#variation partitioning on the FULL model (so NOT on the reduced model)
vpart.Tw01 <- varpart(T_w01[,traits_w01], ~ elevation + heatload2 + fog_freq, ~ trans_rock_soil + soil_depth + PC1s + PC2s + PC3s, data = e)
vpart.Tw01

p5_label <- paste('total~R^2 == ', round(vpart.Tw01$part$fract$Adj.R.squared[3]*100, 1))
p5 <- draw.pairwise.venn(
  round((vpart.Tw01$part$indfract$Adj.R.squared[1] + vpart.Tw01$part$indfract$Adj.R.squared[2])*100/ vpart.Tw01$part$fract$Adj.R.squared[3], digits = 1),
  round((vpart.Tw01$part$indfract$Adj.R.squared[3] + vpart.Tw01$part$indfract$Adj.R.squared[2])*100/ vpart.Tw01$part$fract$Adj.R.squared[3], digits = 1), 
  round(vpart.Tw01$part$indfract$Adj.R.squared[2]*100/ vpart.Tw01$part$fract$Adj.R.squared[3], digits = 1),
  c("climate proxies", "soil"), fill = c(viridis(3)[2], viridis(3)[3]), col = c(viridis(3)[2], viridis(3)[3]), alpha = c(0.6, 0.6), ind = FALSE, cat.pos = c(330, 30),
  scaled = T, cat.dist = -0.08,cat.cex = 1.5, inverted = T)
p5 <- ggdraw(p5) + 
  annotate(geom = 'text', x = 0.65, y = 0.02, label = p5_label, parse = TRUE, size = 5) +
  theme(plot.background = element_rect(fill = NA), plot.margin = margin(15, 15, 15, 15)) + theme(aspect.ratio=1/1)
p5 #!!! need manual fix (switch percentages in graph!!)

#ferns trait RDA

rda.Tft.0 <- rda(T_ft[, traits_ft] ~ 1, e)                                         #intercept model
rda.Tft.1 <- rda(T_ft[, traits_ft] ~ 1 + elevation + heatload2 + fog_freq + 
                   trans_rock_soil + soil_depth + PC1s + PC2s + PC3s, e)      #full model

#the COMBINED model
ordiR2step(rda.Tft.0, scope = formula(rda.Tft.1), direction = "forward", perm.max = 9999, R2scope = T) #automatic stepwise model selection based on R2
rda.Tft <- rda(T_ft[, traits_ft] ~ 1 + elevation + PC2s + PC3s + fog_freq, e)   #final model 1 
summary(rda.Tft)
plot(rda.Tft)

RsquareAdj(rda.Tft) #model adj R2 0.4221794
anova.cca(rda.Tft)
anova.cca(rda.Tft, by = "axis", permutations = 9999)
anova.cca(rda.Tft, by = "margin", permutations = 9999)

#variation partitioning on the FULL model (so NOT on the reduced model)

vpart.Tft <- varpart(T_ft[, traits_ft], ~ elevation + heatload2 + fog_freq, ~ trans_rock_soil + soil_depth + PC1s + PC2s + PC3s, data = e)
vpart.Tft

p4 <- draw.pairwise.venn(
  round((vpart.Tft$part$indfract$Adj.R.squared[1]+vpart.Tft$part$indfract$Adj.R.squared[2])*100/vpart.Tft$part$fract$Adj.R.squared[3], digits = 1),
  round((vpart.Tft$part$indfract$Adj.R.squared[3]+vpart.Tft$part$indfract$Adj.R.squared[2])*100/vpart.Tft$part$fract$Adj.R.squared[3], digits = 1),
  round(vpart.Tft$part$indfract$Adj.R.squared[2]*100/ vpart.Tft$part$fract$Adj.R.squared[3], digits = 1),
  c("topography", "soil"), fill = c(viridis(3)[2], viridis(3)[3]), col = c(viridis(3)[2], viridis(3)[3]), alpha = c(0.6, 0.6), ind = FALSE, cat.pos = c(330, 30),
  scaled = T, cat.dist = -0.08,cat.cex = 1.5, inverted = F)
p4_label <- paste('total~R^2 == ', round(vpart.Tft$part$fract$Adj.R.squared[3]*100, 1))
p4 <- ggdraw(p4) + 
  annotate(geom = 'text', x = 0.65, y = 0.02, label = p4_label, parse = TRUE, size = 5) +
  theme(plot.background = element_rect(fill = NA), plot.margin = margin(15, 15, 15, 15)) + theme(aspect.ratio=1/1)
p4

#explore if woody trait composition has additional effect on fern trait composition (residuals of RDA) ----
#environmental vars as 'conditional' = partialling effect out
summary(rda.Tw01.0) #first 3 axes contain most trait information
sum(rda.Tw01.0$CA$eig)/length(rda.Tw01.0$CA$eig) #need to include all axes with eig > 1  --> PC1 to 3! (this is 81.48% of total variation)

rda.Tft_pr.0 <- rda(T_ft[, traits_ft] ~ 1 + Condition(elevation + heatload2 + fog_freq + trans_rock_soil + soil_depth + PC1s + PC2s + PC3s), e) 
rda.Tft_prt <- rda(T_ft[, traits_ft] ~ 1 + rda.Tw01.0$CA$u[,1] + rda.Tw01.0$CA$u[,2] + rda.Tw01.0$CA$u[,3] + Condition(elevation + heatload2 + fog_freq + trans_rock_soil + soil_depth + PC1s + PC2s + PC3s), e) #NOT including tree density!!

anova.cca(rda.Tft_prt, permutations = 9999) #significant! F=2.4, p=0.0109
RsquareAdj(rda.Tft_prt) #model adj R2 0.04475993
anova.cca(rda.Tft_prt, by = "axis", permutations = 9999)
anova.cca(rda.Tft_prt, by = "margin", permutations = 99)
plot(rda.Tft_prt)

# ordiR2step(rda.Tft_pr.0, scope = formula(rda.Tft_prt), direction = "forward", perm.max = 9999, R2scope = TRUE) #automatic stepwise model selection based on R2
# rda(formula = T_ft[, traits_ft] ~ 1 + Condition(elevation + heatload2 + fog_freq +
#                                                   trans_rock_soil + soil_depth + PC1s + PC2s + PC3s), data = e)


# Combining all venn diagramsto 2 composite figure ----

tiff(file = "venn diagrams.tiff", width = 18, height = 18, units = "cm", res = 300,compression = "lzw")
plot_grid(p2, p1, p5, p4, labels = c('A overstory species', 'B understory species', 'C overstory traits', 'D understory traits'), ncol = 2, rel_heights=c(1,1))
dev.off()

rm(p1, p2, p4, p5)

#exploring species richness & functional diversity link environment ----
#standardize all predictors
e_sc = as.data.frame(scale(e[,2:38], center = T, scale = T))

# woody species (presabs) ----
MS_w <- glm(FD_w01$S ~ 1 + e_sc$elevation + e_sc$heatload2 + e_sc$fog_freq +  e_sc$trans_rock_soil + e_sc$soil_depth + e_sc$PC1s + e_sc$PC2s + e_sc$PC3s, family = poisson(link = "log"), na.action = na.fail) #NOTE: S is count variable, so use poisson model

dredge(MS_w) 
MS_w2 <- glm(FD_w01$S ~ 1 + e_sc$elevation + e_sc$fog_freq + e_sc$heatload2, family = poisson(link = "log"), na.action = na.fail) 
summary(MS_w2)
car::Anova(MS_w2, type = "II", test='F') #Type II because no interaction terms present

rsq(MS_w2, adj = TRUE, type = 'v') #R2 = 0.3649562
#type 'v' is the (default) - variance-function-based (Zhang, 2016)

#Testing of assumptions & goodness of fit
#for Poisson distribution you do not need normal distribution of residuals, but you need to check for overdispersion
DHARMa::testDispersion(MS_w) #this tests if you have significant overdispersion (problem), or not. The p-value is > 0.05, thus NO overdispersion,.
simulationOutput <- simulateResiduals(fittedModel = MS_w, plot = T)
testResiduals(simulationOutput) #all assumption tests OK (non-significant)
plotResiduals(simulationOutput, e$elevation) #OK
plotResiduals(simulationOutput, e$heatload) #OK
plotResiduals(simulationOutput, e$fog_freq) #OK

#variation partitioning. Method for binomial GLM/GAMM: https://rdrr.io/cran/ecospat/man/ecospat.varpart.html#heading-6
MS_w_topo <- glm(FD_w01$S ~ 1 + e_sc$elevation + e_sc$heatload2 + e_sc$fog_freq, family = poisson(link = "log"), na.action = na.fail) #NOTE: S is count variable, so use poisson model
MS_w_soil <- glm(FD_w01$S ~ 1 + e_sc$trans_rock_soil + e_sc$soil_depth + e_sc$PC1s + e_sc$PC2s + e_sc$PC3s, family = poisson(link = "log"), na.action = na.fail) #NOTE: S is count variable, so use poisson model

A_FS <- rsq(MS_w_topo, adj = TRUE, type = 'v')
B_FS <- rsq(MS_w_soil, adj = TRUE, type = 'v')
T_FS <- rsq(MS_w, adj = TRUE, type = 'v')
T_FS <- A_FS

p1_label <- paste('total~R^2 == ', round(T_FS*100, 1))
p1 <- draw.pairwise.venn(
  round((A_FS*100/ T_FS), digits = 1),
  round((B_FS*100/ T_FS), digits = 1),
  round(((A_FS+B_FS-T_FS)*100/ T_FS), digits = 1),
  c("c", "s"), fill = c(viridis(3)[2], viridis(3)[3]), col = c(viridis(3)[2], viridis(3)[3]), alpha = c(0.6, 0.6), ind = F, cat.pos = c(330, 30),
  scaled = T, cat.dist = -0.08,cat.cex = 1.5)
p1 <- ggdraw(p1) + 
  annotate(geom = 'text', x = 0.65, y = 0.02, label = p1_label, parse = TRUE, size = 5) +
  theme(plot.background = element_rect(fill = NA), plot.margin = margin(15, 15, 15, 15))+ theme(aspect.ratio=1/1)
p1

# Functional magnitude/ dispersion (M')

MFmag_w <- glm(FD_w01$M.prime_g ~ 1 + e_sc$elevation + e_sc$heatload2 + e_sc$fog_freq  + e_sc$trans_rock_soil + e_sc$soil_depth + e_sc$PC1s + e_sc$PC2s + e_sc$PC3s, family = Gamma(link = "inverse"), na.action = na.fail) #NOTE: M.prime is continuous positive, so use gamma model

dredge(MFmag_w)
MFmag_w2 <- glm(FD_w01$M.prime_g ~ 1 + e_sc$elevation + e_sc$heatload2, family = Gamma(link = "inverse"), na.action = na.fail)
summary(MFmag_w2)
car::Anova(MFmag_w2, type = "II", test='F') #Type II because no interaction terms present
rsq(MFmag_w2, adj = TRUE, type = 'v') #R2 = 0.6418392
rsq.partial(MFmag_w2, adj = TRUE ,type = 'v')

#Testing of assumptions & goodness of fit
#for gamma distribution you do not need normal distribution of residuals, but you need to check for overdispersion
DHARMa::testDispersion(MFmag_w) #this tests if you have significant overdispersion (problem), or not. The p-value is > 0.05, thus NO overdispersion.
simulationOutput <- simulateResiduals(fittedModel = MFmag_w, plot = T)
testResiduals(simulationOutput) #all assumption tests OK (non-significant)
plotResiduals(simulationOutput, e$elevation) #OK

#variation partitioning. Method for binomial GLM/GAMM: https://rdrr.io/cran/ecospat/man/ecospat.varpart.html#heading-6
MFmag_w_topo <- glm(FD_w01$M.prime_g ~ 1 + e_sc$elevation + e_sc$heatload2 + e_sc$fog_freq, family = Gamma(link = "inverse"), na.action = na.fail) #NOTE: S is count variable, so use poisson model
MFmag_w_soil <- glm(FD_w01$M.prime_g ~ 1 + e_sc$trans_rock_soil + e_sc$soil_depth + e_sc$PC1s + e_sc$PC2s + e_sc$PC3s, family = Gamma(link = "inverse"), na.action = na.fail) #NOTE: S is count variable, so use poisson model

A_Fmag <- rsq(MFmag_w_topo, adj = TRUE, type = 'v')
B_Fmag <- rsq(MFmag_w_soil, adj = TRUE, type = 'v')
T_Fmag <- rsq(MFmag_w, adj = TRUE, type = 'v')
T_Fmag <- A_Fmag

p2_label <- paste('total~R^2 == ', round(T_Fmag*100, 1))
p2 <- draw.pairwise.venn(
  round((A_Fmag*100/ T_Fmag), digits = 1),
  round((B_Fmag*100/ T_Fmag), digits = 1),
  round(((A_Fmag+B_Fmag-T_Fmag)*100/ T_Fmag), digits = 1),
  c("c", "s"), fill = c(viridis(3)[2], viridis(3)[3]), col = c(viridis(3)[2], viridis(3)[3]), alpha = c(0.6, 0.6), ind = F, cat.pos = c(330, 30),
  scaled = T, cat.dist = -0.08,cat.cex = 1.5)
p2 <- ggdraw(p2) + 
  annotate(geom = 'text', x = 0.65, y = 0.02, label = p2_label, parse = TRUE, size = 5) +
  theme(plot.background = element_rect(fill = NA), plot.margin = margin(15, 15, 15, 15))+ theme(aspect.ratio=1/1)
p2

# Functional variability/ equability (qET) 

MFvar_w <- glm(FD_w01$qEt_g ~ 1 + e_sc$elevation + e_sc$heatload2 + e_sc$fog_freq + I(e_sc$fog_freq**2) + e_sc$trans_rock_soil + e_sc$soil_depth + e_sc$PC1s + e_sc$PC2s + e_sc$PC3s, family = Gamma(link = "inverse"), na.action = na.fail)

dredge(MFvar_w)
MFvar_w2 <- glm(FD_w01$qEt_g ~ 1 + e_sc$elevation + e_sc$heatload2 + e_sc$fog_freq + I(e_sc$fog_freq**2) + e_sc$PC1s, family = Gamma(link = "inverse"), na.action = na.fail)
summary(MFvar_w2)
car::Anova(MFvar_w2, type = "II", test='F') #Type II because no interaction terms present
rsq(MFvar_w2, adj = TRUE, type = 'v') #R2 = 0.4847534
rsq.partial(MFvar_w2, adj = TRUE ,type = 'v')

#Testing of assumptions & goodness of fit
#for gamma distribution you do not need normal distribution of residuals, but you need to check for overdispersion
DHARMa::testDispersion(MFvar_w) #this tests if you have significant overdispersion (problem), or not. The p-value is > 0.05, thus NO overdispersion,.
simulationOutput <- simulateResiduals(fittedModel = MFvar_w, plot = T)
testResiduals(simulationOutput) #all assumption tests OK (non-significant)
plotResiduals(simulationOutput, e$elevation) #OK
plotResiduals(simulationOutput, e$heatload2) #more or less OK
plotResiduals(simulationOutput, e$fog_freq) #small deviation

#variation partitioning. Method for binomial GLM/GAMM: https://rdrr.io/cran/ecospat/man/ecospat.varpart.html#heading-6
MFvar_w_topo <- glm(FD_w01$qEt_g ~ 1 + e_sc$elevation + e_sc$heatload2 + e_sc$fog_freq + I(e_sc$fog_freq**2), family = Gamma(link = "inverse"), na.action = na.fail)
MFvar_w_soil <- glm(FD_w01$qEt_g ~ 1 + e_sc$trans_rock_soil + e_sc$soil_depth + e_sc$PC1s + e_sc$PC2s + e_sc$PC3s, family = Gamma(link = "inverse"), na.action = na.fail)

A_Fvar <- rsq(MFvar_w_topo, adj = TRUE, type = 'v')
B_Fvar <- rsq(MFvar_w_soil, adj = TRUE, type = 'v')
T_Fvar <- rsq(MFvar_w, adj = TRUE, type = 'v')

p3_label <- paste('total~R^2 == ', round(T_Fvar*100, 1))
p3 <- draw.pairwise.venn(
  round((A_Fvar*100/ T_Fvar), digits = 1),
  round((B_Fvar*100/ T_Fvar), digits = 1),
  round(((A_Fvar+B_Fvar-T_Fvar)*100/ T_Fvar), digits = 1),
  c("c", "s"), fill = c(viridis(3)[2], viridis(3)[3]), col = c(viridis(3)[2], viridis(3)[3]), alpha = c(0.6, 0.6), ind = F, cat.pos = c(330, 30),
  scaled = T, cat.dist = -0.08,cat.cex = 1.5)
p3 <- ggdraw(p3) + 
  annotate(geom = 'text', x = 0.65, y = 0.02, label = p3_label, parse = TRUE, size = 5) +
  theme(plot.background = element_rect(fill = NA), plot.margin = margin(15, 15, 15, 15))+ theme(aspect.ratio=1/1)
p3

# fern species (presabs) ----

e$PC2s2 <- (e$PC2s + 2)**2 #to satisfy assumptions, transformation needed for S model
e_sc$PC2s2 <- (e_sc$PC2s + 2)**2 #to satisfy assumptions, transformation needed for S model

MS_f <- glm(FD_ft$S ~ 1 + e_sc$elevation + I(e_sc$elevation**2) + e_sc$heatload2 + e_sc$trans_rock_soil +
              e_sc$soil_depth + e_sc$PC1s + e_sc$PC2s2 + e_sc$PC3s + e_sc$fog_freq, family = poisson(link = "log"), na.action = na.fail) #NOTE: S is count variable, so use poisson model

dredge(MS_f)
MS_f2 <- glm(FD_ft$S ~ 1 + e_sc$elevation + I(e_sc$elevation**2) + 
               e_sc$PC1s + e_sc$PC2s2 + e_sc$trans_rock_soil, family = poisson(link = "log"), na.action = na.fail) #NOTE: S is count variable, so use poisson model
# MS_f2 <- glm(FD_ft$S ~ 1 + e_sc$elevation + e_sc$fog_freq + 
#                e_sc$PC1s + e_sc$PC2s2 + e_sc$PC3s, family = poisson(link = "log"), na.action = na.fail) #NOTE: S is count variable, so use poisson model


summary(MS_f2)
car::Anova(MS_f2, type = "II", test='F') #Type II because no interaction terms present
rsq(MS_f2, adj = TRUE, type = 'v') #R2 = 0.6570462
rsq.partial(MS_f2, adj = TRUE ,type = 'v')

#Testing of assumptions & goodness of fit
#for Poisson distribution you do not need normal distribution of residuals, but you need to check for overdispersion
DHARMa::testDispersion(MS_f) #some small overdispersion when using Poisson dist --> solution: use quasi-poisson
simulationOutput <- simulateResiduals(fittedModel = MS_f, plot = T)
testResiduals(simulationOutput) #all assumption tests OK (non-significant)
plotResiduals(simulationOutput, e$elevation) #OK
plotResiduals(simulationOutput, e$PC2s2) #some deviation...
plotResiduals(simulationOutput, e$PC3s) #OK
plotResiduals(simulationOutput, e$PC1s) #OK
plotResiduals(simulationOutput, e$fog_freq) #some deviation...

#variation partitioning. Method for binomial GLM/GAMM: https://rdrr.io/cran/ecospat/man/ecospat.varpart.html#heading-6
MS_f_topo <- glm(FD_ft$S ~ 1 + e_sc$elevation + I(e_sc$elevation**2) + e_sc$heatload2 + e_sc$fog_freq, family = poisson(link = "log"), na.action = na.fail) #NOTE: S is count variable, so use poisson model
MS_f_soil <- glm(FD_ft$S ~ 1 + e_sc$trans_rock_soil + e_sc$soil_depth + e_sc$PC1s + e_sc$PC2s2 + e_sc$PC3s, family = poisson(link = "log"), na.action = na.fail) #NOTE: S is count variable, so use poisson model

A_FSf <- rsq(MS_f_topo, adj = TRUE, type = 'v')
B_FSf <- rsq(MS_f_soil, adj = TRUE, type = 'v')
T_FSf <- rsq(MS_f, adj = TRUE, type = 'v')

p4_label <- paste('total~R^2 == ', round(T_FSf*100, 1))
p4 <- draw.pairwise.venn(
  round((A_FSf*100/ T_FSf), digits = 1),
  round((B_FSf*100/ T_FSf), digits = 1),
  round(((A_FSf+B_FSf-T_FSf)*100/ T_FSf), digits = 1),
  c("c", "s"), fill = c(viridis(3)[2], viridis(3)[3]), col = c(viridis(3)[2], viridis(3)[3]), alpha = c(0.6, 0.6), ind = F, cat.pos = c(330, 30),
  scaled = T, cat.dist = -0.08,cat.cex = 1.5)
p4 <- ggdraw(p4) + 
  annotate(geom = 'text', x = 0.65, y = 0.02, label = p4_label, parse = TRUE, size = 5) +
  theme(plot.background = element_rect(fill = NA), plot.margin = margin(15, 15, 15, 15))+ theme(aspect.ratio=1/1)
p4

# Functional magnitude/ dispersion (M')

FD_ft$M.prime_g <- as.vector(t(naniar::replace_with_na_all(FD_ft$M.prime_g, ~.x < 0.001))) #first remove the 0's (1 species plots, N=3)
FD_ftredux1 <- cbind(e_sc[,c(1:6,8:38)], FD_ft)
FD_ftredux1 <- na.omit(FD_ftredux1)

MFmag_f <- glm(M.prime_g ~ 1 + elevation + heatload2 + trans_rock_soil +
                 soil_depth + PC1s + PC2s + PC3s + fog_freq, data = FD_ftredux1, family = Gamma(link = "inverse"), na.action = na.fail)
dredge(MFmag_f)
MFmag_f2 <- glm(M.prime_g ~ 1 + elevation + fog_freq, data = FD_ftredux1, family = Gamma(link = "inverse"), na.action = na.omit)
summary(MFmag_f2)
car::Anova(MFmag_f2, type="II", test='F') #Type II because no interaction terms present
rsq(MFmag_f2, adj=TRUE, type='v') #0.3166937

#Testing of assumptions & goodness of fit
#for gamma distribution you do not need normal distribution of residuals, but you need to check for overdispersion
DHARMa::testDispersion(MFmag_f) #some small overdispersion
simulationOutput <- simulateResiduals(fittedModel = MFmag_f, plot = T)
testResiduals(simulationOutput) #all assumption tests OK (non-significant)
plotResiduals(simulationOutput, FD_ftredux1$elevation) #OK
plotResiduals(simulationOutput, FD_ftredux1$PC1s) #OK
plotResiduals(simulationOutput, FD_ftredux1$PC2s) #OK

#variation partitioning. Method for binomial GLM/GAMM: https://rdrr.io/cran/ecospat/man/ecospat.varpart.html#heading-6
MFmag_f_topo <- glm(M.prime_g ~ 1 + elevation + heatload2 + fog_freq, data = FD_ftredux1, family = Gamma(link = "inverse"), na.action = na.fail)
MFmag_f_soil <- glm(M.prime_g ~ 1 + trans_rock_soil + soil_depth + PC1s + PC2s + PC3s, data = FD_ftredux1, family = Gamma(link = "inverse"), na.action = na.fail)


A_Fmagf <- rsq(MFmag_f_topo, adj = TRUE, type = 'v')
B_Fmagf <- rsq(MFmag_f_soil, adj = TRUE, type = 'v')
T_Fmagf <- rsq(MFmag_f, adj = TRUE, type = 'v')
T_Fmagf <- A_Fmagf

p5_label <- paste('total~R^2 == ', round(T_Fmagf*100, 1))
p5 <- draw.pairwise.venn(
  round((A_Fmagf*100/ T_Fmagf), digits = 1),
  round((B_Fmagf*100/ T_Fmagf), digits = 1),
  round(((A_Fmagf+B_Fmagf-T_Fmagf)*100/ T_Fmagf), digits = 1),
  c("c", "s"), fill = c(viridis(3)[2], viridis(3)[3]), col = c(viridis(3)[2], viridis(3)[3]), alpha = c(0.6, 0.6), ind = F, cat.pos = c(330, 30),
  scaled = T, cat.dist = -0.08,cat.cex = 1.5)
p5 <- ggdraw(p5) + 
  annotate(geom = 'text', x = 0.65, y = 0.02, label = p5_label, parse = TRUE, size = 5) +
  theme(plot.background = element_rect(fill = NA), plot.margin = margin(15, 15, 15, 15))+ theme(aspect.ratio=1/1)
p5


# Functional variability/ equability (qET)
FD_ft$qEt_g <- as.vector(t(naniar::replace_with_na_all(FD_ft$qEt_g, ~.x > 0.999))) #first remove the 1's (1 species AND 2 species plots, N=5)
FD_ftredux2 <- cbind(e_sc[,c(1:6,8:38)], FD_ft)
FD_ftredux2 <- na.omit(FD_ftredux2)


MFvar_f <- glm(qEt_g ~ 1 + elevation + heatload2 + trans_rock_soil +
                 soil_depth + PC1s + PC2s + PC3s + fog_freq , data = FD_ftredux2, family = Gamma(link = "inverse"), na.action = na.fail)
dredge(MFvar_f) 
MFvar_f2 <- glm(qEt_g ~ 1 + elevation + heatload2 + soil_depth + 
                 PC2s + PC3s, data = FD_ftredux2, family = Gamma(link = "inverse"), na.action = na.omit) 
summary(MFvar_f2)
car::Anova(MFvar_f2, type = "II", test='F') #Type II because no interaction terms present
rsq(MFvar_f2, adj = TRUE, type = 'v') #0.4646476

#Testing of assumptions & goodness of fit
#for gamma distribution you do not need normal distribution of residuals, but you need to check for overdispersion
DHARMa::testDispersion(MFvar_f) #this tests if you have significant overdispersion (problem), or not. The p-value is > 0.05, thus NO overdispersion,.
simulationOutput <- simulateResiduals(fittedModel = MFvar_f, plot = T)
testResiduals(simulationOutput) #all assumption tests OK (non-significant)
plotResiduals(simulationOutput, FD_ftredux2$elevation) #small deviation
plotResiduals(simulationOutput, FD_ftredux2$heatload2) #OK
plotResiduals(simulationOutput, FD_ftredux2$soil_depth) #slight deviation
plotResiduals(simulationOutput, FD_ftredux2$PC2s) #OK
plotResiduals(simulationOutput, FD_ftredux2$PC3s) #OK

#variation partitioning. Method for binomial GLM/GAMM: https://rdrr.io/cran/ecospat/man/ecospat.varpart.html#heading-6

MFvar_f_topo <- glm(qEt_g ~ 1 + elevation + heatload2 + fog_freq, data = FD_ftredux2, family = Gamma(link = "inverse"), na.action = na.fail) 
MFvar_f_soil <- glm(qEt_g ~ 1 + soil_depth + PC1s + PC2s + PC3s + soil_depth + trans_rock_soil, data = FD_ftredux2, family = Gamma(link = "inverse"), na.action = na.fail) 

A_Fvarf <- rsq(MFvar_f_topo, adj = TRUE, type = 'v')
B_Fvarf <- rsq(MFvar_f_soil, adj = TRUE, type = 'v')
T_Fvarf <- rsq(MFvar_f, adj = TRUE, type = 'v')

p6_label <- paste('total~R^2 == ', round(T_Fvarf*100, 1))
p6 <- draw.pairwise.venn(
  round((A_Fvarf*100/ T_Fvarf), digits = 1),
  round((B_Fvarf*100/ T_Fvarf), digits = 1),
  round(((A_Fvarf+B_Fvarf-T_Fvarf)*100/ T_Fvarf), digits = 1),
  c("c", "s"), fill = c(viridis(3)[2], viridis(3)[3]), col = c(viridis(3)[2], viridis(3)[3]), alpha = c(0.6, 0.6), ind = F, cat.pos = c(330, 30),
  scaled = T, cat.dist = -0.08,cat.cex = 1.5)
p6 <- ggdraw(p6) + 
  annotate(geom = 'text', x = 0.65, y = 0.02, label = p6_label, parse = TRUE, size = 5) +
  theme(plot.background = element_rect(fill = NA), plot.margin = margin(15, 15, 15, 15))+ theme(aspect.ratio=1/1)
p6

#combine VENN diagram plots ----
tiff(file = "venn diagrams FD.tiff", width = 27, height = 18, units = "cm", res = 300,compression = "lzw")
plot_grid(p1, p2, p3, p4, p5, p6, labels = c('A woody S', 'B woody M', 'C woody 1E(T)', 'D fern S', 'E fern M', 'F fern 1E(T)'), ncol = 3, rel_heights=c(1,1))
dev.off()

plot(e_sc$fog_freq, FD_w01$qEt_g)
plot(e_sc$elevation, FD_w01$qEt_g)

#make plots diversity - environment links ----

xx <- dplyr::rename(FD_ft, S_f = S, M.prime_gf = M.prime_g, qEt_gf = qEt_g)
trts <- c("S_f", "M.prime_gf", "qEt_gf")
dataplotFD <- cbind(xx[,trts], FD_w01, e)
rm(xx)

#"#66CC00" #light green
#"#4169E1" #light blue
#"#3C7700" #dark green (ferns)
#"#191970" #dark blue (trees)

sFD1 <- ggplot(dataplotFD, aes(x = elevation, y = S_f)) + geom_point(colour = "#3C7700", shape = 17) + ggtitle(label = "") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 15)) + 
  scale_x_continuous("") + scale_y_continuous("species richness", limits = c(1, 22)) + theme(legend.position = "none") +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black")) + 
  geom_point(data = dataplotFD, aes(y = S), colour = "#191970") +
  stat_smooth(aes(y = S_f), method = "lm", formula = y ~ x + I(x^2), size = 1, color = "#3C7700") +
  stat_smooth(aes(y = S), method = "lm", size = 1, color = "#191970")

sFD2 <- ggplot(dataplotFD, aes(x = elevation, y = M.prime_gf)) + geom_point(colour = "#3C7700", shape = 17) + ggtitle(label = "") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 15)) + 
  scale_x_continuous("") + scale_y_continuous("functional dispersion") + theme(legend.position = "none") +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black")) + 
  geom_point(data = dataplotFD, aes(y = M.prime_g), colour = "#191970") +
  stat_smooth(aes(y = M.prime_gf), method = "lm", size = 1, color = "#3C7700") +
  stat_smooth(aes(y = M.prime_g), method = "lm", size = 1, color = "#191970")

sFD3 <- ggplot(dataplotFD, aes(x = elevation, y = qEt_gf)) + geom_point(colour = "#3C7700", shape = 17) + ggtitle(label = "") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 15)) + 
  scale_x_continuous("elevation (m)") + scale_y_continuous("functional equability", breaks = c(0.90, 0.92, 0.94, 0.96, 0.98)) + theme(legend.position = "none") +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black")) + 
  geom_point(data = dataplotFD, aes(y = qEt_g), colour = "#191970") +
  stat_smooth(aes(y = qEt_gf), method = "lm", size = 1, color = "#3C7700") +
  stat_smooth(aes(y = qEt_g), method = "lm", size = 1, color = "#191970")

sFD4 <- ggplot(dataplotFD, aes(x = PC2s, y = S_f)) + geom_point(colour = "#3C7700", shape = 17) + ggtitle(label = "") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 15)) + 
  scale_x_continuous("", breaks = c(-1, -0.5, 0, 0.5, 1, 1.5)) + scale_y_continuous("", limits = c(1, 22)) + theme(legend.position = "none") +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black")) + 
  geom_point(data = dataplotFD, aes(y = S),colour = "#191970") +
  stat_smooth(aes(y = S_f), method = "lm", size = 1, color = "#3C7700")
# + stat_smooth(aes(y = S), method = "lm", linetype = "dashed", se = FALSE, size = 1, color = "#191970")

sFD5 <- ggplot(dataplotFD, aes(x = PC2s, y = M.prime_gf)) + geom_point(colour = "#3C7700", shape = 17) + ggtitle(label = "")+
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 15)) + 
  scale_x_continuous("", breaks = c(-1, -0.5, 0, 0.5, 1, 1.5)) + scale_y_continuous("") + theme(legend.position = "none") +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black")) + 
  geom_point(data = dataplotFD, aes(y = M.prime_g),colour = "#191970") 
#  + stat_smooth(aes(y = M.prime_gf), method = "lm", linetype = "dashed", se = FALSE, size = 1, color = "#3C7700") +
#  stat_smooth(aes(y = M.prime_g), method = "lm", linetype = "dashed", se = FALSE, size = 1, color = "#191970")

sFD6 <- ggplot(dataplotFD, aes(x = PC2s, y = qEt_gf)) + geom_point(colour = "#3C7700", shape = 17) + ggtitle(label = "")+
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 15)) + 
  scale_x_continuous("soil pH (PC 2)", breaks = c(-1, -0.5, 0, 0.5, 1, 1.5)) + scale_y_continuous("", breaks = c(0.90, 0.92, 0.94, 0.96, 0.98)) + theme(legend.position = "none") +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black")) + 
  geom_point(data = dataplotFD, aes(y = qEt_g),colour = "#191970") +
  stat_smooth(aes(y = qEt_gf), method = "lm", size = 1, color = "#3C7700") 
# + stat_smooth(aes(y = qEt_g), method = "lm", linetype = "dashed", se = FALSE, size = 1, color = "#191970")

tiff(file = "div env.tiff", width = 18, height = 25, units = "cm", res = 300,compression = "lzw")
plot_grid(sFD1, sFD4, sFD2, sFD5, sFD3, sFD6, ncol = 2, labels="AUTO")
dev.off()

#explore diversity - diversity correlations (incl plots) ----

# cor.test(dataplotFD$S, dataplotFD$S_f, method="pearson")
# cor.test(dataplotFD$M.prime_g, dataplotFD$M.prime_gf, method="pearson")
# cor.test(dataplotFD$qEt_g, dataplotFD$qEt_gf, method="pearson")

modelS <- glm(S_f ~ 1 + S, data = dataplotFD, na.action = na.omit)
car::Anova(modelS, type = "II", test='F') 
rsq(modelS, adj = TRUE, type = 'v') 

modelM <- glm(M.prime_gf ~ 1 + M.prime_g, data = dataplotFD, na.action = na.omit)
car::Anova(modelM, type = "II", test='F') 
rsq(modelM, adj = TRUE, type = 'v') 

modelEt <- glm(qEt_gf ~ 1 + qEt_g, data = dataplotFD, na.action = na.omit)
car::Anova(modelEt, type = "II", test='F') 
rsq(modelEt, adj = TRUE, type = 'v') 

sDiv1 <- ggplot(dataplotFD, aes(x = S, y = S_f, color = elevation)) + geom_point() + ggtitle(label = "") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradientn(colours = rev(viridis(10))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 13)) + 
  scale_x_continuous("overstory S") + scale_y_continuous("understory S") + theme(legend.position = "none") +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black"))

sDiv2 <- ggplot(dataplotFD, aes(x = M.prime_g, y = M.prime_gf, color = elevation)) + geom_point() + ggtitle(label = "") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradientn(colours = rev(viridis(10))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 13)) + 
  scale_x_continuous("overstory M'") + scale_y_continuous("understory M'") + theme(legend.position = "none") +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black")) + 
  geom_smooth(method = "lm", se = TRUE, color = "grey25")

sDiv3 <- ggplot(dataplotFD, aes(x = qEt_g, y = qEt_gf, color = elevation)) + geom_point() + ggtitle(label = "") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +  
  scale_color_gradientn(colours = rev(viridis(10))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 13)) + 
  scale_x_continuous(bquote("overstory "^1*"E(T)")) + scale_y_continuous(bquote("understory "^1*"E(T)"), breaks = c(0.90, 0.92, 0.94, 0.96, 0.98)) + theme(legend.position = "right") +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black")) + 
  geom_smooth(method = "lm", linetype = "dashed", se = F, color = "grey25")
theme(legend.title = element_text(size = 12), legend.text = element_text(size = 8))

tiff(file = "FD FD2.tiff", width = 25, height = 7.5, units = "cm", res = 300,compression = "lzw")
plot_grid(sDiv1, sDiv2, sDiv3, ncol = 3,  axis = "rlbt", labels="AUTO", rel_widths = c(1,1,1.333))
dev.off()

#Explore trait-trait correlations ----

MSLA <- glm(T_ft$SLA ~ 1 + T_w01$SLA + I(T_w01$SLA^2), na.action = na.fail) #for SLA we have x2 relationship!
summary(MSLA)
car::Anova(MSLA, type = "II", test = 'F') 
rsq(MSLA, adj = TRUE, type = 'v') 

MLDMC <- glm(T_ft$LDMC ~ 1 + T_w01$LDMC, na.action = na.fail) 
car::Anova(MLDMC, type = "II", test = 'F') 
rsq(MLDMC, adj = TRUE, type = 'v') 

Mlog_LA <- glm(T_ft$log_LA ~ 1 + T_w01$log_LA, na.action = na.fail) 
car::Anova(Mlog_LA, type = "II", test = 'F') 
rsq(Mlog_LA, adj = TRUE, type = 'v') 

MLth <- glm(T_ft$Lth ~ 1 + T_w01$Lth, na.action = na.fail) 
car::Anova(MLth, type = "II", test = 'F') 
rsq(MLth, adj = TRUE, type = 'v') 

MChl_SPAD <- glm(T_ft$Chl_SPAD ~ 1 + T_w01$Chl_SPAD, na.action = na.fail) 
car::Anova(MChl_SPAD, type = "II", test = 'F') 
rsq(MChl_SPAD, adj = TRUE, type = 'v') 

MEWT <- glm(T_ft$EWT ~ 1 + T_w01$EWT, na.action = na.fail) 
car::Anova(MEWT, type = "II", test = 'F') 
rsq(MEWT, adj = TRUE, type = 'v') 

Md13C <- glm(T_ft$d13C ~ 1 + T_w01$d13C, na.action = na.fail) 
car::Anova(Md13C, type = "II", test = 'F') 
rsq(Md13C, adj = TRUE, type = 'v') 

#Md15N <- glm(T_ft$log_d15N ~ 1 + T_w01$d15N2, na.action = na.fail) 
#car::Anova(Md15N, type = "II", test = 'F') 
#rsq(Md15N, adj = TRUE, type = 'v') 

Mleaf_N <- glm(T_ft$leaf_N ~ 1 + T_w01$leaf_N, na.action = na.fail) 
car::Anova(Mleaf_N, type = "II", test = 'F') 
rsq(Mleaf_N, adj = TRUE, type = 'v') 

# #for Model II MA regression, you CANNOT use Z-transformed trait data! 
# xx <- dplyr::rename(T_ftb, SLA_f = SLA, LDMC_f = LDMC, log_LA_f = log_LA, Lth_f = Lth, Chl_SPAD_f = Chl_SPAD,
#                     EWT_f = EWT, d13C_f = d13C, log_d15N_f = log_d15N, leaf_N_f = leaf_N, d15N_f = d15N)
# yy <- dplyr::rename(T_w01b, SLA_t = SLA, LDMC_t = LDMC, log_LA_t = log_LA, Lth_t = Lth, Chl_SPAD_t = Chl_SPAD,
#                     EWT_t = EWT, d13C_t = d13C, d15N2_t = d15N2, leaf_N_t = leaf_N, d15N_t = d15N)
# combinedT <- cbind(yy, xx)
# rm(xx, yy)
# 
# lmodel2(formula =  LDMC_t ~ LDMC_f, data = combinedT, range.x = "interval", range.y = "interval", nperm = 999)
# lmodel2(formula =  LDMC_t ~ LDMC_f, data = combinedT, range.x = "relative", range.y = "relative", nperm = 999)
# 
# lmodel2(formula =  Chl_SPAD_t ~ Chl_SPAD_f, data = combinedT, range.x = "relative", range.y = "relative", nperm = 999)
# 
# lmodel2(formula =  SLA_t ~ SLA_f + SLA_f^2, data = combinedT, range.x = "relative", range.y = "relative", nperm = 999)

# graph pairwise trait - trait correlations ----

xx <- dplyr::rename(T_ftgraph, SLA_f = SLA, LDMC_f = LDMC, log_LA_f = log_LA, Lth_f = Lth, Chl_SPAD_f = Chl_SPAD,
                    EWT_f = EWT, d13C_f = d13C, log_d15N_f = log_d15N, leaf_N_f = leaf_N, d15N_f = d15N)
dataplot <- cbind(T_w01graph, xx, e)
rm(xx)

s1 <- ggplot(dataplot, aes(x = SLA, y = SLA_f, color = elevation)) + geom_point() + 
  geom_smooth(method = lm, formula = y ~ x + I(x^2), colour = "grey25") + 
  ggtitle(label = bquote('SLA (mm'^2*'/mg)')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradientn(colours = rev(viridis(10))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 13)) + scale_x_continuous("") + scale_y_continuous("understory CM trait") + 
  theme(legend.position = "none") + theme(aspect.ratio=1) +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black"))
s2 <- ggplot(dataplot, aes(x = LDMC, y = LDMC_f, color = elevation)) + geom_point() + 
  geom_smooth(method = lm, colour = "grey25") + 
  ggtitle(label = bquote('LDMC (mg/g)')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradientn(colours = rev(viridis(10))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 13)) + scale_x_continuous("", breaks = c(375, 400, 425, 450, 470)) + 
  scale_y_continuous("") + 
  theme(legend.position = "none") + theme(aspect.ratio=1) +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black"))
s3 <- ggplot(dataplot, aes(x = (exp(log_LA))/100, y = (exp(log_LA_f))/100, color = elevation)) + geom_point() + 
  ggtitle(label = bquote('leaf area (cm'^2*')')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradientn(colours = rev(viridis(10))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 13)) + scale_x_continuous("", trans = "log10", breaks = c(8.5, 10, 12.5, 15, 20, 25)) + 
  scale_y_continuous("", trans = "log10", breaks = c(150, 200, 300, 500, 800)) + 
  theme(legend.position = "none") + theme(aspect.ratio=1) +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black"))
s4 <- ggplot(dataplot, aes(x = Lth, y = Lth_f, color = elevation)) + geom_point() + 
  geom_smooth(method = lm, colour = "grey25") + 
  ggtitle(label = bquote('leaf thickness (mm)')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradientn(colours = rev(viridis(10))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 13)) + scale_x_continuous("") + scale_y_continuous("understory CM trait") + 
  theme(legend.position = "none") + theme(aspect.ratio=1) +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black"))
s5 <- ggplot(dataplot, aes(x = Chl_SPAD, y = Chl_SPAD_f, color = elevation)) + geom_point() + 
  ggtitle(label = bquote('leaf chlorophyll content')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradientn(colours = rev(viridis(10))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 13)) + scale_x_continuous("") + scale_y_continuous("") + 
  theme(legend.position = "none") + theme(aspect.ratio=1) +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black"))
s6 <- ggplot(dataplot, aes(x = EWT, y = EWT_f, color = elevation)) + geom_point() + 
  ggtitle(label = bquote('EWT (mg/mm'^2*')')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradientn(colours = rev(viridis(10))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 13)) + scale_x_continuous("", breaks = c(0.11, 0.13, 0.15, 0.17)) + 
  scale_y_continuous("") +  theme(aspect.ratio=1) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black"))
s7 <- ggplot(dataplot, aes(x = d13C, y = d13C_f, color = elevation)) + geom_point() + 
  geom_smooth(method = lm, colour = "grey25", linetype = "dashed", se = F) + 
  ggtitle(label = expression(delta^13*'C (‰)')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradientn(colours = rev(viridis(10))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 13)) + scale_x_continuous("overstory CM trait", breaks = c(-31.5, -31, -30.5, -30.1)) + 
  scale_y_continuous("understory CM trait", breaks = c(-32, -31.5, -31, -30.5, -30.1)) + 
  theme(legend.position = "none") + theme(aspect.ratio=1) +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black"))
s8 <- ggplot(dataplot, aes(x = d15N2, y = d15N_f+2.5, color = elevation)) + geom_point() + 
  geom_smooth(method = lm, colour = "grey25") + 
  ggtitle(label = expression(delta^15*'N (‰)')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradientn(colours = rev(viridis(10))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 13)) +  theme(aspect.ratio=1) +
  scale_x_continuous("overstory CM trait", breaks = c(1, 2.25, 4, 6.25, 8.41), labels = c('-3.0', '-2.5', '-2.0', '-1.5', '-1.1')) + 
  scale_y_continuous("", trans = "log10", breaks = c(0.05, 0.1, 0.25, 0.50, 1, 2), labels = c('-2.45', '-2.40', '-2.25', '-2.00', '-1.50', '-0.50')) + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black"))
s9 <- ggplot(dataplot, aes(x = leaf_N, y = leaf_N_f, color = elevation)) + geom_point() + 
  geom_smooth(method = lm, colour = "grey25") + 
  ggtitle(label = bquote('leaf N (mg/g)')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradientn(colours = rev(viridis(10))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 13)) + 
  scale_x_continuous("overstory CM trait") + scale_y_continuous("", breaks = c(20, 22, 24, 26, 28,30)) + 
  theme(legend.position = "none") + theme(aspect.ratio=1) +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black")) +
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 8))
s9b <- ggplot(dataplot, aes(x = leaf_N, y = leaf_N_f, color = elevation)) + geom_point() + 
  geom_smooth(method = lm, colour = "grey25") + 
  ggtitle(label = bquote('leaf N (mg/g)')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradientn(colours = rev(viridis(10))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 13)) + 
  scale_x_continuous("overstory CM trait") + scale_y_continuous("", breaks = c(20, 22, 24, 26, 28,30)) + 
  theme(legend.position = "right") + theme(aspect.ratio=1) +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black")) +
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 8))

tiff(file = "CWM trait trait.tiff", width = 25, height = 25, units = "cm", res = 300,compression = "lzw")
plot_grid(s1, s2, s3, s4, s5, s6,s7, s8, s9, ncol = 3, axis = "rlbt", labels="AUTO")
dev.off()

#Explore impact of deciduous species (trait RDA excluding deciduous species) ----

t_w2 <- merge(t_w, terr_list[,c("species_abbrev","leaf_habit")], by.x = "row.names", by.y = "species_abbrev")
row.names(t_w2) <- t_w2[,1]
t_w2 <- t_w2[,-1]

bluh <- as.data.frame(t(L_w.t))
bluh[bluh>0] <- 1
bluh2 <- merge(bluh, terr_list[,c("species_abbrev","leaf_habit")], by.x = "row.names", by.y = "species_abbrev")
bluh2 <- bluh2[,-1]

evergreen <- bluh2[bluh2$leaf_habit!='DB',]
t_w_ev  <- t_w2[t_w2$leaf_habit!='DB',]
t_w_ev$log_LA <- log(t_w_ev$LA) 

#filter trait dataframes to include only those species present in the L dataframes. 
t_w_ev <- t_w_ev[rownames(t_w_ev) %in% colnames(L_w), ] #select only traits for woody species that occur in the 100m2 plots
rownames(t_w_ev) <- sub("_T", "", rownames(t_w_ev))

#woody species CWMs (presabs)
colnames(L_w) <- sub("_T", "", colnames(L_w))
L_w_ev.t <- L_w[, colnames(L_w) %in% rownames(t_w_ev)] #select only species for which traits are measured (needed for FD function to work)
L_w_ev.t <- L_w_ev.t[, order(colnames(L_w_ev.t))] 
t_w_ev <- t_w_ev[order(rownames(t_w_ev)), ] 

x <- dbFD(t_w_ev[, traits], L_w_ev.t, w.abun = FALSE, stand.x = FALSE, calc.CWM = T, 
          calc.FDiv = F, calc.FGR = F, calc.FRic = F, corr='none')
T_w_ev01 <- x$CWM

#evaluate CWM values
hist(T_w_ev01$LDMC)
hist(T_w_ev01$Lth)
hist(T_w_ev01$Chl_SPAD)
hist(T_w_ev01$log_LA) 
hist(T_w_ev01$SLA)
hist(T_w_ev01$EWT)
hist(T_w_ev01$d13C)
hist((T_w_ev01$d15N+4)**2) #transform!!!
hist(T_w_ev01$leaf_N)

T_w_ev01$d15N2 <- (T_w_ev01$d15N + 4)**2

#T_w_ev01b <- T_w_ev01 #for making PCA of merged trait datasets
T_w_ev01graph <- T_w_ev01 #retain unstanderdized dataset for graphing
T_w_ev01 <- as.data.frame(scale(T_w_ev01)) #standardize all data! (mean=0, sd=1)

rda.Tw_ev01.0 <- rda(T_w_ev01[, traits_w01] ~ 1, e)                                         #intercept model
rda.Tw_ev01.1 <- rda(T_w_ev01[, traits_w01] ~ 1 + elevation + heatload2 + fog_freq + 
                       trans_rock_soil + soil_depth + PC1s + PC2s + PC3s, e)      #full model

#the COMBINED model
ordiR2step(rda.Tw_ev01.0, scope = formula(rda.Tw_ev01.1), direction = "forward", perm.max = 9999, R2scope = T) #automatic stepwise model selection based on R2
rda.Tw_ev01 <- rda(T_w_ev01[, traits_w01] ~ 1 + elevation + PC2s + PC3s + PC1s + fog_freq, e)   #final model 1 
summary(rda.Tw_ev01)
plot(rda.Tw_ev01)

RsquareAdj(rda.Tw_ev01) #model adj R2 0.5092502
anova.cca(rda.Tw_ev01)
anova.cca(rda.Tw_ev01, by = "axis", permutations = 9999)
anova.cca(rda.Tw_ev01, by = "margin", permutations = 9999)

vpart.Tw_ev01 <- varpart(T_w_ev01[,traits_w01], ~ elevation + heatload2 + fog_freq, ~ trans_rock_soil + soil_depth + PC1s + PC2s + PC3s, data = e)
vpart.Tw_ev01

p5_label <- paste('total~R^2 == ', round(vpart.Tw_ev01$part$fract$Adj.R.squared[3]*100, 1))
p5_ev <- draw.pairwise.venn(
  round((vpart.Tw_ev01$part$indfract$Adj.R.squared[1] + vpart.Tw_ev01$part$indfract$Adj.R.squared[2])*100/ vpart.Tw_ev01$part$fract$Adj.R.squared[3], digits = 1),
  round((vpart.Tw_ev01$part$indfract$Adj.R.squared[3] + vpart.Tw_ev01$part$indfract$Adj.R.squared[2])*100/ vpart.Tw_ev01$part$fract$Adj.R.squared[3], digits = 1), 
  round(vpart.Tw_ev01$part$indfract$Adj.R.squared[2]*100/ vpart.Tw_ev01$part$fract$Adj.R.squared[3], digits = 1),
  c("topography", "soil"), fill = c(viridis(3)[2], viridis(3)[3]), col = c(viridis(3)[2], viridis(3)[3]), alpha = c(0.6, 0.6), ind = FALSE, cat.pos = c(330, 30),
  scaled = T, cat.dist = -0.08,cat.cex = 1.5, inverted = F)
p5_ev <- ggdraw(p5_ev) + 
  annotate(geom = 'text', x = 0.65, y = 0.02, label = p5_label, parse = TRUE, size = 5) +
  theme(plot.background = element_rect(fill = NA), plot.margin = margin(15, 15, 15, 15)) + theme(aspect.ratio=1/1)
p5_ev #

#pairwise Pearson correlations ----

MSLA <- glm(T_ft$SLA ~ 1 + T_w_ev01$SLA + I(T_w_ev01$SLA^2), na.action = na.fail) #for SLA we have x2 relationship!
car::Anova(MSLA, type = "II", test = 'F') 
rsq(MSLA, adj = TRUE, type = 'v') 

MLDMC <- glm(T_ft$LDMC ~ 1 + T_w_ev01$LDMC, na.action = na.fail) 
car::Anova(MLDMC, type = "II", test = 'F') 
rsq(MLDMC, adj = TRUE, type = 'v') 

Mlog_LA <- glm(T_ft$log_LA ~ 1 + T_w_ev01$log_LA, na.action = na.fail) 
car::Anova(Mlog_LA, type = "II", test = 'F') 
rsq(Mlog_LA, adj = TRUE, type = 'v') 

MLth <- glm(T_ft$Lth ~ 1 + T_w_ev01$Lth, na.action = na.fail) 
car::Anova(MLth, type = "II", test = 'F') 
rsq(MLth, adj = TRUE, type = 'v') 

MChl_SPAD <- glm(T_ft$Chl_SPAD ~ 1 + T_w_ev01$Chl_SPAD, na.action = na.fail) 
car::Anova(MChl_SPAD, type = "II", test = 'F') 
rsq(MChl_SPAD, adj = TRUE, type = 'v') 

MEWT <- glm(T_ft$EWT ~ 1 + T_w_ev01$EWT, na.action = na.fail) 
car::Anova(MEWT, type = "II", test = 'F') 
rsq(MEWT, adj = TRUE, type = 'v') 

Md13C <- glm(T_ft$d13C ~ 1 + T_w_ev01$d13C, na.action = na.fail) 
car::Anova(Md13C, type = "II", test = 'F') 
rsq(Md13C, adj = TRUE, type = 'v') 

Md15N <- glm(T_ft$log_d15N ~ 1 + T_w_ev01$d15N2, na.action = na.fail) 
car::Anova(Md15N, type = "II", test = 'F') 
rsq(Md15N, adj = TRUE, type = 'v') 

Mleaf_N <- glm(T_ft$leaf_N ~ 1 + T_w_ev01$leaf_N, na.action = na.fail) 
car::Anova(Mleaf_N, type = "II", test = 'F') 
rsq(Mleaf_N, adj = TRUE, type = 'v') 

# graph pairwise trait correlations ----

xx <- dplyr::rename(T_ftgraph, SLA_f = SLA, LDMC_f = LDMC, log_LA_f = log_LA, Lth_f = Lth, Chl_SPAD_f = Chl_SPAD,
                    EWT_f = EWT, d13C_f = d13C, log_d15N_f = log_d15N, leaf_N_f = leaf_N, d15N_f = d15N)
dataplot2 <- cbind(T_w_ev01graph, xx, e)
rm(xx)

s1 <- ggplot(dataplot2, aes(x = SLA, y = SLA_f, color = elevation)) + geom_point() + 
  geom_smooth(method = lm, formula = y ~ x + I(x^2), colour = "grey25") + 
  ggtitle(label = bquote('SLA (mm'^2*'/mg)')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradientn(colours = rev(viridis(10))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 13)) + scale_x_continuous("") + scale_y_continuous("understory CM trait") + 
  theme(legend.position = "none") + theme(aspect.ratio=1) +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black"))
s2 <- ggplot(dataplot2, aes(x = LDMC, y = LDMC_f, color = elevation)) + geom_point() + 
  geom_smooth(method = lm, colour = "grey25") + 
  ggtitle(label = bquote('LDMC (mg/g)')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradientn(colours = rev(viridis(10))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 13)) + scale_x_continuous("", breaks = c(375, 400, 425, 450, 470)) + 
  scale_y_continuous("") + 
  theme(legend.position = "none") + theme(aspect.ratio=1) +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black"))
s3 <- ggplot(dataplot2, aes(x = (exp(log_LA))/100, y = (exp(log_LA_f))/100, color = elevation)) + geom_point() + 
  ggtitle(label = bquote('leaf area (cm'^2*')')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradientn(colours = rev(viridis(10))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 13)) + scale_x_continuous("", trans = "log10", breaks = c(8.5, 10, 12.5, 15, 20, 25)) + 
  scale_y_continuous("", trans = "log10", breaks = c(150, 200, 300, 500, 800)) + 
  theme(legend.position = "none") + theme(aspect.ratio=1) +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black"))
s4 <- ggplot(dataplot2, aes(x = Lth, y = Lth_f, color = elevation)) + geom_point() + 
  geom_smooth(method = lm, colour = "grey25") + 
  ggtitle(label = bquote('leaf thickness (mm)')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradientn(colours = rev(viridis(10))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 13)) + scale_x_continuous("") + scale_y_continuous("understory CM trait") + 
  theme(legend.position = "none") + theme(aspect.ratio=1) +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black"))
s5 <- ggplot(dataplot2, aes(x = Chl_SPAD, y = Chl_SPAD_f, color = elevation)) + geom_point() + 
  ggtitle(label = bquote('leaf chlorophyll content')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradientn(colours = rev(viridis(10))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 13)) + scale_x_continuous("") + scale_y_continuous("") + 
  theme(legend.position = "none") + theme(aspect.ratio=1) +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black"))
s6 <- ggplot(dataplot2, aes(x = EWT, y = EWT_f, color = elevation)) + geom_point() + 
  ggtitle(label = bquote('EWT (mg/mm'^2*')')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradientn(colours = rev(viridis(10))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 13)) + scale_x_continuous("", breaks = c(0.11, 0.13, 0.15, 0.17)) + 
  scale_y_continuous("") +  theme(aspect.ratio=1) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black"))
s7 <- ggplot(dataplot2, aes(x = d13C, y = d13C_f, color = elevation)) + geom_point() + 
  ggtitle(label = expression(delta^13*'C (‰)')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradientn(colours = rev(viridis(10))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 13)) + scale_x_continuous("overstory CM trait", breaks = c(-31.5, -31, -30.5, -30.1)) + 
  scale_y_continuous("understory CM trait", breaks = c(-32, -31.5, -31, -30.5, -30.1)) + 
  theme(legend.position = "none") + theme(aspect.ratio=1) +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black"))
s8 <- ggplot(dataplot2, aes(x = d15N2, y = d15N_f+2.5, color = elevation)) + geom_point() + 
  geom_smooth(method = lm, colour = "grey25") + 
  ggtitle(label = expression(delta^15*'N (‰)')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradientn(colours = rev(viridis(10))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 13)) +  theme(aspect.ratio=1) +
  scale_x_continuous("overstory CM trait", breaks = c(1, 2.25, 4, 6.25, 8.41), labels = c('-3.0', '-2.5', '-2.0', '-1.5', '-1.1')) + 
  scale_y_continuous("", trans = "log10", breaks = c(0.05, 0.1, 0.25, 0.50, 1, 2), labels = c('-2.45', '-2.40', '-2.25', '-2.00', '-1.50', '-0.50')) + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black"))
s9 <- ggplot(dataplot2, aes(x = leaf_N, y = leaf_N_f, color = elevation)) + geom_point() + 
  geom_smooth(method = lm, colour = "grey25") + 
  ggtitle(label = bquote('leaf N (mg/g)')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradientn(colours = rev(viridis(10))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 13)) + 
  scale_x_continuous("overstory CM trait") + scale_y_continuous("", breaks = c(20, 22, 24, 26, 28,30)) + 
  theme(legend.position = "none") + theme(aspect.ratio=1) +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black")) +
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 8))

s9b <- ggplot(dataplot2, aes(x = leaf_N, y = leaf_N_f, color = elevation)) + geom_point() + 
  geom_smooth(method = lm, colour = "grey25") + 
  ggtitle(label = bquote('leaf N (mg/g)')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradientn(colours = rev(viridis(10))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 13)) + 
  scale_x_continuous("overstory CM trait") + scale_y_continuous("", breaks = c(20, 22, 24, 26, 28,30)) + 
  theme(legend.position = "right") + theme(aspect.ratio=1) +
  theme(axis.text.x = element_text(colour = "black")) + theme(axis.text.y = element_text(colour = "black")) +
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 8))

tiff(file = "CWM trait trait_exclDeciduoussp2.tiff", width = 25, height = 25, units = "cm", res = 300,compression = "lzw")
plot_grid(s1, s2, s3, s4, s5, s6,s7, s8, s9, ncol = 3, axis = "rlbt", labels="AUTO")
dev.off()

#RDA (& varpart) on presence/absence L matrices on species matrix of ONLY species with traits ----

#woody species
#covert L_w to presence/absence (L_w01), to match data from understory
L_w01.t <- L_w.t
L_w01.t[L_w01.t>0] <- 1

rda.w01.0 <- rda(L_w01.t ~ 1, e)                                         #intercept model
rda.w01.1 <- rda(L_w01.t ~ 1 + elevation + heatload2 + fog_freq + 
                   trans_rock_soil + soil_depth + PC1s + PC2s + PC3s, e)      #full model

#the COMBINED model
ordiR2step(rda.w01.0, scope = formula(rda.w01.1), direction = "forward", perm.max = 9999, R2scope = TRUE) #automatic stepwise model selection based on R2
rda.w01 <- rda(L_w01.t ~ 1 + elevation + fog_freq + PC2s + PC3s + PC1s + heatload2, e)   #final model 1 
summary(rda.w01)
plot(rda.w01)

RsquareAdj(rda.w01) #model adj R2 0.3111928
anova.cca(rda.w01)
anova.cca(rda.w01, by = "axis", permutations = 9999)
anova.cca(rda.w01, by = "margin", permutations = 9999)

vpart.w01 <- varpart(L_w01.t, ~ elevation + heatload2 + fog_freq, ~ trans_rock_soil + soil_depth + PC1s + PC2s + PC3s, data = e)
vpart.w01

p2x_label <- paste('total~R^2 == ', round(vpart.w01$part$fract$Adj.R.squared[3]*100, 1))
p2x <- draw.pairwise.venn(
  round((vpart.w01$part$indfract$Adj.R.squared[1] + vpart.w01$part$indfract$Adj.R.squared[2])*100/ vpart.w01$part$fract$Adj.R.squared[3], digits=1),
  round((vpart.w01$part$indfract$Adj.R.squared[3] + vpart.w01$part$indfract$Adj.R.squared[2])*100/ vpart.w01$part$fract$Adj.R.squared[3], digits=1), 
  round(vpart.w01$part$indfract$Adj.R.squared[2]*100/ vpart.w01$part$fract$Adj.R.squared[3], digits = 1),
  c("topography", "soil"), fill = c(viridis(3)[2], viridis(3)[3]), col = c(viridis(3)[2], viridis(3)[3]), alpha = c(0.6, 0.6), ind = FALSE, cat.pos = c(330, 30), 
  scaled = T, cat.dist = -0.08,cat.cex = 1.5)
p2x <- ggdraw(p2x) + 
  annotate(geom = 'text', x = 0.65, y = 0.02, label = p2x_label, parse = TRUE, size = 5) +
  theme(plot.background = element_rect(fill = NA), plot.margin = margin(15, 15, 15, 15)) + theme(aspect.ratio=1/1)
p2x

#ferns

rda.ft.0 <- rda(L_ft.t ~ 1, e)                                         #intercept model
rda.ft.1 <- rda(L_ft.t ~ 1 + elevation + heatload2 + fog_freq + 
                  trans_rock_soil + soil_depth + PC1s + PC2s + PC3s, e)      #full model

#the COMBINED model
ordiR2step(rda.ft.0, scope = formula(rda.ft.1), direction = "forward", perm.max = 9999, R2scope = TRUE) #automatic stepwise model selection based on R2
rda.ft <- rda(L_ft.t ~ 1 + elevation + fog_freq + PC2s + heatload2, e)   #final model 1 
summary(rda.ft)
plot(rda.ft)

RsquareAdj(rda.ft) #model adj R2 0.3113454
anova.cca(rda.ft)
anova.cca(rda.ft, by = "axis", permutations = 9999)
anova.cca(rda.ft, by = "margin", permutations = 9999)

vpart.ft <- varpart(L_ft.t, ~ elevation + heatload2 + fog_freq, ~ trans_rock_soil + soil_depth + PC1s + PC2s + PC3s, data = e)
vpart.ft

p1x_label <- paste('total~R^2 == ', round(vpart.ft$part$fract$Adj.R.squared[3]*100, 1))
p1x <- draw.pairwise.venn(
  round((vpart.ft$part$indfract$Adj.R.squared[1] + vpart.ft$part$indfract$Adj.R.squared[2])*100/ vpart.ft$part$fract$Adj.R.squared[3], digits = 1),
  round((vpart.ft$part$indfract$Adj.R.squared[3] + vpart.ft$part$indfract$Adj.R.squared[2])*100/ vpart.ft$part$fract$Adj.R.squared[3], digits = 1),
  round(vpart.ft$part$indfract$Adj.R.squared[2]*100/ vpart.ft$part$fract$Adj.R.squared[3], digits = 1),
  c("topography", "soil"), fill = c(viridis(3)[2], viridis(3)[3]), col = c(viridis(3)[2], viridis(3)[3]), alpha = c(0.6, 0.6), ind = F, cat.pos = c(330, 30),
  scaled = T, cat.dist = -0.08,cat.cex = 1.5)
p1x <- ggdraw(p1x) + 
  annotate(geom = 'text', x = 0.65, y = 0.02, label = p1x_label, parse = TRUE, size = 5) +
  theme(plot.background = element_rect(fill = NA), plot.margin = margin(15, 15, 15, 15))+ theme(aspect.ratio=1/1)
p1x

#combine venn diagram plots for the alternative RDA analyses (traits excl decidous & species excluding species with no traits)
tiff(file = "venn diagrams redux analyses.tiff", width V = 27, height = 9, units = "cm", res = 300,compression = "lzw")
plot_grid(p2x, p1x, p5_ev, labels = c('A overstory species subset', 'B understory species subset', 'C overstory traits subset'), nrow = 1, ncol = 3, rel_heights=c(1,1))
dev.off()

#FDR on p-values for trait-trait & div-div correlations ----

FDR_in <- read_csv('input FDR pvalues.csv')
FDR_out <- cbind(FDR_in, p.adjust(FDR_in$p, "fdr"))
write.csv(FDR_out, 'FDR_out.csv')
