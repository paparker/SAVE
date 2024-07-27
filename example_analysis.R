library(tidycensus)
library(tigris)
library(spdep)
library(ggplot2)
library(dplyr)
library(tidyr)
library(maptools)
library(ggthemes)
library(Matrix)
library(mvtnorm)
library(readr)
library(sampling)
library(mase)
library(truncnorm)
library(statmod)
library(matrixStats)
library(VGAM)
library(rstan)
library(expm)
library(parallel)
library(MASS)
library(spatialreg)
source('R/Models/MCMC.R')


### Prep data
temp <- read_rds("example_data.rds")

sf <- pumas('06', cb=T, year=2019)
W <-as.matrix(as_dgRMatrix_listw(nb2listw(poly2nb(sf), style='B')))


XmatA <- cbind(1, log(temp$TotalPop))


### Fit model
iter <- 3000
burn <- 1000
set.seed(1)
modSpat <- saeICAR2(Y=log(temp$Y), X=XmatA, Psi=diag(length(temp$Y)), S2=temp$Var/temp$Y^2, SS=temp$N, W=W, iter=iter, burn=burn, propSD = 4)


### Plot results
lowSpat <- apply(exp(modSpat$Preds[,-c(1:burn)]), 1, quantile, probs=0.025)
highSpat <- apply(exp(modSpat$Preds[,-c(1:burn)]), 1, quantile, probs=0.975)
predSpat <- rowMeans(exp(modSpat$Preds[,-c(1:burn)]))
sdSpat <- apply(exp(modSpat$Preds[,-c(1:burn)]), 1, sd)

df2 <- data.frame(PUMA=temp$PUMA, Est=predSpat, SD=sdSpat)

sf2 <- pumas('06', cb=T, year=2019)
sf2 <- sf2 %>% left_join(df2, by=c("PUMACE10"="PUMA")) 
p1 <- ggplot(sf2)+
  geom_sf(size=0, aes(fill=(Est)))+
  theme_map()+
  scale_fill_viridis_c(name="Mean Income")

p2 <- ggplot(sf2)+
  geom_sf(size=0, aes(fill=log(SD)))+
  theme_map()+
  scale_fill_viridis_c(name="Log Standard Error")
cowplot::plot_grid(p1, p2)





