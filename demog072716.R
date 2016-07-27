rm(list=ls())
graphics.off()

require(ggplot2)
require(epinet)
require(animation)

set.seed(123)

setwd('~/Hagellopolis/')

source("fastHag.R")
source("simseir.R")

## define global parameters
beta = 2
ke = 17
ki = 22
thetae = 0.65
thetai = 0.4
latencydist = "fixed"
latencyperiod = 1

## birthrate should be fraction of susceptibles
birthrate <- 0.4

#R0 <- (1 - (1/(1+beta*thetai))^ki)

#Create Hagellopolis
n_pop <- 300  #Total population
n_hags <- 5 # number of villages
hagopolis_size = 1e6# size of landscape (square sides)
hag_size = n_pop/n_hags # size of a single village

## how many years?
n_years <- 25

## burn in time?
## need this to get to equilibrium
burn_in <- 4

## how often do you rewire the network

rewire.time <- 365

ptm = proc.time()

# get covariates for each village
ndcv = sim_node_cov(n_pop, n_hags, city_size = hagopolis_size, village_size = hag_size)

births <- S <- rep(0,n_years)

setwd('~/Hagellopolis/gifs/')

saveGIF({
  
for(it in 1:n_years){
  
  print(sprintf('starting year %d',it))
  
  ## rewire things each year
  
  # make up eta vector
  #eta_vec <- -.25 * randomLHS(1, 8)
  eta_vec <- -0.1*randomLHS(1,8)
  #eta_vec[4]<- -.004 * n_pop / n_hags
  #eta_vec[8] <- -.004
  eta_vec <- c(eta_vec)
  
  # record info about villages
  village = data.frame(table(ndcv$village), unique(ndcv$village_xpos), unique(ndcv$village_ypos))
  names(village) <- c("id", "pop", "xpos", "ypos")
  village$index = cumsum(c(1,village$pop))[1:n_hags]
  #village$index = cumsum(c(1,village$pop))[1:nrow(village)]
  
  # simulate edges within villages
  #within_edges = sim_within_village_edges(node_cov = ndcv, eta = eta_vec, village_vec = 1:n_hags)
  within_edges = sim_within_village_edges(node_cov = ndcv, eta = eta_vec, village_vec = 1:nrow(village))
  
  # simulate edges between villages
  #btwn_edges = sim_between_village_edges(village, eta = -0.04, ndcv)
  
  btwn_edges = sim_between_village_edges(village, eta = -0.02, ndcv)
  
  # stick them together
  my_ergm = rbind(within_edges, btwn_edges)
  
  M = my_ergm
  N = n_pop
  
  if(n_pop < 500){
    plot.network(as.network(M))
    title(sprintf('Year %d Network',it))
  }
  
  if(it == 1){
    I0 <- 1
  }else{
    I0 <- max(still.inf,1)
  }
  
  S[it] <- nrow(ndcv)
  
  epidemic = simseir(M = my_ergm, N = nrow(ndcv), I0 = I0, rewire.time = rewire.time,
                     beta = beta, ki = ki, thetai = thetai,
                     latencydist = latencydist, latencyperiod = latencyperiod)
  
  new.epi <- as.data.frame(epidemic)
  
  ## who got measles   
  removed <- new.epi$'Node ID'[which(is.na(new.epi$Rtime) == F)]
  
  ## find how are still infectious
  ## remove the all NA columns
  if(length(removed) == nrow(ndcv)){
    still.inf <- 0
  }else{
    still.inf <- new.epi[-c(which(is.na(new.epi$Parent) == TRUE)[2]:nrow(new.epi)),]
    still.inf <- sum(is.na(still.inf$Rtime))
  }
  ## record ages of who got infected
  R.ages <- ndcv[removed,]
  R.ages <- R.ages$age
  
  if(it == 1){
    ages <- R.ages
  }
  if(it > 1){
    ages <- c(ages,R.ages)
  }
  
  #if(n_pop < 500){
  #  hist(ages,breaks = 20,main=sprintf('Year %d Network',it))#,xlim = c(0,35))
  #  title(sprintf('Year %d Network',it))
  #}
  
  ## remove them from the pop
  new.ndcv <- ndcv[-removed,]
  
  ## age susceptibles by a year
  new.ndcv$age <- new.ndcv$age + 1
  
  ## add births / new susceptibles 
  new.susceptibles <- length(removed)
  new.susceptibles <- round(birthrate*n_pop*runif(1,0.5,1.5))
  
  replacements <- sim_node_cov(N=new.susceptibles, n_village=n_hags, city_size = hagopolis_size, village_size = hag_size)
  replacements$age <- round(runif(new.susceptibles,2,5))
  
  births[it] <- nrow(replacements)
  
  village.opts <- unique(ndcv[,4:6])
  
  replacements[,4:6] <- village.opts[sample(nrow(village.opts), nrow(replacements),replace = T), ]
  new.ndcv <- rbind(new.ndcv,replacements)
  new.ndcv <- new.ndcv[order(new.ndcv$village),]
  
  ##rewire the network
  #new.ndcv$id <- 1:n_pop
  new.ndcv$id <- 1:nrow(new.ndcv)
  rownames(new.ndcv) <- new.ndcv$id
  
  ndcv <- new.ndcv
  
  new.epi$Etime <- new.epi$Etime + 365*(it-1)
  new.epi$Itime <- new.epi$Itime + 365*(it-1)
  new.epi$Rtime <- new.epi$Rtime + 365*(it-1)
  
  if(it == 1){
    epi <- new.epi
  }
  
  if(it > 1){
    epi <- rbind(new.epi,epi)
  }
  
}

},movie.name = "Hagellopolis_network.gif")

Rtime <- epi$Rtime
Rtime <- na.omit(Rtime)
Rtime <- sort(Rtime)

hist(Rtime,breaks = 52*n_years/5)
plot(Rtime)

## throw away burn ins
throwaway <- which(Rtime >= burn_in*365)
burnin.Rtime <- Rtime[throwaway]
burnin.ages <- ages[throwaway]

hist(ages,breaks=length(unique(ages))-1)

roundUp <- function(x,to=7)
{
  to*(x%/%to + as.logical(x%%to))
}


weekvec <- seq(from=0,to=roundUp(max(Rtime)),by=7)
casesvec <- rep(0,length(weekvec))

for(it in 1:length(weekvec)){
  week <- weekvec[it]
  casesvec[it] <- length(which(Rtime >= week & Rtime < week+7))
}


proc.time() - ptm


data <- as.data.frame(cbind('time'=weekvec,'cases'=casesvec))

ggplot(data,aes(time,cases))+geom_line()

na.rm.epi <- na.omit(epi)

# #Plot spread of epidemic
plotepitree(na.rm.epi)

## use the profvis library for examing the speed of code
# library(profvis)
# profvis({
#   within_edges1 = sim_within_village_edges(node_cov = ndcv, eta = eta_vec, village_vec = 1:n_hags)
#   btwn_edges1 = sim_between_village_edges(village, eta = -0.4, ndcv)
# })
# 
# 