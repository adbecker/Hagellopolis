rm(list=ls())
graphics.off()

require(Rmisc)
require(ggplot2)
require(network)
require(epinet)
require(animation)

set.seed(123)

setwd('~/Hagellopolis/')

source("fastHag.R")
source("simseir.R")

## define global parameters
beta = .1
ke = 17
ki = 22
thetae = 0.65
thetai = 0.4
latencydist = "fixed"
latencyperiod = 1

## birthrate should be fraction of susceptibles
birthrate <- 0.25

#R0 <- (1 - (1/(1+beta*thetai))^ki)

#Create Hagellopolis
n_pop <- 499
#Total population
n_hags <- 5 # number of villages
hagopolis_size = 1e6# size of landscape (square sides)
hag_size = n_pop/n_hags # size of a single village

## how many years?
n_years <- 10

## burn in time?
## need this to get to equilibrium
burn_in <- 4

## how often do you rewire the network
## this is in days
rewire.time <- 365

## set timer
ptm = proc.time()

## get covariates for each village
## produces a dataframe of n_pop individuals and assigns them:
## house loc, village loc village ID, house, class, age, sex status
ndcv = sim_node_cov(n_pop, n_hags, city_size = hagopolis_size, village_size = hag_size)

## start with ages such that most are 3 years old, going down until 8 years old
## roughly the school ages
## say 40% 3 year olds, 20% 4, 20% 5, 10% 6, 5% 7 5% 8
init.ages <- c(rep(3,n_pop*.4),rep(4,n_pop*.2),rep(5,n_pop*.2),rep(6,n_pop*.1),rep(7,n_pop*0.05),rep(8,n_pop*0.05))

## if this doesnt generate enough people (nrow(ndcv)) just add more 3 year olds
if(length(init.ages) != nrow(ndcv)){
  init.ages <- c(init.ages, rep(3,nrow(ndcv) - length(init.ages)))
}

## replace the hagelloch ages to be more like an endemic scenario
ndcv$age <- init.ages

## assign empty vector to keep track of yearly births and susceptibles
births <- S <- rep(0,n_years)

#setwd('~/Hagellopolis/gifs/')

#saveGIF({
  
## make up eta vector
## assign random weights for the covariate strengths
eta_vec <- -0.3*randomLHS(1,8)
#eta_vec <- -.25 * randomLHS(1, 8)
#eta_vec[4]<- -.004 * n_pop / n_hags
#eta_vec[8] <- -.004
eta_vec <- c(eta_vec)

for(it in 1:n_years){
  
  print(sprintf('starting year %d',it))
  
  ## record info about villages
  ## at each year this stays the same unless everyone in that village gets infected
  village = data.frame(table(ndcv$village), unique(ndcv$village_xpos), unique(ndcv$village_ypos))
  names(village) <- c("id", "pop", "xpos", "ypos")
  village$index = cumsum(c(1,village$pop))[1:n_hags]
  #village$index = cumsum(c(1,village$pop))[1:nrow(village)]
  
  
  ## simulate edges within villages
  #within_edges = sim_within_village_edges(node_cov = ndcv, eta = eta_vec, village_vec = 1:n_hags)
  within_edges = sim_within_village_edges(node_cov = ndcv, eta = eta_vec, village_vec = 1:nrow(village))
  
  ## simulate edges between villages
  #btwn_edges = sim_between_village_edges(village, eta = -0.04, ndcv)
  btwn_edges = sim_between_village_edges(village, eta = mean(eta_vec), ndcv)
  
  ## stick them together
  my_ergm = rbind(within_edges, btwn_edges)
  
  ## this rewires the network each year
  ## should come up with a way to add the connections without really changing the existing ones
  
  M = my_ergm
  N = n_pop
  
  ## plot the network if it isn't too big
  if(n_pop < 500){
    plot.network(as.network(M))
    title(sprintf('Year %d Network',it))
  }
  
  ## assign initial conditions to start the outbreak
  ## at the first time, set I0 to be ceiling(1e-5*pop)
  ## this is ~~ about in line with historical UK
  ## we can draw parms from He et al (2010) Plug and Play inference or
  ## Grenfell and Finkenstadt (2000) dynamics of tsir model
  
  if(it == 1){
    I0 <- ceiling(1e-5*n_pop)
  }else{
    I0 <- max(still.inf,ceiling(1e-5*n_pop))
  }
  
  ## number of people susceptible in the pop is number of people in the ndcv df
  ## since this is yearly its not super important to keep track of 
  S[it] <- nrow(ndcv)
  
  ## now run the epidemic for that year
  ## this is a slightly modified simseir function with following changes:
  ## 1) instead of one infected this allows for I0 initial infected
  ## 2) added a rewire time that if the next recovering is occuring after rewire.time, end the epidemic
  epidemic = simseir(M = my_ergm, N = nrow(ndcv), I0 = I0, rewire.time = rewire.time,
                     beta = beta, ki = ki, thetai = thetai,
                     latencydist = latencydist, latencyperiod = latencyperiod)
  
  ## store the results in a new data frame 
  new.epi <- as.data.frame(epidemic)
  
  ## see who got infected, i.e. who recovered 
  ## pull their indices
  removed <- new.epi$'Node ID'[which(is.na(new.epi$Rtime) == F)]
  
  ## find how are still infectious
  ## remove the all NA columns
  ## if number of people removed is the whole population, then no one is infected to carry over into next year
  if(length(removed) == nrow(ndcv)){
    still.inf <- 0
  }else{
    still.inf <- new.epi[-c(which(is.na(new.epi$Parent) == TRUE)[2]:nrow(new.epi)),]
    still.inf <- sum(is.na(still.inf$Rtime))
  }
  
  ## record ages of who got infected
  R.ages <- ndcv[removed,]
  R.ages <- R.ages$age
  
  ## store all of these into a vector called ages as the loop goes on
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
  new.susceptibles <- round(birthrate*nrow(ndcv)*runif(1,0.5,1.5))
  
  ## make replacement covariates
  ## make all replacements have age 3
  ## consistent with UK entry to school
  replacements <- sim_node_cov(N=new.susceptibles, n_village=n_hags, city_size = hagopolis_size, village_size = hag_size)
  replacements$age <- 3
  
  ## the new replacements are births, so store them
  births[it] <- nrow(replacements)
  
  ## these are the village options based on ID, and xpos and ypos
  village.opts <- unique(ndcv[,4:6])
  
  ## sample the village options and then assign this random sample to the replacement 'births'
  replacements[,4:6] <- village.opts[sample(nrow(village.opts), nrow(replacements),replace = T), ]
  
  ## bind the ndcv wthat has removed infecteds with the replacement births
  new.ndcv <- rbind(new.ndcv,replacements)
  
  ## order it by village
  new.ndcv <- new.ndcv[order(new.ndcv$village),]
  
  ##rewire the network
  ## replace IDs as the should be the null contact covariate  
  new.ndcv$id <- 1:nrow(new.ndcv)
  rownames(new.ndcv) <- new.ndcv$id
  
  ## replace the ndcv with the new updated yearly one
  ndcv <- new.ndcv
  
  ## add 365 to the epi times based on the year to keep track of date of infection
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

#},movie.name = "Hagellopolis_network.gif")

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
#plotepitree(na.rm.epi)

## use the profvis library for examing the speed of code
# library(profvis)
# profvis({
#   within_edges1 = sim_within_village_edges(node_cov = ndcv, eta = eta_vec, village_vec = 1:n_hags)
#   btwn_edges1 = sim_between_village_edges(village, eta = -0.4, ndcv)
# })
# 
# 