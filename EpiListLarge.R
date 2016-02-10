#Load Packages
#install.packages(epinet)
#install.packages(network)
#install.packages(lhs)

require(epinet)
require(network)
require(lhs)

#Function to generate Dyadic Cov Matrix of connected mini-Hagellochs
createNewHagelloch <- function(N, id_list, Nhags,town_btw=1000, town_in=100){
  #generate Class vectors
  Class1_vec = sample(0:1,N,replace = T)
  Class2_vec =abs(1-Class1_vec)
  
  #generate sex vectors
  Male_vec = sample(0:1,N,replace=T)
  Female_vec = abs(1-Male_vec)
  
  #generate house positions
  #1) get pop size of each mini hagelloch
  popSize_list <-sort(sample(1:N,(Nhags-1)))
  popSize_list <-c(0,popSize_list,N)
  popSizes <-diff(popSize_list)
  
  #2) get town x,y coordinates
  town_xposvec <- lapply(1:Nhags, function(x) rep( town_btw*runif(1),popSizes[x]))
  town_xposvec <- unlist(town_xposvec)
  town_yposvec <- lapply(1:Nhags, function(x) rep( town_btw*runif(1),popSizes[x]))
  town_yposvec <- unlist(town_yposvec)
  house_xposvec <- town_in*runif(N)+town_xposvec
  house_yposvec <- town_in*runif(N)+town_yposvec
  
  #set up data frame with spatial positions of N individuals (random x and y) and class, gender, age
  mycov <- data.frame(id =id_list, house_xpos = house_xposvec, house_ypos = house_yposvec,HouseHold = rep(0,N),
                      Class1 = Class1_vec, Class12 =rep(0,N),
                      Male = Male_vec, FemaleM = rep(0,N),
                      age = sample(1:13,N, replace = T), 
                      town_xpos=town_xposvec, town_ypos = town_yposvec)
  #generate dyadic cov matrix
  dyadCov <- BuildDyadicCovMatrix(mycov, unaryCol = c(4,9), unaryFunc = c("absdiff","absdiff"), 
                                  binaryCol = list(c(2,3), c(5,6),c(7,8),c(10,11)), binaryFunc = c("euclidean","euclidean","euclidean","euclidean"))
  
  #create columns with correct matching indicators (reverse from distance)
  dyadCov <- cbind(dyadCov, Class2 = abs(rep(1,dim(dyadCov)[1])-dyadCov[,'Class1.Class12.L2Dist']), Female = abs(rep(1,dim(dyadCov)[1])-dyadCov[,'Male.FemaleM.L2Dist']))
  dyadCov[,'Class1.Class12.L2Dist'] <- dyadCov[,'Class2']
  dyadCov[,'Male.FemaleM.L2Dist'] <- dyadCov[,'Female']
  
  #replace values with 0s for match
  dyadCov[Class2_vec[dyadCov[,'node.1']] == 0,'Class2'] <-0                
  dyadCov[Class1_vec[dyadCov[,'node.1']] == 0,'Class1.Class12.L2Dist'] <-0                
  dyadCov[Female_vec[dyadCov[,'node.1']] == 0,'Female'] <-0                
  dyadCov[Male_vec[dyadCov[,'node.1']] == 0,'Male.FemaleM.L2Dist'] <-0                
  
  #Reorder columns to match format of original Hagelloch data matrix
  dyadCov <- dyadCov[,c(1:4,7,10, 6, 8, 11,5,9)]
  colnames(dyadCov) <- c(colnames(HagellochDyadCov), "Town Distance")
  colnames(dyadCov)[11] <- 'Town Distance'
  dyadCov[dyadCov[,'House Distance']<3,'House Distance'] <-0
  dyadCov[dyadCov[,'House Distance'] <3, 'Household'] <-1
  dyadCov[dyadCov[,'Town Distance']>0,c('Classroom 2', 'Classroom 1', 'Male Match',
                                        'Female Match', 'Household','Age Diff')]<-0
  
  
  newHagelloch <-dyadCov
  return(newHagelloch)
}

#generate SEIR function
SEIR.simulatorMod <- function (dyadcov, N, beta, ki, thetai, ke = ki, thetae = thetai, latencydist = "fixed", 
                               latencyperiod = 0, eta_vec= NULL, regenFreq = .05) 
{
  require(epinet)
  # N=1000
  # beta=2
  # ki = 5
  # thetai = 5
  # ke = 2
  # latencydist = "gamma"
  # thetae =thetai
  
  
  
  #if (!is.matrix(M)) 
  #stop("Input error: Network M must be an edgelist matrix.")
  #if ((length(dim(M)) != 2) || (dim(M)[2] != 2)) 
  #  stop("Input error: Network M must an a 2-dimensional edgelist matrix.")
  #get initial infected case
  init <- sample(1:N, 1)
  #initialize arrays of transmission time and recovery time
  t.time = array(dim = N)
  r.time = array(dim = N)
  #set transmission time for initial case
  t.time[init] <- ifelse(latencydist == "fixed", latencyperiod, 
                         rgamma(1, ke, scale = thetae))
  #set recovery time for initial case
  r.time[init] = rgamma(1, ki, scale = thetai) + t.time[init]
  #set id of next individual to recover (?)
  nextrec = init
  #initialize list of infecteds
  inf.list <- matrix(c(init, NA, 0, t.time[init], NA), nrow = 1)
  #set time to start time of initial infection
  time <- cm.time <- t.time[init]
  #set next transmittor to init ID
  nexttrans = init
  #reset transmission time for initial case to Inf
  t.time[init] <- Inf
  #initialize susceptibles, exposed and infected
  s.list <- (1:N)[-init]
  e.list <- NULL
  i.list <- init
  inf <- list(i.list)
  susc <- list(s.list)
  expo <- list(e.list)
  #generate network
  M = BuildDyadicLinearERGM(N=N, dyadiccovmat = dyadcov, eta = eta_vec)
  for (i in 2:(N * 3)) {
    #regenerate network with some probability
    if (runif(1) < regenFreq)
      M = BuildDyadicLinearERGM(N=N, dyadiccovmat = dyadcov, eta = eta_vec)
    #list of IDs of susceptibles
    s.list <- array(susc[[i - 1]])
    i.list <- array(inf[[i - 1]])
    #put together boolean vector of network connections between an infected and susceptible
    si.ex <- ((M[, 1] %in% i.list) & (M[, 2] %in% s.list)) | 
      ((M[, 2] %in% i.list) & (M[, 1] %in% s.list))
    #add up the number of connections
    n.si <- sum(si.ex)
    
    #get duration of infection for next one to recover
    dwt <- ifelse(length(inf[[i - 1]]) > 0, r.time[nextrec] - 
                    cm.time, Inf)
    #if there are infection events, get time in susceptible class
    bwt <- ifelse(n.si != 0, rexp(1, beta * n.si), Inf)
    #get  transition time
    twt <- t.time[nexttrans] - cm.time
    #choose which event to do next
    ewt <- min(bwt, dwt, twt, na.rm = TRUE)
    #increment time vector
    time <- c(time, ewt)
    #move up new start time
    cm.time <- cm.time + ewt
    #check which event to do
    if (ewt == bwt) 
      test <- "Infect"
    else if (ewt == dwt) 
      test <- "removal"
    else test <- "transition"
    #if event is infection
    if (test == "Infect") {
      #get pairs of connections btw infected and susceptible
      is.pairs <- which(si.ex == 1)
      #choose one of the pairs
      smp.ind <- ifelse(length(is.pairs) == 1, is.pairs, 
                        sample(is.pairs, 1))
      #get the index of the parent (i.e. infected) -- which column of the network
      #matrix it is in
      parentindex <- which(M[smp.ind, ] %in% i.list)
      # the newly infected susceptible is whichever ID is in the other column
      new.inf <- M[smp.ind, 3 - parentindex]
      #parent is the one in the original column
      parent <- M[smp.ind, parentindex]
      #generate latency time
      lat <- ifelse(latencydist == "fixed", latencyperiod, 
                    rgamma(1, ke, scale = thetae))
      #time of new infection infectiousness is time of infection plus latency
      t.time[new.inf] <- cm.time + lat
      #append updated list of susceptibles with new infection removed to original
      susc <- append(susc, list(susc[[i - 1]][-which(susc[[i - 
                                                             1]] == new.inf)]))
      #add new infection to list of exposed
      expo <- append(expo, list(c(expo[[i - 1]], new.inf)))
      #repeat list of infecteds and append to list from previous timestep
      inf <- append(inf, list(inf[[i - 1]]))
      #add row to inf list with infection timing of new infected
      inf.list <- rbind(inf.list, c(new.inf, parent, cm.time, 
                                    NA, NA))
      nexttrans <- which(t.time == min(t.time, na.rm = TRUE))
    }
    else if (test == "removal") {
      if (i == 2) {
        inf.list[1, 5] <- cm.time
        break
      }
      new.rec <- nextrec
      susc <- append(susc, list(susc[[i - 1]]))
      expo <- append(expo, list(expo[[i - 1]]))
      inf <- append(inf, list(inf[[i - 1]][-which(inf[[i - 
                                                         1]] == new.rec)]))
      inf.list[which(inf.list[, 1] == new.rec), 5] <- cm.time
      r.time[nextrec] <- NA
      if (length(inf[[i]]) > 0) 
        nextrec <- which(r.time == min(r.time, na.rm = TRUE))
      else if (length(expo[[i]]) > 0) {
        nextrec <- which(t.time == min(t.time, na.rm = TRUE))
        r.time[nextrec] <- Inf
      }
    }
    else {
      new.trans <- nexttrans
      susc <- append(susc, list(susc[[i - 1]]))
      expo <- append(expo, list(expo[[i - 1]][-which(expo[[i - 
                                                             1]] == new.trans)]))
      inf <- append(inf, list(c(inf[[i - 1]], new.trans)))
      inf.list[which(inf.list[, 1] == new.trans), 4] <- cm.time
      t.time[nexttrans] <- NA
      nexttrans <- which(t.time == min(t.time, na.rm = TRUE))
      r.time[new.trans] <- cm.time + rgamma(1, ki, scale = thetai)
      if (r.time[new.trans] < r.time[nextrec]) 
        nextrec <- new.trans
    }
    if (length(inf[[i]]) + length(expo[[i]]) == 0) {
      break
    }
  }
  inf.list[, 3:5] <- inf.list[, 3:5] - min(inf.list[, 5])
  if (length(s.list) > 0) {
    for (i in 1:length(s.list)) inf.list <- rbind(inf.list, 
                                                  c(s.list[i], NA, NA, NA, NA))
  }
  colnames(inf.list) <- c("Node ID", "Parent", "Etime", "Itime", 
                          "Rtime")
  return(inf.list)
}


#Create Hagellopolis
N<-10000  #Total population
id_list <- 1:N  #ids
Nhags <- 20
newHag <- createNewHagelloch(N, id_list, Nhags,town_btw=1000, town_in=100 )
newHag<-newHag[,-11]
#Set up etas for network
eta_vec <- -.25*randomLHS(1,7)
eta_vec[4]<- -.0004*N/Nhags
#eta_vec[8] <- -.00004
eta_vec <- c(0,eta_vec)

#create vector of probabilities of regenerating network
regenVec <-c(rep(.001,3),rep(.005,3),rep(.01,3),rep(.02,3),rep(.05,3))#c(rep(.001,5),rep(.005,5),rep(.01,5),rep(.02,5),rep(.05,5))
eta_distance_vec <-c(-.00004*N/Nhags,-.0004*N/Nhags,-.004*N/Nhags )
#create empty list
epilist <- list()

for (i in regenVec) { 
  for (j in eta_distance_vec) {
    eta_vec[4] <- j
    epi <-SEIR.simulatorMod(newHag, N=N, beta=2,ki = 5, thetai = 5, ke = 2, 
                            latencydist = "gamma", eta_vec = eta_vec,  
                            regenFreq = i)
    epilist <-append(epilist, list(epi))
    save.image("~/EpiSpaceLarge.RData")
  } 
}

plotepitreeMod <-function (epi, regenFreq, lwd = 1, leaf.labs = TRUE, leaf.cex = 0.75, zero.at.start = FALSE, 
                           main = "Spread of Epidemic", xlab = "Time", ylab = "", e.col = "black", 
                           i.col = "red", lty.transmission = 3, marktransitions = TRUE, 
                           label.trans = "|", cex.trans = 0.5, ...) 
{
  ninf <- nrow(epi)
  if (ninf < 2) 
    stop("Plot error: need at least two infecteds to plot the epidemic.")
  susc <- NULL
  for (i in 1:ninf) if (sum(is.na(epi[i, 3:5])) == 3) 
    susc <- c(susc, i)
  if (length(susc) > 0) 
    epi <- epi[-susc, ]
  ninf <- nrow(epi)
  if (ninf < 2) 
    stop("Plot error: need at least two infecteds to plot the epidemic.")
  if (sum(is.na(epi[, 3:5])) > 0) 
    stop("Plot error: need full data to plot the epidemic -- run again with all times known.")
  if (sum(is.na(epi[, 1])) > 0) 
    stop("Plot error: missing entries in first (Node ID) column.")
  if (sum(is.na(epi[, 2])) > 1) 
    stop("Plot error: too many NA values in second (Parent) column -- should only have one initial infected node.")
  if (zero.at.start) 
    epi[, 3:5] = epi[, 3:5] - min(epi[, 3:5])
  epi <- epi[order(epi[, 3]), ]
  nchild = array(0, max(epi[, 1]))
  for (i in ninf:2) nchild[epi[i, 2]] = nchild[epi[i, 2]] + 
    nchild[epi[i, 1]] + 1
  ypos = array(0, max(epi[, 1]))
  ypos[epi[1, 1]] = 1
  for (i in 2:ninf) {
    ypos[epi[i, 1]] = ypos[epi[i, 2]] + nchild[epi[i, 2]] - 
      nchild[epi[i, 1]]
    nchild[epi[i, 2]] = nchild[epi[i, 2]] - nchild[epi[i, 
                                                       1]] - 1
  }
  plot(x = c(min(epi[, 3]), max(epi[, 5]) + 2), y = c(0.9, 
                                                      ninf), pch = " ", xlab = xlab, ylab = ylab, main = c(main,"Regeneration Frequency =", 100*regenFreq,"%"), 
       yaxt = "n", bty = "n", ...)
  for (i in 1:ninf) {
    lines(x = c(epi[i, 3], epi[i, 4]), y = c(ypos[epi[i, 
                                                      1]], ypos[epi[i, 1]]), lwd = lwd, col = e.col)
    lines(x = c(epi[i, 4], epi[i, 5]), y = c(ypos[epi[i, 
                                                      1]], ypos[epi[i, 1]]), lwd = lwd, col = i.col)
    lines(x = c(epi[i, 3], epi[i, 3]), y = c(ypos[epi[i, 
                                                      2]], ypos[epi[i, 1]]), lwd = lwd, lty = lty.transmission)
  }
  if (leaf.labs) 
    text(epi[, 5], ypos[epi[, 1]], labels = epi[, 1], pos = 4, 
         offset = 0.25, cex = leaf.cex)
  if (marktransitions) 
    for (i in 1:ninf) text(epi[i, 4], ypos[epi[i, 1]], labels = label.trans, 
                           cex = cex.trans)
}

# par(mfrow=c(3,3))
# x<- c(1:9)
# plot1 <-function(x) plotepitreeMod(epilist[[x]],regenVec[x])
# lapply(x, plot1)
# 
# par(mfrow=c(3,3))
# x<- c(10:18)
# plot1 <-function(x) plotepitreeMod(epilist[[x]],regenVec[x])
# lapply(x, plot1)
# 
# par(mfrow=c(3,2))
# x<- c(19:25)
# plot1 <-function(x) plotepitreeMod(epilist[[x]],regenVec[x])
# lapply(x, plot1)


save.image("~/EpiSpaceLarge.RData")