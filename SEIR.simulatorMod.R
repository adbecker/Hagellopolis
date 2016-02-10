SEIR.simulatorMod <- function (dyadcov, N, beta, ki, thetai, ke = ki, thetae = thetai, latencydist = "fixed", 
          latencyperiod = 0, eta_vec= NULL, regenFreq = .05) 
{

# N=1000
# beta=2
# ki = 5
# thetai = 5
# ke = 2
# latencydist = "gamma"
# thetae =thetai



  if (!is.matrix(M)) 
    stop("Input error: Network M must be an edgelist matrix.")
  if ((length(dim(M)) != 2) || (dim(M)[2] != 2)) 
    stop("Input error: Network M must an a 2-dimensional edgelist matrix.")
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

