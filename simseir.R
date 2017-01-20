simseir <- function(M, N, beta, I0, ki, thetai, ke = ki, thetae = thetai, latencydist = "fixed", 
                    latencyperiod = 0, rewire.time){
  
  # First do some input format checking
  if (!is.null(M))
  {
    # M should be an edgelist matrix
    if(!is.matrix(M))
      stop("Input error: Network M must be an edgelist matrix.")
    
    # Make sure M has the correct dimensions
    if ( (length(dim(M)) != 2) || (dim(M)[2] != 2) )
      stop("Input error: Network M must an a 2-dimensional edgelist matrix.")
  }
  
  # convert to network data structure
  net = network.initialize(N, directed = FALSE)
  add.edges(net, head = M[,1], tail = M[,2])
  
  # make another empty network of same size to keep track of potential infections  
  trans_net = network.initialize(N, directed = FALSE)
  # double
  t_net = array(FALSE, network.edgecount(net))
  
  # lists of current states of state
  is_infectious = is_exposed = is_removed = array(FALSE, N)
  is_susceptible = array(TRUE, N)
  
  # get location of first infections
  ## pull I0 samples from N and those are the individuals infected to start
  init <- sample(N,I0)  # Inital infected individual is chosen at random
  if(N == 1){
    init <- N
  }
  # update state vectors
  is_infectious[init] = TRUE
  is_susceptible[init] = FALSE
  
  # update live transmission graph
  # get all neighbours
  neighbour = get.neighborhood(net, init)
  # double
  neighbour_id = get.edgeIDs(net,init)
  # find susceptoble ones
  susc_neighbour = neighbour[is_susceptible[neighbour]]
  # double
  susc_neighbour_id = neighbour_id[is_susceptible[neighbour]]
  
  # add them
  if (length(susc_neighbour) > 0) add.edges(trans_net, tail = rep(init, length(susc_neighbour)), head = susc_neighbour) #, edge.check = TRUE)
  # double
  t_net[susc_neighbour_id] = TRUE
  
  # Keep a list of all upcoming transition and recovery times (t.time[i] = r.time = NA if i is susceptible)
  t.time = array(dim = N)
  r.time = array(dim = N)
  
  # Generate a transition time and recovery time for initial infected
  t.time[init] <- ifelse( latencydist=="fixed", latencyperiod, rgamma(1, ke, scale = thetae) )  
  r.time[init] = rgamma(1, ki, scale = thetai) + t.time[init]
  
  nextrec = init  	# Keep track of who, of current infecteds, is next to recover
  
  inf.list <- matrix(c(init, NA, 0, t.time[init], NA), nrow = 1)   # Keep track of initial infection
  
  time <- cm.time <- t.time[init]
  
  nexttrans = init # Temporary value
  t.time[init] <- Inf # Temporary value
  
  # keep track of number exposed/infectious
  count_i = 1
  count_e = 0
  
  for( i in 2:(N*3) )	# Maximum number of iterations is 3*N (1 infection ,1 transition from  exposed to infective, and 1 recovery for every node)
  {
    # count number of SI pairs
    n.si = network.edgecount(trans_net)
    n.si_id = sum(t_net)
    
    #     print(paste("n.si", n.si))
    #     print(inf.list)
    #     plot.network(network(trans_net, directed = FALSE), displaylabels = TRUE, edge.label = TRUE, edge.label.col = "blue")
    
    
    # Draw waiting times for the next removal, infection, and transition
    dwt <- ifelse( count_i > 0, r.time[nextrec] - cm.time, Inf ) 
    bwt <- ifelse( n.si != 0, rexp(1,beta*n.si), Inf)
    twt <- t.time[nexttrans] - cm.time
    
    ewt <- min(bwt, dwt, twt, na.rm = TRUE)
    time <- c(time, ewt)			# Increment time
    cm.time <- cm.time + ewt		# Increment cumulative time
    
    if (ewt == bwt) { 
      test <- "Infect"
    } else if (ewt == dwt) {
      test <- "removal" 
    } else {
      test <- "transition"
    }
    
    # print(test)
    
    if (test == "Infect"){ 
      # Event is an infection	
      
      # choose the edge where transmission occured
      trans_edge = ifelse(n.si > 1, sample(valid.eids(trans_net), 1), valid.eids(trans_net))
      trans_edge_id = ifelse(n.si_id > 1, sample(which(t_net), 1), which(t_net))
      
      # print(trans_edge)
      
      # see which end is new infection and which is the parent
      if ( is_susceptible[trans_net$mel[[trans_edge]]$inl] ){
        new.inf = trans_net$mel[[trans_edge]]$inl
        parent = trans_net$mel[[trans_edge]]$outl
      } else {
        new.inf = trans_net$mel[[trans_edge]]$outl
        parent = trans_net$mel[[trans_edge]]$inl
      }
      
      # Generate a transition time for new infected
      lat <- ifelse( latencydist == "fixed", latencyperiod, rgamma(1, ke, scale = thetae) )
      
      t.time[new.inf] <- cm.time + lat
      
      # Update state vectors
      is_exposed[new.inf] = TRUE
      is_susceptible[new.inf] = FALSE
      # update number of exposed 
      count_e = count_e + 1
      
      # update live transmission graph --- remove edges adjacent to newly exposed 
      delete.edges(trans_net, eid = get.edgeIDs(trans_net, new.inf))
      # double
      t_net[get.edgeIDs(net, new.inf)] = FALSE
      
      inf.list <- rbind(inf.list,c(new.inf, parent, cm.time, NA, NA))
      
      nexttrans <- which(t.time == min(t.time, na.rm = TRUE))
      
    } else if (test == "removal")	# Event is a removal
    {								
      if(i==2)	# Only infected dies out
      {
        inf.list[1,5] <- cm.time
        break
      } 	
      
      new.rec <- nextrec
      
      # Update lists of infectious, removed
      is_infectious[new.rec] = FALSE
      is_removed[new.rec] = TRUE
      # update number of infectious
      count_i = count_i - 1
      
      # update live transmission graph --- remove edges adjacent to newly recovered 
      delete.edges(trans_net, eid = get.edgeIDs(trans_net, new.rec))
      # double
      t_net[get.edgeIDs(net, new.rec)] = FALSE
      
      # record removal time
      inf.list[which(inf.list[,1]==new.rec), 5] <- cm.time
      
      
      # Update recovery times and next infected to recover
      r.time[nextrec] <- NA
      if(count_i > 0) 
        nextrec <- which( r.time == min(r.time, na.rm = TRUE) ) 
      else if (count_e > 0)
      {
        nextrec <- which( t.time == min(t.time, na.rm = TRUE) )
        r.time[nextrec] <- Inf
      }
      
    } else 		# Event is a transition
    { 	 
      new.trans <- nexttrans
      
      # Update lists of susceptible, exposed and infecteds	 
      is_exposed[new.trans] = FALSE
      is_infectious[new.trans] = TRUE
      
      # update count of infecteds/exposeds
      count_i = count_i + 1
      count_e = count_e - 1
      
      # update live transmission graph 
      # get all neighbours of newly infectious
      neighbour = get.neighborhood(net, new.trans)
      # double
      neighbour_id = get.edgeIDs(net,new.trans)
      # find susceptoble oness
      susc_neighbour = neighbour[is_susceptible[neighbour]]
      # double
      susc_neighbour_id = neighbour_id[is_susceptible[neighbour]]
      # add them
      if (length(susc_neighbour) > 0) add.edges(trans_net, tail = rep(new.trans, length(susc_neighbour)), head = susc_neighbour)
      # double
      t_net[susc_neighbour_id] = TRUE
      
      # record time of transition
      inf.list[which(inf.list[,1]==new.trans),4] <- cm.time
      
      # Update transition times, recovery times, and assign a recovery time to the new transition
      t.time[nexttrans] <- NA      
      nexttrans <-  which(t.time == min(t.time,na.rm = TRUE))
      r.time[new.trans] <- cm.time + rgamma(1,ki, scale = thetai)
      if (r.time[new.trans] < r.time[nextrec]) nextrec <- new.trans
      
      
    }
    #print(r.time[nextrec])
    
    #if( is.na(r.time[nextrec]) == F && r.time[nextrec] > 365){print('updating to next year'); break}
    if( is.na(r.time[nextrec]) == F && r.time[nextrec] > rewire.time && r.time[nextrec] < Inf){print('breaking'); break}
    
    
    
    if (count_e + count_i == 0) {break}	# No more infectious or exposed members, so epidemic ends
    
    
    
    
  }		
  # reset time 0 to be time of first R
  inf.list[,3:5] <- inf.list[,3:5] - min(inf.list[,5],na.rm = T)
  
  # add any never infecteds to list
  if (any(is_susceptible)) inf.list <- rbind(inf.list, cbind(which(is_susceptible), NA, NA, NA, NA) )
  
  # fix up column names
  colnames(inf.list) <- c("Node ID", "Parent", "Etime", "Itime", "Rtime")
  
  return(inf.list)
}
