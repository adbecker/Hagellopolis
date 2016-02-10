#Load Packages
require(epinet)
require(network)
require(lhs)

rand_no_overlap <- function(n, x_min = 0, x_max = 1000, y_min = x_min, y_max = x_max, x_size = 10, y_size = x_size){
  # input checking
  if (n <= 0) stop("n must be an integer greater than 0")
  if (x_min > x_max) stop("x_min needs to be less than x_max")
  if (y_min > y_max) stop("y_min needs to be less than y_max")
  if (x_max - x_min < x_size) stop("x_size needs to be smaller than x_max - x_min")
  if (y_max - y_min < y_size) stop("y_size needs to be smaller than y_max - y_min")
  
  # initialise return values
  x = y = rep(0, n)
  
  # draw first value
  x[1] = runif(1, x_min, x_max - x_size)
  y[1] = runif(1, y_min, y_max - y_size)
  
  # if only want one draw, stop here
  if (n == 1) return(list(x = x, y = y))
  
  for (i in 2:n){
    drawn = FALSE 
    try = 0
    # keep making draws until we have a non-overlapping patch
    while(!drawn){
      # make a draw
      x[i] = runif(1, x_min, x_max - x_size)
      y[i] = runif(1, y_min, y_max - y_size)
      # check it doesn't overlap
      if ( any((x[i] + x_size > x[1:i-1]) & (x[i] < (x[1:i-1] + x_size)) & (y[i] + y_size > y[1:i-1]) & (y[i] < (y[1:i-1] + y_size))) ){
        # have overlap
        try = try + 1
      } else {
        # no overlap
        drawn = TRUE
      }
      if (try > 50) stop(paste("unable to find", n,"non-overlapping patches after 50 attempts. Got up to", i,"patches. Try making patches smaller or landscape larger"))
    }
  }
  return(list(x = x, y = y))
}

sim_node_cov <- function(N, n_village, city_size, village_size){

  # decide on size of villages
  pop_size_boundary <-c(0,sort(sample(1:(N-1), (n_village - 1), replace = FALSE)),N)
  pop_size <-diff(pop_size_boundary)
  
  # do village
  village = rep(1:n_village, times = pop_size)
  
  # get village positions
  village_pos = rand_no_overlap(n = n_village, x_max = city_size, x_size = village_size)
  village_x = rep(village_pos$x, times = pop_size)
  village_y = rep(village_pos$y, times = pop_size)
  
  # do classrooms
  classroom = sample(0:2, N, prob = c(.2, .35, .45), replace = T)
  class1 = as.numeric(classroom == 1)
  class2 = as.numeric(classroom == 2)
  
  # do sex
  is_male = sample(0:1, N, replace=T)
  is_female = 1 - is_male
  #do something with households - get household distribution
  #position houses instead of people
  
  # make covariate data frame
  node_cov <- data.frame(id = 1:N, 
                         house_xpos = runif(N, max = village_size) + village_x, 
                         house_ypos = runif(N, max = village_size) + village_y,
                         village_xpos = village_x,
                         village_ypos = village_y,
                         village = village,
                         household = rep(0, N), 
                         class1 = class1, 
                         class2 = class2, 
                         male = is_male, 
                         female = is_female, 
                         age = sample(1:13, N, replace = T))
  
  return(node_cov)
}

comb2 <- function(n){
  ## returns an an array with n choose 2 rows and 2 columns giving all possible pairs of 1:N
  k = n * (n-1) / 2
  x = array(rep(1:(n-1), (n-1):1), dim = c(k,2))
  for (i in 1:(n-1)) {
    z = k - ((n-i+1)*(n-i)/2)
    x[(z + 1):(z + n - i),2] =  (i+1):n
  }
  return(x)
}

build_village_dcm <- function(node_cov, village){
  
  # pull out indidvuals from stated village 
  in_village = which(node_cov[,"village"] == village)
  n_node = length(in_village)
  dyad_cov = matrix(nrow = n_node*(n_node - 1) / 2, ncol = 10)
  
  if (n_node == 1) return(dyad_cov)
  
  # check ids match row numbers for selected group and that they are contiguous
  if (any(node_cov[in_village,1] != in_village[1]:(in_village[1]+length(in_village)-1) )) stop("ids of nodes need to match row number and all nodes in village must be contiguous in larger matrix")
  
  # pull out relevant rows of matrix, and adjust ids
  first_index = in_village[1]
  node_cov = node_cov[in_village,]
  node_cov[,1] = node_cov[,1] - (first_index - 1)
  
  # generate dyadic cov matrix (by hand)
  colnames(dyad_cov) <-  c("node1", "node2", "intercept", "household", "class1", "class2" , "housedist", "male", "female", "agediff")
  # node ids
  dyad_cov[,1:2] = comb2(n_node)
  # intercept
  dyad_cov[,"intercept"] = 1
  # household
  dyad_cov[,"household"] = 0
  # do male, female, class1 and class2 all in same way
  for(i in c("male", "female", "class1", "class2")){
    dyad_cov[,i] = as.numeric(node_cov[dyad_cov[,1],i] & node_cov[dyad_cov[,2],i])
  }  
  # house distance
  dyad_cov[,"housedist"] = sqrt((node_cov[dyad_cov[,1],"house_xpos"] - node_cov[dyad_cov[,2],"house_xpos"])^2 + (node_cov[dyad_cov[,1],"house_ypos"] - node_cov[dyad_cov[,2],"house_ypos"])^2)  
  # age
  dyad_cov[,"agediff"] = abs(node_cov[dyad_cov[,1], "age"] - node_cov[dyad_cov[,2],"age"])

  # fix up node ids
  dyad_cov[,1:2] =  dyad_cov[,1:2] + (first_index - 1)
  
#   # generate dyadic cov matrix (using BuildX function which is actually pretty useless)
#   cov_name = colnames(node_cov)
#   dyad_cov = BuildX(node_cov, 
#                     unaryCol = c(which(cov_name == "household"), which(cov_name == "age")), 
#                     unaryFunc = c("absdiff", "absdiff"), 
#                     binaryCol = list(c(which(cov_name == "house_xpos"), which(cov_name == "house_ypos"))), 
#                     binaryFunc = c("euclidean") )
#   
#   # add extra columns to dyadic covariate matrix (need this since in BuildX "match" doesn't work in unary func when column is logical)
#   for(i in c("male", "female", "class1", "class2")){
#     dyad_cov = cbind(dyad_cov, as.numeric(node_cov[dyad_cov[,1],i] & node_cov[dyad_cov[,2],i]))
#     colnames(dyad_cov)[ncol(dyad_cov)] <- i
#   }  
#   
#   # fix up node ids
#   dyad_cov[,1:2] =  dyad_cov[,1:2] + (first_index - 1)
#   
#   # dyad_cov now has columns named: 
#   # node.1,    node.2,    (Intercept), household.diff, age.diff,    house_xpos.house_ypos.L2Dist, male,       female,       class1,    class2
#   # want columns to follow Hagelooch example: 
#   # Node ID 1, Node ID 2, Household,   Classroom 1,    Classroom 2, House Distance,               Male Match, Female Match, Age Diff
#   # so reorder columns
#   dyad_cov <- dyad_cov[,c(1:4, 9:10, 5, 7:8, 6)]
#   
  return(dyad_cov)
  
}

invlogit <- function(y) return(exp(y)/(1 + exp(y)))

fast_sim_ergm <- function(dyadiccovmat, eta)  {
  n_par <- length(eta)
  # n_dyad <- N * (N - 1)/2
  # if (dim(dyadiccovmat)[1] != n_dyad) 
  #  stop("Invalid Input. Number of rows of dyadic covariate matrix should match number of possible dyads")
  if (dim(dyadiccovmat)[2] != (n_par) + 2) 
    stop("Invalid Input. Need exactly one eta parameter for each non-id column of  dyadic covariate matrix")
  
  # get log odds for each dyad
  logodds = dyadiccovmat[,2+(1:n_par)] %*% eta
  edge_prob = invlogit(logodds)
  
  # decide which dyads have edges present and return edge list matrix
  return(dyadiccovmat[edge_prob > runif(length(edge_prob)),1:2])
}

sim_within_village_edges <- function(node_cov, eta, village_vec){
  edge_list = matrix(0, nrow = 0, ncol = 2)
  for (v in village_vec){
    dcm = build_village_dcm(node_cov, village = v)
    if (nrow(dcm) > 0) edge_list = rbind(edge_list, fast_sim_ergm(dcm, eta))
  }
  return(edge_list)  
}

sim_between_village_edges <- function(village, eta = 0, node_cov){
  n_village = nrow(village)
  
  # get all pairs of villages and distances between them
  comb = comb2(n_village)
  dist = sqrt((village[comb[,1],"xpos"] - village[comb[,2],"xpos"])^2 + (village[comb[,1],"ypos"] - village[comb[,2],"ypos"])^2)  
  
  # generate random number of edges 
  n = rbinom(n_village * (n_village - 1) / 2, size = village[comb[,1],"pop"]*village[comb[,2],"pop"], prob = invlogit(dist * eta))
  
  have_edge = which(n > 0)
  
  edge_list = matrix(nrow = sum(n), ncol = 2)
  j = 1
  for (i in have_edge){
    # want sample of edges from bipartite graph with k1*k2 poss edges
    index = sample(village[comb[i,1],"pop"]*village[comb[i,2],"pop"], size = n[i], replace = FALSE) 
    # use integer division to convert to edge indices
    
    # get v1 pop, call it z. then v1 indices are 
    # index - 1 %/% z + 1
    # and v2 indices are
    # index - 1 %% z + 1
    z = village[comb[i,1],"pop"]
    edge_list[j:(j + n[i] - 1),1] = (index - 1) %% z + village[comb[i,1],"index"]
    edge_list[j:(j + n[i] - 1),2] = (index - 1) %/% z + village[comb[i,2],"index"]
    j = j + n[i]
  }
  return(edge_list)  
}