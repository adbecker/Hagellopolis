source("fastHag.R")
source("simseir.R")

#Create Hagellopolis
n_pop <- 10000  #Total population
n_hags <- 10 # number of villages
hagopolis_size = 1000 # size of landscape (square sides)
hag_size = 10 # size of a single village

# get covariates for each village
ndcv = sim_node_cov(n_pop, n_hags, city_size = hagopolis_size, village_size = hag_size)

# record info about villages
village = data.frame(table(ndcv$village), unique(ndcv$village_xpos), unique(ndcv$village_ypos))
names(village) <- c("id", "pop", "xpos", "ypos")
village$index = cumsum(c(1,village$pop))[1:n_hags]

# make up eta vector
eta_vec <- -.25 * randomLHS(1, 8)
eta_vec[4]<- -.004 * n_pop / n_hags
eta_vec[8] <- -.004
eta_vec <- c(eta_vec)

# simulate edges within villages
ptm = proc.time()
within_edges = sim_within_village_edges(node_cov = ndcv, eta = eta_vec, village_vec = 1:n_hags)
proc.time() - ptm

# simulate edges between villages
ptm1 = proc.time()
btwn_edges = sim_between_village_edges(village, eta = -0.04, ndcv)
proc.time() - ptm1

# stick them together
my_ergm = rbind(within_edges, btwn_edges)

M = my_ergm
N = n_pop
beta = 0.3
ke = 2
ki = 2
thetae = 5
thetai = 5
latencydist = "fixed"
latencyperiod = 0

#plot(network(M, directed = FALSE))
pt = proc.time()
epidemic = simseir(M = my_ergm, N = n_pop, beta = 0.3, ki = 1, thetai = 5, latencydist = "fixed", latencyperiod = 1)
proc.time() - pt
#Plot Network
# plot.network(as.network(my_ergm))
# 


# 
# #Plot spread of epidemic
# plotepitree(epi)
# 
## use the profvis library for examing the speed of code
# library(profvis)
# profvis({
#   within_edges1 = sim_within_village_edges(node_cov = ndcv, eta = eta_vec, village_vec = 1:n_hags)
#   btwn_edges1 = sim_between_village_edges(village, eta = -0.4, ndcv)
# })
# 
# 
