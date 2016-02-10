#California network simulations

require(epinet)
require(ggplot2)
require(network)

#setwd("~/Dropbox/Princeton/1904schoolmeasles")
data <- read.csv('immunization-levels-california.csv')

test <- subset(data,Enrollment > 100)
test <- subset(test, MMR.immunization.rate < 0.70)
LA <- subset(test, District == "LOS ANGELES UNIFIED")
keep <- c('County','District','Enrollment','MMR.immunization.rate', 'lat','lon')

LA <-LA[keep]

ggplot(LA, aes(lon,lat,size=Enrollment))+geom_point(aes(colour=MMR.immunization.rate))+
  scale_size_continuous(range = c(5,15))
LA$scale_lat <- 69*(LA$lat - mean(LA$lat)) #*(max(LA$lat)-min(LA$lat))
LA$scale_lon <- 58*(LA$lon - mean(LA$lon)) #(max(LA$lon)-min(LA$lon))*

ggplot(LA, aes(scale_lon,scale_lat,size=Enrollment))+geom_point(aes(colour=MMR.immunization.rate))+
  scale_size_continuous(range = c(5,15))



LA_susc <- LA$Enrollment*(1-LA$MMR.immunization.rate)
LA_susc_pop <-sum(LA_susc)
LA_idlist <- 1:sum(LA_susc)
last_id <-0
old_mat <- rep(NULL, 4)
for (i in 1:dim(LA)[1]){
  #get school info
  temp_row <- LA[i,]
  #calculate number of susceptibles
  temp_susc <- temp_row['Enrollment']*(1-temp_row['MMR.immunization.rate'])
  #generate child IDs
  temp_ids <- last_id + 1:temp_susc[,'Enrollment']
  #update last for next time
  last_id <-tail(temp_ids,1)
  #get latitude
  temp_lat <- jitter(rep(temp_row[,'scale_lat'],temp_susc[,'Enrollment']),temp_susc[,'Enrollment'])
  #get longitude
  temp_lon <- jitter(rep(temp_row[,'scale_lon'],temp_susc[,'Enrollment']),temp_susc[,'Enrollment'])
  #get school ID 
  school_id <- rep(i,temp_susc[,'Enrollment'])
  # school_id,
  temp_mat <- cbind(temp_ids, temp_lat, temp_lon,school_id)
  old_mat <- rbind(old_mat, temp_mat)  
}
mycov <- old_mat
#generate dyadic cov matrix unaryCol = c(2), unaryFunc = c("match"), 
dyadCov <- BuildDyadicCovMatrix(mycov, unaryCol = c(4), unaryFunc = c("match"), 
                                binaryCol = list(c(2,3)), binaryFunc = c("euclidean"))


N = last_id
eta_vec = c(0,3.5,-5)
test_network1 <- BuildDyadicLinearERGM(N=N, dyadiccovmat = dyadCov, eta = eta_vec)

pdf(file = 'test_network1.pdf')
plot.network(as.network(test_network1))
dev.off()

epi10 <- SEIR.simulator(M = test_network1, N=N, beta=10/13,ki = 13, thetai = 1, latencyperiod = 0)
epi40 <- SEIR.simulator(M = test_network1, N=N, beta=40/13,ki = 13, thetai = 1, latencyperiod = 0)
epi99 <- SEIR.simulator(M = test_network1, N=N, beta=99/13,ki = 13, thetai = 1, latencyperiod = 0)

pdf(file ='epi10.pdf')
plotepitree(epi10)
dev.off() 
pdf(file ='epi40.pdf')
plotepitree(epi40)
dev.off()
pdf(file ='epi99.pdf')
plotepitree(epi99)
dev.off()
# eta_vec = c(0,-.25,-4)
# test_network4 <- BuildDyadicLinearERGM(N=N, dyadiccovmat = dyadCov, eta = eta_vec)
# 
# epi410 <- SEIR.simulator(M = test_network4, N=N, beta=10/13,ki = 13, thetai = 1, ke = 0, latencydist = "gamma")
# epi440 <- SEIR.simulator(M = test_network4, N=N, beta=40/13,ki = 13, thetai = 1, ke = 0, latencydist = "gamma")
# epi499 <- SEIR.simulator(M = test_network4, N=N, beta=99/13,ki = 13, thetai = 1, ke = 0, latencydist = "gamma")
# 
# eta_vec = c(0,-.25,-10)
# test_network10 <- BuildDyadicLinearERGM(N=N, dyadiccovmat = dyadCov, eta = eta_vec)
# 
# epi1010 <- SEIR.simulator(M = test_network10, N=N, beta=10/13,ki = 13, thetai = 1, ke = 0, latencydist = "gamma")
# epi1040 <- SEIR.simulator(M = test_network10, N=N, beta=40/13,ki = 13, thetai = 1, ke = 0, latencydist = "gamma")
# epi1099 <- SEIR.simulator(M = test_network10, N=N, beta=99/13,ki = 13, thetai = 1, ke = 0, latencydist = "gamma")
# 
# eta_vec = c(0,-.25,-50)
# test_network50 <- BuildDyadicLinearERGM(N=N, dyadiccovmat = dyadCov, eta = eta_vec)
# 
# epi5010 <- SEIR.simulator(M = test_network50, N=N, beta=10/13,ki = 13, thetai = 1, ke = 0, latencydist = "gamma")
# epi5040 <- SEIR.simulator(M = test_network50, N=N, beta=40/13,ki = 13, thetai = 1, ke = 0, latencydist = "gamma")
# epi5099 <- SEIR.simulator(M = test_network50, N=N, beta=99/13,ki = 13, thetai = 1, ke = 0, latencydist = "gamma")
# 


save.image("~/CaliSpace_ALL.RData")

