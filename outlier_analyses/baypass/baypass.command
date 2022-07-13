baypass -gfile bft.kinless.baypass -nthreads 20 -nval 10000 -burnin 10000 -pilotlength 1000

# running baypass on the simulations
baypass -gfile G.bft.BP.sim -nthreads 20 -nval 10000 -burnin 10000 -pilotlength 1000

# population names/numbers
1 <- BRZ
2 <- KEY
3 <- MRT
4 <- PNS
5 <- PR
6 <- SCA
7 <- TX
8 <- VZ
