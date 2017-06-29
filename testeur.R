rm(list=ls())
set.seed(123)
require(ppseg)
require(parallel)
require(msm)
data("BritishCoalMiningDisasterData")
x <- matrix(x$count, nrow = 1)

# donnees <- x
# g <- 2
# res <- ppsegestim(x, test_group = 2, S = 3, nb_tests = 3)
#resEM <- selection_EM(x, 2, nb_tests = 20)

load("~/Bureau/tmpJulien.rda")
donnees <- x
g <- 2
Integratedlikelihoodcompleted(S=5000, donnees, resEM)
