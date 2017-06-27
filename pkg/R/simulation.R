
# ************************************************************************************** #
#                                                                                        # 
#                                 Simulation                                             #
#                                                                                        #
# ************************************************************************************** #

# source("outils_poids.R")

# ------------------------------------------------- #
#       simulation d'un seul individu               #
#                                                   #

simulation_individu <- function(beta,lambda,TT){
  return(sapply(1:TT,function(k) rpois(1,lambda[sample(1:dim(beta)[1],1,prob=poids_intervalle(beta,k,TT))]/TT)))
}

# ------------------------------------------------- #
#       simulation de plusieurs individus           #
#                                                   #

simulation <- function(n,beta,lambda,TT){
  # execution de plusieurs fois la fonction "simulation_individu"
  return(t(sapply(1:n,function(i) simulation_individu(beta,lambda,TT))))
}  


