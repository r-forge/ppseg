
# ************************************************************************************** #
#                                                                                        # 
#             Calcul des poids d'un intervalle & Calcul ensemble des poids               #
#                                                                                        #
# ************************************************************************************** #


# ------------------------------------------------- #
#           poids d'un intervalle k                 #
#                                                   #

poids_intervalle <- function(beta,k,TT){
  m <- exp(rowSums(cbind(beta[,1], beta[,2]*k/TT)))
  return(m /sum(m))
}

poids_intervalle_2 <- function(betaVec,k,TT,g){
  tmp <- c(0,sapply(2:g, function(j) betaVec[j-1]+betaVec[(g-1)+(j-1)]*(k/TT)))
  m <- max(tmp)
  return(sapply(1:g,function(j) exp(tmp[j]-m-log(sum(sapply(1:g,function(jj) exp(tmp[jj]-m)))))))
}




# ------------------------------------------------- #
#           poids de chaque intervalle              #
#                                                   #

ensemble_poids <- function(beta,TT){
  # execution de plusieurs fois la fonction "poids_intervalle"
  return(t(sapply(1:TT, poids_intervalle, beta=beta, TT=TT)))
}

ensemble_poids_2 <- function(betaVec,TT,g){
  # execution de plusieurs fois la fonction "poids_intervalle"
  return(t(sapply(1:TT, poids_intervalle_2, betaVec=betaVec, TT=TT, g=g)))
}

