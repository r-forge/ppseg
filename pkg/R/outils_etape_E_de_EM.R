
# ************************************************************************************** #
#                                                                                        # 
#                                 Etape (E) de EM                                        #
#                                                                                        #
# ************************************************************************************** #


# ------------------------------------------------- #
#   Calcul du tableau H[i,k,j] = E( Z[i,k,j] | X )  #
#                                                   #
matrice_H <- function(donnees,lambda,poids){
  TT <- length(donnees[1,])
  H <- array(0,c(n,TT,g))
  
  for (i in 1:nrow(donnees)){
    tmp <- sapply(1:length(lambda), function(j) poids[,j] * dpois(donnees[i,], lambda[j]/TT))
    H[i,,] <- sweep(tmp, 1, rowSums(tmp), "/")
    
  }
  return(H)
}
# ------------------------------------------------- #



#  w2 <- rowSums(sapply(1:nrow(donnees),function(i){
#    tmp <- sapply(1:length(lambda), function(j) poids[,j] * dpois(donnees[i,], lambda[j]/T))
#    as.numeric(t(sweep(tmp, 1, rowSums(tmp), "/")))
# })) 
# w2 - w
