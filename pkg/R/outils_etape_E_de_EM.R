
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
  n <- length(donnees[,1])
  g <- length(lambda)
  H <- array(0,c(n,TT,g))

  for (i in 1:nrow(donnees)){
    tmp <- sapply(1:length(lambda), function(j) poids[,j] * dpois(donnees[i,], lambda[j]/TT))
    H[i,,] <- sweep(tmp, 1, rowSums(tmp), "/")

  }
  return(H)
}
# ------------------------------------------------- #

matrice_Hreg <- function(donnees,cova,alpha,poids,reg=TRUE){
  TT <- length(donnees[1,])
  n <- length(donnees[,1])
  g <- length(alpha[,1])

  H <- array(0,c(n,TT,g))

  for (i in 1:nrow(donnees)){
    tmp <- sapply(1:g, function(k){
            if(reg==TRUE){
              sapply(1:TT, function(j) poids[j,k] * dpois(donnees[i,j], exp(sum(alpha[k,]*c(1,cova[i,j,])))/TT) )
            }else{
              sapply(1:TT, function(j) poids[j,k] * dpois(donnees[i,j], exp(alpha[k,1])/TT) )
            }
            })
    H[i,,] <- sweep(tmp, 1, rowSums(tmp), "/")
  }
  return(H)
}

