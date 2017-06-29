# trouve la partition
MAP <- function(donnees,g,lambda,beta){
  TT <- length(donnees[1,])
  n <- length(donnees[,1])
  
  poids <- ensemble_poids(beta,TT)
  h <- matrice_H(donnees,lambda,poids)
  # Reinitialisation de notre partition z 
  z <- array(0,c(n,TT,g))
  for (i in 1:n){
    for(k in 1:TT){
      z[i,k,which.max(h[i,k,])] = 1
    }
  }
  return(z)
}

# Astuce
logsum <- function(vec){
  max(vec) + log(sum(exp(vec - max(vec))))
}
 


# calcul de la vraisemblance
likelihood <- function(donnees,betaVec,lambda,n,TT,g,log=TRUE){
  ret <- 0
  for(k in 1:TT){
    e <- c(0,sapply(2:g,function(jj) betaVec[jj-1]+betaVec[g-1+jj-1]*(k/TT)))
    ss2 <- logsum(e)
    for(i in 1:n){
      dp <- sapply(1:g, function(j) dpois(donnees[i,k],lambda[j]/TT, log=TRUE))
      ss1 <- logsum(e+dp)
      ret <- ret + ss1 - ss2
    }
  }
  if(!log){
    ret <- exp(ret)
  }
  return(ret)
}

# calcul de la vraisemblance complete pour la partie des poids logistic (utilise dans p(x,z|m))
likelihoodcompleted_weight <- function(donnees,betaVec,z,n,TT,g,log=TRUE){
  z_weight <- rowSums(sapply(1:(dim(z)[1]), function(i) as.numeric(t(z[i,,]))))
  ret <- 0
  for(k in 1:TT){
    e <- c(0,sapply(2:g,function(jj) betaVec[jj-1]+betaVec[g-1+jj-1]*(k/TT)))
    ss2 <- logsum(e)
    ret <- ret + sum(sapply(1:g,function(j)  z_weight[g*(k-1)+j]*(e[j] - ss2) ))
  }
  if(!log){
    ret <- exp(ret)
  }
  return(ret)
}




