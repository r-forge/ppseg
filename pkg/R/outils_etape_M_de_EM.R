 
# ************************************************************************************** #
#                                                                                        # 
#                                 Etape (M) de EM                                        #
#                                                                                        #
# ************************************************************************************** #


# ------------------------------------------------------------ #
#                Calcul des lambda a l etape m                 #
#                                                              #
estimation_max_lambda <- function(n,g,TT,donnees,H,w){
  lambda <- sapply(1:g,function(j) sum(sapply(1:TT,function(k) H[1:n,k,j]*donnees[1:n,k]*TT))/sum(sapply(1:TT,function(k) w[g*(k-1)+j])) )
  return(lambda)
}

# ------------------------------------------------------------- #



# ------------------------------------------------------------ #
#    Utilisation de la fonction optim() pour calcul de beta    #
#                                                              #


log_vraisemblance_poids <- function(betaVec,w,TT,g){
  # calcul de la log vraisemblance des poids sur chaque intervalle I_k
  res <- sum(sapply(1:TT,function(k){
    # Calcul du denominateur
    tmp <- c(1,sapply(2:g, function(jj) exp(betaVec[jj-1]+betaVec[(g-1)+(jj-1)]*(k/TT))))
    s <- sum(tmp)
    
    sum(w[g*(k-1)+c(1:g)] * log(tmp/s))
  }))
  return(res)
}


# Version optimisee
log_vraisemblance_poids_2 <- function(betaVec,w,TT,g){
  # calcul de la log vraisemblance des poids sur chaque intervalle I_k  -------> avec Astuce <-------
  res <- sum(sapply(1:TT,function(k){
    tmp_2 <- c(0,sapply(2:g, function(j) betaVec[j-1]+betaVec[(g-1)+(j-1)]*(k/TT)))
    m <- max(tmp_2)
    s1 <- sum(sapply(1:g,function(j) w[g*(k-1)+j]))
    sum(sapply(2:g,function(j) (betaVec[j-1]+betaVec[(g-1)+(j-1)]*(k/TT))*w[g*(k-1)+j])) - m*s1 -  s1*log(sum(sapply(1:g,function(j) exp(tmp_2[j]-m))))
  }))
  return(res)
}


grad_log_vraisemblance_poids <- function(betaVec,w,TT,g){
  # Calcul des derivee partiel des 
  #         \beta_0 sur la 1 ere ligne
  #         \beta_1 sur la 2 eme ligne
  m_grad <- sapply(2:g,function(j){
    rowSums(sapply(1:TT,function(k){
      d1 <- sum(sapply(1:g,function(jj) w[g*(k-1)+jj])) * exp(betaVec[j-1]+betaVec[(g-1)+(j-1)]*(k/TT))
      d2 <- sum(c(1,sapply(2:g,function(jj) exp(betaVec[jj-1]+betaVec[(g-1)+(jj-1)]*(k/TT)))))
      c(w[g*(k-1)+j] - d1/d2,(k/TT)*(w[g*(k-1)+j] - d1/d2)) }))
  })
  return(c(m_grad[1,],m_grad[2,]))
}

# Version optimisee
grad_log_vraisemblance_poids_2 <- function(betaVec,w,TT,g){
  tmp_3 <- t(sapply(1:TT,function(k) c(0,sapply(2:g, function(jj) betaVec[jj-1]+betaVec[(g-1)+(jj-1)]*(k/TT)))))
  m <- sapply(1:TT,function(k) max(tmp_3[k,]))
  m_grad <- sapply(2:g,function(j){
    rowSums(sapply(1:TT,function(k){
      s1 <- tmp_3[k,j] - m[k] - log(sum(sapply(1:g,function(jj) exp(tmp_3[k,jj] - m[k]))))
      s2 <- sum(sapply(1:g,function(jj) w[g*(k-1)+jj]))
      c(w[g*(k-1)+j] - exp(s1)*s2,(k/TT)*(w[g*(k-1)+j] - exp(s1)*s2))
    }))
  }) 
  return(c(m_grad[1,],m_grad[2,]))
}


# ------------------------------------------------------------- #


log_vraisemblance_composante <- function(alphaVec,donnees,cova,H,n,TT,gg,pp){
  res <- sum(sapply(1:n, function(i){
            sum(sapply(1:TT, function(j){
              sum(sapply(1:gg, function(k){
                # a <- 7 # (alphaVec[((k-1)*p+1):(k*p)]%*%cova[k,j,])[1]
                H[i,j,k] * dpois(donnees[i,j],exp(sum(alphaVec[((k-1)*pp+1):(k*pp)]*cova[i,j,])),log=TRUE)
              }))
            })) 
          }))
  return(res)
}

grad_log_vraisemblance_composante <- function(alphaVec,donnees,cova,H,n,TT,gg,pp){
#   m_grad <- sapply(1:p, function(l){
#               lapply(1:g, function(k){
#                 sum(sapply(1:n, function(i){
#                   sum(sapply(1:TT, function(j){
#                     a <- 7 # (alphaVec[((k-1)*p+1):(k*p)]%*%cova[k,j,])[1]
#                     H[i,j,k] * ( donnees[i,j]*cova[k,j,l] - cova[k,j,l]*exp(a)/TT ) 
#                   }))
#                 }))  
#               })
#             })
  m_grad <- matrix(0,nrow = gg,ncol = pp)
#   for(l in 1:p){
#     m_grad[,l] <- sapply(1:g, function(k){
#                     sum(sapply(1:TT, function(j){
#                       sum(sapply(1:n, function(i){
#                         a <- 7 # (alphaVec[((k-1)*p+1):(k*p)]%*%cova[k,j,])[1]
#                         H[i,j,k] * ( donnees[i,j]*cova[k,j,l] - cova[k,j,l]*exp(a)/TT ) 
#                       }))
#                     }))  
#                   })
#   }
  m_grad[,1] <- sapply(1:gg, function(k){
    sum(sapply(1:TT, function(j){
      sum(sapply(1:n, function(i){
        a <- sum(alphaVec[((k-1)*pp+1):(k*pp)]*cova[i,j,])
        H[i,j,k] * ( donnees[i,j]*cova[i,j,l] - cova[i,j,l]*exp(a)/TT ) 
      }))
    }))  
  })
  return(c(m_grad))
}

















