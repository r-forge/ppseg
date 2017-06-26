 
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

# log_vraisemblance_poids_2 <- function(betaVec,w,TT,g){
#   res <- 0
#   tmp_2 <- rep(0,g)
#   for(k in 1:TT){
#     p1 <- 0
#     tmp_2[1] <- 0 
#     for(j in 2:g){
#       p1 <- p1 + (betaVec[j-1]+betaVec[(g-1)+(j-1)]*(k/TT))*w[g*(k-1)+j]
#       tmp_2[j] <- betaVec[j-1]+betaVec[(g-1)+(j-1)]*(k/TT)
#     }
#     m <- max(tmp_2)
#     s1 <- 0
#     for(j in 1:g){
#       s1 <- s1 +  w[g*(k-1)+j]
#     }
#     p2 <- 0
#     p2 <- m*s1
#     s2 <- 0
#     for(j in 1:g){
#       s2 <- s2 + exp(tmp_2[j]-m)
#     }
#     p2 <- p2 + s1*log(s2)
#     res <- res + p1 - p2
#   }
#   return(res)
# }

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

# grad_log_vraisemblance_poids_2 <- function(betaVec,w,TT,g){
#   res <- rep(0,2*(g-1))
#   tmp_3 <- matrix(0,TT,g)
#   for(k in 1:TT){
#     for(j in 2:g){
#       tmp_3[k,j] <- betaVec[j-1]+betaVec[(g-1)+(j-1)]*(k/TT) 
#     }
#     m[k] <- max(tmp_3[k,])
#   }
#   for(j in 2:g){
#     for(k in 1:TT){
#       s1 <- 0
#       ss1 <- 0
#       ss2 <- 0
#       for(jj in 1:g){
#         ss1 <- ss1 + exp(tmp_3[k,jj] - m[k])
#         ss2 <- ss2 + w[g*(k-1)+jj]
#       } 
#       s1 <- tmp_3[k,j] - m[k] - log(ss1)
#       res[j-1] <- res[j-1] + w[g*(k-1)+j] - exp(s1)*ss2
#       res[g-1+j-1] <- res[g-1+j-1] + (k/TT)*(w[g*(k-1)+j] - exp(s1)*ss2)
#     }
#   } 
#   return(res)
# }

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


#  -> Sauvegarde d'une fonction grad qui fonctionne <-
#
# grad_log_vraisemblance_poids <- function(betaVec,H){
#   # --------------------
#   n <- dim(H)[1]
#   TT <- dim(H)[2]
#   g <- dim(H)[3]
#   # --------------------
#   # somme sur les i de H[i,k,j]
#   w <- rowSums(sapply(1:(dim(H)[1]), function(i) as.numeric(t(H[i,,])))) 
#   # --------------------
#   grad <- rep(0,(2*g-2))
#   for(j in 2:g){ # j=2
#     for(k in 1:TT){ # k = 1
#       s1 <-  w[g*(k-1)+1]
#       s2 <- 1
#       for(jj in 2:g){ # jj=2
#         # s1 <- s1 + w[g*(k-1)+jj]*exp(betaVec[jj-1]+betaVec[(g-1)+(jj-1)]*(k/TT))
#         s1 <- s1 + w[g*(k-1)+jj]
#         s2 <- s2 + exp(betaVec[jj-1]+betaVec[(g-1)+(jj-1)]*(k/TT))
#       }
#       s1 <- s1 * exp(betaVec[j-1]+betaVec[(g-1)+(j-1)]*(k/TT))
#       grad[j-1] <- grad[j-1] + w[g*(k-1)+j] - s1/s2
#       grad[(g-1)+j-1] <- grad[(g-1)+j-1] + (k/TT)*(w[g*(k-1)+j] - s1/s2)
#     }
#   }
#   return(grad)
# }




