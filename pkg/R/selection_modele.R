
# ************************************************************************************** #
#                                                                                        # 
#                       Selection du modele (nombre de groupe)                           #
#                                                                                        #
# ************************************************************************************** #


# ************************************************************************************** # 
# ************************************************************************************** # 
#                         FONCTIONS SUPPLEMENTAIRES DU  28/06                            #
# ************************************************************************************** # 
# ************************************************************************************** # 

# ==================================================== #
# A FAIRE FOCNTION
# Fixe hyper parameters

# p(x,htheta|m) = p(x|m,htheta)p(htheta|m)
# ==================================================== #

hyperparameters <- function(donnees,n){
  moyenne <- mean(donnees)
  variance <- sd(donnees)**2
  return(list(logistic = list(mu = 0 , sigma = 10^2),
              poisson = list(a = moyenne**2 / variance , b = moyenne/variance)
              ))
}

# Calcul densite priors pour parameter poisson p(lambda|m)
dpoissonprior <- function(vec,a,b,log=TRUE){
  ret <- sum(dgamma(vec,a,scale=b,log=TRUE))
  if(!log){
    ret <- exp(ret)
  }
  return(ret)
}

# Calcul densite priors pour parameter logistic p(beta|m)
dlogisticprior <- function(vec,mu,sigma,log=TRUE){
  ret <- sum(dnorm(vec,mu,sigma,log=TRUE))
  if(!log){
    ret <- exp(ret)
  }
  return(ret)
}

# rImportancepx
rImportancepx <- function(hbetaVec,hlambda,n){
  return(list(beta = rnorm(length(hbetaVec),hbetaVec,(10/sqrt(n))),
              lambda = rtnorm(length(hlambda),hlambda,(10/sqrt(n)),lower=0,upper=Inf) ))
}


# dImportancepx
dImportancepx <- function(betaVec,lambda,hbetaVec,hlambda,n,log=TRUE){
  ret <- sum(dnorm(betaVec,hbetaVec,(10/sqrt(n)),log=TRUE)) + sum(dtnorm(lambda,hlambda,(10/sqrt(n)),log=TRUE,lower=0,upper=Inf))
  if(!log){
    ret <- exp(ret)
  }
  return(ret)
}

# rImportancepxz
rImportancepxz <- function(hbetaVec,n){
  return( list( beta = rnorm(length(hbetaVec),hbetaVec,(10/sqrt(n))) ) )
}

# dImportancepxz
dImportancepxz <- function(betaVec,hbetaVec,n,log=TRUE){
  ret <- sum(dnorm(betaVec,hbetaVec,(10/sqrt(n)),log=TRUE)) 
  if(!log){
    ret <- exp(ret)
  }
  return(ret)
}

# rImportancepz 
rImportancepz <- function(hpoids,n,TT,g){
  sim_z <- array(0,c(n,TT,g))
  for(i in 1:n){
    for(k in 1:TT){
      sim_z[i,k,sample(g,1,prob=hpoids[k,])] <- 1
    }
  }
  return(sim_z)
}

# dImportancepz
dImportancepz <- function(z,hpoids,n,TT,g,log=TRUE){
  ret <- 0
  for(i in 1:n){
    for(k in 1:TT){
      ret <- ret + log(hpoids[k,which(z[i,k,]==1)])  
    }
  }
  if(!log){
    ret <- exp(ret)
  }
  return(ret)
}

# Calcul de p(x,theta|m)
partOnepx <- function(donnees,betaVec,lambda,hyparameters,n,TT,g,log=TRUE){
  ret <- likelihood(donnees,betaVec,lambda,n,TT,g,log=TRUE) + dlogisticprior(betaVec,hyparameters$logistic$mu,hyparameters$logistic$sigma,log=TRUE) + dpoissonprior(lambda,hyparameters$poisson$a,hyparameters$poisson$b,log=TRUE) 
  if(!log){
    ret <- exp(ret)
  }
  return(ret) 
}


# Calcul de p(x,z,theta|m)
partOnepxz <- function(donnees,betaVec,zMAP,hyparameters,n,TT,g,log=TRUE){
  # partie sur les poids 
  ret <- likelihoodcompleted_weight(donnees,betaVec,zMAP,n,TT,g,log=TRUE) + dlogisticprior(betaVec,hyparameters$logistic$mu,hyparameters$logistic$sigma,log=TRUE) 
  if(!log){
    ret <- exp(ret)
  }
  return(ret) 
}

# OneIterImportancepx
OneIterImportancepx <- function(donnees,betaVec,lambda,hbetaVec,hlambda,hyparameters,n,TT,g,log=TRUE){
  ret <- partOnepx(donnees,betaVec,lambda,hyparameters,n,TT,g,log=TRUE) - dImportancepx(betaVec,lambda,hbetaVec,hlambda,n,log=TRUE)
  if(!log){
    ret <- exp(ret)
  }
  return(ret) 
}

# OneIterImportancepxz
OneIterImportancepxz <- function(donnees,betaVec,hbetaVec,zMAP,hyparameters,n,TT,g,log=TRUE){
  ret <- partOnepxz(donnees,betaVec,zMAP,hyparameters,n,TT,g,log=TRUE) - dImportancepxz(betaVec,hbetaVec,n,log=TRUE)
  if(!log){
    ret <- exp(ret)
  }
  return(ret) 
}

# --------------------------------------------------------- #
#                     p(x,z|m)                              #
# --------------------------------------------------------- #

# Methode 1 : IS + partie exacte poisson
Integratedlikelihoodcompleted <- function(S=100,donnees,hbeta,zMAP,n,TT,g){
  # partie sur les poids 
  hbetaVec <- as.numeric(hbeta[-1,])
  hyparameters <- hyperparameters(donnees,n)
  stock <- sapply(1:S, function(u){
    s <- rImportancepxz(hbetaVec,n)
    OneIterImportancepxz(donnees,s$beta,hbetaVec,zMAP,hyparameters,n,TT,g,log=TRUE)
  })
  p1 <- logsum(stock) - log(S)
  # partie sur les composantes
  # p2 <- estimation_LVC_modele_composante(donnees,zMAP,n,TT,g)
  
  a <- hyparameters$poisson$a
  b <- hyparameters$poisson$b
  p2 <- sum(sapply(1:g, function(j){
    s1 <- sum(sapply(1:n, function(i) sum(sapply(1:TT, function(k) donnees[i,k]*zMAP[i,k,j] )) ))
    s2 <- sum(sapply(1:n, function(i) sum(sapply(1:TT, function(k) zMAP[i,k,j] )) ))
    s3 <- sum(sapply(1:n, function(i) sum(sapply(1:TT, function(k) zMAP[i,k,j]*lgamma(donnees[i,k]+1) )) ))
    a*log(b) + lgamma(s1+a) -  lgamma(a) - (s1 + a)*log(b+s2) - s3
  }))
  
  return(p1 + p2)
}

# Methode 2 : BIClogistic + partie exacte poisson
Integratedlikelihoodcompleted_BIC <- function(donnees,hbeta,zMAP,n,TT,g){
  # partie sur les poids 
  hbetaVec <- as.numeric(hbeta[-1,])
  hyparameters <- hyperparameters(donnees,n)

  p1 <- likelihoodcompleted_weight(donnees,hbetaVec,zMAP,n,TT,g,log=TRUE) - (g-1)*log(n*TT)
  # partie sur les composantes
  # p2 <- estimation_LVC_modele_composante(donnees,zMAP,n,TT,g)
  
  a <- hyparameters$poisson$a
  b <- hyparameters$poisson$b
  p2 <- sum(sapply(1:g, function(j){
    s1 <- sum(sapply(1:n, function(i) sum(sapply(1:TT, function(k) donnees[i,k]*zMAP[i,k,j] )) ))
    s2 <- sum(sapply(1:n, function(i) sum(sapply(1:TT, function(k) zMAP[i,k,j] )) ))
    s3 <- sum(sapply(1:n, function(i) sum(sapply(1:TT, function(k) zMAP[i,k,j]*lgamma(donnees[i,k]+1) )) ))
    a*log(b) + lgamma(s1+a) -  lgamma(a) - (s1 + a)*log(b+s2) - s3
  }))
  
  return(p1 + p2)
}


# --------------------------------------------------------- #
#                     p(x|m)                                #
# --------------------------------------------------------- #

# Methode 1 : IS
Integratedlikelihood <- function(S=100,donnees,hbeta,hlambda,n,TT,g){
  hbetaVec <- as.numeric(hbeta[-1,])
  hyparameters <- hyperparameters(donnees,n)
  stock <- sapply(1:S, function(u){
    s <- rImportancepx(hbetaVec,hlambda,n)
    OneIterImportancepx(donnees,s$beta,s$lambda,hbetaVec,hlambda,hyparameters,n,TT,g,log=TRUE)
  })
  return(logsum(stock) - log(S))
}

# Methode 2 : BIC
BIC <- function(donnees,hbeta,hlambda,n,TT,g){
  hbetaVec <- as.numeric(hbeta[-1,])
  return(likelihood(donnees,hbetaVec,hlambda,n,TT,g,log=TRUE) - ((g-1)*2+g)*log(n*TT)/2) 
}

# Methode 3 : IS sur z avec methode 1 de p(x,z|m)
Integratedlikelihood_z1 <- function(S=100,donnees,hbeta,hpoids,n,TT,g){
  hbetaVec <- as.numeric(hbeta[-1,])
  hyparameters <- hyperparameters(donnees,n)
  stock <- sapply(1:S, function(u){
    s_z <- rImportancepz(hpoids,n,TT,g)
    Integratedlikelihoodcompleted(100,donnees,hbeta,s_z,n,TT,g) - dImportancepz(s_z,hpoids,n,TT,g,log=TRUE)
  })
  return(logsum(stock) - log(S))
}

# Methode 3 : IS sur z avec methode 2 de p(x,z|m)
Integratedlikelihood_z2 <- function(S=100,donnees,hbeta,hpoids,n,TT,g){
  hbetaVec <- as.numeric(hbeta[-1,])
  hyparameters <- hyperparameters(donnees,n)
  stock <- sapply(1:S, function(u){
    s_z <- rImportancepz(hpoids,n,TT,g)
    Integratedlikelihoodcompleted_BIC(donnees,hbeta,s_z,n,TT,g) - dImportancepz(s_z,hpoids,n,TT,g,log=TRUE)
  })
  return(logsum(stock) - log(S))
}



# ************************************************************************************** #        
# ************************************************************************************** # 
#                        ANCIENNES FONCTIONS                                             #
# ************************************************************************************** # 
# ************************************************************************************** # 












# **************************************************************************** # 
#                         p( x , z | m )
# **************************************************************************** # 
# Integratedlikelihood <- function(S=100,donnees,hbeta,hlambda,n,TT,g){
#   hbetaVec <- as.numeric(hbeta[-1,])
#   hyparameters <- hyperparameters(donnees,n)
#   stock <- sapply(1:S, function(u){
#                                   s <- rImportancepx(hbetaVec,hlambda,n)
#                                   OneIterImportancepx(donnees,s$beta,s$lambda,hbetaVec,hlambda,hyparameters,n,TT,g,log=TRUE)
#   })
#   return(logsum(stock) - log(S))
# }
# 
# Integratedlikelihoodcompleted <- function(S=100,donnees,hbeta,zMAP,n,TT,g){
#   # partie sur les poids 
#   hbetaVec <- as.numeric(hbeta[-1,])
#   hyparameters <- hyperparameters(donnees,n)
#   stock <- sapply(1:S, function(u){
#                                   s <- rImportancepxz(hbetaVec,n)
#                                   OneIterImportancepxz(donnees,s$beta,hbetaVec,zMAP,hyparameters,n,TT,g,log=TRUE)
#                                   })
#   p1 <- logsum(stock) - log(S)
#   # partie sur les composantes
#   # p2 <- estimation_LVC_modele_composante(donnees,zMAP,n,TT,g)
#   
#   a <- hyparameters$poisson$a
#   b <- hyparameters$poisson$b
#   p2 <- sum(sapply(1:g, function(j){
#     s1 <- sum(sapply(1:n, function(i) sum(sapply(1:TT, function(k) donnees[i,k]*zMAP[i,k,j] )) ))
#     s2 <- sum(sapply(1:n, function(i) sum(sapply(1:TT, function(k) zMAP[i,k,j] )) ))
#     s3 <- sum(sapply(1:n, function(i) sum(sapply(1:TT, function(k) zMAP[i,k,j]*lgamma(donnees[i,k]+1) )) ))
#     a*log(b) + lgamma(s1+a) -  lgamma(a) - (s1 + a)*log(b+s2) - s3
#   }))
#   
#   return(p1 + p2)
# }

estimation_LVC_modele_poids <- function(S=100,hbeta,z,n=1,TT,g){
  hbetaVec <- as.numeric(hbeta[-1,]) 
  w_2 <- rowSums(sapply(1:(dim(z)[1]), function(i) as.numeric(t(z[i,,]))))
  
  stock <- rep(0,S)
  for(u in 1:S){
    # --- Generation du \theta_s
    b <- rnorm(length(hbetaVec),hbetaVec,(10/sqrt(n)))
    # --- I(\theta_s) ---
    db <- 0
    for(jj in 1:length(hbetaVec)){
      db <- db + log(dnorm(b[jj],hbetaVec[jj],(10/sqrt(n))))
    }
    # -------------------
    # LVC de \theta_s
    
    a1 <- sum(sapply(1:TT,function(k){
      tmp_2 <- c(0,sapply(2:g, function(j) b[j-1]+b[(g-1)+(j-1)]*(k/TT)))
      m <- max(tmp_2)
      s1 <- sum(sapply(1:g,function(j) w_2[g*(k-1)+j]))
      sum(sapply(2:g,function(j) (b[j-1]+b[(g-1)+(j-1)]*(k/TT))*w_2[g*(k-1)+j])) - m*s1 -  s1*log(sum(sapply(1:g,function(j) exp(tmp_2[j]-m))))
    }))
    
#     a1 <- sum(sapply(1:TT,function(k){
#       # Calcul du denominateur
#       tmp <- c(1,sapply(2:g, function(jj) exp(hbetaVec[jj-1]+hbetaVec[(g-1)+(jj-1)]*(k/TT))))
#       s <- sum(tmp)
#       
#       sum(w_2[g*(k-1)+c(1:g)] * log(tmp/s))
#     }))
    
    
    a2 <- 0
    for(jj in 1:length(hbetaVec)){
      a2 <- a2 + log(dnorm(b[jj],0,10^2))
    }
    
    lvb <- a1 + a2
    # --------------------
    
    stock[u] <- lvb - db
  }
  # Astuce
  R <- logsum(stock) - log(S)
  # log ( estimation de notre vraisemblance complété du modèle )
  return(R)
}

estimation_LVC_modele_composante <- function(donnees,z,n,TT,g){
  moyenne <- mean(donnees)
  variance <- sd(donnees)**2
  a <- moyenne**2 / variance
  b <- moyenne / variance
  res <- 0
  for(j in 1:g){
    s1 <- 0
    for(i in 1:n){
      for(k in 1:TT){
        s1 <- s1 + donnees[i,k]*z[i,k,j]
      }
    }
    s2 <- 0
    for(i in 1:n){
      for(k in 1:TT){
        s2 <- s2 + z[i,k,j]
      }
    }
    
    c1 <- a*log(b)
    c2 <- lgamma(s1+a)
    c3 <- - lgamma(a)
    c4 <- - (s1 + a)*log(b+s2)
    c5 <- 0
    for(i in 1:n){
      for(k in 1:TT){
        c5 <- c5 - z[i,k,j]*lgamma(donnees[i,k]+1)
      }
    }
    
    res <- res + c1 + c2 + c3 + c4 + c5  
  }
  return(res)
}

# --- Simplification d'utilisation
log_critere_ICL <- function(S=100,hbeta,donnees,z,n,TT,g){
  return( estimation_LVC_modele_poids(S,hbeta,z,n,TT,g) + estimation_LVC_modele_composante(donnees,z,n,TT,g) )
}
# **************************************************************************** # 
#                         quand g = 1
# **************************************************************************** # 

log_critere_g1 <- function(donnees){
  TT <- length(donnees[1,])          # T est le nombre d intervalles dans [0;1] 
  n <- length(donnees[,1])          # n est le nombre d individus 
  
  moyenne <- mean(donnees)
  variance <- sd(donnees)**2
  a <- moyenne**2 / variance
  b <- moyenne / variance
  c2 <- 0
  for(i in 1:n){
    for(k in 1:TT){
      c2 <- c2 - lgamma(donnees[i,k]+1)
    }
  }
  c <- sum(donnees)
  return(- c*log(TT) + c2 + a*log(b) - (c+a)*log(n+b) + lgamma(c+a) - lgamma(a))
}
  
# **************************************************************************** # 
#                         p( x | m )
# **************************************************************************** # 

log_critere_BIC <- function(S=100,donnees,hbeta,hlambda,g){
  TT <- length(donnees[1,])          # T est le  d intervalles dans [0;1] 
  n <- length(donnees[,1])          # n est le nombre d individus 
  
  # --- Pour prior de lambda ---
  moyenne <- mean(donnees)
  variance <- sd(donnees)**2
  h <- moyenne**2 / variance
  l <- moyenne / variance
  # -------- 
  
  
  hbetaVec <- as.numeric(hbeta[-1,]) 
  
  stock <- rep(0,S)
  for(u in 1:S){
    # --- Generation du \theta_s
    b <- rnorm(length(hbetaVec),hbetaVec,(10/sqrt(n)))
    c <- rtnorm(length(hlambda),hlambda,(10/sqrt(n)),lower=0,upper=Inf) 
    
    # --- log( I(\theta_s) ) ---
#     db <- sum(dnorm(b,hbetaVec,(10/sqrt(n)),log=TRUE))
#     for(jj in 1:length(hbetaVec)){
#       db <- db + dnorm(b[jj],hbetaVec[jj],(10/sqrt(n)),log=TRUE)
#     }
#     dc <- 0
#     for(jj in 1:length(hlambda)){
#       dc <- dc + dtnorm(c[jj],hlambda[jj],(10/sqrt(n)),log=TRUE,lower=0,upper=Inf)
#     }
    d_I <- sum(dnorm(b,hbetaVec,(10/sqrt(n)),log=TRUE)) + sum(dtnorm(c,hlambda,(10/sqrt(n)),log=TRUE,lower=0,upper=Inf))
    
    # -------------------
    # LV de \theta_s
    
    a1 <- 0
    for(k in 1:TT){
      e <- c(0,sapply(2:g,function(jj) b[jj-1]+b[g-1+jj-1]*(k/TT)))
      m2 <- max(e)
      ss2 <- 0
      for(jj in 1:g){
        ss2 <- ss2 + exp(e[jj] - m2)
      }
      ss2 <- log(ss2)
      
      dp <- rep(0,g)
      m1 <- 0
      for(i in 1:n){
        for(j in 1:g){
          dp[j] <- dpois(donnees[i,k],c[j]/TT, log=TRUE)
        }
        m1 <- max(e+dp)
        
        ss1 <- 0
        for(j in 1:g){
          ss1 <- ss1 + exp((e+dp)[j]-m1)
        }
        ss1 <- log(ss1)
        
        a1 <- a1 + m1 - m2 + ss1 - ss2
      }
    }
    # -------------------
    
    a2 <- 0
    for(jj in 1:length(hbetaVec)){
      a2 <- a2 + dnorm(b[jj],0,10^2,log=TRUE)
    }
    for(jj in 1:length(hlambda)){
      a2 <- a2 + dgamma(c[jj],h,scale=l,log=TRUE)
    }
    
    lvb <- a1 + a2
    # --------------------
    
    stock[u] <- lvb - d_I
  }
  # Astuce
  R <- max(stock)  
  somme <- 0
  for(u in 1:S){
    somme <- somme + exp(stock[u]-R)
  }
  ret <- R + log(somme) - log(S)
  # log ( estimation de notre vraisemblance du modèle )
  return(ret)
}  




# **************************************************************************** # 



  # bouge
# **************************** # 
#    Hessienne                 #

  
hessienne <- function(beta,w,TT,g){
  betaVec <- as.numeric(beta[-1,])
  hess <- matrix(0,2*g,2*g)
  for(k in 1:TT){
    w_bis <- sum(w[g*(k-1)+c(1:g)])
    element_exp <- rep(0,g)
    element_exp[1] <- 1
    for(jj in 2:g){
      element_exp[jj] <- exp(betaVec[jj-1]+betaVec[g-1+jj-1]*(k/TT))
    } 
    s <- sum(element_exp)
    s2 <- sum(element_exp)**2
    
    for(j in 1:g){
      # elements diagonaux
      hess[j,j] <- hess[j,j] - w_bis * element_exp[j] * (s - element_exp[j]) / s2 
      hess[j,g+j] <- hess[j,g+j] - (k/TT) * w_bis * element_exp[j] * (s - element_exp[j]) / s2 
      hess[g+j,g+j] <- hess[g+j,g+j] - (k/TT)**2 * w_bis * element_exp[j] * (s - element_exp[j]) / s2 
      # ------------------
      # autres elements
      if(j+1!=g+1){
        for(l in (j+1):g){
          hess[j,l] <- hess[j,l] + w_bis * element_exp[j] * element_exp[l] / s2
          hess[j,g+l] <- hess[j,g+l] + (k/TT) * w_bis * element_exp[j] * element_exp[l] / s2
          hess[g+j,g+l] <- hess[g+j,g+l] - (k/TT)**2 * w_bis * element_exp[j] * element_exp[l] / s2
        }
      }
      # ----------------- 
    }
  }
  # Comme hess est symetrique
  for(j in 1:g){
    hess[g+j,j] <- hess[j,g+j]
    # ------------------
    # autres elements
    if(j+1!=g){
      for(l in (j+1):g){
        hess[l,j] <- hess[j,l]
        hess[l,g+j] <- hess[j,g+l]
        hess[g+l,j] <- hess[j,g+l]
        hess[g+j,l] <- hess[j,g+l]
        hess[g+j,g+l] <- hess[g+l,g+j]
      }
    }
    # ----------------- 
  }
  return(hess)
} 
  


