
# ************************************************************************************** #
#                                                                                        # 
#                       Selection du modele (nombre de groupe)                           #
#                                                                                        #
# ************************************************************************************** #

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
  R <- max(stock) 
  somme <- 0
  for(u in 1:S){
    somme <- somme + exp(stock[u]-R)
  }
  R <- R + log(somme) - log(S)
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
  

log_critere_BIC <- function(S=100,donnees,hbeta,hlambda,g){
  TT <- length(donnees[1,])          # T est le nombre d intervalles dans [0;1] 
  n <- length(donnees[,1])          # n est le nombre d individus 
  
  # --- Pour prior de lambda ---
  moyenne <- mean(donnees)
  variance <- sd(donnees)**2
  h <- moyenne**2 / variance
  k <- moyenne / variance
  # -------- 
  
  
  hbetaVec <- as.numeric(hbeta[-1,]) 
  
  stock <- rep(0,S)
  for(u in 1:S){
    # --- Generation du \theta_s
    b <- rnorm(length(hbetaVec),hbetaVec,(10/sqrt(n)))
    c <- rtnorm(length(hlambda),hlambda,(10/sqrt(n)),lower=0,upper=Inf) 
    
    # --- log( I(\theta_s) ) ---
    db <- 0
    for(jj in 1:length(hbetaVec)){
      db <- db + dnorm(b[jj],hbetaVec[jj],(10/sqrt(n)),log=TRUE)
    }
    dc <- 0
    for(jj in 1:length(hlambda)){
      dc <- dc + dtnorm(c[jj],hlambda[jj],(10/sqrt(n)),log=TRUE,lower=0,upper=Inf)
    }
    d_I <- db + dc
    
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
      a2 <- a2 + dgamma(c[jj],h,scale=k,log=TRUE)
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
  


