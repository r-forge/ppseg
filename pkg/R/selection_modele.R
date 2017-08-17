
# ************************************************************************************** #
#                                                                                        #
#                       Selection du modele (nombre de groupe)                           #
#                                                                                        #
# ************************************************************************************** #

# Bloc pour matrice covariance
Bloc_diag <- function(bloc,t,betaVec,n,g){
  e <- c(0,sapply(1:(g-1), function(k) betaVec[2*(k-1)+1]+betaVec[2*(k-1)+2]*t))

  B <- matrix(c(1,t,t,t**2),2,2)
  B <- n*exp( e[bloc+1] + logsum(e[-(bloc+1)]) - 2*logsum(e) )*B
  return(B)
}


Bloc_non_diag <- function(bloc,t,betaVec,n,g){
  e <- c(0,sapply(1:(g-1), function(k) betaVec[2*(k-1)+1]+betaVec[2*(k-1)+2]*t))

  B <- matrix(c(1,t,t,t**2),2,2)
  B <- - n*exp( e[bloc[1]+1] + e[bloc[2]+1] - 2*logsum(e) )*B
  return(B)
}

Assemblage_Matcov <-function(t,betaVec,n,g){
  M2 <- matrix(0,2*(g-1),2*(g-1))
  for(k in 1:(g-1)){
    for(j in 1:(g-1)){
      if(j==k){
        M2[ (2*(k-1)+1) : (2*(k-1)+2) , (2*(j-1)+1) : (2*(j-1)+2) ] <- Bloc_diag( bloc=k ,t,betaVec,n,g)
      }else{
        M2[ (2*(k-1)+1) : (2*(k-1)+2) , (2*(j-1)+1) : (2*(j-1)+2) ] <- Bloc_non_diag( bloc=c(k,j) ,t,betaVec,n,g)
      }
    }
  }
  return(M2)
}

Matcovbeta <- function(betaVec,n,TT,g){
  M <- lapply(1:TT,function(t) Assemblage_Matcov(t/TT,betaVec,n,g))
  ret <- M[[1]]
  for(j in 2:TT){
    ret <- ret + M[[j]]
  }
  return(ret)
}

# *******************************************************************************************
# *******************************************************************************************
# *******************************************************************************************
# dImportancepxz avec estimation matrice cov
dImportancepxz_bis <- function(betaVec,hbetaVec,n,TT,g,log=TRUE){
  Mcov <- Matcovbeta(hbetaVec,n,TT,g)
  ret <- dmvnorm(betaVec,mean=hbetaVec,sigma=ginv(Mcov),log=TRUE)
  if(!log){
    ret <- exp(ret)
  }
  return(ret)
}

# rImportancepz avec estimation matrice cov
rImportancepxz_bis <- function(hbetaVec,n,TT,g){
  Mcov <- Matcovbeta(hbetaVec,n,TT,g)
  return( list( beta = rmvnorm(1,mean=hbetaVec,sigma=ginv(Mcov)) ) )
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

# OneIterImportancepxz
OneIterImportancepxz_bis <- function(donnees,betaVec,hbetaVec,zMAP,hyparameters,n,TT,g,log=TRUE){
  ret <- partOnepxz(donnees,betaVec,zMAP,hyparameters,n,TT,g,log=TRUE) - dImportancepxz_bis(betaVec,hbetaVec,n,TT,g,log=TRUE)
  if(!log){
    ret <- exp(ret)
  }
  return(ret)
}
# *******************************************************************************************
# *******************************************************************************************
# *******************************************************************************************



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

likelihoodcompleted_component <- function(donnees,lambda,z,n,TT,g,log=TRUE){
  ret <- sum(sapply(1:n,function(i){
    sum(sapply(1:TT, function(j){
      sum(sapply(1:g,function(k)  z[i,j,k]*dpois(donnees[i,j],lambda[k],log=TRUE) ))
    }))
  }))
  if(!log){
    ret <- exp(ret)
  }
  return(ret)
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

# rImportancepxz IS
rImportancepxz_IS <- function(hlambda,n){
  return( list( lambda = rtnorm(length(hlambda),hlambda,(1000/sqrt(n)),lower=0,upper=Inf) ) )
}

# dImportancepxz IS
dImportancepxz_IS <- function(lambda,hlambda,n,log=TRUE){
  ret <- sum(rtnorm(lambda,hlambda,(1000/sqrt(n)),lower=0,upper=Inf))
  if(!log){
    ret <- exp(ret)
  }
  return(ret)
}



# rImportancepz
rImportancepz <- function(donnees,hlambda,hpoids,n,TT,g){
  sim_z <- array(0,c(n,TT,g))
  h <- matrice_H(donnees,hlambda,hpoids)
  for(i in 1:n){
    for(k in 1:TT){
      sim_z[i,k,sample(g,1,prob=h[i,k,])] <- 1
    }
  }
  return(sim_z)
}


# dImportancepz
dImportancepz <- function(donnees,z,hlambda,hpoids,n,TT,g,log=TRUE){
  ret <- 0
  h <- matrice_H(donnees,hlambda,hpoids)
  for(i in 1:n){
    for(k in 1:TT){
      ret <- ret + log(h[i,k,which(z[i,k,]==1)])
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

# Calcul de p(x,z,lambda|m)
partOnepxz_IS <- function(donnees,lambda,zMAP,hyparameters,n,TT,g,log=TRUE){
  # partie sur les poids
  ret <- likelihoodcompleted_component(donnees,lambda,zMAP,n,TT,g,log=TRUE) + dpoissonprior(lambda,hyparameters$poisson$a,hyparameters$poisson$b,log=TRUE)
  if(!log){
    ret <- exp(ret)
  }
  return(ret)
}

# OneIterImportancepxz IS
OneIterImportancepxz_IS <- function(donnees,lambda,hlambda,zMAP,hyparameters,n,TT,g,log=TRUE){
  ret <- partOnepxz_IS(donnees,lambda,zMAP,hyparameters,n,TT,g,log=TRUE) - dImportancepxz_IS(lambda,hlambda,n,log=TRUE)
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

  return(p1+p2)
}

# Methode 1 : IS + IS
Integratedlikelihoodcompleted_IS <- function(S=100,donnees,hbeta,hlambda,zMAP,n,TT,g){
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
  stock <- sapply(1:S, function(u){
    s <- rImportancepxz_IS(hlambda,n)
    OneIterImportancepxz_IS(donnees,s$lambda,hlambda,zMAP,hyparameters,n,TT,g,log=TRUE)
  })
  p2 <- logsum(stock) - log(S)
  
  return(p1+p2)
}

# Methode 1 BIS : IS + partie exacte poisson
Integratedlikelihoodcompleted_bis <- function(S=100,donnees,hbeta,zMAP,n,TT,g){
  # partie sur les poids
  hbetaVec <- as.numeric(hbeta[-1,])
  hyparameters <- hyperparameters(donnees,n)
  stock <- sapply(1:S, function(u){
    s <- rImportancepxz_bis(hbetaVec,n,TT,g)
    OneIterImportancepxz_bis(donnees,s$beta,hbetaVec,zMAP,hyparameters,n,TT,g,log=TRUE)
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
BIC_select <- function(donnees,hbeta,hlambda,n,TT,g){
  hbetaVec <- as.numeric(hbeta[-1,])
  return(likelihood(donnees,hbetaVec,hlambda,n,TT,g,log=TRUE) - ((g-1)*2+g)*log(n*TT)/2)
}


# Methode 3 : IS sur z avec methode 1 de p(x,z|m)
Integratedlikelihood_z1 <- function(S=100,donnees,hbeta,hpoids,n,TT,g){
  hbetaVec <- as.numeric(hbeta[-1,])
  hyparameters <- hyperparameters(donnees,n)
  stock <- sapply(1:S, function(u){
    s_z <- rImportancepz(hpoids,n,TT,g)
    Integratedlikelihoodcompleted(100,donnees,hbeta,s_z,n,TT,g) - dImportancepz(donnees,s_z,hpoids,n,TT,g,log=TRUE)
  })
  return(logsum(stock) - log(S))
}

# Methode 4 : IS sur z avec methode 2 de p(x,z|m)
Integratedlikelihood_z2 <- function(S=100,donnees,hbeta,hpoids,n,TT,g){
  hbetaVec <- as.numeric(hbeta[-1,])
  hyparameters <- hyperparameters(donnees,n)
  stock <- sapply(1:S, function(u){
    s_z <- rImportancepz(hpoids,n,TT,g)
    Integratedlikelihoodcompleted_BIC(donnees,hbeta,s_z,n,TT,g) - dImportancepz(donnees,s_z,hpoids,n,TT,g,log=TRUE)
  })
  return(logsum(stock) - log(S))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# modification du la methode pour ne pas faire des calcul deja fait
# Methode 3 : IS sur z avec methode 1 de p(x,z|m)
Integratedlikelihood_z1 <- function(S=100,donnees,hbeta,hlambda,hpoids,zMAP,n,TT,g){
  hbetaVec <- as.numeric(hbeta[-1,])
  hyparameters <- hyperparameters(donnees,n)

  elem <- vector("list",S+1)
  elem[[1]]$z <- zMAP
  elem[[1]]$val <- Integratedlikelihoodcompleted(100,donnees,hbeta,zMAP,n,TT,g) - dImportancepz(donnees,zMAP,hlambda,hpoids,n,TT,g,log=TRUE)
  dif <- 1
  new <- TRUE

  stock <- rep(0,S)
  for(u in 1:S){
    s_z <- rImportancepz(donnees,hlambda,hpoids,n,TT,g)
    new <- TRUE
    for(ind in 1:dif){
      if(sum(s_z==elem[[ind]]$z)==TT*g){
        new <- FALSE
        stock[u] = elem[[ind]]$val
      }
    }
    if(new){
      dif <- dif + 1
      elem[[dif]]$z <- s_z
      elem[[dif]]$val <- Integratedlikelihoodcompleted(100,donnees,hbeta,s_z,n,TT,g) - dImportancepz(donnees,s_z,hlambda,hpoids,n,TT,g,log=TRUE)
      stock[u] = elem[[dif]]$val
    }
  }

  return(logsum(stock) - log(S))
}

# Methode 4 : IS sur z avec methode 2 de p(x,z|m)
Integratedlikelihood_z2 <- function(S=100,donnees,hbeta,hlambda,hpoids,zMAP,n,TT,g){
  hbetaVec <- as.numeric(hbeta[-1,])
  hyparameters <- hyperparameters(donnees,n)

  elem <- vector("list",S+1)
  elem[[1]]$z <- zMAP
  elem[[1]]$val <- Integratedlikelihoodcompleted_BIC(donnees,hbeta,zMAP,n,TT,g) - dImportancepz(donnees,zMAP,hlambda,hpoids,n,TT,g,log=TRUE)
  dif <- 1
  new <- TRUE

  stock <- rep(0,S)
  for(u in 1:S){
    s_z <- rImportancepz(donnees,hlambda,hpoids,n,TT,g)
    new <- TRUE
    for(ind in 1:dif){
      if(sum(s_z==elem[[ind]]$z)==TT*g){
        new <- FALSE
        stock[u] = elem[[ind]]$val
      }
    }
    if(new){
      dif <- dif + 1
      elem[[dif]]$z <- s_z
      elem[[dif]]$val <- Integratedlikelihoodcompleted_BIC(donnees,hbeta,s_z,n,TT,g) - dImportancepz(donnees,s_z,hlambda,hpoids,n,TT,g,log=TRUE)
      stock[u] = elem[[dif]]$val
    }
  }

  return(logsum(stock) - log(S))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
# **************************************************************************** #
# **************************************************************************** #
#
#        LASSO
#
# **************************************************************************** #
# **************************************************************************** #
# **************************************************************************** #

Fonction_Lasso <- function(lambda_Lasso,donneesVec,coef){
  ret <- dpois(donneesVec[1],lambda_Lasso[1])
  for(t in 2:length(donneesVec)){
    ret <- ret + dpois(donneesVec[t],lambda_Lasso[t]) - coef*(lambda_Lasso[t]-lambda_Lasso[t-1])**2
  }
  return(-ret)
}
Grad_Lasso <- function(lambda_Lasso,donneesVec,coef){
  TT <- length(donneesVec)
  ret <- rep(0,TT)
  ret[1] <- ((1 - lambda_Lasso[1]/donneesVec[1]) * exp(-lambda_Lasso[1]) * lambda_Lasso[1]**(donneesVec[1]-1)/factorial(donneesVec[1]-1)) + 2*(lambda_Lasso[2]- lambda_Lasso[1])
  for(t in 2:(TT-1)){
    ret[t] <- ((1 - lambda_Lasso[t]/donneesVec[t]) * exp(-lambda_Lasso[t]) * lambda_Lasso[t]**(donneesVec[t]-1)/factorial(donneesVec[t]-1)) + 2*(lambda_Lasso[t+1] + lambda_Lasso[t-1] - 2*lambda_Lasso[t])
  }
  ret[TT] <- ((1 - lambda_Lasso[TT]/donneesVec[TT]) * exp(-lambda_Lasso[TT]) * lambda_Lasso[TT]**(donneesVec[TT]-1)/factorial(donneesVec[TT]-1)) + 2*(lambda_Lasso[TT-1]- lambda_Lasso[TT])
  return(ret)
}


# optim(par=c(rep(42,50),rep(160,50)),fn=Fonction_Lasso,gr = Grad_Lasso,donneesVec=donnees,coef=20,lower = 0,method = "L-BFGS-B")
# Fonction_Lasso(c(rep(25,50),rep(100,50)),donnees,coef=100)
# Grad_Lasso(c(rep(25,50),rep(100,50)),donnees,coef=100)
# Fonction_Lasso(optim(par=rep(50,100),fn=Fonction_Lasso,donneesVec=donnees,coef=1000)$par,donnees,coef=50)

# ******************************
#       Avec Cova
# *****************************


hyperparameters_reg <- function(){
  return(list(logistic = list(mu = 0 , sigma = 10^2),
              poisson = list(mu = 0 , sigma = 10^2)
              ))
}


# *************************************************************
#          Critere de selection avec IS pour ppsegreg
# *************************************************************

likelihoodcompleted_component_reg <- function(donnees,cova,alphaVec,z,n,TT,g,log=TRUE){
  l <- length(alphaVec)/g
  ret <- sum(sapply(1:n,function(i){
            sum(sapply(1:TT, function(j){
              sum(sapply(1:g,function(k)  z[i,j,k]*dpois(donnees[i,j],exp(sum(alphaVec[(l*(k-1)+1):(l*k)]*c(1,cova[i,j,]))),log=TRUE) ))
            }))
          }))
  if(!log){
    ret <- exp(ret)
  }
  return(ret)
}

rImportancepxz_reg <- function(halphaVec,n){
  return( list( alpha = rnorm(length(halphaVec),halphaVec,(10/sqrt(n))) ) )
}

dImportancepxz_reg <- function(alphaVec,halphaVec,n,log=TRUE){
  ret <- sum(dnorm(alphaVec,halphaVec,(10/sqrt(n)),log=TRUE))
  if(!log){
    ret <- exp(ret)
  }
  return(ret)
}

dpoissonprior_reg <- function(vec,mu,sigma,log=TRUE){
  ret <- sum(dnorm(vec,mu,sigma,log=TRUE))
  if(!log){
    ret <- exp(ret)
  }
  return(ret)
}


partOnepxz_reg <- function(donnees,cova,alphaVec,zMAP,hyparameters,n,TT,g,log=TRUE){
  ret <- likelihoodcompleted_component_reg(donnees,cova,alphaVec,zMAP,n,TT,g,log=TRUE) + dpoissonprior_reg(alphaVec,hyparameters$poisson$mu,hyparameters$poisson$sigma,log=TRUE)
  if(!log){
    ret <- exp(ret)
  }
  return(ret)
}

OneIterImportancepxz_reg <- function(donnees,cova,alphaVec,halphaVec,zMAP,hyparameters,n,TT,g,log=TRUE){
  ret <- partOnepxz_reg(donnees,cova,alphaVec,zMAP,hyparameters,n,TT,g,log=TRUE) - dImportancepxz_reg(alphaVec,halphaVec,n,log=TRUE)
  if(!log){
    ret <- exp(ret)
  }
  return(ret)
}

# Methode 1 : IS
Integratedlikelihoodcompleted_reg <- function(S=1000,donnees,cova,hbeta,halpha,zMAP,n,TT,g){
  hyparameters <- hyperparameters_reg()
  # partie sur les poids
  hbetaVec <- as.numeric(hbeta[-1,])
  stock <- sapply(1:S, function(u){
    s <- rImportancepxz(hbetaVec,n)
    OneIterImportancepxz(donnees,s$beta,hbetaVec,zMAP,hyparameters,n,TT,g,log=TRUE)
  })
  p1 <- logsum(stock) - log(S)
  # partie sur les composantes
  halphaVec <- as.numeric(t(halpha))
  stock <- sapply(1:S, function(u){
    s <- rImportancepxz_reg(halphaVec,n)
    OneIterImportancepxz_reg(donnees,cova,s$alpha,halphaVec,zMAP,hyparameters,n,TT,g,log=TRUE)
  })
  p2 <- logsum(stock) - log(S)
  
  
  return(p1 + p2)
}

# ***************************
#          g = 1
# ***************************

log_critere_g1_reg <- function(S=1000,donnees,cova,halpha,n,TT){
  hyparameters <- hyperparameters_reg()
  halphaVec <- as.numeric(t(halpha))
  zMAP <- array(1,dim = c(n,TT,1))
  stock <- sapply(1:S, function(u){
    s <- rImportancepxz_reg(halphaVec,n)
    OneIterImportancepxz_reg(donnees,cova,s$alpha,halphaVec,zMAP,hyparameters,n,TT,1,log=TRUE)
  })
  return(logsum(stock) - log(S))
}




