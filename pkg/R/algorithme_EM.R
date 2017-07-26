# ************************************************************************************** #
#                                                                                        # 
#                           EM  avec la fonction optim()                                 #
#                                                                                        #
# ************************************************************************************** #

# 
# setwd("~/Documents/programmation/EM_15_juin")
# source("outils_poids.R")
# source("outils_etape_E_de_EM.R")             # Fonction pour (E)
# source("outils_etape_M_de_EM.R")             # Fonction pour (M)
# 


# ************************************************* #
#                    EM                             #
# ************************************************* #

EM <- function(donnees,g,nb_iteration_max=1000,eps=2*10^-8){

  # *************************************************************** #
  #                       Initialisation                            #
  
  # Necessite de "donnees" et "g" comme parametre et le fichier "outils_poids.R"
  
  # --- Si  g = 1 => maximun de vraisamblance --- #
  if(g==1){
    TT <- length(donnees[1,])          # T est le nombre d intervalles dans [0;1] 
    n <- length(donnees[,1])          # n est le nombre d individus 
    poids <- matrix(1,nrow=TT,ncol=1)
    beta <- matrix(0,nrow=1,ncol=2)
    lambda <- sum(donnees)/n
    croissance <- TRUE
    m <- 1
    Log_vraisemblance <- sum(sapply(1:n,function(i) sum(sapply(1:TT,function(k) log(dpois(donnees[i,k],lambda/TT))))))
    return(list(poids=poids,beta=beta,lambda=lambda,croissance_algo=croissance,nb_iteration=m,log_vraisemblance=Log_vraisemblance))
  }else{
  
    TT <- length(donnees[1,])          # T est le nombre d intervalles dans [0;1] 
    n <- length(donnees[,1])          # n est le nombre d individus 
    
    # ---------------------------------------------------------------
    
    # S est le tableau utilise lors de l estimation de beta par la fonction "glm" 
    S <- data.frame(compo=as.factor(rep(1:g,TT)),time=rep((1:TT)/TT,each=g))
    
    # ---------------------------------------------------------------
    
    # Initialisation des parametres \theta_0
    lambda <- 0.01 + sample((TT*min(donnees)):(TT*max(donnees)),g,replace=TRUE)
    beta <- rbind(c(0,0),matrix(runif(2*(g-1),-10,10),nrow=g-1,ncol=2))
    betaVec <- as.numeric(beta[-1,])
    poids <- ensemble_poids_2(betaVec,TT,g)     # fonction presente dans fichier "outils_poids.R"
    
    # ---------------------------------------------------------------
    
    # parametres pour la boucle
    Log_vraisemblance <- -Inf
    Log_vraisemblance_precedent <- -Inf
    # La log_vraisemblance doit croitre a chaque etape
    
    m <- 0                            # nombre d iteration 
    
    # parametres pour le critere d arret de la boucle "while"
    difference <- Inf
    croissance <- TRUE
    
    # *************************************************************** #
    #                         Boucle EM                               #
    while((m<nb_iteration_max)&(abs(difference)>eps)){ # &(croissance)
      # ------------------------------------------------- #
      m <- m + 1
      Log_vraisemblance_precedent <- Log_vraisemblance
      
      # ------------------------------------------------- #
      #                    Etape (E)                      #
      H <- matrice_H(donnees,lambda,poids)
      #    somme sur les i de H[i,k,j]
      w <- rowSums(sapply(1:(dim(H)[1]), function(i) as.numeric(t(H[i,,]))))
  
      # ------------------------------------------------- #
      #                    Etape (M)                      #
      lambda <- estimation_max_lambda(n,g,TT,donnees,H,w)
  
      betaVec <- as.numeric(beta[-1,]) 
      tmp <- optim(betaVec, fn=log_vraisemblance_poids_2, gr=grad_log_vraisemblance_poids_2, w=w,TT=TT,g=g, control=list(fnscale=-1),method="L-BFGS-B" ) 
      beta <- rbind(c(0,0),matrix(tmp$par, ncol=2))
      betaVec <- as.numeric(beta[-1,]) 
      # print(beta)
      
      # ------------------------------------------------- #
      #            Calcul de la log vraisemblance         #
      
      poids <- ensemble_poids_2(betaVec,TT,g)
      Log_vraisemblance <- sum(sapply(1:n,function(i) sum(sapply(1:TT,function(k) log(sum(sapply(1:g,function(j) poids[k,j]*dpois(donnees[i,k],lambda[j]/TT))))))))
      
      # ---------------------------------------------------------------------- #
      #          Verification de la croissante de la log vraisemblance         #
      
      difference <- Log_vraisemblance - Log_vraisemblance_precedent
      if(difference<(-eps)){
        croissance <- FALSE
        print("l'algorithme EM n'a pas fait croitre la vraisemblance a l'etape : ")
        print(m)
      }
  
      # ---------------------------------------------------------------------- #
    }
    # *************************************************************** #
    H <- matrice_H(donnees,lambda,poids)
    return(list(poids=poids,beta=beta,lambda=lambda,H=H,croissance_algo=croissance,nb_iteration=m,log_vraisemblance=Log_vraisemblance))
  }
}


# ************************************************* #
#                    EM reg                         #
# ************************************************* #

EMreg <- function(donnees,cova=array(1,dim=c(n,length(donnees[1,]),1)),g,nb_iteration_max=1000,eps=2*10^-8){
  
  # *************************************************************** #
  #                       Initialisation                            #
  
  # Necessite de "donnees" et "g" comme parametre et le fichier "outils_poids.R"
  
  # --- Si  g = 1 => maximun de vraisamblance --- #
  if(g==1){
    TT <- length(donnees[1,])          # T est le nombre d intervalles dans [0;1] 
    n <- length(donnees[,1])          # n est le nombre d individus 
    p <- dim(cova)[3]
    
    poids <- matrix(1,nrow=TT,ncol=1)
    beta <- matrix(0,nrow=1,ncol=2)
    lambda <- sum(donnees)/n
    croissance <- TRUE
    m <- 1
    Log_vraisemblance <- sum(sapply(1:n,function(i) sum(sapply(1:TT,function(k) log(dpois(donnees[i,k],lambda/TT))))))
    return(list(poids=poids,beta=beta,lambda=lambda,croissance_algo=croissance,nb_iteration=m,log_vraisemblance=Log_vraisemblance))
  }else{
    
    TT <- length(donnees[1,])          # T est le nombre d intervalles dans [0;1] 
    n <- length(donnees[,1])          # n est le nombre d individus 
    p <- dim(cova)[3]
    
    # Structure pour glm des composantes 
    slist <- vector("list",p)
    for(l in 1:p){
      slist[[l]] <- as.numeric(cova)[((l-1)*(n*TT)+1):(l*n*TT)]
    }
    Scova <- data.frame(slist)
    for(l in 1:p){
      names(Scova)[l] <- paste("cova",l,sep="_")
    }
    S <- data.frame(X = as.numeric(donnees))
    S <- cbind(S*TT,Scova)
    #cbind(as.numeric(donnees),as.numeric(cova[,,1]))
    # ---------------------------------------------------------------
    
    # Initialisation des parametres \theta_0
    # alpha <- log(matrix(rep(0.01 + sample((TT*min(donnees)):(TT*max(donnees)),g,replace=TRUE),p),g,p))
    
    alpha <- matrix(0,nrow=g,ncol=p+1)
    alpha[,1] <- log(rep(0.01 + sample((TT*min(donnees)):(TT*max(donnees)),g,replace=TRUE)),g)
    # lambda <- 0.01 + sample((TT*min(donnees)):(TT*max(donnees)),g,replace=TRUE)
    beta <- rbind(c(0,0),matrix(runif(2*(g-1),-10,10),nrow=g-1,ncol=2))
    betaVec <- as.numeric(beta[-1,])
    poids <- ensemble_poids_2(betaVec,TT,g)     # fonction presente dans fichier "outils_poids.R"
    
    # ---------------------------------------------------------------
    
    # parametres pour la boucle
    Log_vraisemblance <- -Inf
    Log_vraisemblance_precedent <- -Inf
    # La log_vraisemblance doit croitre a chaque etape
    
    m <- 0                            # nombre d iteration 
    
    # parametres pour le critere d arret de la boucle "while"
    difference <- Inf
    croissance <- TRUE
    
    # *************************************************************** #
    #                         Boucle EM                               #
    while((m<nb_iteration_max)&(abs(difference)>eps)){ # &(croissance)
      # ------------------------------------------------- #
      m <- m + 1
      Log_vraisemblance_precedent <- Log_vraisemblance
      
      # ------------------------------------------------- #
      #                    Etape (E)                      #
      H <- matrice_Hreg(donnees,cova,alpha,poids)
      #    somme sur les i de H[i,k,j]
      w <- rowSums(sapply(1:(dim(H)[1]), function(i) as.numeric(t(H[i,,]))))
      
      # ------------------------------------------------- #
      #                    Etape (M)                      #
      
      for(k in 1:g){
        alpha[k,] <- glm(X~.,data=S,weights=as.numeric(H[,,k]),family=poisson())$coefficients
      }
      
#       alphaVec <- c(alpha)
#       alpha <- matrix(data = optim(alphaVec,
#                                    fn=log_vraisemblance_composante,
#                                    gr=grad_log_vraisemblance_composante,
#                                    donnees=donnees,cova=cova,H=H,n=n,TT=TT,gg=g,pp=p,control=list(fnscale=-1))$par,nrow = g,ncol = p)
#       
      betaVec <- as.numeric(beta[-1,]) 
      tmp <- optim(betaVec, fn=log_vraisemblance_poids_2, gr=grad_log_vraisemblance_poids_2, w=w,TT=TT,g=g, control=list(fnscale=-1),method="L-BFGS-B" ) 
      beta <- rbind(c(0,0),matrix(tmp$par, ncol=2))
      betaVec <- as.numeric(beta[-1,]) 
      
      # ------------------------------------------------- #
      #            Calcul de la log vraisemblance         #
      
      poids <- ensemble_poids_2(betaVec,TT,g)
      Log_vraisemblance <- sum(sapply(1:n,function(i) sum(sapply(1:TT,function(j) log(sum(sapply(1:g,function(k) poids[j,k]*dpois(donnees[i,j],exp(sum(alpha[k,]*c(1,cova[i,j,])))/TT))))))))
      
      # ---------------------------------------------------------------------- #
      #          Verification de la croissante de la log vraisemblance         #
      
      difference <- Log_vraisemblance - Log_vraisemblance_precedent
      if(difference<(-eps)){
        croissance <- FALSE
        print("l'algorithme EM n'a pas fait croitre la vraisemblance a l'etape : ")
        print(m)
      }
      
      # ---------------------------------------------------------------------- #
    }
    # *************************************************************** #
    H <- matrice_Hreg(donnees,cova,alpha,poids)
    return(list(poids=poids,alpha=alpha,lambda=lambda,H=H,croissance_algo=croissance,nb_iteration=m,log_vraisemblance=Log_vraisemblance))
  }
}

