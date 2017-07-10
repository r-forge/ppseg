

selection_EM <- function(donnees, g, nb_tests=4, nbcpu=3){
  # ------------------------------------------------- #
  #         premiere utilisation de EM                #
  res_EM <- mclapply(1:nb_tests, function(i) EM(donnees,g), mc.silent=TRUE, mc.cores=nbcpu)
  # print(res_EM)
  res_EM <- res_EM[[which.max(sapply(1:nb_tests, function(i) res_EM[[i]]$log_vraisemblance))]]
  return(res_EM)
}


ppsegestim <- function(donnees, test_group, S=1000, nb_tests=4, nbcpu=3)
  sapply(test_group, ppsegestim, donnees=donnees, S=S, nb_tests, nbcpu=nbcpu)



ppsegestimgfixed <- function(donnees, g, S=1000, nb_tests=4, nbcpu=3){
  TT <- length(donnees[1,])
  n <- length(donnees[,1])
  
  print("A")
  ppseg <- new("PPseg")
  ppseg@nb_group <- g
  resEM <- selection_EM(donnees, g, nb_tests, nbcpu)
  ppseg@segmentation$estimated <- resEM$poids
  ppseg@parameter$beta <- resEM$beta
  ppseg@parameter$lambda <- resEM$lambda
  ppseg@partition$estimated <- resEM$H
  ppseg@log_likelihood <- resEM$log_vraisemblance
  
  print("B")
  if(g==1){
    ppseg@segmentation$MAP <- ppseg@partition$estimated
    ppseg@partition$MAP <- ppseg@partition$estimated
    # criteria
    ppseg@criteria$log_likelihoodinteg_IS <- log_critere_g1(donnees) 
    ppseg@criteria$log_likelihoodinteg_BIC <- ppseg@criteria$log_likelihoodinteg_IS
    ppseg@criteria$log_likelihoodinteg_IS_z_IS <- ppseg@criteria$log_likelihoodinteg_IS
    ppseg@criteria$log_likelihoodinteg_IS_z_BIClogi <- ppseg@criteria$log_likelihoodinteg_IS
    ppseg@criteria$log_likelihoodcompletedinteg_IS <- ppseg@criteria$log_likelihoodinteg_IS
    ppseg@criteria$log_likelihoodcompletedinteg_IS_bis <- ppseg@criteria$log_likelihoodinteg_IS
    ppseg@criteria$log_likelihoodcompletedinteg_BIClogi <- ppseg@criteria$log_likelihoodinteg_IS
  }else{
    zMAP <- MAP(donnees,g,lambda=resEM$lambda,beta=resEM$beta)
    ppseg@segmentation$MAP <- MAP_seg_Vec(resEM$beta,TT)
    ppseg@partition$MAP <- MAP_mat(zMAP)
    # criteria
    # ppseg@criteria$log_likelihoodinteg_IS <- Integratedlikelihood(S,donnees,hbeta=resEM$beta,hlambda=resEM$lambda,n,TT,g)
    # ppseg@criteria$log_likelihoodinteg_BIC <- BIC_select(donnees,hbeta=resEM$beta,hlambda=resEM$lambda,n,TT,g)
    # ppseg@criteria$log_likelihoodinteg_IS_z_IS <- Integratedlikelihood_z1(S,donnees,hbeta=resEM$beta,hlambda=resEM$lambda,hpoids=resEM$poids,zMAP=zMAP,n,TT,g)
    # ppseg@criteria$log_likelihoodinteg_IS_z_BIClogi <- Integratedlikelihood_z2(S,donnees,hbeta=resEM$beta,hlambda=resEM$lambda,hpoids=resEM$poids,zMAP=zMAP,n,TT,g)
    ppseg@criteria$log_likelihoodcompletedinteg_IS <- Integratedlikelihoodcompleted(S=1000,donnees,hbeta=resEM$beta,zMAP=zMAP,n,TT,g)
    ppseg@criteria$log_likelihoodcompletedinteg_IS_bis <- Integratedlikelihoodcompleted_bis(S=1000,donnees,hbeta=resEM$beta,zMAP=zMAP,n,TT,g)
    # ppseg@criteria$log_likelihoodcompletedinteg_BIClogi <- Integratedlikelihoodcompleted_BIC(donnees,hbeta=resEM$beta,zMAP=zMAP,n,TT,g)
  }
  return(ppseg)
}

ppseg <- function(donnees,gtests,S=1000,nb_tests=4,nbcpu=3){
  out <- PPsegOutput(data=PPsegData(data=donnees),results=vector("list",length(gtests)))
  names(out@results) <- paste("g =",gtests)
  
  for(k in gtests){
    print(k)
    out@results[[k]] <- ppsegestimgfixed(donnees,g=k,S,nb_tests,nbcpu)
  }
  return(out)
}

PPsegData <- function(data, range=c(0,1)){
  if (is.numeric(data)) data <- matrix(data, nrow = 1)
  if (!(is.matrix(data))) stop("Input data must be an instance of matrix or numeric")
  new("PPsegData", data=data, n=nrow(data), TT=ncol(data), range=range)
}



