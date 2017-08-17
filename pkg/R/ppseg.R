

# **********************************************************************
#' @title Selection EM
#'
#' @description selection_EM is a function that performs nb_tests algorithms EM
#'  and select the results with the minimum likelihood. Algorithm EM is based on a Poisson process mixing model with logistic weights
#'
#' @author Julien Gheysens (email : julien.gheysens.59@gmail.com)
#'
#' @usage selection_EM(donnees, g, nb_tests = 4, nbcpu = 3)
#'
#' @param donnees is a matrix with the individuals in rows and the time step in column
#' @param g is the number of groups
#' @param nb_tests is the number of EM algorithms performed
#' @param nbcpu is the number of cores used with a function coming from "parallel"
#'
#' @return a list
#'
#' @examples
#' set.seed(123)
#' x <- c(rpois(30,5),rpois(10,10),rpois(30,15))
#' out  <- selection_EM(donnees = x, g = 3, nb_tests = 6, nbcpu = 3)
#'
#'
#'
selection_EM <- function(donnees, g, nb_tests=4, nbcpu=1){
  # ------------------------------------------------- #
  #         premiere utilisation de EM                #
  res_EM <- mclapply(1:nb_tests, function(i) EM(donnees,g), mc.silent=TRUE, mc.cores=nbcpu)
  # print(res_EM)
  res_EM <- res_EM[[which.max(sapply(1:nb_tests, function(i) res_EM[[i]]$log_vraisemblance))]]
  return(res_EM)
}




# **********************************************************************
#' @title Selection EMreg
#'
#' @description selection_EMreg is a function that performs nb_tests algorithms EM
#'  and select the results with the minimum likelihood.
#'  Algorithm EM is based on a Poisson process mixing model with logistic weights and Poisson parameters dependent on covariates
#'
#' @author Julien Gheysens (email : julien.gheysens.59@gmail.com)
#'
#' @usage selection_EMreg(donnees, cova = array(0), g, nb_tests = 4, nbcpu = 3)
#'
#' @param donnees is a matrix with the individuals in rows and the time step in column
#' @param cova is a array individuals x time step x covariates
#' @param g is the number of groups
#' @param nb_tests is the number of EM algorithms performed
#' @param nbcpu is the number of cores used with a function coming from "parallel"
#'
#' @return a list
#'
#'

selection_EMreg <- function(donnees, cova=array(0), g, nb_tests=4, nbcpu=1){
  # ------------------------------------------------- #
  #         premiere utilisation de EM                #
  res_EM <- mclapply(1:nb_tests, function(i) EMreg(donnees,cova,g), mc.silent=TRUE, mc.cores=nbcpu)
  # print(res_EM)
  res_EM <- res_EM[[which.max(sapply(1:nb_tests, function(i) res_EM[[i]]$log_vraisemblance))]]
  return(res_EM)
}



# **********************************************************************




ppsegestim <- function(donnees, test_group, S=1000, nb_tests=4, nbcpu=1)
  sapply(test_group, ppsegestimgfixed, donnees=donnees, S=S, nb_tests, nbcpu=nbcpu)


ppsegestimgfixed <- function(donnees, g, S=1000, nb_tests=4, nbcpu=1,
                             likelihoodinteg_IS = FALSE,
                            likelihoodinteg_BIC = TRUE,
                             likelihoodinteg_IS_z_IS = FALSE,
                             likelihoodinteg_IS_z_BIClogi = FALSE,
                             likelihoodcompletedinteg_IS = FALSE,
                             likelihoodcompletedinteg_IS_bis = FALSE,
                             likelihoodcompletedinteg_BIClogi = FALSE){
  TT <- length(donnees[1,])
  n <- length(donnees[,1])

  ppseg <- new("PPseg")
  ppseg@nb_group <- g
  resEM <- selection_EM(donnees, g, nb_tests, nbcpu)
  ppseg@segmentation$estimated <- resEM$poids
  ppseg@parameter$beta <- resEM$beta
  ppseg@parameter$lambda <- resEM$lambda
  ppseg@partition$estimated <- resEM$H
  ppseg@log_likelihood <- resEM$log_vraisemblance

   if(g==1){
    ppseg@segmentation$MAP <- ppseg@partition$estimated
    ppseg@partition$MAP <- ppseg@partition$estimated
    val <- log_critere_g1(donnees)
    # criteria
    if(likelihoodinteg_IS) ppseg@criteria$log_likelihoodinteg_IS <- val
    if(likelihoodinteg_BIC) ppseg@criteria$log_likelihoodinteg_BIC <- val
    if(likelihoodinteg_IS_z_IS) ppseg@criteria$log_likelihoodinteg_IS_z_IS <- val
    if(likelihoodinteg_IS_z_BIClogi) ppseg@criteria$log_likelihoodinteg_IS_z_BIClogi <- val
    if(likelihoodcompletedinteg_IS) ppseg@criteria$log_likelihoodcompletedinteg_IS <- val
    if(likelihoodcompletedinteg_IS_bis) ppseg@criteria$log_likelihoodcompletedinteg_IS_bis <- val
    if(likelihoodcompletedinteg_BIClogi) ppseg@criteria$log_likelihoodcompletedinteg_BIClogi <- val
   }else{
    zMAP <- MAP(donnees,g,lambda=resEM$lambda,beta=resEM$beta)
    ppseg@segmentation$MAP <- MAP_seg_Vec(resEM$beta,TT)
    ppseg@partition$MAP <- MAP_mat(zMAP,n,TT)
    # criteria
    if(likelihoodinteg_IS) ppseg@criteria$log_likelihoodinteg_IS <- Integratedlikelihood(S,donnees,hbeta=resEM$beta,hlambda=resEM$lambda,n,TT,g)
    if(likelihoodinteg_BIC) ppseg@criteria$log_likelihoodinteg_BIC <- BIC_select(donnees,hbeta=resEM$beta,hlambda=resEM$lambda,n,TT,g)
    if(likelihoodinteg_IS_z_IS) ppseg@criteria$log_likelihoodinteg_IS_z_IS <- Integratedlikelihood_z1(S,donnees,hbeta=resEM$beta,hlambda=resEM$lambda,hpoids=resEM$poids,zMAP=zMAP,n,TT,g)
    if(likelihoodinteg_IS_z_BIClogi) ppseg@criteria$log_likelihoodinteg_IS_z_BIClogi <- Integratedlikelihood_z2(S,donnees,hbeta=resEM$beta,hlambda=resEM$lambda,hpoids=resEM$poids,zMAP=zMAP,n,TT,g)
    if(likelihoodcompletedinteg_IS) ppseg@criteria$log_likelihoodcompletedinteg_IS <- Integratedlikelihoodcompleted(S,donnees,hbeta=resEM$beta,zMAP=zMAP,n,TT,g)
    # if(likelihoodcompletedinteg_IS_bis) ppseg@criteria$log_likelihoodcompletedinteg_IS_bis <- Integratedlikelihoodcompleted_bis(S,donnees,hbeta=resEM$beta,zMAP=zMAP,n,TT,g)
    if(likelihoodcompletedinteg_BIClogi) ppseg@criteria$log_likelihoodcompletedinteg_BIClogi <- Integratedlikelihoodcompleted_BIC(donnees,hbeta=resEM$beta,zMAP=zMAP,n,TT,g)
   }
   return(ppseg)
}





# ********************************************************
#                   Avec Covariable
# ********************************************************


ppsegregestimgfixed <- function(donnees, cova, g, S=1000, nb_tests=4, nbcpu=1,
                               likelihoodcompletedinteg_reg_IS = TRUE){
  TT <- length(donnees[1,])
  n <- length(donnees[,1])

  ppseg <- new("PPsegreg")
  ppseg@nb_group <- g
  resEM <- selection_EMreg(donnees, cova, g, nb_tests, nbcpu)
  ppseg@segmentation$estimated <- resEM$poids
  ppseg@parameter$beta <- resEM$beta
  ppseg@parameter$alpha <- resEM$alpha
  ppseg@partition$estimated <- resEM$H
  ppseg@log_likelihood <- resEM$log_vraisemblance

  if(g==1){
    ppseg@segmentation$MAP <- ppseg@partition$estimated
    ppseg@partition$MAP <- ppseg@partition$estimated
    # criteria
    val <- log_critere_g1_reg(S,donnees,cova,halpha=resEM$alpha,n,TT)
    if(likelihoodcompletedinteg_reg_IS) ppseg@criteria$log_likelihoodcompletedinteg_reg_IS <- val
 }else{
   zMAP <- MAPreg(donnees,cova,g,alpha=resEM$alpha,beta=resEM$beta)
   ppseg@segmentation$MAP <- MAP_seg_Vec_reg(as.numeric(resEM$beta[-1,]),TT,g)
   ppseg@partition$MAP <- MAP_mat(zMAP,n,TT)
   # criteria
   if(likelihoodcompletedinteg_reg_IS) ppseg@criteria$log_likelihoodcompletedinteg_reg_IS <- Integratedlikelihoodcompleted_reg(S,donnees,cova,hbeta=resEM$beta,halpha=resEM$alpha,zMAP=zMAP,n,TT,g)
 }
 return(ppseg)
}



# **********************************************************************
#' @title ppseg
#'
#' @description ppseg is a function that performs nb_tests algorithms EM
#'  and select the results with the minimum likelihood. Then, calculates the different criteria.
#'  This operation is performed for each number of groups in the vector gtests.
#'  Algorithm EM is based on a Poisson process mixing model with logistic weights.
#'
#' @author Julien Gheysens (email : julien.gheysens.59@gmail.com)
#'
#' @usage ppseg(donnees, gtests, S = 1000, nb_tests = 4, nbcpu = 1)
#'
#' @param donnees is a matrix with the individuals in rows and the time step in column
#' @param cova is a array individuals x time step x covariates
#' @param gtests is the vector of groups will be used
#' @param S is the number of iterations used in the different criteria that used a importance sampling
#' @param nb_tests is the number of EM algorithms performed
#' @param nbcpu is the number of cores used with a function coming from "parallel"
#'
#' @return a PPsegOutputReg object or a PPSegOutput object if cova is empty
#'
#'@examples
#' set.seed(123)
#' x <- c(rpois(30,5),rpois(10,10),rpois(30,15))
#' out  <- ppseg(donnees = x, g = c(1,2,3,4,5), S = 1000, nb_tests = 6, nbcpu = 3)
#'
#'
ppseg <- function(donnees,cova=array(0),gtests,S=1000,nb_tests=4,nbcpu=1,
                    likelihoodinteg_IS = FALSE,
                    likelihoodinteg_BIC = TRUE,
                    likelihoodinteg_IS_z_IS = FALSE,
                    likelihoodinteg_IS_z_BIClogi = FALSE,
                    likelihoodcompletedinteg_IS = FALSE,
                    likelihoodcompletedinteg_IS_bis = FALSE,
                    likelihoodcompletedinteg_BIClogi = FALSE,
                  # -------------------------------------------------
                    likelihoodcompletedinteg_reg_IS = TRUE){

   if(length(cova)==1){
   out <- PPsegOutput(data=PPsegData(data=donnees),results=vector("list",length(gtests)))
   names(out@results) <- paste("g = ",gtests)
   cpt <- 1
   for(k in gtests){
     out@results[[cpt]] <- ppsegestimgfixed(donnees,g=k,S,nb_tests,nbcpu,
                                            likelihoodinteg_IS,
                                            likelihoodinteg_BIC,
                                            likelihoodinteg_IS_z_IS,
                                            likelihoodinteg_IS_z_BIClogi,
                                            likelihoodcompletedinteg_IS,
                                            likelihoodcompletedinteg_IS_bis,
                                            likelihoodcompletedinteg_BIClogi)
     cpt <- cpt + 1
   }
   return(out)
 }else{
   out <- PPsegOutputReg(data=PPsegDataReg(data=donnees,y=cova),results=vector("list",length(gtests)))
   names(out@results) <- paste("g = ",gtests)
   cpt <- 1
   for(k in gtests){
     out@results[[cpt]] <- ppsegregestimgfixed(donnees,cova,g=k,S,nb_tests,nbcpu,
                                               likelihoodcompletedinteg_reg_IS)
     cpt <- cpt + 1
   }
   return(out)
 }
}





