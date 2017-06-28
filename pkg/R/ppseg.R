

selection_EM <- function(donnees, g, nb_tests=4, nbcpu=3){
  # ------------------------------------------------- #
  #         premiere utilisation de EM                #
  res_EM <- mclapply(1:nb_tests, function(i) EM(donnees,g), mc.silent=TRUE, mc.cores=nbcpu)
  # print(res_EM)
  res_EM <- res_EM[[which.max(sapply(1:nb_tests, function(i) res_EM[[i]]$log_vraisemblance))]]
  return(res_EM)
}

# **************************************************************************** # 
ppsegestimgfixed <- function(donnees, g, S=1000, nb_tests=4, nbcpu=3){
  resEM <- selection_EM(donnees, g, nb_tests)
  resEM$zMAP <- MAP(donnees, g, resEM$lambda, resEM$beta)
  resEM$integratedlikelihood <- log_critere_BIC(S, donnees, resEM$beta, resEM$lambda, g)
  resEM$integratedcompletedlikelihood <- log_critere_ICL(S, 
                                                         resEM$beta,
                                                         donnees,
                                                         resEM$zMAP,
                                                         nrow(donnees), 
                                                         ncol(donnees),
                                                         g)
  resEM
  # Ici creer le S4
}

ppsegestim <- function(donnees, test_group, S=1000, nb_tests=4, nbcpu=3)
  sapply(test_group, ppsegestim, donnees=donnees, S=S, nb_tests, nbcpu=nbcpu)

