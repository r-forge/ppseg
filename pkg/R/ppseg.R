

selection_EM <- function(donnees,g,nb_tests=4){
  # ------------------------------------------------- #
  #         premiere utilisation de EM                #
  res_EM <- mclapply(1:nb_tests, function(i) EM(donnees,g), mc.silent=TRUE, mc.cores=2)
  res_EM <- res_EM[[which.max(sapply(1:nb_tests, function(i) res_EM[[i]]$log_vraisemblance))]]
  #   res_EM_selection <- EM(donnees,g)
  #   # ------------------------------------------------- #
  #   #         recherche d'un meilleur resultat          #
  #   for (test in 1:nb_tests-1){
  #     res_EM <- EM(donnees,g)
  #     if(res_EM$croissance_algo){
  #       if(res_EM$log_vraisemblance>res_EM_selection$log_vraisemblance){
  #         res_EM_selection <- res_EM
  #       }
  #     }else{ # croissance_algo == FALSE
  #       print("L'une des executions de l algorithme EM ne fait pas croitre le vraisemblance")
  #       return(0)
  #     }
  #   }
  # ------------------------------------------------- #
  #          sortie du meilleur resultat              #
  return(res_EM)
}



