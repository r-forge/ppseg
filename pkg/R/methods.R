# ************************************************
#'@title Plot with PPseg 
#'
#'@description 
#'The function plot make a graph data as a function of a discretization of the interval "range" in the PPsegOutput object.
#'Moreover, there are on the graph different lines representing the value of the Poisson parameters of each group with their starting time of a section.
#'
#'@param 
#' x : a PPsegOutput object or a PPsegOutputReg object
#' y : the number of groups displayed ( the values must be in PPsegOutput )
#'

setMethod(
  f="plot",
  signature = c("PPsegOutput", "numeric"),
  definition = function(x, y){
    plot(seq(x@data@range[1], x@data@range[2], length.out = x@data@TT),
         x@data@data,
         pch=20,
         ylab="Time",
         xlab="Count")
    gestim <- sapply(x@results, function(j) nrow(j@parameter$beta))
    if (!(y %in% gestim)) stop(paste0("Model with", y, "components is not estimated"))
    y <- which(y == gestim)
    segmentation <- cbind(seq(x@data@range[1], x@data@range[2], length.out = x@data@TT),
                          apply(x@results[[y]]@segmentation$estimated, 1, which.max))
    if (any(segmentation[-1,2] != segmentation[-x@data@TT,2])){
      change <- which(segmentation[-1,2] != segmentation[-x@data@TT,2])
      lev <- list(list(x=c(x@data@range[1], segmentation[change[1]-1,1]), y=segmentation[change[1]-1,2]))
      if (length(change)>1)    for (j in 1:(length(change)-1)) lev <- c(lev, list(list(x=c(segmentation[change[j],1], segmentation[change[j+1]-1,1]), y=segmentation[change[j+1],2])))
      lev <- c(lev, list(list(x=c(segmentation[change[length(change)],1], x@data@range[2]), y=segmentation[x@data@TT,2])))
    }else{
      lev <- list(x=x@data@range, y=segmentation[1,2])
    }
    for (j in 1:length(lev)){
      lines(lev[[j]]$x, rep(x@results[[y]]@parameter$lambda[lev[[j]]$y]/x@data@TT,2), col=lev[[j]]$y+1, lwd=3)
      points(lev[[j]]$x[1],
             x@results[[y]]@parameter$lambda[lev[[j]]$y]/x@data@TT,
             pch = 1, cex = 2)
      text(lev[[j]]$x[1],
           x@results[[y]]@parameter$lambda[lev[[j]]$y]/x@data@TT + min(5,mean(x@data@data)/4),
           paste(lev[[j]]$x[1]))
    } 
  }
)

setMethod(
  f="plot",
  signature = c("PPsegOutputReg", "numeric"),
  definition = function(x, y){
    plot(seq(x@data@range[1], x@data@range[2], length.out = x@data@TT),
         x@data@data,
         pch=20,
         ylab="Time",
         xlab="Count")
    gestim <- sapply(x@results, function(j) nrow(j@parameter$beta))
    if (!(y %in% gestim)) stop(paste0("Model with", y, "components is not estimated"))
    y <- which(y == gestim)
    segmentation <- cbind(seq(x@data@range[1], x@data@range[2], length.out = x@data@TT),
                          apply(x@results[[y]]@segmentation$estimated, 1, which.max))
    if (any(segmentation[-1,2] != segmentation[-x@data@TT,2])){
      change <- which(segmentation[-1,2] != segmentation[-x@data@TT,2])
      lev <- list(list(x=c(x@data@range[1], segmentation[change[1]-1,1]), y=segmentation[change[1]-1,2]))
      if (length(change)>1)    for (j in 1:(length(change)-1)) lev <- c(lev, list(list(x=c(segmentation[change[j],1], segmentation[change[j+1]-1,1]), y=segmentation[change[j+1],2])))
      lev <- c(lev, list(list(x=c(segmentation[change[length(change)],1], x@data@range[2]), y=segmentation[x@data@TT,2])))
      
      # points(lev[[1]]$x[1],
      #        exp(sum(x@results[[y]]@parameter$alpha[lev[[1]]$y,]*c(1,x@data@y[which(segmentation[,1]==lev[[1]]$x[1])])))/x@data@TT,
      #        pch = 1, cex = 2)
      # text(lev[[1]]$x[1],
      #      exp(sum(x@results[[y]]@parameter$alpha[lev[[1]]$y,]*c(1,x@data@y[which(segmentation[,1]==lev[[1]]$x[1])])))/x@data@TT + 5,
      #      paste(lev[[1]]$x[1]))
      for (j in 1:length(lev)){
        lines(lev[[j]]$x, 
              c(exp(sum(x@results[[y]]@parameter$alpha[lev[[j]]$y,]*c(1,x@data@y[which(segmentation[,1]==lev[[j]]$x[1])])))/x@data@TT,
                exp(sum(x@results[[y]]@parameter$alpha[lev[[j]]$y,]*c(1,x@data@y[which(segmentation[,1]==lev[[j]]$x[2])])))/x@data@TT),
              col=lev[[j]]$y+1, lwd=3)
        points(lev[[j]]$x[1],
               exp(sum(x@results[[y]]@parameter$alpha[lev[[j]]$y,]*c(1,x@data@y[which(segmentation[,1]==lev[[j]]$x[1])])))/x@data@TT,
               pch = 1, cex = 2)
        text(lev[[j]]$x[1],
             exp(sum(x@results[[y]]@parameter$alpha[lev[[j]]$y,]*c(1,x@data@y[which(segmentation[,1]==lev[[j]]$x[1])])))/x@data@TT + min(5,mean(x@data@data)/4),
             paste(lev[[j]]$x[1]))
      }
    }else{
      lev <- list(x=x@data@range, y=segmentation[1,2])
      lines(lev$x, 
            c(exp(sum(x@results[[y]]@parameter$alpha[lev$y,]*c(1,x@data@y[which(segmentation[,1]==lev$x[1])])))/x@data@TT,
              exp(sum(x@results[[y]]@parameter$alpha[lev$y,]*c(1,x@data@y[which(segmentation[,1]==lev$x[2])])))/x@data@TT),
            col=lev$y+1, lwd=3)
      
    }
  }
)


# ************************************************
#'@title Summary with PPseg 
#'
#'@description 
#'The function summary show the number of groups selected by the different criteria with the maximum values used in a PPsegOutput object or PPsegOutput object.
#'
#'@param 
#' object : a PPsegOutput object or a PPsegOutputReg object
#' 
#'

setMethod(
  f="summary",
  signature = c("PPsegOutputReg"),
  definition = function(object){
    x <- object
    cat("*** Data used ***\n")
    cat("n =",x@data@n," & T =",x@data@TT,"\n")
    cat("*** Result ***\n")
    gestim <- sapply(x@results, function(j) nrow(j@parameter$beta))
    crit1 <- sapply(x@results, function(j) j@criteria$log_likelihoodcompletedinteg_reg_IS)
    if(length(crit1[[1]])!=0){
      cat("With the likelihoodcompletedinteg_reg_IS criteria \n")
      cat(names(which.max(crit1)),"with criteria = ",max(crit1),"\n")
      cat("**************\n")
    }
  }
)
 
setMethod(
  f="summary",
  signature = c("PPsegOutput"),
  definition = function(object){
    x <- object
    cat("*** Data used ***\n")
    cat("n =",x@data@n," & T =",x@data@TT,"\n")
    cat("*** Result ***\n")
    gestim <- sapply(x@results, function(j) nrow(j@parameter$beta))
    crit1 <- sapply(x@results, function(j) j@criteria$log_likelihoodinteg_IS)
    crit2 <- sapply(x@results, function(j) j@criteria$log_likelihoodinteg_BIC)
    crit3 <- sapply(x@results, function(j) j@criteria$log_likelihoodinteg_IS_z_IS)
    crit4 <- sapply(x@results, function(j) j@criteria$log_likelihoodinteg_IS_z_BIClogi)
    crit5 <- sapply(x@results, function(j) j@criteria$log_likelihoodcompletedinteg_IS)
    crit6 <- sapply(x@results, function(j) j@criteria$log_likelihoodcompletedinteg_IS_bis)
    crit7 <- sapply(x@results, function(j) j@criteria$log_likelihoodcompletedinteg_BIClogi)
    if(length(crit1[[1]])!=0){
      cat("With the likelihoodinteg_IS criteria \n")
      cat(names(which.max(crit1)),"with criteria = ",max(crit1),"\n")
      cat("**************\n")
    }
    if(length(crit2[[1]])!=0){
      cat("With the likelihoodinteg_BIC criteria \n")
      cat(names(which.max(crit2)),"with criteria = ",max(crit2),"\n")
      cat("**************\n")
    }
    if(length(crit3[[1]])!=0){
      cat("With the likelihoodinteg_IS_z_IS criteria \n")
      cat(names(which.max(crit3)),"with criteria = ",max(crit3),"\n")
      cat("**************\n")
    }
    if(length(crit4[[1]])!=0){
      cat("With the likelihoodinteg_IS_z_BIClogi criteria \n")
      cat(names(which.max(crit4)),"with criteria = ",max(crit4),"\n")
      cat("**************\n")
    }
    if(length(crit5[[1]])!=0){
      cat("With the likelihoodcompletedinteg_IS criteria \n")
      cat(names(which.max(crit5)),"with criteria = ",max(crit5),"\n")
      cat("**************\n")
    }
    if(length(crit6[[1]])!=0){
      cat("With the likelihoodcompletedinteg_IS_bis criteria \n")
      cat(names(which.max(crit6)),"with criteria = ",max(crit6),"\n")
      cat("**************\n")
    }
    if(length(crit7[[1]])!=0){
      cat("With the likelihoodcompletedinteg_BIClogi criteria \n")
      cat(names(which.max(crit7)),"with criteria = ",max(crit7),"\n")
      cat("**************\n")
    }    
  }
)






