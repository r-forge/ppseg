
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
    for (j in 1:length(lev)) lines(lev[[j]]$x, rep(x@results[[y]]@parameter$lambda[lev[[j]]$y]/x@data@TT,2), col=lev[[j]]$y+1, lwd=3)
  }
)

# setMethod(
#   f="summary",
#   signature = c("PPsegOutput"),
#   definition = function(x){
#     cat("*** Data used ***\n")
#     cat("n =",x@data@n," & T =",x@data@TT,"\n")
#     cat("*** Result ***\n")
#     gestim <- sapply(x@results, function(j) nrow(j@parameter$beta))
#     crit <- sapply(x@results, function(j) j@criteria$log_likelihoodinteg_BIC)
#     cat("With the BIC criteria \n")
#     cat(names(which.max(crit)),"with criteria = ",max(crit),"\n")
#     cat("**************\n")        
#   }
# )
# 
# 
# 







