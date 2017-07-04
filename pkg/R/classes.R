setClass(Class = "PPsegData", 
         representation = representation(data="matrix", n="numeric", TT="numeric", range="numeric"), 
         prototype = prototype(data=matrix(), n=numeric(), TT=numeric(), range=numeric())
)

PPsegData <- function(data, range=c(0,1)){
  if (is.numeric(data)) data <- matrix(data, nrow = 1)
  if (!(is.matrix(data))) stop("Input data must be an instance of matrix or numeric")
  new("PPsegData", data=data, n=nrow(data), TT=ncol(data), range=range)
}

setClass(Class = "PPsegOutput", 
         representation = representation(data="PPsegData", results="list"), 
         prototype = prototype(data=new("PPsegData"), results=list())
)

PPsegOutput <- function(data, results)
  new("PPsegOutput", data=data, results=results)

