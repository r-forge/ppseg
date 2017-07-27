

# removeClass("PPseg")

setClass(
  Class = "PPseg",
  representation=representation(
    nb_group = "numeric",
    log_likelihood = "numeric",
    parameter = "list",
    partition = "list",
    segmentation = "list",
    criteria = "list"
  ),
  prototype=prototype(
    nb_group = numeric(),
    log_likelihood = numeric(),
    parameter = list(beta = matrix(nrow=0,ncol=0),
                     lambda = numeric()),
    partition = list(estimated = matrix(nrow=0,ncol=0),
                     MAP = matrix(nrow=0,ncol=0)),                          # t_ik + MAP 
    segmentation = list(estimated = matrix(nrow=0,ncol=0),
                        MAP = numeric()),                                   # \pi_k + MAP 
    criteria = list(log_likelihoodinteg_IS = numeric(),
                    log_likelihoodinteg_BIC = numeric(),
                    log_likelihoodinteg_IS_z_IS = numeric(),
                    log_likelihoodinteg_IS_z_BIClogi = numeric(),
                    log_likelihoodcompletedinteg_IS = numeric(),
                    log_likelihoodcompletedinteg_IS_bis = numeric(),
                    log_likelihoodcompletedinteg_BIClogi = numeric())
  )  
)

setClass(
  Class = "PPsegreg",
  representation=representation(
    nb_group = "numeric",
    log_likelihood = "numeric",
    parameter = "list",
    partition = "list",
    segmentation = "list",
    criteria = "list"
  ),
  prototype=prototype(
    nb_group = numeric(),
    log_likelihood = numeric(),
    parameter = list(beta = matrix(nrow=0,ncol=0),
                     alpha = matrix(nrow=0,ncol=0)),
    partition = list(estimated = matrix(nrow=0,ncol=0),
                     MAP = matrix(nrow=0,ncol=0)),                          # t_ik + MAP 
    segmentation = list(estimated = matrix(nrow=0,ncol=0),
                        MAP = numeric()),                                   # \pi_k + MAP 
    criteria = list(log_likelihoodinteg_IS = numeric(),
                    log_likelihoodinteg_BIC = numeric(),
                    log_likelihoodinteg_IS_z_IS = numeric(),
                    log_likelihoodinteg_IS_z_BIClogi = numeric(),
                    log_likelihoodcompletedinteg_IS = numeric(),
                    log_likelihoodcompletedinteg_IS_bis = numeric(),
                    log_likelihoodcompletedinteg_BIClogi = numeric())
  )  
)

setClass(Class = "PPsegData", 
         representation = representation(data="matrix", n="numeric", TT="numeric", range="numeric"), 
         prototype = prototype(data=matrix(), n=numeric(), TT=numeric(), range=numeric())
)

PPsegData <- function(data, range=c(0,1)){
  if (is.numeric(data)) data <- matrix(data, nrow = 1)
  if (!(is.matrix(data))) stop("Input data must be an instance of matrix or numeric")
  new("PPsegData", data=data, n=nrow(data), TT=ncol(data), range=range)
}

setClass(Class = "PPsegDataReg", 
         representation = representation(data="matrix",y="array", n="numeric", TT="numeric", range="numeric"), 
         prototype = prototype(data=matrix(),y=array(), n=numeric(), TT=numeric(), range=numeric())
)

PPsegDataReg <- function(data, y, range=c(0,1)){
  if (is.numeric(data)) data <- matrix(data, nrow = 1)
  if (!(is.matrix(data))) stop("Input data must be an instance of matrix or numeric")
  if (!(is.array(data))) stop("Input y must be an instance of array")
  new("PPsegDataReg", data=data, y=y, n=nrow(data), TT=ncol(data), range=range)
}

setClass(Class = "PPsegOutput", 
         representation = representation(data="PPsegData", results="list"), 
         prototype = prototype(data=new("PPsegData"), results=list())
)

PPsegOutput <- function(data, results)
  new("PPsegOutput", data=data, results=results)

setClass(Class = "PPsegOutputReg", 
         representation = representation(data="PPsegDataReg", results="list"), 
         prototype = prototype(data=new("PPsegDataReg"), results=list())
)

PPsegOutputReg <- function(data, results)
  new("PPsegOutputReg", data=data, results=results)

