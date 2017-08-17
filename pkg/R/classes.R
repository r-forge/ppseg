

# ************************************************
#'@title PPseg Class
#'
#'@description 
#'PPseg Class groups all the output elements of EM and the selection criteria values for a certain value of nb_group.
#'
#'@return 
#'    Output elements of EM 
#' -> nb_group : number of groups
#' -> log_likelihood
#' -> parameter : a list with beta and lambda parameters
#' -> partition : a list with partition estimated by EM ( last step (E) ) and MAP the maximum posteriori of the partition estimated
#' -> segmentation : a list with segmentation estimated by the logistic part with beta estiamted by EM and MAP the maximum posteriori of the segmentation estimated
#' 
#'    Output elements of criteria : list of different criterias perform 
#'

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

# ************************************************
#'@title PPsegreg Class
#'
#'@description 
#'PPsegreg Class groups all the output elements of EM and the selection criteria values for a certain value of nb_group.
#'
#'@return 
#'    Output elements of EM 
#' -> nb_group : number of groups
#' -> log_likelihood
#' -> parameter : a list with beta and alpha parameters
#' -> partition : a list with partition estimated by EM ( last step (E) ) and MAP the maximum posteriori of the partition estimated
#' -> segmentation : a list with segmentation estimated by the logistic part with beta estiamted by EM and MAP the maximum posteriori of the segmentation estimated
#' 
#'    Output elements of criteria : list of different criterias perform 
#'

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
    criteria = list(log_likelihoodcompletedinteg_reg_IS = numeric()
                    )  
  )
)
# ************************************************
#'@title PPsegData Class
#'
#'@description 
#'PPsegData Class groups all input elements of EM.
#'
#'@usage 
#'PPsegData(data,range=(0,1))
#'
#'@return 
#'    Input elements of EM 
#' -> Data
#' -> n : Number of individuals
#' -> TT : Number of steps of time
#' -> range : data's range 
#'
#'
setClass(Class = "PPsegData", 
         representation = representation(data="matrix", n="numeric", TT="numeric", range="numeric"), 
         prototype = prototype(data=matrix(), n=numeric(), TT=numeric(), range=numeric())
)

PPsegData <- function(data, range=c(0,1)){
  if (is.numeric(data)) data <- matrix(data, nrow = 1)
  if (!(is.matrix(data))) stop("Input data must be an instance of matrix or numeric")
  new("PPsegData", data=data, n=nrow(data), TT=ncol(data), range=range)
}

# ************************************************
#'@title PPsegDataReg Class
#'
#'@description 
#'PPsegDataReg Class groups all input elements of EM.
#'
#'@usage 
#'PPsegDataReg(data = data,y = covariates,range = (0,1))
#'
#'@return 
#'    Input elements of EM 
#' -> Data
#' -> y : covariates
#' -> n : Number of individuals
#' -> TT : Number of steps of time
#' -> range : data's range 
#'
#'
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


# ************************************************
#'@title PPsegOutput Class
#'
#'@description 
#'PPsegOutput Class groups all input and output elements of EM.
#'
#'@return 
#'    Input elements of EM : a PPsegData object
#'    Output elementsof EM : a list with PPseg objects
#'
#'
setClass(Class = "PPsegOutput", 
         representation = representation(data="PPsegData", results="list"), 
         prototype = prototype(data=new("PPsegData"), results=list())
)

PPsegOutput <- function(data, results)
  new("PPsegOutput", data=data, results=results)


# ************************************************
#'@title PPsegOutputReg Class
#'
#'@description 
#'PPsegOutput Class groups all input and output elements of EM.
#'
#'@return 
#'    Input elements of EM : a PPsegDataReg object
#'    Output elementsof EM : a list with PPsegreg objects
#'
#'
setClass(Class = "PPsegOutputReg", 
         representation = representation(data="PPsegDataReg", results="list"), 
         prototype = prototype(data=new("PPsegDataReg"), results=list())
)

PPsegOutputReg <- function(data, results)
  new("PPsegOutputReg", data=data, results=results)

