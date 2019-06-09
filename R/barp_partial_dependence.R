#' barp_partial_dependence
#' 
#' This function calculates the partial dependence for selected covariates based on a \code{BARP} object. The user can specify up to three covariates, resulting in a three-way interaction. 
#' @param barp A \code{\link{barp}} object.
#' @param vars The variable names or column indices of interest. #' @param vars The variable names or column indices of interest. If a character vector, values at which to evaluate partial dependences based on quantiles c(0.05,seq(.1,.9,by = .1),.95).
#' @param prop_data The percentage of the original data to use for estimation. Larger values are more time consuming.
#' @param credible_interval The credible interval for the estimates.
#' @param setSeed Seed to control random number generation.
#' @param ... Additional arguments to be passed to bartMachine.
#' @return Returns an object of class "barpd", containing a list of the following components:
#' \item{summary}{A \code{data.frame} listing the mean and upper and lower bounds of the predicted outcome at each level of the covariate(s).}
#' \item{raw}{A \code{data.frame} where each row is a posterior draw and each column is a level for the covariate(s).}
#' @keywords MRP BARP Mister P multilevel regression poststratification
#' @seealso \code{\link{barp}} which generates the prerequisite \code{barp} object.
#' @examples 
#' data("gaymar")
#' barp.obj <- barp(y = "supp_gaymar",
#'              x = c("pvote","religcon","age","educ","gXr","stateid","region"),
#'              dat = svy,
#'              census = census06,
#'              geo.unit = "stateid",
#'              proportion = "n")
#' bpd <- barp_partial_dependence(barp = barp.obj,
#'                                vars = list(region = 1:4,educ = 1:4),
#'                                prop_data = .2,
#'                                credible_interval = c(0.025,0.975),
#'                                setSeed = 1021)
#' @rdname barp_partial_dependence
#' @export

barp_partial_dependence <- function (barp.obj,vars = NULL, prop_data = .1,
                                     credible_interval = c(0.025,.975),setSeed = NULL)
  
{
  if(is.null(vars)) {
    stop("You must supply at least one variable for analysis.")
  }
  if(is.null(barp.obj$trees)) {
    stop(paste0("This BARP object was built using ",paste(barp.obj$methods,collapse = ", "),
                ". To evaluate partial dependencies, please build with method=bartMachine."))
  } else {
    bart_machine <- barp.obj$trees
    bartMachine:::check_serialization(bart_machine)
    
    factors <- names(barp.obj$factors)
    
    if(class(vars) == "character") {
      tmpvars <- list()
      for(v in vars) {
        tmpvars[[v]] <- unique(bart_machine$X[[v]])
      }
      vars <- tmpvars
    }
    
    for(v in names(vars)) {
      if (!(v %in% colnames(bart_machine$X)) & !(v %in% factors)) {
        stop(paste0("var ",v," must be the name of one of the training features (see \"<bart_model>$training_data_features\")"))
      }
    }
    
    n_samp = round(bart_machine$n * prop_data)
    sum.pred <- sum.prob <- as.data.frame(expand.grid(vars))
    sum.pred <- as.data.frame(cbind(sum.pred,matrix(NA,nrow = nrow(sum.pred),ncol = 3)))
    colnames(sum.pred) <- c(names(vars),"pred","lb","ub")
    raw.pred <- raw.prob <- as.data.frame(matrix(NA,nrow = bart_machine$num_iterations_after_burn_in,
                                                 ncol = nrow(sum.pred)))
    cat(paste0("\nCalculating partial dependencies for ",
               paste(names(vars),collapse = " by ")," using ",prop_data*100,"% of the data.\n"))
    pb <- txtProgressBar(min = 0, max = nrow(sum.pred), style = 3)
    set.seed(setSeed)
    for(n in 1:nrow(sum.pred)) {
      indices = sample(1:bart_machine$n, n_samp)
      test_data = bart_machine$X[indices, ]
      for(col in 1:length(vars)) {
        test_data[, names(vars)[col]] = rep(sum.pred[n,col], n_samp)
      }
      
      tmp <- extract_predictions(bart_machine,test_data,credible_interval)
      sum.pred[n,c("pred","lb","ub")] <- tmp$summary.pred
      
      raw.pred[,n] <- tmp$raw.pred
      colnames(raw.pred)[n] <- paste(paste0(names(vars),sum.pred[n,1:length(vars)]),collapse = "_")
      setTxtProgressBar(pb, n)
    }
    close(pb)
    res <- list(summary = sum.pred,raw = raw.pred)
    class(res) <- "bpd"
    return(res)
  }
}
