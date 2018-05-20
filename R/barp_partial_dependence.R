#' barp_partial_dependence
#' 
#' This function calculates the partial dependence for selected covariates based on a \code{BARP} object. The user can specify up to three covariates, resulting in a three-way interaction. 
#' @param BARP A \code{BARP} object returned by using the \code{BARP} function.
#' @param vars The variable names or column indices of interest.
#' @param prop_data The percentage of the original data to use for estimation. Larger values are more time consuming.
#' @param levs A list of values at which to estimate the partial dependence for each variable. If NULL (the default), quantiles c(0.05,seq(.1,.9,by = .1),.95) are used.
#' @param credible_interval The credible interval for the estimates.
#' @param ... Additional arguments to be passed to bartMachine.
#' @return Returns an object of class "barpd", containing a list of the following components:
#' \item{summary}{A \code{data.frame} listing the mean and upper and lower bounds of the predicted outcome at each level of the covariate(s).}
#' \item{raw}{A \code{data.frame} where each row is a posterior draw and each column is a level for the covariate(s).}
#' @keywords MRP BARP Mister P multilevel regression poststratification
#' @seealso \code{\link{BARP}} which generates the prerequisite \code{BARP} object.
#' @examples 
#' data("gaymar")
#' BARP <- barp(y = "supp_gaymar",
#'              x = c("pvote","religcon","age","educ","gXr","stateid","region"),
#'              dat = svy,
#'              census = census06,
#'              geo.unit = "stateid",
#'              proportion = "n")
#' bpd <- barp_partial_dependence(BARP = BARP,
#'                                vars = c("region","educ"),
#'                                prop_data = .2,
#'                                levs = list(c(1:4),c(1:4)),
#'                                credible_interval = c(0.025,0.975))
#' @export barp

barp_partial_dependence <- function (BARP, vars = c("region","educ"), prop_data = .1,
                                     levs = NULL,credible_interval = c(0.025,.975))
  
{
  bart_machine <- BARP$trees
  bartMachine:::check_serialization(bart_machine)
  for(v in 1:length(vars)) {
    if (class(vars[v]) == "integer") {
      vars[v] = as.numeric(vars[v])
    }
    if (class(vars[v]) == "numeric" && (vars[v] < 1 || vars[v] > bart_machine$p)) {
      stop(paste("You must set var ",v, " to a number between 1 and p =", 
                 bart_machine$p))
    } else if (class(vars[v]) == "character" && !(vars[v] %in% bart_machine$training_data_features)) {
      stop(paste0("var ",v," must be the name of one of the training features (see \"<bart_model>$training_data_features\")"))
    } else if (!(class(vars[v]) == "numeric" || class(vars[v]) == "character")) {
      stop(paste0("var ",v," must be a column number or column name"))
    }
  }
  
  if(is.null(levs)) {
    quants <- c(.05,seq(.1,.9,by = .1),.95)
    levs <- list(NA)
    for(v in 1:length(vars)) {
      levs[[v]] <- unique(quantile(bart_machine$model_matrix_training_data[,vars[v]], quants, na.rm = TRUE))
      if (length(unique(levs[[v]])) <= 1) {
        warning(paste0("There must be more than one unique value in",vars[v],"."))
        return()
      }
    }
  } else if(class(levs) == "list" && length(levs) != length(vars)) {
    stop("Not enough levs for all provided vars.")
  } else if(class(levs) != "list") {
    stop("The levs argument should be a list of values for each provided variable.")
  }
  
  n_samp = round(bart_machine$n * prop_data)
  
  sum.pred <- sum.prob <- as.matrix(expand.grid(levs))
  sum.pred <- as.data.frame(cbind(sum.pred,matrix(NA,nrow = nrow(sum.pred),ncol = 3)))
  # sum.prob <- as.data.frame(cbind(sum.prob,matrix(NA,nrow = nrow(sum.prob),ncol = 3)))
  colnames(sum.pred) <- c(vars,"pred","lb","ub")
  # colnames(sum.prob) <- c(vars,"pred","lb","ub")
  raw.pred <- raw.prob <- as.data.frame(matrix(NA,nrow = bart_machine$num_iterations_after_burn_in,ncol = nrow(sum.pred)))
  cat(paste0("\nCalculating partial dependencies for ",paste(vars,collapse = " by ")," using ",prop_data*100,"% of the data.\n"))
  pb <- txtProgressBar(min = 0, max = nrow(sum.pred), style = 3)
  for(n in 1:nrow(sum.pred)) {
    indices = sample(1:bart_machine$n, n_samp)
    test_data = bart_machine$X[indices, ]
    for(col in 1:length(vars)) {
      test_data[, vars[col]] = rep(sum.pred[n,col], n_samp)
    }
    tmp <- extract_predictions(bart_machine,test_data,credible_interval)
    sum.pred[n,c("pred","lb","ub")] <- tmp$summary.pred
    
    raw.pred[,n] <- tmp$raw.pred
    # raw.prob[,n] <- tmp$raw.prob
    # sum.prob[n,c("pred","lb","ub")] <- tmp$summary.prob
    colnames(raw.pred)[n] <- paste(paste0(vars,sum.pred[n,1:length(vars)]),collapse = "_")
    setTxtProgressBar(pb, n)
  }
  close(pb)
  res <- list(summary = sum.pred,raw = raw.pred) #,probits = list(summary = sum.prob,raw = raw.prob))
  class(res) <- "barpd"
  return(res)
}