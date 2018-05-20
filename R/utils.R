.javaMem <- function(cmd) {
  return(tryCatch({
    cmd
  },warning = function(w) {
    cat(paste0(w))
  }, error = function(e) {
    if(grepl("OutOfMemoryError",e)) {
      cat(paste0("Not enough memory allocated to Java.\n",
                 "Please restart R and run the command:\n",
                 "    options(java.parameters = '-Xmx[NUM]g')\n",
                 "before loading BARP where [NUM] is the amount of RAM in gb.\n",
                 "I.e. to allocate 3gb, run:\n",
                 "    options(java.parameters = '-Xmx3g')\n",
                 "    library(BARP)"))
    } else {
      cat(paste0(e))
    }
  }))
}

extract_predictions <- function(bart_machine,test_data,credible_interval) {
  tmp <- bart_machine_get_posterior(bart_machine, test_data)$y_hat_posterior_samples
  tmp.gibbs <- apply(tmp,2,mean)
  tmp.mean <- mean(tmp.gibbs)
  tmp.lb <- quantile(tmp.gibbs,credible_interval[1])
  tmp.ub <- quantile(tmp.gibbs,credible_interval[2])
  
  # if(bart_machine$pred_type == "classification") {
  #   tmp.prob <- qnorm(tmp)
  #   tmp.prob.gibbs <- apply(tmp.prob,2,mean)
  #   tmp.prob.mean <- mean(tmp.prob.gibbs)
  #   tmp.prob.lb <- quantile(tmp.prob.gibbs,lower_ci)
  #   tmp.prob.ub <- quantile(tmp.prob.gibbs,upper_ci)
  #   
  #   res <- list("summary.pred" = data.frame('pred' = tmp.mean,'lb' = tmp.lb,'ub' = tmp.ub,row.names = c()),
  #               "summary.prob" = data.frame('pred' = tmp.prob.mean,'lb' = tmp.prob.lb,'ub' = tmp.prob.ub,row.names = c()),
  #               "raw.pred" = tmp.gibbs,
  #               "raw.prob" = tmp.prob.gibbs)
  # } else {
  res <- list("summary.pred" = data.frame('pred' = tmp.mean,'lb' = tmp.lb,'ub' = tmp.ub,row.names = c()),
              "raw.pred" = tmp.gibbs)
  return(res)
}

permute_test <- function(perm_mat,val) {
  pctile <- ecdf(quantile(perm_mat,seq(0,1,by = .01)))
  return(1 - pctile(val))
}