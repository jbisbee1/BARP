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
  res <- list("summary.pred" = data.frame('pred' = tmp.mean,'lb' = tmp.lb,'ub' = tmp.ub,row.names = c()),
              "raw.pred" = tmp.gibbs)
  return(res)
}

permute_test <- function(perm_mat,val) {
  pctile <- ecdf(quantile(perm_mat,seq(0,1,by = .01)))
  return(1 - pctile(val))
}


combine.factors <- function(mat,var,inter = F) {
  inds <- which(substr(colnames(mat),1,nchar(var)) == var)
  tmp.col <- apply(mat[,inds],1,sum)
  tmp.mat <- cbind(mat[,-inds],tmp.col)
  if(inter) {
    inds <- which(substr(rownames(mat),1,nchar(var)) == var)
    tmp.row <- apply(tmp.mat[inds,],2,sum)
    tmp.mat <- rbind(tmp.mat[-inds,],tmp.row)
    colnames(tmp.mat)[which(colnames(tmp.mat) == "tmp.col")] <- var
    rownames(tmp.mat)[which(rownames(tmp.mat) == "tmp.row")] <- var 
  } else {
    colnames(tmp.mat)[which(colnames(tmp.mat) == "tmp.col")] <- var
  }
  return(tmp.mat)
}


prep.dat <- function(dat,census,x,y) {
  classes <- unlist(lapply(dat,class))
  factors <- classes[which(classes %in% c("character","factor","integer"))]
  x.train <- dat[,which(colnames(dat) %in% x)]
  x.train[,which(colnames(x.train) %in% factors)] <- lapply(x.train[,which(colnames(x.train) %in% factors)],factor)
  y.train <- dat[,which(colnames(dat) %in% y)]
  x.test <- census[,which(colnames(census) %in% x)]
  x.test[,which(colnames(x.test) %in% factors)] <- lapply(x.test[,which(colnames(x.test) %in% factors)],factor)
  x_train <- as.data.frame(model.frame( ~ .-1, x.train,na.action = NULL))
  x_test <- as.data.frame(model.frame( ~ .-1,x.test,na.action  =NULL))
  missings <- colnames(x_test)[which(!colnames(x_test) %in% colnames(x_train))]
  if(length(missings) > 0) {
    warning(paste("Variables",paste(missings,collapse = ", ")," are not in the survey data.\n"),sep = "")
    for(miss in missings) {
      x_train$tmp <- 0
      colnames(x_train)[which(colnames(x_train) == "tmp")] <- miss
    }
  }
  missings <- colnames(x_train)[which(!colnames(x_train) %in% colnames(x_test))]
  if(length(missings) > 0) {
    warning(paste("Variables",paste(missings,collapse = ", ")," are not in the census data.\n"),sep = "")
    for(miss in missings) {
      x_test$tmp <- 0
      colnames(x_test)[which(colnames(x_test) == "tmp")] <- miss
    }
  }
  
  x_test <- x_test[,order(colnames(x_test))]
  x_train <- x_train[,order(colnames(x_train))]
  return(list(x_train = x_train,x_test = x_test,y.train = y.train))
}


quiet <- function(x) {
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 


barp.SL <- function (Y, X, newX = NULL, family = gaussian(), SL.library, 
                     method = "method.NNLS", id = NULL, verbose = FALSE, control = list(), 
                     cvControl = list(), obsWeights = NULL, env = parent.frame(),BSSD = FALSE) 
{
  
  time_start = proc.time()
  if (is.character(method)) {
    if (exists(method, mode = "list")) {
      method <- get(method, mode = "list")
    }
    else if (exists(method, mode = "function")) {
      method <- get(method, mode = "function")()
    }
  }
  else if (is.function(method)) {
    method <- method()
  }
  if (!is.list(method)) {
    stop("method is not in the appropriate format. Check out help('method.template')")
  }
  if (!is.null(method$require)) {
    sapply(method$require, function(x) require(force(x), 
                                               character.only = TRUE))
  }
  control <- do.call("SuperLearner.control", control)
  cvControl <- do.call("SuperLearner.CV.control", cvControl)
  library <- SuperLearner:::.createLibrary(SL.library)
  SuperLearner:::.check.SL.library(library = c(unique(library$library$predAlgorithm), 
                                               library$screenAlgorithm))
  call <- match.call(expand.dots = TRUE)
  if (!inherits(X, "data.frame")) 
    message("X is not a data frame. Check the algorithms in SL.library to make sure they are compatible with non data.frame inputs")
  varNames <- colnames(X)
  N <- dim(X)[1L]
  p <- dim(X)[2L]
  k <- nrow(library$library)
  kScreen <- length(library$screenAlgorithm)
  Z <- matrix(NA, N, k)
  libraryNames <- paste(library$library$predAlgorithm, library$screenAlgorithm[library$library$rowScreen], 
                        sep = "_")
  if (p < 2 & !identical(library$screenAlgorithm, "All")) {
    warning("Screening algorithms specified in combination with single-column X.")
  }
  fitLibEnv <- new.env()
  assign("fitLibrary", vector("list", length = k), envir = fitLibEnv)
  assign("libraryNames", libraryNames, envir = fitLibEnv)
  evalq(names(fitLibrary) <- libraryNames, envir = fitLibEnv)
  errorsInCVLibrary <- rep(0, k)
  errorsInLibrary <- rep(0, k)
  if (is.null(newX)) {
    newX <- X
  }
  if (!identical(colnames(X), colnames(newX))) {
    stop("The variable names and order in newX must be identical to the variable names and order in X")
  }
  if ((sum(is.na(X)) > 0 | sum(is.na(newX)) > 0 | sum(is.na(Y)) > 0) & !any(grepl("bartMachine",libraryNames))) {
    stop("missing data is currently not supported. Check Y, X, and newX for missing values")
  }
  if (!is.numeric(Y) & !any(grepl("SL.bartMachine",SL.library))) {
    stop("the outcome Y must be a numeric vector unless running with bartMachine.")
  }
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if (family$family != "binomial" & isTRUE("cvAUC" %in% method$require)) {
    stop("'method.AUC' is designed for the 'binomial' family only")
  }
  validRows <- CVFolds(N = N, id = id, Y = Y, cvControl = cvControl)
  if (is.null(id)) {
    id <- seq(N)
  }
  if (!identical(length(id), N)) {
    stop("id vector must have the same dimension as Y")
  }
  if (is.null(obsWeights)) {
    obsWeights <- rep(1, N)
  }
  if (!identical(length(obsWeights), N)) {
    stop("obsWeights vector must have the same dimension as Y")
  }
  .crossValFUN <- function(valid, dataY, dataX, id, obsWeights, 
                           library, kScreen, k, p, libraryNames) {
    
    tempLearn <- dataX[-valid, , drop = FALSE]
    tempOutcome <- dataY[-valid]
    tempValid <- dataX[valid, , drop = FALSE]
    tempWhichScreen <- matrix(NA, nrow = kScreen, ncol = p)
    tempId <- id[-valid]
    tempObsWeights <- obsWeights[-valid]
    
    for (s in seq(kScreen)) {
      screen_fn = get(library$screenAlgorithm[s], envir = env)
      testScreen <- try(do.call(screen_fn, list(Y = tempOutcome, 
                                                X = tempLearn, family = family, id = tempId, 
                                                obsWeights = tempObsWeights)))
      if (inherits(testScreen, "try-error")) {
        warning(paste("replacing failed screening algorithm,", 
                      library$screenAlgorithm[s], ", with All()", 
                      "\n "))
        tempWhichScreen[s, ] <- TRUE
      }
      else {
        tempWhichScreen[s, ] <- testScreen
      }
      if (verbose) {
        message(paste("Number of covariates in ", library$screenAlgorithm[s], 
                      " is: ", sum(tempWhichScreen[s, ]), sep = ""))
      }
    }
    out <- matrix(NA, nrow = nrow(tempValid), ncol = k)
    for (s in seq(k)) {
      if(!grepl("bartMachine|BARP",library$library$predAlgorithm[s])) {
        tmpX <- subset(tempLearn, select = tempWhichScreen[library$library$rowScreen[s], ], drop = FALSE)
        classes <- sapply(tmpX,class)
        factors <- classes[which(classes %in% c("character","factor","integer"))]
        tmpX[,which(colnames(tmpX) %in% factors)] <- lapply(tmpX[,which(colnames(tmpX) %in% factors)],factor)
        tmpX <- as.data.frame(model.matrix( ~ .-1, tmpX))
        
        tmpNewX <- subset(tempValid, select = tempWhichScreen[library$library$rowScreen[s], ], drop = FALSE)
        classes <- sapply(tmpNewX,class)
        factors <- classes[which(classes %in% c("character","factor","integer"))]
        tmpNewX[,which(colnames(tmpNewX) %in% factors)] <- lapply(tmpNewX[,which(colnames(tmpNewX) %in% factors)],factor)
        tmpNewX <- as.data.frame(model.matrix( ~ .-1, tmpNewX))
        missings <- colnames(tmpNewX)[which(!colnames(tmpNewX) %in% colnames(tmpX))]
        if(length(missings) > 0) {
          for(miss in missings) {
            tmpX$tmp <- 0
            colnames(tmpX)[which(colnames(tmpX) == "tmp")] <- miss
          }
        }
        missings <- colnames(tmpX)[which(!colnames(tmpX) %in% colnames(tmpNewX))]
        if(length(missings) > 0) {
          for(miss in missings) {
            tmpNewX$tmp <- 0
            colnames(tmpNewX)[which(colnames(tmpNewX) == "tmp")] <- miss
          }
        }
      } else {
        tmpX <- subset(tempLearn, select = tempWhichScreen[library$library$rowScreen[s], ], drop = FALSE)
        tmpNewX <- subset(tempValid, select = tempWhichScreen[library$library$rowScreen[s], ], drop = FALSE)
      }
      tmpX <- tmpX[,order(colnames(tmpX))]
      tmpNewX <- tmpNewX[,order(colnames(tmpNewX))]
      
      pred_fn = get(library$library$predAlgorithm[s], envir = env)
      
      testAlg <- try(do.call(pred_fn, list(Y = tempOutcome, 
                                           X = tmpX, 
                                           newX = tmpNewX, 
                                           family = family, id = tempId, 
                                           obsWeights = tempObsWeights)))
      if (inherits(testAlg, "try-error")) {
        warning(paste("Error in algorithm", library$library$predAlgorithm[s], 
                      "\n  The Algorithm will be removed from the Super Learner (i.e. given weight 0) \n"))
      }
      else {
        out[, s] <- testAlg$pred
      }
      if (verbose) 
        message(paste("CV", libraryNames[s]))
    }
    invisible(out)
  }
  
  # Beginning training
  time_train_start = proc.time()
  if(k == 1) {
    if(!BSSD) {
      cat("Only one algorithm. Skipping training.\n")
    }
    coef <- 1
    names(coef) <- libraryNames
    getCoef <- list(cvRisk = NA)
  } else {
    if(!BSSD) {
      cat("Starting training...\n")
    }
    Z[unlist(validRows, use.names = FALSE), ] <- do.call("rbind", 
                                                         lapply(validRows, FUN = .crossValFUN, dataY = Y, dataX = X, 
                                                                id = id, obsWeights = obsWeights, library = library, 
                                                                kScreen = kScreen, k = k, p = p, 
                                                                libraryNames = libraryNames))
    errorsInCVLibrary <- apply(Z, 2, function(x) anyNA(x))
    if (sum(errorsInCVLibrary) > 0) {
      Z[, as.logical(errorsInCVLibrary)] <- 0
    }
    if (all(Z == 0)) {
      stop("All algorithms dropped from library")
    }
    getCoef <- method$computeCoef(Z = Z, Y = Y, libraryNames = libraryNames, 
                                  obsWeights = obsWeights, control = control, verbose = verbose, 
                                  errorsInLibrary = errorsInCVLibrary)
    coef <- getCoef$coef
    names(coef) <- libraryNames
  }
  time_train = proc.time() - time_train_start
  if (!("optimizer" %in% names(getCoef))) {
    getCoef["optimizer"] <- NA
  }
  
  # Beginning prediction
  m <- dim(newX)[1L]
  predY <- matrix(NA, nrow = m, ncol = k)
  .screenFun <- function(fun, list) {
    screen_fn = get(fun, envir = env)
    testScreen <- try(do.call(screen_fn, list))
    if (inherits(testScreen, "try-error")) {
      warning(paste("replacing failed screening algorithm,", 
                    fun, ", with All() in full data", "\n "))
      out <- rep(TRUE, ncol(list$X))
    }
    else {
      out <- testScreen
    }
    return(out)
  }
  time_predict_start = proc.time()
  whichScreen <- sapply(library$screenAlgorithm, FUN = .screenFun, 
                        list = list(Y = Y, X = X, family = family, id = id, obsWeights = obsWeights), 
                        simplify = FALSE)
  whichScreen <- do.call(rbind, whichScreen)
  .predFun <- function(index, lib, dataY, dataX, newX, whichScreen, 
                       family, id, obsWeights, verbose, control, libraryNames) {
    pred_fn = get(lib$predAlgorithm[index], envir = env)
    if(!grepl("bartMachine",lib$predAlgorithm[index])) {
      tmpX <- subset(dataX, select = whichScreen[lib$rowScreen[index], ], drop = FALSE)
      classes <- sapply(tmpX,class)
      factors <- classes[which(classes %in% c("character","factor","integer"))]
      tmpX[,which(colnames(tmpX) %in% factors)] <- lapply(tmpX[,which(colnames(tmpX) %in% factors)],factor)
      tmpX <- as.data.frame(model.matrix( ~ .-1, tmpX))
      
      tmpNewX <- subset(newX, select = whichScreen[lib$rowScreen[index], ], drop = FALSE)
      classes <- sapply(tmpNewX,class)
      factors <- classes[which(classes %in% c("character","factor","integer"))]
      tmpNewX[,which(colnames(tmpNewX) %in% factors)] <- lapply(tmpNewX[,which(colnames(tmpNewX) %in% factors)],factor)
      tmpNewX <- as.data.frame(model.matrix( ~ .-1, tmpNewX))
      missings <- colnames(tmpNewX)[which(!colnames(tmpNewX) %in% colnames(tmpX))]
      if(length(missings) > 0) {
        for(miss in missings) {
          tmpX$tmp <- 0
          colnames(tmpX)[which(colnames(tmpX) == "tmp")] <- miss
        }
      }
      missings <- colnames(tmpX)[which(!colnames(tmpX) %in% colnames(tmpNewX))]
      if(length(missings) > 0) {
        for(miss in missings) {
          tmpNewX$tmp <- 0
          colnames(tmpNewX)[which(colnames(tmpNewX) == "tmp")] <- miss
        }
      }
    } else {
      tmpX <- subset(dataX, select = whichScreen[lib$rowScreen[index], ], drop = FALSE)
      tmpNewX <- subset(newX, select = whichScreen[lib$rowScreen[index], ], drop = FALSE)
    }
    tmpX <- tmpX[,order(colnames(tmpX))]
    tmpNewX <- tmpNewX[,order(colnames(tmpNewX))]
    testAlg <- try(do.call(pred_fn, list(Y = dataY, X = tmpX, 
                                         newX = tmpNewX, family = family, id = id, obsWeights = obsWeights)))
    if (inherits(testAlg, "try-error")) {
      warning(paste("Error in algorithm", lib$predAlgorithm[index], 
                    " on full data", "\n  The Algorithm will be removed from the Super Learner (i.e. given weight 0) \n"))
      out <- rep.int(NA, times = nrow(newX))
    }
    else {
      out <- testAlg$pred
      if (control$saveFitLibrary) {
        eval(bquote(fitLibrary[[.(index)]] <- .(testAlg$fit)), 
             envir = fitLibEnv)
      }
    }
    if (verbose) {
      message(paste("full", libraryNames[index]))
    }
    invisible(out)
  }
  if(!BSSD) {
    cat("Starting prediction...\n")
  }
  predY <- do.call("cbind", lapply(seq(k), FUN = .predFun, 
                                   lib = library$library, dataY = Y, dataX = X, newX = newX, 
                                   whichScreen = whichScreen, family = family, id = id, 
                                   obsWeights = obsWeights, verbose = verbose, control = control, 
                                   libraryNames = libraryNames))
  errorsInLibrary <- apply(predY, 2, function(algorithm) anyNA(algorithm))
  if (sum(errorsInLibrary) > 0) {
    if (sum(coef[as.logical(errorsInLibrary)]) > 0) {
      warning(paste0("Re-running estimation of coefficients removing failed algorithm(s)\n", 
                     "Original coefficients are: \n", paste(coef, 
                                                            collapse = ", "), "\n"))
      Z[, as.logical(errorsInLibrary)] <- 0
      if (all(Z == 0)) {
        stop("All algorithms dropped from library")
      }
      getCoef <- method$computeCoef(Z = Z, Y = Y, libraryNames = libraryNames, 
                                    obsWeights = obsWeights, control = control, verbose = verbose, 
                                    errorsInLibrary = errorsInLibrary)
      coef <- getCoef$coef
      names(coef) <- libraryNames
    }
    else {
      warning("Coefficients already 0 for all failed algorithm(s)")
    }
  }
  getPred <- method$computePred(predY = predY, coef = coef, 
                                control = control)
  time_predict = proc.time() - time_predict_start
  colnames(predY) <- libraryNames
  if (sum(errorsInCVLibrary) > 0) {
    getCoef$cvRisk[as.logical(errorsInCVLibrary)] <- NA
  }
  time_end = proc.time()
  times = list(everything = time_end - time_start, train = time_train, 
               predict = time_predict)
  out <- list(call = call, libraryNames = libraryNames, SL.library = library, 
              SL.predict = getPred, coef = coef, library.predict = predY, 
              Z = Z, cvRisk = getCoef$cvRisk, family = family, 
              fitLibrary = get("fitLibrary", envir = fitLibEnv), 
              varNames = varNames, validRows = validRows, 
              method = method, whichScreen = whichScreen, control = control, 
              cvControl = cvControl, errorsInCVLibrary = errorsInCVLibrary, 
              errorsInLibrary = errorsInLibrary, metaOptimizer = getCoef$optimizer, 
              env = env, times = times)
  class(out) <- c("SuperLearner")
  return(out)
}


SL.bartMachine <- function (Y, X, newX, family, obsWeights, id, num_trees = 50,
                            num_burn_in = 250, verbose = F, alpha = 0.95, beta = 2, k = 2,
                            q = 0.9, nu = 3, num_iterations_after_burn_in = 1000,
                            prob_rule_class = 0.5, mh_prob_steps = c(2.5, 2.5, 4)/9,
                            debug_log = FALSE, run_in_sample = TRUE, s_sq_y = "mse",
                            sig_sq_est = NULL, cov_prior_vec = NULL, use_missing_data = FALSE,
                            covariates_to_permute = NULL, num_rand_samps_in_library = 10000,
                            use_missing_data_dummies_as_covars = FALSE, replace_missing_data_with_x_j_bar = FALSE,
                            impute_missingness_with_rf_impute = FALSE, impute_missingness_with_x_j_bar_for_lm = TRUE,
                            mem_cache_for_speed = TRUE, serialize = FALSE, seed = NULL, ...)
{
  SuperLearner:::.SL.require("bartMachine")
  model = bartMachine::bartMachine(X, Y, num_trees = num_trees,
                                   num_burn_in = num_burn_in, verbose = verbose, alpha = alpha,
                                   beta = beta, k = k, q = q, nu = nu,
                                   num_iterations_after_burn_in = num_iterations_after_burn_in,
                                   prob_rule_class = prob_rule_class, mh_prob_steps = mh_prob_steps,
                                   debug_log = debug_log, run_in_sample = run_in_sample, s_sq_y = s_sq_y,
                                   sig_sq_est = sig_sq_est, cov_prior_vec = cov_prior_vec, use_missing_data = use_missing_data,
                                   covariates_to_permute = covariates_to_permute, num_rand_samps_in_library = num_rand_samps_in_library,
                                   use_missing_data_dummies_as_covars = use_missing_data_dummies_as_covars,
                                   replace_missing_data_with_x_j_bar = replace_missing_data_with_x_j_bar,
                                   impute_missingness_with_rf_impute = impute_missingness_with_rf_impute,
                                   impute_missingness_with_x_j_bar_for_lm = impute_missingness_with_x_j_bar_for_lm,
                                   mem_cache_for_speed = mem_cache_for_speed, serialize = serialize, seed = seed)
  pred <- predict(model, newX)
  fit <- list(object = model)
  class(fit) <- c("SL.bartMachine")
  out <- list(pred = pred, fit = fit)
  return(out)
}

