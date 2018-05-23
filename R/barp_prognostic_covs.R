#' barp_prognostic_covs
#' 
#' This function extracts the relative prognostic power of the covariates used in the BARP model. The user can choose to evaluate the statistical significance of these metrics by setting perm_test to TRUE. Doing so re-esetimates variable inclusion proportions (in either trees or splits) after randomly permuting the outcome variable. These permuted inclusion proportions are used as a null distribution against which the average observed proportions are compared.
#' @param BARP A \code{BARP} object.
#' @param interactions A logical statement for whether to evaluate individual variables or their pairwise interactions. Defaults to \code{FALSE}.
#' @param perm_test A logical statement for whether to evaluate covariate significance via permutation tests. Defaults to \code{FALSE}.
#' @param num_reps The number of reps used to estimate the average inclusion proportions and standard deviation.
#' @param num_trees The number of trees to be used.
#' @param num_permute The number of permutation simulations.
#' @param type The context in which to evaluate the proportion of variable inclusion, either in terms of \code{trees} or in terms of \code{splits}.
#' @return Returns an object of class "barpcov", containing a list of the following components:
#' \item{covariate_importance}{A \code{matrix} containing the raw estimated inclusion proportions for the covariates.}
#' \item{type}{An indicator for which \code{type} was used to estimate variable importance.}
#' \item{permutation_test}{A \code{matrix} containing the raw estimated inclusion proportions for the covariates estimated under the permuted null relationship. Only included if \code{perm_test = TRUE}.}
#' \item{p_vals}{A \code{vector} containing the p-values for each variable from the permutation test. These specifically correspond to the percentabe of the null proportions that are less than or equal to the observed average proportion. Only included if \code{perm_test = TRUE}.}
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
#' barpcov <- barp_prognostic_covs(BARP,
#'                                 interactions = F,
#'                                 perm_test = T,
#'                                 num_reps = 30,
#'                                 num_trees = 20,
#'                                 num_permute = 30,
#'                                 type = 'splits')
#' @export

barp_prognostic_covs <- function(BARP,
                                 interactions = F,
                                 perm_test = T,
                                 num_reps = 10,
                                 num_trees = 20,
                                 num_permute = 10,
                                 type = "splits",...)
{
  bart_machine <- BARP$trees
  factors <- names(bart_machine$X)[which(sapply(bart_machine$X, class) == "factor")]
  
  if(!type %in% c("trees","splits")) {
    stop("'type' must be one of either 'trees' or 'splits'.")
  }
  if(interactions) {
    cat(paste0("Calculating covariate importance with interactions based on ",type,".\n"))
    interaction_counts = array(NA, c(bart_machine$p, bart_machine$p,num_reps))
    pb <- txtProgressBar(min = 0, max = num_reps, style = 3)
    for (r in 1:num_reps) {
      if (r == 1 & num_trees == bart_machine$num_trees) {
        interaction_counts[, , r] = sapply(.jcall(bart_machine$java_bart_machine, "[[I", "getInteractionCounts"), .jevalArray)
      }
      else {
        bart_machine_dup = bartMachine:::bart_machine_duplicate(bart_machine, num_trees = num_trees)
        interaction_counts[, , r] = sapply(.jcall(bart_machine_dup$java_bart_machine, "[[I", "getInteractionCounts"), .jevalArray)
      }
      setTxtProgressBar(pb, r)
    }
    close(pb)
    rownames(interaction_counts) = colnames(interaction_counts) = bart_machine$training_data_features
    
    if(length(factors) > 0) {
      interaction_counts2 <- array(NA,c(ncol(bart_machine$X),ncol(bart_machine$X),num_reps))
      for(r in 1:num_reps) {
        tmp <- combine.factors(interaction_counts[,,r],factors[1],inter = T)
        for(f in 2:length(factors)) {
          tmp <- combine.factors(tmp,factors[f],inter = T)
        }
        interaction_counts2[,,r] <- tmp
      }
      
      rownames(interaction_counts2) = colnames(interaction_counts2)  = rownames(tmp)
      interaction_counts <- interaction_counts2
    }
    num_total_interactions = nrow(interaction_counts) * (nrow(interaction_counts) + 1)/2
    inter_counts <- matrix(NA,nrow = num_reps,ncol = num_total_interactions)
    colnames(inter_counts) <- rep("tmp",num_total_interactions)
    for(r in 1:num_reps) {
      c <- 1
      for(i in 1:nrow(interaction_counts)) {
        for(j in 1:ncol(interaction_counts)) {
          if(j <= i) {
            inter_counts[r,c] <- interaction_counts[i,j,r]
            colnames(inter_counts)[c] <- paste(rownames(interaction_counts[,,r])[i], 
                                               "x",colnames(interaction_counts[,,r])[j])
            c <- c + 1
          }
        }
      }
    }
    inter_counts <- inter_counts / apply(inter_counts,1,sum)
    if(perm_test) {
      cat(paste0("Permutation test for significance.\n"))
      pb <- txtProgressBar(min = 0, max = num_permute, style = 3)
      permute_counts = array(NA, c(bart_machine$p, bart_machine$p,num_replicates_for_avg=num_permute))
      for (r in 1:num_permute) {
        y_permuted = sample(bart_machine$y, replace = FALSE)
        bart_machine_with_permuted_y = build_bart_machine(bart_machine$X, 
                                                          y_permuted, 
                                                          num_trees = as.numeric(num_trees), 
                                                          num_burn_in = bart_machine$num_burn_in, 
                                                          num_iterations_after_burn_in = bart_machine$num_iterations_after_burn_in, 
                                                          run_in_sample = FALSE, use_missing_data = bart_machine$use_missing_data, 
                                                          use_missing_data_dummies_as_covars = bart_machine$use_missing_data_dummies_as_covars, 
                                                          num_rand_samps_in_library = bart_machine$num_rand_samps_in_library, 
                                                          replace_missing_data_with_x_j_bar = bart_machine$replace_missing_data_with_x_j_bar, 
                                                          impute_missingness_with_rf_impute = bart_machine$impute_missingness_with_rf_impute, 
                                                          impute_missingness_with_x_j_bar_for_lm = bart_machine$impute_missingness_with_x_j_bar_for_lm, 
                                                          verbose = FALSE)
        
        permute_counts[, , r] = sapply(.jcall(bart_machine_with_permuted_y$java_bart_machine, "[[I", "getInteractionCounts"), .jevalArray)
        setTxtProgressBar(pb, r)
      }
      close(pb)
      rownames(permute_counts) = colnames(permute_counts) = bart_machine$training_data_features
      
      if(length(factors) > 0) {
        permute_counts2 <- array(NA,c(ncol(bart_machine$X),ncol(bart_machine$X),num_permute))
        for(r in 1:num_permute) {
          tmp <- combine.factors(permute_counts[,,r],factors[1],inter = T)
          for(f in 2:length(factors)) {
            tmp <- combine.factors(tmp,factors[f],inter = T)
          }
          permute_counts2[,,r] <- tmp
        }
        
        rownames(permute_counts2) = colnames(permute_counts2)  = rownames(tmp)
        permute_counts <- permute_counts2
      }
      
      perm_counts <- matrix(NA,nrow = num_permute,ncol = num_total_interactions)
      colnames(perm_counts) <- rep("tmp",num_total_interactions)
      for(r in 1:num_permute) {
        c <- 1
        for(i in 1:nrow(permute_counts)) {
          for(j in 1:nrow(permute_counts)) {
            if(j <= i) {
              perm_counts[r,c] <- permute_counts[i,j,r]
              colnames(perm_counts)[c] <- paste(rownames(permute_counts[,,r])[i], 
                                                "x",colnames(permute_counts[,,r])[j])
              c <- c + 1
            }
          }
        }
      }
      perm_counts <- perm_counts / apply(perm_counts,1,sum)
      p.vals <- sapply(1:ncol(perm_counts), function(x) permute_test(perm_counts[,x],mean(inter_counts[,x])))
      names(p.vals) <- colnames(perm_counts)
      res <- list(covariate_importance = inter_counts,permutation_test = perm_counts,p_vals = p.vals)
    } else {
      res <- list(covariate_importance = inter_counts)
    }
  } else {
    cat(paste0("Calculating covariate importance based on ",type,".\n"))
    var_props = array(0, c(num_reps, bart_machine$p))
    pb <- txtProgressBar(min = 0, max = num_reps, style = 3)
    for (i in 1:num_reps) {
      bart_machine_dup = bartMachine:::bart_machine_duplicate(bart_machine, num_trees = num_trees)#,...)
      var_props[i,] = .jcall(bart_machine_dup$java_bart_machine, "[D", "getAttributeProps", type)
      setTxtProgressBar(pb, i)
    }
    close(pb)
    colnames(var_props) = bart_machine$training_data_features_with_missing_features
    if(length(factors) > 0) {
      tmp <- combine.factors(var_props,factors[1],inter = F)
      for(f in 2:length(factors)) {
        tmp <- combine.factors(tmp,factors[f],inter = F)
      }
      var_props <- tmp
    }
    
    
    
    if(perm_test) {
      cat(paste0("Permutation test for significance.\n"))
      permute_mat = matrix(NA, nrow = num_permute, ncol = bart_machine$p)
      colnames(permute_mat) = bart_machine$training_data_features_with_missing_features
      pb <- txtProgressBar(min = 0, max = num_permute, style = 3)
      for (p in 1:num_permute) {
        y_permuted = sample(bart_machine$y, replace = FALSE)
        bart_machine_with_permuted_y = build_bart_machine(bart_machine$X, 
                                                          y_permuted, 
                                                          num_trees = as.numeric(num_trees), 
                                                          num_burn_in = bart_machine$num_burn_in, 
                                                          num_iterations_after_burn_in = bart_machine$num_iterations_after_burn_in, 
                                                          run_in_sample = FALSE, use_missing_data = bart_machine$use_missing_data, 
                                                          use_missing_data_dummies_as_covars = bart_machine$use_missing_data_dummies_as_covars, 
                                                          num_rand_samps_in_library = bart_machine$num_rand_samps_in_library, 
                                                          replace_missing_data_with_x_j_bar = bart_machine$replace_missing_data_with_x_j_bar, 
                                                          impute_missingness_with_rf_impute = bart_machine$impute_missingness_with_rf_impute, 
                                                          impute_missingness_with_x_j_bar_for_lm = bart_machine$impute_missingness_with_x_j_bar_for_lm, 
                                                          verbose = FALSE)
        permute_mat[p,] = .jcall(bart_machine_with_permuted_y$java_bart_machine, "[D", "getAttributeProps", type)
        setTxtProgressBar(pb, p)
      }
      close(pb)
      if(length(factors) > 0) {
        tmp <- combine.factors(permute_mat,factors[1],inter = F)
        for(f in 2:length(factors)) {
          tmp <- combine.factors(tmp,factors[f],inter = F)
        }
        permute_mat <- tmp
      }
      p.vals <- sapply(1:ncol(permute_mat), function(x) permute_test(permute_mat[,x],mean(var_props[,x])))
      names(p.vals) <- colnames(permute_mat)
      res <- list(covariate_importance = var_props,permutation_test = permute_mat,p_vals = p.vals)
    } else {
      res <- list(covariate_importance = var_props)
    }
  }
  res$type <- type
  class(res) <- "barpcov"
  return(res)
}
