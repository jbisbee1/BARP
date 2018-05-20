#' plot.barpcov
#' 
#' This function plots the covariate importance. 
#' @param barpcov A \code{BARP} covariate importance object.
#' @param var_names A vector of variable names. If \code{NULL} (the default) the variable names are taken from the \code{BARP} covariate importance object.
#' @param sig_level The significance level at which to evaluate covariate significance. Only applicable for BARP covariate importance objects that included a permutation test. Defaults to \code{0.05}.
#' @keywords MRP BARP Mister P multilevel regression poststratification
#' @seealso \code{\link{barp_prognostic_covs}} which creates the prerequisite \code{barpcov} object.
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
#' plot.barpcov(barpcov,
#'              var_names = c("Presidential Vote","Religious-Conservative",
#'                            "Age","Education","Gender X Race","State","Region"),
#'              sig_level = 0.05)
#' @export barp

plot.barpcov <- function(barpcov,
                         var_names = NULL,
                         sig_level = NULL) 
{
  ords <- order(apply(barpcov$covariate_importance,2,mean),decreasing = T)
  if(is.null(var_names)) {
    var_names <- colnames(barpcov$covariate_importance)
  }
  var_names <- var_names[ords]
  barpcov$covariate_importance <- barpcov$covariate_importance[,ords]
  if(!is.null(sig_level) & !any(grepl("permutation_test",names(barpcov)))) {
    warning("You have indicated a significance threshold for evaluating covariate importance without running permutation analysis. As such, significance is not assessable.")
  }
  
  df <- data.frame(apply(barpcov$covariate_importance,2,mean))
  colnames(df)[1] <- "mean"
  df$var <- factor(var_names,levels = rev(var_names))
  df$sd <- apply(barpcov$covariate_importance,2,sd)
  if(any(grepl("permutation_test",names(barpcov)))) {
    barpcov$permutation_test <- barpcov$permutation_test[,ords]
    if(is.null(sig_level)) {
      sig_level <- .1
    }
    df$perm <- apply(barpcov$permutation_test,2,quantile,1-sig_level)
    df$pvals <- barpcov$p_vals[ords]
    df$cols <- ifelse(df$pvals < sig_level,rgb(0,0,0,.6),rgb(0,0,0,.3))
    title <- paste0("Covariate Significance at ",round((sig_level)*100,0),"%")
    subtitle <- paste0("Based on ",paste0(toupper(substr(barpcov$type,1,1)),
                                          substr(barpcov$type,2,10))," with Permutation")
  } else {
    df$perm <- 0
    df$cols <- rgb(0,0,0,.3)
    title <- "Covariate Importance"
    subtitle <- paste0("Based on ",paste0(toupper(substr(barpcov$type,1,1)),
                                          substr(barpcov$type,2,10)))
  }
  
  ggplot(df, aes(var,mean))+ 
    geom_bar(stat = "identity",fill = rev(df$cols)) +
    coord_flip() + 
    labs(title = title,
         subtitle = subtitle,
         y = paste0("Share of ",paste0(toupper(substr(barpcov$type,1,1)),
                                       substr(barpcov$type,2,10))),
         x = "") + 
    geom_errorbar(aes(ymin = mean - sd,ymax = mean+sd),width = .2,col = rgb(0,0,0,.4)) + 
    geom_errorbar(aes(ymin = perm,ymax = perm),col = "black") +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
}
