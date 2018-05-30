#' plot.bpd
#' 
#' This function plots the partial dependence results estimated returned by running barp_partial_dependence.
#' @param bpd A \code{\link{barp_partial_dependence} partial dependence object of class "bpd".
#' @param var_names A vector of variable names. If \code{NULL} (the default), the variable names from the training data columns are used.
#' @param var_labs A list of variable labels corresponding to the levels at which the partial dependence was calculated. If \code{NULL} (the default), the raw values are used.
#' @param is_categorical A vector of logicals (T,F) corresponding to which variables are categorical and which are continuous.If \code{NULL} (the default), the function determines whether each variable is categorical or continuous based on the data class.
#' @param ... Additional arguments to be passed to plot.
#' @keywords MRP BARP Mister P multilevel regression poststratification
#' @seealso \code{\link{barp_partial_dependence}} which is a prerequisite to plotting.
#' @examples 
#' data("gaymar")
#' barp.obj <- barp(y = "supp_gaymar",
#'              x = c("pvote","religcon","age","educ","gXr","stateid","region"),
#'              dat = svy,
#'              census = census06,
#'              geo.unit = "stateid",
#'              proportion = "n")
#' bpd <- barp_partial_dependence(BARP = barp.obj,
#'                                vars = c("age","educ"),
#'                                prop_data = .2,
#'                                levs = list(c(1:4),c(1:4),
#'                                credible_interval = c(0.025,0.975)))
#' plot(bpd,
#'      var_names = c("Age","Education"),
#'      var_labs = list(c("18-30","31-50","51-65","65+"),
#'                      c("LTHS","HS","Some Coll","Coll+")),
#'      is_categorical = c(T,T))
#' @rdname plot.bpd
#' @export

plot.bpd <- function(bpd,
                     var_names = NULL,
                     var_labs = NULL,
                     is_categorical = NULL,
                     cols = c('#f2f0f7','#cbc9e2','#9e9ac8','#6a51a3'),...)
{
  #while (!is.null(dev.list()))  dev.off()
  pd.plot <- bpd$summary
  
  vars <- colnames(pd.plot)[which(!(colnames(pd.plot) %in% c("pred","prob","lb","ub")))]
  if(is.null(var_names)) {
    var_names <- vars
  }
  
  levs <- list()
  for(v in 1:length(vars)) {
    levs[[v]] <- unique(pd.plot[,vars[v]])
  }
  if(is.null(var_labs)) {
    var_labs <- levs
  }
  
  ylab_name = "Partial Effect"
  
  if(is.null(is_categorical)) {
    for(v in vars) {
      if(class(pd.plot[[v]]) == "factor" || (all(as.integer(pd.plot[[v]]) == pd.plot[[v]]) && length(unique(pd.plot[[v]])) <= 50)) {
        is_categorical <- c(is_categorical,T)
      } else {
        is_categorical <- c(is_categorical,F)
      }
    }
  }
  
  for(v in 1:length(vars)) {
    pd.plot[[vars[v]]] <- factor(pd.plot[[vars[v]]],levs[[v]],labels = var_labs[[v]])
  }
  
  if(length(vars) == 1) {
    if(is_categorical[1]) {
      ggplot(pd.plot,aes(y = pred,x = get(vars[1]))) +
      geom_point() + 
      geom_errorbar(aes(ymin = lb,ymax = ub),width = .2,col = rgb(0,0,0,.4)) + 
      theme_classic() + 
      labs(title = "Predicted Values and Credible Intervals",
           y = "Predicted Value",
           x = var_names[1])
    } else {
      ggplot(pd.plot,aes(y = pred,x = as.numeric(get(vars[1])))) +
        geom_line() + 
        geom_ribbon(aes(ymin = lb,ymax = ub),alpha = .2) + 
        theme_classic() + 
        labs(title = "Predicted Values and Credible Intervals",
             subtitle = paste0(var_names,collapse = " by "),
             y = "Predicted Value",
             x = var_names[1])
    }
  } else if(length(vars) == 2) {
    if(all(is_categorical)) {
      ggplot(pd.plot,aes(y = pred,x = get(vars[2]),color = factor(get(vars[1])))) +
        geom_point() + 
        scale_color_manual(values = cols) + 
        geom_errorbar(aes(ymin = lb,ymax = ub),width = .2,col = rgb(0,0,0,.4)) + 
        theme_classic() + 
        labs(title = "Predicted Values and Credible Intervals",
             subtitle = paste0(var_names,collapse = " by "),
             y = "Predicted Value",
             x = var_names[2]) +
        guides(color=guide_legend(title = var_names[1]))
    } else if(any(is_categorical)) {
      cat_ind <- which(is_categorical)
      cont_ind <- which(!is_categorical)

      ggplot(pd.plot,aes(y = pred,x = as.numeric(get(vars[cont_ind])))) +
        geom_line() + 
        geom_ribbon(aes(ymin = lb,ymax = ub),alpha = .2) + 
        facet_wrap(~ get(vars[cat_ind])) + 
        theme_classic() + 
        labs(title = "Predicted Values and Credible Intervals",
             subtitle = paste0(var_names,collapse = " by "),
             y = "Predicted Value",
             x = var_names[cont_ind]) +
        theme(legend.position="bottom")
    } else {
      ggplot(pd.plot, aes(get(vars[1]), get(vars[2]) )) +
        geom_tile(aes(fill = pred), color = "white") +
        scale_fill_gradient(low = "white", high = "steelblue") +
        ylab(var_names[2]) +
        xlab(var_names[1]) +
        labs(fill = "Predicted Value") +
        theme_minimal() +
        theme(plot.title = element_text(size=16),
              axis.title=element_text(size=14,face="bold"),
              axis.text.x = element_text(angle = 90, hjust = 1),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 9)) +
        guides(fill = guide_colorbar(title.position = "top",title.hjust = .5)) +
        ggtitle(paste0("Partial Dependence Plot:\n",paste(var_names,collapse = " by "))) +
        geom_text(aes(label = paste0(round(pred,2),"\n(",round(lb,2),", ",round(ub,2),")")))
    } 
  } else if(length(vars) == 3) {
    ncols <- length(unique(pd.plot[[vars[3]]]))
    ggplot(pd.plot, aes(get(vars[1]), get(vars[2]) )) +
      geom_tile(aes(fill = pred), color = "white") +
      scale_fill_gradient(low = "white", high = "steelblue") +
      ylab(var_names[2]) +
      xlab(var_names[1]) +
      labs(fill = "Predicted Value") +
      theme_minimal() +
      theme(plot.title = element_text(size=16),
            axis.title=element_text(size=14,face="bold"),
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 9)) +
      guides(fill = guide_colorbar(title.position = "top",title.hjust = .5)) +
      facet_wrap(~ get(vars[3]),ncol = floor(sqrt(ncols))) +
      ggtitle(paste0("Partial Dependence Plot:\n",paste(var_names,collapse = " by ")))
  } else {
    cat("More than 3-way interaction. Not plottable by default methods.")
  }
}
