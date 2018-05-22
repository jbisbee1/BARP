#' plot.bpd
#' 
#' This function plots the partial dependence results estimated returned by running barp_partial_dependence.
#' @param bpd A \code{BARP} partial dependence object.
#' @param var_names A vector of variable names. If \code{NULL} (the default), the variable names from the training data columns are used.
#' @param var_labs A list of variable labels corresponding to the levels at which the partial dependence was calculated. If \code{NULL} (the default), the raw values are used.
#' @param is_categorical A vector of logicals (T,F) corresponding to which variables are categorical and which are continuous.If \code{NULL} (the default), the function determines whether each variable is categorical or continuous based on the data class.
#' @param ... Additional arguments to be passed to plot.
#' @keywords MRP BARP Mister P multilevel regression poststratification
#' @seealso \code{\link{barp_partial_dependence}} which is a prerequisite to plotting.
#' @examples 
#' data("gaymar")
#' BARP <- barp(y = "supp_gaymar",
#'              x = c("pvote","religcon","age","educ","gXr","stateid","region"),
#'              dat = svy,
#'              census = census06,
#'              geo.unit = "stateid",
#'              proportion = "n")
#' bpd <- barp_partial_dependence(BARP = BARP,
#'                                vars = c("age","educ"),
#'                                prop_data = .2,
#'                                levs = list(c(1:4),c(1:4),
#'                                credible_interval = c(0.0)))
#' plot.bpd(bpd,
#'      var_names = c("Age","Education"),
#'      var_labs = list(c("18-30","31-50","51-65","65+"),
#'                      c("LTHS","HS","Some Coll","Coll+")),
#'      is_categorical = c(T,T))
#' @export barp

plot.bpd <- function(bpd,
                     var_names = NULL,
                     var_labs = NULL,
                     is_categorical = NULL,
                     cols = c('#f2f0f7','#cbc9e2','#9e9ac8','#6a51a3'),...)
{
  while (!is.null(dev.list()))  dev.off()
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
    pd.plot[[vars[v]]] <- factor(pd.plot[[vars[v]]],rev(levs[[v]]),labels = var_labs[[v]])
  }
  
  min.1 <- min(pd.plot$lb,pd.plot$ub)
  max.1 <- max(pd.plot$lb,pd.plot$ub)
  if(length(vars) == 1) {
    plot(1:length(levs[[1]]),rep(0,length(levs[[1]])),main = paste0("Partial Dependence Plot:\n",paste(var_names,collapse = " by ")),
         ylim = c(min.1,max.1),xlim = c(1,length(levs[[1]])),type = 'n',las = 1,xaxt = 'n',bty = 'n',xlab = var_names[1],ylab = ylab_name)
    axis(1, at = levs[[1]],labels = var_labs[[1]])
    if(is_categorical[1]) {
      segments(levs[[1]],pd.plot$lb,levs[[1]],pd.plot$ub)
      points(pd.plot$pred,pch = 21,bg = "white")
    } else {
      polygon(c(levs[[1]], rev(levs[[1]])), c(pd.plot$ub, rev(pd.plot$lb)), col = rgb(0,0,0,.2), border = NA)
      lines(levs[[1]], pd.plot$lb, type = "o", col = "blue", lty = 2)
      lines(levs[[1]], pd.plot$ub, type = "o", col = "blue", lty = 2)
      lines(levs[[1]], pd.plot$pred, type = "o", lwd = 2)
    }
  } else if(length(vars) == 2) {
    if(all(is_categorical)) {
      par(oma = c(0,0,0,4))
      plot(1:length(levs[[1]]),rep(0,length(levs[[1]])),main = paste0("Partial Dependence Plot:\n",paste(var_names,collapse = " by ")),
           ylim = c(min.1,max.1),xlim = c(1,length(levs[[1]])),type = 'n',las = 1,xaxt = 'n',bty = 'n',xlab = var_names[1],ylab = ylab_name)
      axis(1, at = levs[[1]],labels = var_labs[[1]])
      sapply(levs[[1]],function(x) segments(levs[[2]],pd.plot[which(as.numeric(pd.plot[[vars[1]]]) == x),"lb"],levs[[2]],pd.plot[which(as.numeric(pd.plot[[vars[1]]]) == x),"ub"],col = cols[x]))
      sapply(levs[[1]],function(x) points(pd.plot[which(as.numeric(pd.plot[[vars[1]]]) == x),"pred"],pch = 21,bg = cols[x]))
      legend(par('usr')[ 2 ]*1.025, mean(par('usr')[c(3,4)]),legend = var_labs[[2]],pch = 21,pt.bg = cols,xpd = NA,cex = .8,title = var_names[2])
    } else if(any(is_categorical)) {
      cat_ind <- which(is_categorical)
      cont_ind <- which(!is_categorical)
      cats <- length(levs[[cat_ind]])
      rows <- ceiling(sqrt(cats))
      colms <- floor(sqrt(cats))
      
      mat.layout <- matrix(NA,nrow = rows + 2,ncol = colms + 1)
      mat.layout[1,2:ncol(mat.layout)] <- 2
      mat.layout[nrow(mat.layout),2:ncol(mat.layout)] <- cats+3
      mat.layout[2:(rows+1),2:(colms+1)] <- matrix(3:(cats+2),nrow = rows,ncol = colms,byrow = T)
      mat.layout[,1] <- 1
      nf <- layout(mat.layout,heights = c(1,rep(5,rows),1),widths = c(1,rep(5,colms)))
      # layout.show(nf)
      par(mar = c(0,0,0,0))
      plot(0,0,type = 'n',xaxt = 'n',yaxt = 'n',bty = 'n',ylab = "",xlab = "")
      mtext(text = ylab_name,side = 2,srt = 90,cex = 1.3,line = -2)
      plot(0,0,type = 'n',xaxt = 'n',yaxt = 'n',bty = 'n',xlab = "",ylab = "")
      mtext(text = paste0("Partial Dependence Plot:\n",var_names[1]," X ",var_names[2]),side = 3,cex = 1.3,line = -4)
      
      par(mar = c(2,2,2,2))
      for(j in levs[[cat_ind]]) {
        plot(pd.plot[which(as.numeric(pd.plot[[vars[cat_ind]]]) == j),"pred"],type = 'l',
             ylim = c(min.1,max.1),las = 1,xlab = "",ylab = "",yaxt = 'n',xaxt = 'n')
        text(par('usr')[ 2 ]*.985, par('usr')[4],adj = c(1.05,1.5),
             labels = paste0(var_names[cont_ind],":\n ",var_labs[[cont_ind]][j]),cex = .9)
        if(j == 1) {
          axis(2,las = 2,tck = -.02,cex.axis = .8)
        }
        if(j == cats) {
          axis(1,tck = -.02,cex.axis = .8)
        }
        
        polygon(c(levs[[cont_ind]], rev(levs[[cont_ind]])),
                c(pd.plot[which(as.numeric(pd.plot[[vars[cat_ind]]]) == j),"ub"],
                  rev(pd.plot[which(as.numeric(pd.plot[[vars[cat_ind]]]) == j),"lb"])), col = rgb(0,0,0,.2), border = NA)
        lines(levs[[cont_ind]], pd.plot[which(as.numeric(pd.plot[[vars[cat_ind]]]) == j),"ub"], col = "blue", lty = 2)
        lines(levs[[cont_ind]], pd.plot[which(as.numeric(pd.plot[[vars[cat_ind]]]) == j),"lb"], col = "blue", lty = 2)
      }
      par(mar = c(0,0,0,0))
      plot(0,0,type = 'n',xaxt = 'n',yaxt = 'n',bty = 'n',xlab = "",ylab = "")
      mtext(text = var_names[cat_ind],side = 1,cex = 1.3,line = -2)
    } else {
      heat.dat <- pd.plot
      ggplot(heat.dat, aes(get(vars[1]), get(vars[2]) )) +
        geom_tile(aes(fill = pred), color = "white") +
        scale_fill_gradient(low = "white", high = "steelblue") +
        ylab(var_names[2]) +
        xlab(var_names[1]) +
        theme(legend.title = element_text(size = 10),
              legend.text = element_text(size = 12),
              plot.title = element_text(size=16),
              axis.title=element_text(size=14,face="bold"),
              axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(fill = "Predicted Support") +
        theme_minimal() +
        ggtitle(paste0("Partial Dependence Plot:\n",paste(var_names,collapse = " by "))) +
        geom_text(aes(label = paste0(round(pred,2),"\n(",round(lb,2),", ",round(ub,2),")")))
    } 
  } else if(length(vars) == 3) {
    heat.dat <- pd.plot
    ggplot(heat.dat, aes(get(vars[1]), get(vars[2]) )) +
      geom_tile(aes(fill = pred), color = "white") +
      scale_fill_gradient(low = "white", high = "steelblue") +
      ylab(var_names[2]) +
      xlab(var_names[1]) +
      theme(legend.title = element_text(size = 10),
            legend.text = element_text(size = 12),
            plot.title = element_text(size=16),
            axis.title=element_text(size=14,face="bold"),
            axis.text.x = element_text(angle = 90, hjust = 1),
            strip.text.x = element_text(size=8, angle=75)) +
      labs(fill = "Predicted Support") +
      theme_minimal() +
      facet_wrap(~ get(vars[3]),ncol = 2) +
      ggtitle(paste0("Partial Dependence Plot:\n",paste(var_names,collapse = " by ")))
  } else {
    cat("More than 3-way interaction. Not plottable by default methods.")
  }
}
