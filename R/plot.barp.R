#' plot.barp
#' 
#' This function plots the predictions by geographic unit \code{geo.unit}. 
#' @param barp.obj A \code{\link{barp}} object.
#' @param evaluate_model Plot model convergence diagnostics? If left to the default \code{FALSE}, return a plot of predictions and credible intervals.
#' @param algorithm Select which algorithm's predictions should be plotted. If left to the default \code{"all"}, plot all predictions without credible intervals.
#' @keywords MRP BARP Mister P multilevel regression poststratification
#' @seealso \code{\link{barp}} which creates the prerequisite \code{barp} object.
#' @examples 
#' data("gaymar")
#' barp.obj <- barp(y = "supp_gaymar",
#'              x = c("pvote","religcon","age","educ","gXr","stateid","region"),
#'              dat = svy,
#'              census = census06,
#'              geo.unit = "stateid",
#'              proportion = "n")
#' plot(barp.obj,
#'      evaluate_model = F,
#'      algorithm = "BARP")
#' @rdname plot.barp
#' @export

plot.barp <- function(barp.obj,
                      evaluate_model = F,
                      algorithm = "all")
{
  if(evaluate_model) {
    if(is.null(barp.obj$trees)) {
      stop("barp object not built with bartMachine. Diagnostics unavailable.")
    }
    plot_convergence_diagnostics(barp.obj$trees)
  } else {
    geo.unit <- colnames(barp.obj$pred.opn)[1]
    if(barp.obj$BSSD) {
      range <- max(barp.obj$pred.opn$opn.ub) - min(barp.obj$pred.opn$opn.lb)
      pad <- c(min(barp.obj$pred.opn$opn.lb)-(range*.1),
               max(barp.obj$pred.opn$opn.ub)+(range*.1))
      ggplot(barp.obj$pred.opn,aes(y = pred.opn,x = get(geo.unit))) +
        geom_point(aes(reorder(get(geo.unit),pred.opn))) + 
        geom_text(aes(reorder(get(geo.unit),pred.opn),label = get(geo.unit),y = opn.lb),hjust = 1.5,size = 2.75) + 
        geom_errorbar(aes(reorder(get(geo.unit),pred.opn),ymin = opn.lb,ymax = opn.ub),width = .2,col = rgb(0,0,0,.4)) + 
        coord_flip() + 
        ylim(pad) + 
        expand_limits(x = 0) + 
        theme_classic() + 
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank()) + 
        labs(title = "Predicted Values and Credible Intervals",
             subtitle = paste0("Bootstrapped Confidence Intervals Across Ensemble Predictions"),
             y = "Predicted Value") 
      } else if(algorithm == "all") {
        methods <- colnames(barp.obj$pred.opn)
        methods <- methods[-which(methods %in% c(geo.unit,"opn.lb","opn.ub"))]
        # And only BARP
        if(length(methods) == 1) {
          if(grepl("bartMachine",methods)) {
            algorithm <- methods[which(grepl("bartMachine",methods))][1]
            range <- max(barp.obj$pred.opn$opn.ub) - min(barp.obj$pred.opn$opn.lb)
            pad <- c(min(barp.obj$pred.opn$opn.lb)-(range*.1),
                     max(barp.obj$pred.opn$opn.ub)+(range*.1))
            ggplot(barp.obj$pred.opn,aes(y = get(algorithm),x = get(geo.unit))) +
              geom_point(aes(reorder(get(geo.unit),get(algorithm)))) + 
              geom_text(aes(reorder(get(geo.unit),get(algorithm)),label = get(geo.unit),y = opn.lb),hjust = 1.5,size = 2.75) + 
              geom_errorbar(aes(reorder(get(geo.unit),get(algorithm)),ymin = opn.lb,ymax = opn.ub),width = .2,col = rgb(0,0,0,.4)) + 
              coord_flip() + 
              ylim(pad) + 
              expand_limits(x = 0) + 
              theme_classic() + 
              theme(axis.title.y = element_blank(),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.line.y = element_blank()) + 
              labs(title = "Predicted Values and Credible Intervals",
                   subtitle = paste0("Algorithm: ",algorithm),
                   y = "Predicted Value") 
          } else {
            # And only one method (not BARP)
            algorithm <- methods[1]
            range <- max(barp.obj$pred.opn[[algorithm]]) - min(barp.obj$pred.opn[[algorithm]])
            pad <- c(min(barp.obj$pred.opn[[algorithm]])-(range*.1),
                     max(barp.obj$pred.opn[[algorithm]])+(range*.1))
            ggplot(barp.obj$pred.opn,aes(y = get(algorithm),x = get(geo.unit))) +
              geom_point(aes(reorder(get(geo.unit),get(algorithm)))) + 
              geom_text(aes(reorder(get(geo.unit),get(algorithm)),label = get(geo.unit),y = get(algorithm)),hjust = 1.5,size = 2.75) + 
              # geom_errorbar(aes(reorder(get(geo.unit),get(algorithm)),ymin = opn.lb,ymax = opn.ub),width = .2,col = rgb(0,0,0,.4)) + 
              coord_flip() + 
              ylim(pad) + 
              expand_limits(x = 0) + 
              theme_classic() + 
              theme(axis.title.y = element_blank(),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.line.y = element_blank()) + 
              labs(title = "Predicted Values and Credible Intervals",
                   subtitle = paste0("Algorithm: ",algorithm),
                   y = "Predicted Value")
          }
        } else {
          # Multiple methods
          toplot <- barp.obj$pred.opn %>% dplyr::select(methods,geo.unit) %>%
            dplyr::mutate(max.ind = get(methods[1])) %>% 
            tidyr::gather(algorithm,pred.opn,methods)
          range <- max(toplot$pred.opn) - min(toplot$pred.opn)
          pad <- c(min(toplot$pred.opn)-(range*.1),
                   max(toplot$pred.opn)+(range*.1))
          toplot %>%
            ggplot(aes(y = pred.opn,x = get(geo.unit),color = algorithm)) +
            geom_point(aes(reorder(get(geo.unit),pred.opn))) + 
            geom_text(aes(reorder(get(geo.unit),pred.opn),label = get(geo.unit),y = max.ind),hjust = 1.5,size = 2.75) + 
            # geom_errorbar(aes(reorder(get(geo.unit),pred.opn),ymin = opn.lb,ymax = opn.ub),width = .2,col = rgb(0,0,0,.4)) + 
            coord_flip() + 
            ylim(pad) + 
            expand_limits(x = 0) + 
            theme_classic() + 
            theme(axis.title.y = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.line.y = element_blank()) + 
            labs(title = "Predicted Values by Algorithm",
                 subtitle = "All Algorithms & Ensemble",
                 y = "Predicted Value")
        }
    } else {
      # A specific chosen method
      algorithm <- paste0("SL.",gsub("SL.","",gsub("barp|BARP","bartMachine",algorithm)))
      methods <- colnames(barp.obj$pred.opn)
      methods <- methods[-which(methods %in% c(geo.unit,"opn.lb","opn.ub"))]
      algorithm <- methods[which(grepl(algorithm,methods))]
      if(grepl("bartMachine",algorithm)) {
        # If it is BARP, add error bars
        range <- max(barp.obj$pred.opn$opn.ub) - min(barp.obj$pred.opn$opn.lb)
        pad <- c(min(barp.obj$pred.opn$opn.lb)-(range*.1),
                 max(barp.obj$pred.opn$opn.ub)+(range*.1))
        ggplot(barp.obj$pred.opn,aes(y = get(algorithm),x = get(geo.unit))) +
          geom_point(aes(reorder(get(geo.unit),get(algorithm)))) + 
          geom_text(aes(reorder(get(geo.unit),get(algorithm)),label = get(geo.unit),y = opn.lb),hjust = 1.5,size = 2.75) + 
          geom_errorbar(aes(reorder(get(geo.unit),get(algorithm)),ymin = opn.lb,ymax = opn.ub),width = .2,col = rgb(0,0,0,.4)) + 
          coord_flip() + 
          ylim(pad) + 
          expand_limits(x = 0) + 
          theme_classic() + 
          theme(axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.y = element_blank()) + 
          labs(title = "Predicted Values and Credible Intervals",
               subtitle = paste0("Algorithm: ",algorithm),
               y = "Predicted Value") 
      } else {
        # Otherwise plot without
        range <- max(barp.obj$pred.opn[[algorithm]]) - min(barp.obj$pred.opn[[algorithm]])
        pad <- c(min(barp.obj$pred.opn[[algorithm]])-(range*.1),
                 max(barp.obj$pred.opn[[algorithm]])+(range*.1))
        ggplot(barp.obj$pred.opn,aes(y = get(algorithm),x = get(geo.unit))) +
          geom_point(aes(reorder(get(geo.unit),get(algorithm)))) + 
          geom_text(aes(reorder(get(geo.unit),get(algorithm)),label = get(geo.unit),y = get(algorithm)),hjust = 1.5,size = 2.75) + 
          # geom_errorbar(aes(reorder(get(geo.unit),get(algorithm)),ymin = opn.lb,ymax = opn.ub),width = .2,col = rgb(0,0,0,.4)) + 
          ylim(pad) + 
          expand_limits(x = 0) + 
          coord_flip() + 
          theme_classic() + 
          theme(axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.y = element_blank()) + 
          labs(title = "Predicted Values and Credible Intervals",
               subtitle = paste0("Algorithm: ",algorithm),
               y = "Predicted Value")
      }
    }
  }
}
