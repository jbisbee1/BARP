#' plot.brp
#' 
#' This function plots the predictions by geographic unit \code{geo.unit}. 
#' @param barp A \code{\link{barp}} object.
#' @param evaluate_model Plot model convergence diagnostics? If left to the default \code{FALSE}, return a plot of predictions and credible intervals.
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
#'      evaluate_model = F)
#' @rdname plot.barp
#' @export

plot.brp <- function(barp,
                      evaluate_model = F)
{
  if(evaluate_model) {
    plot_convergence_diagnostics(barp$trees)
  } else {
    geo.unit <- colnames(barp$pred.opn)[1]
    ggplot(barp$pred.opn,aes(y = pred.opn,x = get(geo.unit))) +
      geom_point(aes(reorder(get(geo.unit),pred.opn))) + 
      geom_text(aes(reorder(get(geo.unit),pred.opn),label = get(geo.unit),y = opn.lb),hjust = 1.5) + 
      geom_errorbar(aes(reorder(get(geo.unit),pred.opn),ymin = opn.lb,ymax = opn.ub),width = .2,col = rgb(0,0,0,.4)) + 
      coord_flip() + 
      theme_classic() + 
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank()) + 
      labs(title = "Predicted Values and Credible Intervals",
           subtitle = "Sorted by Value",
           x = "Predicted Value")
  }
}
