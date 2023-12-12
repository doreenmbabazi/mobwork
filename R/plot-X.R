
#' Basic plotting: default plot for models of human infection and immunity
#'
#' @param model a [list] specifying the model
#' @param clrs a [list] of colors
#' @param stable [logical]
#' @param add [logical]
#'
#' @export
plot_X = function(model, clrs="black", stable=FALSE, add=FALSE){
  UseMethod("plot_X", model$Xpar)
}

#' Basic plotting: add outputs to a plot for models of human infection and immunity
#'
#' @param XH a [list] from exDE::parse_deout
#' @param model a [list] specifying the model
#' @param clrs a [list] of colors
#'
#' @export
lines_X = function(XH, model, clrs="black"){
  UseMethod("lines_X", model$Xpar)
}

