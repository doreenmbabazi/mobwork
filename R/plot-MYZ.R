
#' Basic plotting: plot mosquito population density
#'
#' @param model a [list] specifying the model
#' @param clrs a [list] of colors
#' @param stable [logical]
#' @param add [logical]
#'
#' @export
plot_M = function(model, clrs, stable, add){
  UseMethod("plot_M", model$MYZpar)
}

#' Basic plotting: add to a plot for mosquito population density
#'
#' @param MYZ a [list] the outputs of exDE::parse_deout
#' @param model a [list] specifying the model
#' @param clrs a [list] of colors
#'
#' @export
lines_M = function(MYZ, model, clrs){
  UseMethod("lines_M", model$MYZpar)
}

#' Basic plotting: plot the density of infected and infective mosquitoes
#'
#' @param model a [list] specifying the model
#' @param Yclrs a [list] of colors
#' @param Zclrs a [list] of colors
#' @param stable [logical]
#' @param add [logical]
#'
#' @export
plot_YZ = function(model, Yclrs, Zclrs, stable, add){
  UseMethod("plot_YZ", model$MYZpar)
}

#' Basic plotting: add lines for the density of infected and infective mosquitoes
#'
#' @param MYZ a [list] the output of exDE::parse_deout
#' @param model a [list] specifying the model
#' @param Yclrs a [list] of colors
#' @param Zclrs a [list] of colors
#'
#' @export
lines_YZ = function(MYZ, model, Yclrs, Zclrs){
  UseMethod("lines_YZ", model$MYZpar)
}

#' Basic plotting: plot the proportion of mosquitoes that are infected and infective.
#'
#' @param model a [list] specifying the model
#' @param Yclrs a [list] of colors
#' @param Zclrs a [list] of colors
#' @param stable [logical]
#' @param add [logical]
#'
#' @export
plot_YZ_fracs = function(model, Yclrs, Zclrs, stable, add){
  UseMethod("plot_YZ_fracs", model$MYZpar)
}

#' Add lines for the fraction of mosquitoes infected and infective
#'
#' @param MYZ a [list] the output of exDE::parse_deout
#' @param model a [list] specifying the model
#' @param Yclrs a [list] of colors
#' @param Zclrs a [list] of colors
#'
#' @export
lines_YZ_fracs = function(MYZ, model, Yclrs, Zclrs){
  UseMethod("lines_YZ_fracs", model$MYZpar)
}
