
#' Split a stratum into two strata, assigning
#'
#' @param pars a [list] defining a model
#' @param i the host species index
#' @param j the stratum to split
#' @param p the fraction in the higher exposure stratum
#' @param fac the factor increase
#'
#' @return pars a list
#' @export
split_stratum_by_biting = function(pars, i, j, p, fac){
  UseMethod("split_stratum_by_biting", pars$Xpar[[i]])
}


