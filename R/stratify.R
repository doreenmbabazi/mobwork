
#' Split a stratum into two strata, assigning
#'
#' @param pars a [list] defining a model
#' @param i the stratum to split
#' @param p the fraction in the higher exposure stratum
#' @param fac the factor increase
#'
#' @return pars a list
#' @export
split_stratum_by_biting = function(pars, i, p, fac){
  UseMethod("split_stratum_by_biting", pars$Xpar)
}


