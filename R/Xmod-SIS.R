#' Split a stratum into two strata, assigning a fraction `p` a new biting weight that is multiplied by a factor `fac` compared with the old one. The biting weight for the remaining `1-p` gets a new factor `1/fac`
#'
#' @inheritParams split_stratum_by_biting
#'
#' @return pars a list
#' @param i the stratum to split
#' @param p the fraction that gets multiplied by `fac`
#' @param fac a factor
#' @return a [list]
#' @export
split_stratum_by_biting.SIS = function(pars, i, p, fac){
  stopifnot(i <= pars$nStrata)
  pars$nStrata = pars$nStrata + 1

  Xpar = pars$Xpar[[1]]
  Xpar$b <- c(Xpar$b, Xpar$b[i])
  Xpar$r <- c(Xpar$r, Xpar$r[i])
  Xpar$c <- c(Xpar$c, Xpar$c[i])

  Xinits = pars$inits[[1]]
  Xinits$X0 <- c(Xinits$X0, Xinits$X0[i])
  Xinits$X0[i] <- Xinits$X0[i]*p
  Xinits$X0[pars$nStrata] <- Xinits$X0[i]*(1-p)

  pars$Hpar$residence = c(pars$Hpar$residence, pars$Hpar$residence[i])
  pars$Hpar$wts_f = c(pars$Hpar$wts_f, pars$Hpar$wts_f[i])
  pars$Hpar$wts_f[pars$nStrata] = pars$Hpar$wts_f[i]*fac
  pars$Hpar$wts_f[i] = pars$Hpar$wts_f[i]/fac
  pars$Hpar$rbr = with(pars$Hpar, wts_f*sum(H)/sum(wts_f*H))
  pars$Hpar$H = c(pars$Hpar$H, pars$Hpar$H[i])
  pars$Hpar$H[pars$nStrata] = pars$Hpar$H[i]*p
  pars$Hpar$H[i] = pars$Hpar$H[i]*(1-p)
  pars$Hpar$TaR = cbind(pars$Hpar$TaR, pars$Hpar$TaR[,i])

  pars <- exDE::make_indices(pars)

  return(pars)
}
