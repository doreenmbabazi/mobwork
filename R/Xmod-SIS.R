#' Basic plotting: plot the density of infected humans for an SIS model
#'
#' @inheritParams plot_X
#' @export
plot_X.SIS = function(model, clrs="black", stable=FALSE, add=FALSE){
  vars = if(stable==TRUE){
    model$outputs$stable_orbits
  }else{
    model$outputs$orbits
  }

  if(add==FALSE) with(vars$XH,
      plot(time, 0*time, type = "n", ylim = c(0, max(H)),
           ylab = "# Infected", xlab = "Time"))

  lines_X(vars$XH, model, clrs)
}

#' Add the orbits for the SIS model to a plot for models of human infection and immunity
#'
#' @inheritParams lines_X
#'
#' @export
lines_X.SIS = function(XH, model, clrs="black"){
  with(XH,{
    if(model$nStrata==1) lines(time, X, col=clrs)
    if(model$nStrata>1){
      if (length(clrs)==1) clrs=rep(clrs, model$nStrata)
      for(i in 1:model$nStrata){
        lines(time, X[,i], col=clrs[i])
      }
    }
})}

#' Split a stratum into two strata, assigning a fraction `p` a new biting weight that is multiplied by a factor `fac` compared with the old one. The biting weight for the remaining `1-p` gets a new factor `1/fac`
#'
#' @inheritParams split_stratum_by_biting
#'
#' @return pars a list
#' @return i [integer] -- the stratum to split
#' @return p [numeric] -- the fraction that gets multiplied by `fac`
#' @return fac [numeric] -- the factor
#' @export
split_stratum_by_biting.SIS = function(pars, i, p, fac){
  stopifnot(i <= pars$nStrata)
  pars$nStrata = pars$nStrata + 1

  pars$Xpar$b <- c(pars$Xpar$b, pars$Xpar$b[i])
  pars$Xpar$r <- c(pars$Xpar$r, pars$Xpar$r[i])
  pars$Xpar$c <- c(pars$Xpar$c, pars$Xpar$c[i])

  pars$Xinits$X0 <- c(pars$Xinits$X0, pars$Xinits$X0[i])
  pars$Xinits$X0[i] <- pars$Xinits$X0[i]*p
  pars$Xinits$X0[pars$nStrata] <- pars$Xinits$X0[i]*(1-p)

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
