## ----eval=F-------------------------------------------------------------------
#  dXdt.SIS <- function(t, y, pars, FoI) {
#    with(pars$Xpar, {
#      X <- y[X_ix]
#      H <- F_H(t, y, pars)
#  
#      dX <- FoI*(H - X) - r*X
#  
#      return(c(dX))
#    })
#  }

## ----eval=F-------------------------------------------------------------------
#  dXdt.SISdH <- function(t, y, pars, FoI) {
#    with(pars$Xpar, {
#  
#      H <- F_H(t, y, pars)
#      X <- y[X_ix]
#  
#      dX <- FoI*(H - X) - r*X + dHdt(t, X, pars)
#      dH <- Births(t, H, pars) + dHdt(t, H, pars)
#  
#      return(c(dX, dH))
#    })
#  }

## ----eval=F-------------------------------------------------------------------
#  setup_X.SIS = function(pars, Xname, Xopts=list()){
#    pars$Xname = "SIS"
#    pars = make_Xpar_SIS(pars, Xopts)
#    pars = make_Xinits_SIS(pars, Xopts)
#  
#    return(pars)
#  }

## ----eval=F-------------------------------------------------------------------
#  make_Xpar_SIS = function(pars, Xopts=list(),
#                           b=0.55, r=1/180, c=0.15){
#    with(Xopts,{
#      Xpar = list()
#      class(Xpar) <- c("SIS", "SISdX")
#  
#      Xpar$b = checkIt(b, pars$nStrata)
#      Xpar$c = checkIt(c, pars$nStrata)
#      Xpar$r = checkIt(r, pars$nStrata)
#  
#      pars$Xpar = Xpar
#      return(pars)
#  })}

## ----eval=F-------------------------------------------------------------------
#  make_Xinits_SIS = function(pars, Xopts = list(), X0=1){with(Xopts,{
#    inits = list()
#    inits$X0 = checkIt(X0, pars$nStrata)
#    pars$Xinits = inits
#    return(pars)
#  })}
#  

## ----eval=F-------------------------------------------------------------------
#  #' @title Parse the output of deSolve and return variables for the SIS model
#  #' @description Implements [parse_deout_X] for the SIS model
#  #' @inheritParams parse_deout_X
#  #' @return none
#  #' @export
#  parse_deout_X.SIS <- function(deout, pars) {
#    time = deout[,1]
#    Hlist <- parse_deout_H(deout, pars)
#    with(Hlist,{
#      X = deout[,pars$Xpar$X_ix+1]
#      return(list(time=time, X=X, H=H))
#  })}

## ----eval=F-------------------------------------------------------------------
#  #' @title Compute the HTC for the SIS model
#  #' @description Implements [HTC] for the SIS model with demography.
#  #' @inheritParams HTC
#  #' @return a [numeric] vector
#  #' @export
#  HTC.SIS <- function(pars) {
#    with(pars$Xpar,
#      return(c/r)
#    )
#  }

## -----------------------------------------------------------------------------
#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SIS model.
#' @inheritParams make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.SIS <- function(pars) {
  pars$Xpar$X_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$Xpar$X_ix, 1)
  return(pars)
}

## ----eval=F-------------------------------------------------------------------
#  #' @title Update inits for the SIS human model from a vector of states
#  #' @param pars a [list]
#  #' @param y0 a vector of initial values
#  #' @return none
#  #' @export
#  update_inits_X.SIS <- function(pars, y0) {
#    X0 = y0[pars$Xpar$X_ix]
#    pars = make_inits_X_SIS(pars, X0)
#    return(pars)
#  }

## ----eval=F-------------------------------------------------------------------
#  #' @title Return initial values as a vector
#  #' @description This method dispatches on the type of `pars$Xpar`.
#  #' @param pars a [list]
#  #' @return none
#  #' @export
#  get_inits_X.SIS <- function(pars){
#    pars$Xinits$X0
#  }

## ----eval=F-------------------------------------------------------------------
#  #' Plot the density of infected individuals for the SIS model
#  #'
#  #' @inheritParams xde_plot_X
#  #' @export
#  xde_plot_X.SIS = function(pars, clrs="black", llty=1, stable=FALSE, add_axes=TRUE){
#    vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})
#  
#    if(add_axes==TRUE)
#      with(vars$XH,
#           plot(time, 0*time, type = "n", ylim = c(0, max(H)),
#                ylab = "# Infected", xlab = "Time"))
#  
#    xde_lines_X(vars$XH, pars, clrs, llty)
#  }

## ----eval=F-------------------------------------------------------------------
#  #' Add lines for the density of infected individuals for the SIS model
#  #'
#  #' @inheritParams xde_lines_X
#  #'
#  #' @export
#  xde_lines_X.SIS = function(XH, pars, clrs="black", llty=1){
#    with(XH,{
#      if(pars$nStrata==1) lines(time, X, col=clrs[1], lty = llty[1])
#      if(pars$nStrata>1){
#        if (length(clrs)==1) clrs=rep(clrs, pars$nStrata)
#        if (length(llty)==1) llty=rep(llty, pars$nStrata)
#        for(i in 1:pars$nStrata){
#          lines(time, X[,i], col=clrs[i], lty = llty[i])
#        }
#      }
#    })}

