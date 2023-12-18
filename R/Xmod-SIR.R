# specialized methods for the human SIR model

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SIR model.
#' @inheritParams exDE::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SIR <- function(t, y, pars) {
  X <- y[pars$ix$X$X_ix]
  with(pars$Xpar, X*c)
}

#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SIR model.
#' @inheritParams exDE::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SIR <- function(varslist, pars) {
  pr = with(varslist$XH, X/H)
  return(pr)
}


#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SIR model.
#' @inheritParams exDE::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SIR <- function(y, pars) {
  with(pars$Xpar, b)
}

#' @title Derivatives for human population
#' @description Implements [dXdt] for the SIR model, no demography.
#' @inheritParams exDE::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SIRdX <- function(t, y, pars, FoI) {

  X <- y[pars$ix$X$X_ix]
  R <- y[pars$ix$X$R_ix]
  H <- F_H(t, y, pars)

  with(pars$Xpar, {

    dX <- FoI*(H - X) - r*X
    dR <- r*X

    return(c(dX, dR))
  })
}

#' @title Derivatives for human population
#' @description Implements [dXdt] for the SIR model with demography.
#' @inheritParams exDE::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SIRdXdH <- function(t, y, pars, FoI) {

  X <- y[pars$ix$X$X_ix]
  R <- y[pars$ix$X$R_ix]
  H <- F_H(t, y, pars)

  with(pars$Xpar, {

    dX <- FoI*(H - X) - r*X + dHdt(t, X, pars)
    dR <- r*X + dHdt(t, X, pars)
    dH <- Births(t, H, pars) + dHdt(t, H, pars)

    return(c(dX, dR, dH))
  })
}

#' @title Setup Xpar.SIR
#'
#' @description Implements [setup_X] for the SIR model
#' @inheritParams exDE::setup_X
#' @return a [list] vector
#'
#' @export
setup_X.SIR = function(pars, Xname, Xopts=list()){

  pars$Xname = "SIR"
  pars = make_Xpar_SIR(pars, Xopts)
  pars = make_Xinits_SIR(pars, Xopts)

  return(pars)
}

#' @title Make parameters for SIR human model, with defaults
#'
#' @param pars a [list]
#' @param Xopts a [list] that could overwrite defaults
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @return a [list]
#'
#' @export
make_Xpar_SIR = function(pars, Xopts=list(),
                         b=0.55, r=1/180, c=0.15){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("SIR", "SIRdX")

    Xpar$b = checkIt(b, pars$nStrata)
    Xpar$c = checkIt(c, pars$nStrata)
    Xpar$r = checkIt(r, pars$nStrata)

    pars$Xpar = Xpar
    return(pars)
  })}

#' @title Make initial values for the SIR human model, with defaults
#' @param pars a [list]
#' @param Xopts a [list] to overwrite defaults
#' @param X0 the initial values of the parameter X
#' @param R0 the initial values of the parameter R
#' @return a [list]
#' @export
make_Xinits_SIR = function(pars, Xopts = list(), X0=1, R0=1){with(Xopts,{
  inits = list()
  inits$X0 = checkIt(X0, pars$nStrata)
  inits$R0 = checkIt(R0, pars$nStrata)
  pars$Xinits = inits
  return(pars)
})}

#' @title Parse the output of deSolve and return variables for the SIR model
#' @description Implements [parse_deout_X] for the SIR model
#' @inheritParams exDE::parse_deout_X
#' @return none
#' @export
parse_deout_X.SIR <- function(deout, pars) {
  time = deout[,1]
  Hlist <- parse_deout_H(deout, pars)
  with(Hlist,{
    X = deout[,pars$ix$X$X_ix+1]
    R = deout[,pars$ix$X$R_ix+1]
    return(list(time=time, X=X, R=R, H=H))
  })}

#' @title Compute the HTC for the SIR model
#' @description Implements [HTC] for the SIR model with demography.
#' @inheritParams exDE::HTC
#' @return a [numeric] vector
#' @export
HTC.SIR <- function(pars) {
  with(pars$Xpar,
       return(c/r)
  )
}

#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SIR model.
#' @inheritParams exDE::make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.SIR <- function(pars) {
  pars$ix$X$X_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$ix$X$X_ix, 1)

  pars$ix$X$R_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$ix$X$R_ix, 1)

  return(pars)
}

#' @title Update inits for the SIR human model from a vector of states
#' @param pars a [list]
#' @param y0 a vector of initial values
#' @return none
#' @export
update_inits_X.SIR <- function(pars, y0) {
  X0 = y0[pars$ix$X$X_ix]
  R0 = y0[pars$ix$X$R_ix]
  pars = make_Xinits_SIR(pars, list(), X0, R0)
  return(pars)
}


#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param pars a [list]
#' @return none
#' @export
get_inits_X.SIR <- function(pars){
  with(pars$Xinits, c(X0, R0))
}

#' Plot the density of infected individuals for the SIR model
#'
#' @inheritParams exDE::xde_plot_X
#' @export
xde_plot_X.SIR = function(pars, clrs=c("darkred", "darkblue"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})

  if(add_axes==TRUE)
    with(vars$XH,
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "# Infected", xlab = "Time"))

  xde_lines_X(vars$XH, pars, clrs, llty)
}


#' Add lines for the density of infected individuals for the SIR model
#'
#' @inheritParams exDE::xde_lines_X
#'
#' @export
#'
xde_lines_X.SIR = function(XH, pars, clrs=c("darkred", "darkblue"), llty=1){
  with(XH,{
    if(pars$nStrata==1) {
      lines(time, X, col=clrs[1], lty = llty[1])
      lines(time, R, col=clrs[2], lty = llty[1])
    }
    if(pars$nStrata>1){
      if (length(clrs)==1) clrs=matrix(clrs, 2, pars$nStrata)
      if (length(llty)==1) llty=rep(llty, pars$nStrata)
      for(i in 1:pars$nStrata){
        lines(time, X[,i], col=clrs[1,i], lty = llty[i])
        lines(time, R[,i], col=clrs[2,i], lty = llty[i])
      }
    }
  })}
