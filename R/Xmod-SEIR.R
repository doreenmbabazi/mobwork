# specialized methods for the human SEIR model

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the SEIR model.
#' @inheritParams exDE::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.SEIR <- function(t, y, pars) {
  I <- y[pars$ix$X$I_ix]
  with(pars$Xpar, I*c)
}

#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the SEIR model.
#' @inheritParams exDE::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.SEIR <- function(varslist, pars) {
  pr = with(varslist$XH, I/H)
  return(pr)
}


#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the SEIR model.
#' @inheritParams exDE::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.SEIR <- function(y, pars) {
  with(pars$Xpar, b)
}

#' @title Derivatives for human population
#' @description Implements [dXdt] for the SEIR model, no demography.
#' @inheritParams exDE::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SEIRdX <- function(t, y, pars, FoI) {

  E <- y[pars$ix$X$E_ix]
  I <- y[pars$ix$X$I_ix]
  R <- y[pars$ix$X$R_ix]
  H <- F_H(t, y, pars)

  with(pars$Xpar, {

    dE <- FoI*(H-E-I-R) - tau*E
    dI <- tau*E - r*I
    dR <- r*I

    return(c(dE, dI, dR))
  })
}

#' @title Derivatives for human population
#' @description Implements [dXdt] for the SEIR model with demography.
#' @inheritParams exDE::dXdt
#' @return a [numeric] vector
#' @export
dXdt.SEIRdXdH <- function(t, y, pars, FoI) {

  E <- y[pars$ix$X$E_ix]
  I <- y[pars$ix$X$I_ix]
  R <- y[pars$ix$X$R_ix]
  H <- F_H(t, y, pars)

  with(pars$Xpar, {

    dE <- FoI*(H-E-I-R) - tau*E + dHdt(t, E, pars)
    dI <- tau*E - r*I + dHdt(t, I, pars)
    dR <- r*I + dHdt(t, R, pars)
    dH <- Births(t, H, pars) + dHdt(t, H, pars)

    return(c(dE, dI, dR, dH))
  })
}

#' @title Setup Xpar.SEIR
#'
#' @description Implements [setup_X] for the SEIR model
#' @inheritParams exDE::setup_X
#' @return a [list] vector
#'
#' @export
setup_X.SEIR = function(pars, Xname, Xopts=list()){

  pars$Xname = "SEIR"
  pars = make_Xpar_SEIR(pars, Xopts)
  pars = make_Xinits_SEIR(pars, Xopts)

  return(pars)
}

#' @title Make parameters for SEIR human model, with defaults
#'
#' @param pars a [list]
#' @param Xopts a [list] that could overwrite defaults
#' @param b transmission probability (efficiency) from mosquito to human
#' @param c transmission probability (efficiency) from human to mosquito
#' @param r recovery rate
#' @param tau the latent period is 1/tau
#' @return a [list]
#'
#' @export
make_Xpar_SEIR = function(pars, Xopts=list(),
                         b=0.55, r=1/180, c=0.15, tau=1/5){
  with(Xopts,{
    Xpar = list()
    class(Xpar) <- c("SEIR", "SEIRdX")

    Xpar$b = checkIt(b, pars$nStrata)
    Xpar$c = checkIt(c, pars$nStrata)
    Xpar$r = checkIt(r, pars$nStrata)
    Xpar$tau = checkIt(tau, pars$nStrata)

    pars$Xpar = Xpar
    return(pars)
  })}

#' @title Make initial values for the SEIR human model, with defaults
#' @param pars a [list]
#' @param Xopts a [list] to overwrite defaults
#' @param E0 the initial values of the parameter E
#' @param I0 the initial values of the parameter I
#' @param R0 the initial values of the parameter R
#' @return a [list]
#' @export
make_Xinits_SEIR = function(pars, Xopts = list(), E0=1, I0=0, R0=1){with(Xopts,{
  inits = list()
  inits$E0 = checkEt(E0, pars$nStrata)
  inits$I0 = checkIt(I0, pars$nStrata)
  inits$R0 = checkIt(R0, pars$nStrata)
  pars$Xinits = inits
  return(pars)
})}

#' @title Parse the output of deSolve and return variables for the SEIR model
#' @description Implements [parse_deout_X] for the SEIR model
#' @inheritParams exDE::parse_deout_X
#' @return none
#' @export
parse_deout_X.SEIR <- function(deout, pars) {
  time = deout[,1]
  Hlist <- parse_deout_H(deout, pars)
  with(Hlist,{
    E = deout[,pars$ix$X$E_ix+1]
    I = deout[,pars$ix$X$I_ix+1]
    R = deout[,pars$ix$X$R_ix+1]
    return(list(time=time, E=E, I=I, R=R, H=H))
  })}

#' @title Compute the HTC for the SEIR model
#' @description Implements [HTC] for the SEIR model with demography.
#' @inheritParams exDE::HTC
#' @return a [numeric] vector
#' @export
HTC.SEIR <- function(pars) {
  with(pars$Xpar,
       return(c/r)
  )
}

#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the SEIR model.
#' @inheritParams exDE::make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.SEIR <- function(pars) {

  pars$ix$X$E_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$ix$X$E_ix, 1)

  pars$ix$X$I_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$ix$X$I_ix, 1)

  pars$ix$X$R_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$ix$X$R_ix, 1)

  return(pars)
}

#' @title Update inits for the SEIR human model from a vector of states
#' @param pars a [list]
#' @param y0 a vector of initial values
#' @return none
#' @export
update_inits_X.SEIR <- function(pars, y0) {
  E0 = y0[pars$ix$X$E_ix]
  I0 = y0[pars$ix$X$I_ix]
  R0 = y0[pars$ix$X$R_ix]
  pars = make_Xinits_SEIR(pars, list(), E0, I0, R0)
  return(pars)
}


#' @title Return initial values as a vector
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param pars a [list]
#' @return none
#' @export
get_inits_X.SEIR <- function(pars){
  with(pars$Xinits, c(I0, R0))
}

#' Plot the density of infected individuals for the SEIR model
#'
#' @inheritParams exDE::xde_plot_X
#' @export
xde_plot_X.SEIR = function(pars, clrs=c("darkgreen", "darkred", "darkblue"), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})

  if(add_axes==TRUE)
    with(vars$XH,
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "# Infected", xlab = "Time"))

  xde_lines_X(vars$XH, pars, clrs, llty)
}


#' Add lines for the density of infected individuals for the SEIR model
#'
#' @inheritParams exDE::xde_lines_X
#'
#' @export
#'
xde_lines_X.SEIR = function(XH, pars, clrs=c("darkgreen", "darkred", "darkblue"), llty=1){
  with(XH,{
    if(pars$nStrata==1) {
      lines(time, E, col=clrs[1], lty = llty[1])
      lines(time, I, col=clrs[2], lty = llty[1])
      lines(time, R, col=clrs[3], lty = llty[1])
    }
    if(pars$nStrata>1){
      if (length(clrs)==1) clrs=matrix(clrs, 2, pars$nStrata)
      if (length(llty)==1) llty=rep(llty, pars$nStrata)
      for(i in 1:pars$nStrata){
        lines(time, E[,i], col=clrs[1,i], lty = llty[i])
        lines(time, I[,i], col=clrs[2,i], lty = llty[i])
        lines(time, R[,i], col=clrs[3,i], lty = llty[i])
      }
    }
  })}
