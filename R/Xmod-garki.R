# specialized methods for the human Garki model

#' @title Size of effective infectious human population
#' @description Implements [F_X] for the Garki model.
#' @inheritParams exDE::F_X
#' @return a [numeric] vector of length `nStrata`
#' @export
F_X.garki <- function(t, y, pars){
  y1 <- y[pars$Xpar$y1_ix]
  X = with(pars$Xpar, y1)
  return(X)
}

#' @title Infection blocking pre-erythrocytic immunity
#' @description Implements [F_b] for the Garki model.
#' @inheritParams exDE::F_b
#' @return a [numeric] vector of length `nStrata`
#' @export
F_b.garki <- function(y, pars) {
  with(pars$Xpar, b)
}

#' @title Compute the "true" prevalence of infection / parasite rate
#' @description Implements [F_pr] for the Garki model.
#' @inheritParams exDE::F_pr
#' @return a [numeric] vector of length `nStrata`
#' @export
F_pr.garki <- function(varslist, pars) {
#  pr = with(pars$Xpar,with(varslist$XH, (q1*y1+q2*y2+q3*y3)/H))
  pr = with(pars$Xpar,with(varslist$XH, (y1+y2+y3)/H))
  return(pr)
}
#' @title Derivatives for human population
#' @description Implements [dXdt] for the Garki model.
#' @inheritParams exDE::dXdt
#' @return a [numeric] vector
#' @export
dXdt.garki = function(t, y, pars, FoI){
  x1 <- y[pars$ix$X$x1_ix]
  x2 <- y[pars$ix$X$x2_ix]
  x3 <- y[pars$ix$X$x3_ix]
  x4 <- y[pars$ix$X$x4_ix]
  y1 <- y[pars$ix$X$y1_ix]
  y2 <- y[pars$ix$X$y2_ix]
  y3 <- y[pars$ix$X$y3_ix]
  H  <- y[pars$H_ix]

  with(pars$Xpar,{
    R1 = FoI/(exp(FoI/r1) - 1)
    R2 = FoI/(exp(FoI/r2) - 1)

    #dx1 = Births(t, x1, pars) - FoI*x1 + R1*y2 + dHdt(t, x1, pars)
    dx1 = - FoI*x1 + R1*y2 + dHdt(t, x1, pars)
    dx2 = FoI*x1 - nu*x2 + dHdt(t, x2, pars)
    dy1 = nu*x2 - alpha1*y1  + dHdt(t, y1, pars)
    dy2 = alpha1*y1 - R1*y2 - alpha2*y2 + dHdt(t, y2, pars)
    dy3 = alpha2*y2 + nu*x4 - R2*y3 + dHdt(t, y3, pars)
    dx3 = R2*y3 - FoI*x3 + dHdt(t, x3, pars)
    dx4 = FoI*x3 - nu*x4 + dHdt(t, x4, pars)
    dH =  dHdt(t, H, pars)

    return(c(dx1, dx2, dy1, dy2, dy3, dx3, dx4, dH))
  })}

#' @title Setup Xpar.Garki
#' @description Implements [setup_X] for the Garki model
#' @inheritParams exDE::setup_X
#' @return a [list] vector
#' @export
setup_X.garki = function(pars, Xname, Xopts=list()){

  pars$Xname = "garki"
  pars = make_Xpar_garki(pars, Xopts)
  pars = make_Xinits_garki(pars, Xopts)

  return(pars)
}

#' @title Add indices for human population to parameter list
#' @description Implements [make_indices_X] for the Garki model.
#' @inheritParams exDE::make_indices_X
#' @return none
#' @importFrom utils tail
#' @export
make_indices_X.garki <- function(pars) {
  pars$ix$X$x1_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$ix$X$x1_ix, 1)

  pars$ix$X$x2_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$ix$X$x2_ix, 1)

  pars$ix$X$y1_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$ix$X$y1_ix, 1)

  pars$ix$X$y2_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$ix$X$y2_ix, 1)

  pars$ix$X$y3_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$ix$X$y3_ix, 1)

  pars$ix$X$x3_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$ix$X$x3_ix, 1)

  pars$ix$X$x4_ix <- seq(from = pars$max_ix+1, length.out = pars$nStrata)
  pars$max_ix <- tail(pars$ix$X$x4_ix, 1)
  return(pars)
}

#' @title Make parameters for Garki human model
#' @param pars an [list]
#' @param Xopts an [list]
#' @param b transmission probability (efficiency) from mosquito to human
#' @param r1 a [numeric] recovery rate for non-immunes
#' @param r2 a [numeric]recovery rate for partially immmunes
#' @param nu a [numeric]incubation period
#' @param alpha1 a [numeric] rate of losing infectivity
#' @param alpha2 a [numeric] rate of acquiring immunity
#' @param q1 a [numeric] detection of y1
#' @param q2 a [numeric] detection of y2
#' @param q3 a [numeric] detection of y3
#' @param mu a [numeric] the death rate
#' @return a [list]
#' @export
make_Xpar_garki = function(pars, Xopts=list(), b=0.55,
                           r1=.0023, r2=.023, nu=1/15,
                           alpha1=.002, alpha2=.00019,
                           q1=.7, q2=0.5, q3=0.3, mu=1/65/365){
  with(Xopts,{
    xde <- 'ode'
    class(xde) <- 'ode'

    garki = list()
    class(garki) <- "garki"
    garki$xde <- xde
    garki$b=checkIt(b, pars$nStrata)
    garki$r1=checkIt(r1, pars$nStrata)
    garki$r2=checkIt(r2, pars$nStrata)
    garki$nu=checkIt(nu, pars$nStrata)
    garki$alpha1=checkIt(alpha1, pars$nStrata)
    garki$alpha2=checkIt(alpha2, pars$nStrata)
    garki$mu=checkIt(mu, pars$nStrata)
    garki$q1=checkIt(q1, pars$nStrata)
    garki$q2=checkIt(q2, pars$nStrata)
    garki$q3=checkIt(q3, pars$nStrata)
    pars$Xpar = garki

    return(pars)
  })}

#' @title Make inits for Garki human model. Note that the variables should sum up to H, so the initial value of x1 is not set. The values are passed in the same order as they are presented in the original paper.
#' @param pars an [environment]
#' @param Xopts a [list] with values to override default values
#' @param x1 a [numeric] initial value for the variable x1
#' @param x2 a [numeric] initial value for the variable x2
#' @param y1 a [numeric] initial value for the variable y1
#' @param y2 a [numeric] initial value for the variable y2
#' @param y3 a [numeric] initial value for the variable y3
#' @param x3 a [numeric] initial value for the variable x3
#' @param x4 a [numeric] initial value for the variable x4
#' @return none
#' @export
make_Xinits_garki <- function(pars, Xopts = list(), x1=NULL, x2=0, y1=0, y2=0, y3=0, x3=0, x4=0) {
  stopifnot(is.numeric(x2))
  stopifnot(is.numeric(y1))
  stopifnot(is.numeric(y2))
  stopifnot(is.numeric(y3))
  stopifnot(is.numeric(x3))
  stopifnot(is.numeric(x4))
  if(is.null(x1)) x1 = pars$Hpar$H - x2 - y1 - y2 - y3 - x3 - x4
  stopifnot(x1>0)

  inits = list()
  inits$x1 = checkIt(x1, pars$nStrata)
  inits$x2 = checkIt(x2, pars$nStrata)
  inits$y1 = checkIt(y1, pars$nStrata)
  inits$y2 = checkIt(y2, pars$nStrata)
  inits$y3 = checkIt(y3, pars$nStrata)
  inits$x3 = checkIt(x3, pars$nStrata)
  inits$x4 = checkIt(x4, pars$nStrata)
  pars$Xinits <- inits
  return(pars)
}

#' @title Return initial values as a vector for the Garki model
#' @description This method dispatches on the type of `pars$Xpar`.
#' @param pars a named [list]
#' @return a named [list]
#' @export
get_inits_X.garki <- function(pars){
  with(pars$Xinits, c(x1, x2, y1, y2, y3, x3, x4))
}

#' @title Update Xinits for the Garki model
#' @param pars a [list]
#' @param y0 a vector of initial values
#' @return none
#' @export
update_inits_X.garki <- function(pars, y0){
  with(pars$ix$X,{
    x1 <- y0[x1_ix]
    x2 <- y0[x2_ix]
    x3 <- y0[x3_ix]
    x4 <- y0[x4_ix]
    y1 <- y0[y1_ix]
    y2 <- y0[y2_ix]
    y3 <- y0[y3_ix]
    pars = make_Xinits_garki(pars, x1=x1, x2=x2, y1=y1, y2=y2, y3=y3, x3=x3, x4=x4)
  return(pars)
})}

#' Plot the density of infected individuals for the Garki model
#'
#' @inheritParams exDE::xde_plot_X
#' @export
xde_plot_X.garki = function(pars, clrs=viridisLite::turbo(7), llty=1, stable=FALSE, add_axes=TRUE){
  vars=with(pars$outputs,if(stable==TRUE){stable_orbits}else{orbits})

  if(add_axes==TRUE)
    with(vars$XH,
         plot(time, 0*time, type = "n", ylim = c(0, max(H)),
              ylab = "# Infected", xlab = "Time"))

  xde_lines_X(vars$XH, pars, clrs, llty)
}

#' Add lines for the density of infected individuals for the Garki model
#'
#' @inheritParams exDE::xde_lines_X
#'
#' @export
xde_lines_X.garki= function(XH, pars, clrs=viridisLite::turbo(7), llty=1){
  with(XH,{
    if(pars$nStrata==1){
      lines(time, x1, col=clrs[1], lty = llty[1])
      lines(time, x2, col=clrs[2], lty = llty[1])
      lines(time, y1, col=clrs[3], lty = llty[1])
      lines(time, y2, col=clrs[4], lty = llty[1])
      lines(time, y3, col=clrs[5], lty = llty[1])
      lines(time, x3, col=clrs[6], lty = llty[1])
      lines(time, x4, col=clrs[7], lty = llty[1])
    }
    if(pars$nStrata>1){
      if (length(clrs)==1) clrs=rep(clrs, pars$nStrata)
      if (length(llty)==1) llty=rep(llty, pars$nStrata)
      for(i in 1:pars$nStrata){
        lines(time, x1[,i], col=clrs[1], lty = llty[i])
        lines(time, x2[,i], col=clrs[2], lty = llty[i])
        lines(time, y1[,i], col=clrs[3], lty = llty[i])
        lines(time, y2[,i], col=clrs[4], lty = llty[i])
        lines(time, y3[,i], col=clrs[5], lty = llty[i])
        lines(time, x3[,i], col=clrs[6], lty = llty[i])
        lines(time, x4[,i], col=clrs[7], lty = llty[i])
      }
    }
  })}

#' @title Parse the output of deSolve and return variables for the Garki model
#' @description Implements [parse_deout_X] for the Garki model
#' @inheritParams exDE::parse_deout_X
#' @return none
#' @export
parse_deout_X.garki <- function(deout, pars) {
  time = deout[,1]
  Hlist <- parse_deout_H(deout, pars)
  with(Hlist,{
    x1 = deout[,pars$ix$X$x1_ix+1]
    x2 = deout[,pars$ix$X$x2_ix+1]
    y1 = deout[,pars$ix$X$y1_ix+1]
    y2 = deout[,pars$ix$X$y2_ix+1]
    y3 = deout[,pars$ix$X$y3_ix+1]
    x3 = deout[,pars$ix$X$x3_ix+1]
    x4 = deout[,pars$ix$X$x4_ix+1]
    return(list(time=time, x1=x1, x2=x2, y1=y1, y2=y2, y3=y3, x3=x3, x4=x4, H=H))
  })}

#' @title Compute the HTC for the garki model
#' @description Implements [HTC] for the garki model
#' @inheritParams exDE::HTC
#' @return a [numeric] vector
#' @export
HTC.garki <- function(pars) {
  with(pars$Xpar,
       return(1/r1)
  )
}
