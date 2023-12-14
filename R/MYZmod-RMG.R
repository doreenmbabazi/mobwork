# specialized methods for the adult mosquito RMG model

#' @title Reset bloodfeeding and mortality rates to baseline
#' @description Implements [MBionomics] for the RMG model
#' @inheritParams exDE::MBionomics
#' @return a named [list]
#' @export
MBionomics.RMG <- function(t, y, pars) {

  pars$MYZpar$f <- pars$MYZpar$f0
  pars$MYZpar$q <- pars$MYZpar$q0
  pars$MYZpar$g <- pars$MYZpar$g0
  pars$MYZpar$phi <- 1/pars$MYZpar$phi0
  pars$MYZpar$sigma <- pars$MYZpar$sigma0
  pars$MYZpar$nu <- pars$MYZpar$nu0

  return(pars)
}

#' @title Blood feeding rate of the infective mosquito population
#' @description Implements [F_fqZ] for the RMG model.
#' @inheritParams exDE::F_fqZ
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqZ.RMG <- function(t, y, pars) {
  with(pars$MYZpar, f*q)*y[pars$MYZpar$Z_ix]
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the RMG model.
#' @inheritParams exDE::F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs.RMG <- function(t, y, pars) {
  with(pars$MYZpar, {
    G <- y[Gu_ix] + y[Gy_ix] + y[Gz_ix]
    return(G*nu*eggsPerBatch)
  })
}

#' @title Derivatives for adult mosquitoes
#' @description Implements [dMYZdt] for the RMG ODE model.
#' @inheritParams exDE::dMYZdt
#' @return a [numeric] vector
#' @export
dMYZdt.RMG <- function(t, y, pars, Lambda, kappa) {

  nPatches <- pars$nPatches

  with(pars$MYZpar,{

    U <- y[U_ix]
    Gu <- y[Gu_ix]
    Y <- y[Y_ix]
    Gy <- y[Gy_ix]
    Z <- y[Z_ix]
    Gz <- y[Gz_ix]

    Omega <- make_Omega(g, sigma, calK, nPatches)

    dUdt <- Lambda + nu*Gu - f*U - (Omega %*% U)
    dGudt <- f*(1-q*kappa)*U - nu*Gu - (Omega %*% Gu)
    dYdt <- nu*Gy - f*Y - phi*Y - (Omega %*% Y)
    dGydt <- f*q*kappa*U + f*Y - nu*Gy - phi*Gy - (Omega %*% Gy)
    dZdt <- phi*Y + nu*Gz - f*Z - (Omega %*% Z)
    dGzdt <- phi*Gy + f*Z - nu*Gz - (Omega %*% Gz)

    return(c(dUdt, dGudt, dYdt, dGydt, dZdt, dGzdt))
  })
}

#' @title Setup the RMG model
#' @description Implements [setup_MYZ] for the RMG model
#' @inheritParams exDE::setup_MYZ
#' @return a [list] vector
#' @export
setup_MYZ.RMG = function(pars, MYZname, nPatches=1, MYZopts=list(), calK=diag(1)){

  pars$MYZname = "RMG"
  pars$nPatches = nPatches

  pars = make_MYZpar_RMG(pars, MYZopts, calK)
  pars = make_MYZinits_RMG(pars, MYZopts)
  Omega <- with(pars$MYZpar, make_Omega(g, sigma, calK, nPatches))
  return(pars)
}

#' @title Make parameters for RMG ODE adult mosquito model
#' @param pars a [list]
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param calK mosquito dispersal matrix of dimensions `nPatches` by `nPatches`
#' @param solve_as is either `ode` to solve as an ode or `dde` to solve as a dde
#' @param g mosquito mortality rate
#' @param sigma emigration rate
#' @param f feeding rate
#' @param q human blood fraction
#' @param nu oviposition rate, per mosquito
#' @param eggsPerBatch eggs laid per oviposition
#' @param eip length of extrinsic incubation period
#' @return a [list]
#' @export
make_MYZpar_RMG = function(pars, MYZopts=list(), calK,
                           solve_as = "dde",
                           g=1/12, sigma=1/8,
                           f=0.5, q=0.95, eip=11,
                           nu=1, eggsPerBatch=60){

  stopifnot(is.matrix(calK))
  stopifnot(dim(calK) == c(pars$nPatches, pars$nPatches))

  with(MYZopts,{
    MYZpar <- list()

    MYZpar$xde <- solve_as
    class(MYZpar$xde) <- solve_as
    if(solve_as == 'dde') class(MYZpar) <- c('RMG', 'RMG_dde')
    if(solve_as == 'ode') class(MYZpar) <- c('RMG', 'RMG')

    MYZpar$g0      <- checkIt(g, pars$nPatches)
    MYZpar$sigma0  <- checkIt(sigma, pars$nPatches)
    MYZpar$f0      <- checkIt(f, pars$nPatches)
    MYZpar$q0      <- checkIt(q, pars$nPatches)
    MYZpar$nu0     <- checkIt(nu, pars$nPatches)
    MYZpar$eggsPerBatch <- eggsPerBatch
    MYZpar$phi0 <- 1/eip
    MYZpar$calK <- calK

    pars$MYZpar = MYZpar
    pars = MBionomics.RMG(0, 0, pars)

    return(pars)
  })}

#' @title Make inits for RMG adult mosquito model
#' @param pars a [list]
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param U0 total uninfected, not gravid mosquito density at each patch
#' @param Gu0 total uninfected, gravid mosquito density at each patch
#' @param Y0 infected, not gravid mosquito density at each patch
#' @param Gy0 infected, gravid mosquito density at each patch
#' @param Z0 infectious, not gravid mosquito density at each patch
#' @param Gz0 infectious, gravid mosquito density at each patch
#' @return a [list]
#' @export
make_MYZinits_RMG = function(pars, MYZopts = list(),
                             U0=5, Gu0=1, Y0=1, Gy0=1, Z0=1, Gz0=1){
  with(MYZopts,{
    inits = list()
    inits$U0 = checkIt(U0, pars$nPatches)
    inits$Gu0 = checkIt(Gu0, pars$nPatches)
    inits$Y0 = checkIt(Y0, pars$nPatches)
    inits$Gy0 = checkIt(Gy0, pars$nPatches)
    inits$Z0 = checkIt(Z0, pars$nPatches)
    inits$Gz0 = checkIt(Gz0, pars$nPatches)

    pars$MYZinits = inits
    return(pars)
  })}

#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [make_indices_MYZ] for the RMG model.
#' @inheritParams exDE::make_indices_MYZ
#' @return none
#' @importFrom utils tail
#' @export
make_indices_MYZ.RMG <- function(pars) {

  pars$MYZpar$U_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$U_ix, 1)

  pars$MYZpar$Gu_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$Gu_ix, 1)

  pars$MYZpar$Y_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$Y_ix, 1)

  pars$MYZpar$Gy_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$Gy_ix, 1)

  pars$MYZpar$Z_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$Z_ix, 1)

  pars$MYZpar$Gz_ix <- seq(from = pars$max_ix+1, length.out = pars$nPatches)
  pars$max_ix <- tail(pars$MYZpar$Gz_ix, 1)

  return(pars)
}

#' @title Make parameters for RMG ODE adult mosquito model
#' @param pars a [list]
#' @param g mosquito mortality rate
#' @param sigma emigration rate
#' @param calK mosquito dispersal matrix of dimensions `nPatches` by `nPatches`
#' @param f feeding rate
#' @param q human blood fraction
#' @param nu oviposition rate, per mosquito
#' @param eggsPerBatch eggs laid per oviposition
#' @param eip length of extrinsic incubation period
#' @return none
#' @export
make_parameters_MYZ_RMG <- function(pars, g, sigma, f, q, nu, eggsPerBatch, eip, calK) {
  stopifnot(is.numeric(g), is.numeric(sigma), is.numeric(f), is.numeric(q), is.numeric(nu), is.numeric(eggsPerBatch))

  MYZpar <- list()
  MYZpar$xde = 'ode'
  class(MYZpar$xde) <- 'ode'

  MYZpar$g0 <- g
  MYZpar$sigma0 <- sigma
  MYZpar$f0 <- f
  MYZpar$q0 <- q
  MYZpar$nu0 <- nu
  MYZpar$eggsPerBatch <- eggsPerBatch
  MYZpar$eip <- eip
  MYZpar$calK <- calK

  pars$MYZpar <- MYZpar
  pars = MBionomics.RMG(0, 0, pars)
  return(pars)
}

#' @title Parse the output of deSolve and return variables for the RMG model
#' @description Implements [parse_deout_MYZ] for the RMG model
#' @inheritParams exDE::parse_deout_MYZ
#' @return none
#' @export
parse_deout_MYZ.RMG <- function(deout, pars) {
  time = deout[,1]
  U = deout[,pars$MYZpar$U_ix+1]
  Gu = deout[,pars$MYZpar$Gu_ix+1]
  Y = deout[,pars$MYZpar$Y_ix+1]
  Gy = deout[,pars$MYZpar$Gy_ix+1]
  Z = deout[,pars$MYZpar$Z_ix+1]
  Gz = deout[,pars$MYZpar$Gz_ix+1]
  M =  U+Y+Z+Gu+Gy+Gz
  y = (Y + Gy)/M
  z = (Z + Gz)/M
  return(list(time=time, U=U, Gu=Gu, Y=Y, Gy=Gy, Z=Z, Gz=Gz, M=M, y=y, z=z))
}

#' @title Make inits for RMG adult mosquito model
#' @param pars a [list]
#' @param y0 a vector of initial values
#' @return none
#' @export
update_inits_MYZ.RMG <- function(pars, y0) {
  U0 = y0[pars$MYZpar$U_ix]
  Gu0 = y0[pars$MYZpar$Gu_ix]
  Y0 = y0[pars$MYZpar$Y_ix]
  Gy0 = y0[pars$MYZpar$Gy_ix]
  Z0 = y0[pars$MYZpar$Z_ix]
  Gz0 = y0[pars$MYZpar$Gz_ix]
  pars = make_inits_MYZ_RMG(pars, U0, Gu0, Y0, Gy0, Z0, Gz0)
  return(pars)
}

#' @title Make inits for RMG adult mosquito model
#' @param pars a [list]
#' @param U0 total mosquito density at each patch
#' @param Gu0 total gravid uninfected mosquito density at each patch
#' @param Y0 infected mosquito density at each patch
#' @param Gy0 total gravid infected mosquito density at each patch
#' @param Z0 infectious mosquito density at each patch
#' @param Gz0 total gravid infectious mosquito density at each patch
#' @return none
#' @export
make_inits_MYZ_RMG <- function(pars, U0, Gu0, Y0, Gy0, Z0, Gz0) {
  pars$MYZinits = list(U0=U0, Gu0=Gu0, Y0=Y0, Gy0=Gy0, Z0=Z0, Gz0=Gz0)
  return(pars)
}

#' @title Return initial values as a vector
#' @description Implements [get_inits_MYZ] for the RMG model.
#' @inheritParams exDE::get_inits_MYZ
#' @return none
#' @export
get_inits_MYZ.RMG <- function(pars) {with(pars$MYZinits,{
  c(U0, Gu0, Y0, Gy0, Z0, Gz0)
})}

