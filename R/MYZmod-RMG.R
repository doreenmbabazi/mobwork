# specialized methods for the adult mosquito RMG model

#' @title Reset bloodfeeding and mortality rates to baseline
#' @description Implements [MBionomics] for the RMG model
#' @inheritParams exDE::MBionomics
#' @return a named [list]
#' @export
MBionomics.RMG <- function(t, y, pars, s) {
  with(pars$MYZpar[[s]],{
    pars$MYZpar[[s]]$f <- f0
    pars$MYZpar[[s]]$q <- q0
    pars$MYZpar[[s]]$g <- g0
    pars$MYZpar[[s]]$phi <- 1/phi0
    pars$MYZpar[[s]]$sigma <- sigma0
    pars$MYZpar[[s]]$nu <- nu0

  return(pars)
})}

#' @title Blood feeding rate of the infective mosquito population
#' @description Implements [F_fqZ] for the RMG model.
#' @inheritParams exDE::F_fqZ
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqZ.RMG <- function(t, y, pars, s) {
  fqZ = with(pars$MYZpar[[s]], f*q)*y[pars$ix$MYZ[[s]]$Z_ix]
  return(fqZ)
}

#' @title Blood feeding rate of the infective mosquito population
#' @description Implements [F_fqM] for the RMG model.
#' @inheritParams exDE::F_fqM
#' @return a [numeric] vector of length `nPatches`
#' @export
F_fqM.RMG <- function(t, y, pars, s) {
  M = with(pars$ix$MYZ[[s]], y[U_ix] + y[Y_ix] + y[Z_ix])
  fqM = with(pars$MYZpar[[s]], f*q)*M
  return(fqM)
}

#' @title Number of eggs laid by adult mosquitoes
#' @description Implements [F_eggs] for the RMG model.
#' @inheritParams exDE::F_eggs
#' @return a [numeric] vector of length `nPatches`
#' @export
F_eggs.RMG <- function(t, y, pars, s) {

  G <- with(pars$ix$MYZ[[s]], y[Gu_ix] + y[Gy_ix] + y[Gz_ix])
  eggs = with(pars$MYZpar[[s]], {
    return(G*nu*eggsPerBatch)
  })
  return(eggs)
}

#' @title Derivatives for adult mosquitoes
#' @description Implements [dMYZdt] for the RMG ODE model.
#' @inheritParams exDE::dMYZdt
#' @return a [numeric] vector
#' @export
dMYZdt.RMG <- function(t, y, pars, s){

  Lambda = pars$Lambda[[s]]
  kappa = pars$kappa[[s]]

  with(pars$ix$MYZ[[s]],{
    U <- y[U_ix]
    Gu <- y[Gu_ix]
    Y <- y[Y_ix]
    Gy <- y[Gy_ix]
    Z <- y[Z_ix]
    Gz <- z[Gz_ix]

    with(pars$MYZpar[[s]],{

      Omega <- make_Omega(g, sigma, calK, nPatches)

      dUdt <- Lambda + nu*Gu - f*U - (Omega %*% U)
      dGudt <- f*(1-q*kappa)*U - nu*Gu - (Omega %*% Gu)
      dYdt <- nu*Gy - f*Y - phi*Y - (Omega %*% Y)
      dGydt <- f*q*kappa*U + f*Y - nu*Gy - phi*Gy - (Omega %*% Gy)
      dZdt <- phi*Y + nu*Gz - f*Z - (Omega %*% Z)
      dGzdt <- phi*Gy + f*Z - nu*Gz - (Omega %*% Gz)

      return(c(dUdt, dGudt, dYdt, dGydt, dZdt, dGzdt))
    })
  })
}

#' @title Setup MYZpar for the RMG model
#' @description Implements [setup_MYZpar] for the RM model
#' @inheritParams exDE::setup_MYZpar
#' @return a [list] vector
#' @export
setup_MYZpar.RMG = function(MYZname, pars, s, MYZopts=list(), EIPmod, calK){
  pars$MYZpar[[s]] = make_MYZpar_RMG(pars$nPatches, MYZopts, EIPmod, calK)
  return(pars)
}

#' @title Make parameters for RM ODE adult mosquito model
#' @param nPatches is the number of patches, an integer
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param EIPmod a [list] that defines the EIP model
#' @param calK a mosquito dispersal matrix of dimensions `nPatches` by `nPatches`
#' @param g mosquito mortality rate
#' @param sigma emigration rate
#' @param f feeding rate
#' @param q human blood fraction
#' @param nu oviposition rate, per mosquito
#' @param eggsPerBatch eggs laid per oviposition
#' @return a [list]
#' @export
make_MYZpar_RMG = function(nPatches, MYZopts=list(), EIPmod, calK,
                          g=1/12, sigma=1/8, f=0.3, q=0.95,
                          nu=1, eggsPerBatch=60){

  stopifnot(is.matrix(calK))
  stopifnot(dim(calK) == c(nPatches, nPatches))

  with(MYZopts,{
    MYZpar <- list()
    MYZpar$xde <- "ode"
    class(MYZpar$xde) <- "ode"
    class(MYZpar) <- "RMG"

    MYZpar$nPatches <- nPatches

    MYZpar$g      <- checkIt(g, nPatches)
    MYZpar$sigma  <- checkIt(sigma, nPatches)
    MYZpar$f      <- checkIt(f, nPatches)
    MYZpar$q      <- checkIt(q, nPatches)
    MYZpar$nu     <- checkIt(nu, nPatches)
    MYZpar$eggsPerBatch <- eggsPerBatch

    # Store as baseline values
    MYZpar$g0      <- MYZpar$g
    MYZpar$sigma0  <- MYZpar$sigma
    MYZpar$f0      <- MYZpar$f
    MYZpar$q0      <- MYZpar$q
    MYZpar$nu0     <- MYZpar$nu

    # The EIP model and the eip
    MYZpar$EIPmod <- EIPmod
    MYZpar$eip <- EIP(0, EIPmod)

    MYZpar$calK <- calK

    MYZpar$Omega <- make_Omega(g, sigma, calK, nPatches)
    MYZpar$Upsilon <- with(MYZpar, expm::expm(-Omega*eip))

    return(MYZpar)
})}

#' @title Setup initial values for the RMG model
#' @description Implements [setup_MYZinits] for the RM model
#' @inheritParams exDE::setup_MYZinits
#' @return a [list]
#' @export
setup_MYZinits.RMG = function(pars, s, MYZopts=list()){
  pars$MYZinits[[s]] = with(pars$MYZpar[[s]], make_MYZinits_RMG(nPatches, Upsilon, MYZopts))
  return(pars)
}

#' @title Make inits for RMG adult mosquito model
#' @param nPatches the number of patches in the model
#' @param MYZopts a [list] of values that overwrites the defaults
#' @param U0 total uninfected, not gravid mosquito density at each patch
#' @param Gu0 total uninfected, gravid mosquito density at each patch
#' @param Y0 infected, not gravid mosquito density at each patch
#' @param Gy0 infected, gravid mosquito density at each patch
#' @param Z0 infectious, not gravid mosquito density at each patch
#' @param Gz0 infectious, gravid mosquito density at each patch
#' @return a [list]
#' @export
make_MYZinits_RMG = function(nPatches, MYZopts = list(),
                             U0=5, Gu0=1, Y0=1, Gy0=1, Z0=1, Gz0=1){
  with(MYZopts,{
    U = checkIt(U0, nPatches)
    Gu = checkIt(Gu0, nPatches)
    Y = checkIt(Y0, nPatches)
    Gy = checkIt(Gy0, nPatches)
    Z = checkIt(Z0, nPatches)
    Gz = checkIt(Gz0, nPatches)

    return(list(U=U, Gu=Gu, Y=Y, Gy=Gy, Z=Z, Gz=Gz))
  })}

#' @title Add indices for adult mosquitoes to parameter list
#' @description Implements [make_indices_MYZ] for the RMG model.
#' @inheritParams exDE::make_indices_MYZ
#' @return none
#' @importFrom utils tail
#' @export
make_indices_MYZ.RMG <- function(pars, s) {with(pars,{

  U_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(U_ix, 1)

  Gu_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Gu_ix, 1)

  Y_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Y_ix, 1)

  Gy_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Gy_ix, 1)

  Z_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Z_ix, 1)

  Gz_ix <- seq(from = max_ix+1, length.out=nPatches)
  max_ix <- tail(Gz_ix, 1)

  pars$max_ix = max_ix
  pars$ix$MYZ[[s]] = list(U_ix=U_ix, Gu_ix=Gu_ix,
                          Y_ix=Y_ix, Gy_ix=Gy_ix,
                          Z_ix=Z_ix, Gz_ix=Gz_ix)
  return(pars)
})}

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

   nPatches = pars$nPatches

  MYZpar <- list()
  MYZpar$xde = 'ode'
  class(MYZpar$xde) <- 'ode'
  class(MYZpar) <- "RMG"

  MYZpar$g      <- checkIt(g, nPatches)
  MYZpar$sigma  <- checkIt(sigma, nPatches)
  MYZpar$f      <- checkIt(f, nPatches)
  MYZpar$q      <- checkIt(q, nPatches)
  MYZpar$nu     <- checkIt(nu, nPatches)
  MYZpar$eggsPerBatch <- eggsPerBatch

  # Store as baseline values
  MYZpar$g0      <- MYZpar$g
  MYZpar$sigma0  <- MYZpar$sigma
  MYZpar$f0      <- MYZpar$f
  MYZpar$q0      <- MYZpar$q
  MYZpar$nu0     <- MYZpar$nu

  pars$MYZpar[[1]] <- MYZpar

  return(pars)
}

#' @title Parse the output of deSolve and return variables for the RMG model
#' @description Implements [parse_deout_MYZ] for the RMG model
#' @inheritParams exDE::parse_deout_MYZ
#' @return none
#' @export
parse_deout_MYZ.RMG <- function(deout, pars, s) {
  time = deout[,1]
  with(pars$ix$MYZ[[s]],{
    U = deout[,U_ix+1]
    Gu = deout[,Gu_ix+1]
    Y = deout[,Y_ix+1]
    Gy = deout[,Gy_ix+1]
    Z = deout[,Z_ix+1]
    Gz = deout[,Gz_ix+1]
    M =  U+Y+Z+Gu+Gy+Gz
    y = (Y + Gy)/M
    z = (Z + Gz)/M
  return(list(time=time, U=U, Gu=Gu, Y=Y, Gy=Gy, Z=Z, Gz=Gz, M=M, y=y, z=z))
})}

#' @title Make inits for RMG adult mosquito model
#' @inheritParams exDE::update_inits_MYZ
#' @return none
#' @export
update_inits_MYZ.RMG <- function(pars, y0, s) {
  with(pars$ix$MYZ[[s]],{
    U = y[U_ix]
    Gu = y[Gu_ix]
    Y = y[Y_ix]
    Gy = y[Gy_ix]
    Z = y[Z_ix]
    Gz = y[Gz_ix]
    pars = make_MYZinits_RMG(pars$nPatches, U0=U, Gu0=Gu, Y0=Y, Gy0=Gy, Z0=Z, Gz0=Gz)
    return(pars)
})}

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
  pars$MYZinits[[1]] = list(U=U0, Gu=Gu0, Y=Y0, Gy=Gy0, Z=Z0, Gz=Gz0)
  return(pars)
}

#' @title Return initial values as a vector
#' @description Implements [get_inits_MYZ] for the RMG model.
#' @inheritParams exDE::get_inits_MYZ
#' @return none
#' @export
get_inits_MYZ.RMG <- function(pars, s) {with(pars$MYZinits[[s]],{
  c(U, Gu, Y, Gy, Z, Gz)
})}

