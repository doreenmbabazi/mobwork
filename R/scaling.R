#' @title Compute eir-pr scaling relationships
#' @description This function computes average annual values for the eir, the pr, and other
#' interesting terms and returns a table. It is computed for a model of class "cohort"
#' @param model a [list]
#' @param N an [integer]
#' @export
xde_scaling_eir = function(model, N=25){
  stopifnot(model$nPatches == 1)
  stopifnot(class(model$xde) == "cohort")

  # make sure that eir & ni are output
  model <- exDE::make_indices(model)

  # Get the
  F_eir_base = model$F_eir
  scl_a = with(model$Hpar,sum(wts_f*H)/sum(H))
  scl_b = stats::integrate(F_eir_base, 0, 365, pars=model)$value
  model$scl = 1/scl_a/scl_b
  model$F_eir = function(t, model){with(model,aeir*scl*F_eir_base(t, model))}

  eir = 10^seq(-1,3, length.out=N)
  pr = 0*eir
  ni = 0*eir
  scaling = list()

  for(i in 1:N){
    model$aeir = eir[i]
    xde_tmp <- exDE::xde_stable_orbit(model)
    tmp <- xde_tmp$outputs$stable_orbit
    H = tmp$XH$H
    tot_pr <- rowSums(tmp$terms$pr*H)/rowSums(H)
    wts_f = model$Hpar$wts_f
    mean_ni <- with(tmp$terms, rowSums(ni*wts_f*H)/rowSums(wts_f*H))
    scaling[[i]] = with(tmp$terms, list(aeir=365*eir, eir=eir, pr=tot_pr, ni=mean_ni, pr_t = pr, ni_t = ni))
    pr[i] = mean(tot_pr)
    ni[i] = mean(mean_ni)
  }

  model$outputs$eirpr <- list(aeir=eir, eir=eir/365, pr=pr, ni=ni, scaling=scaling)
  model$F_eir = F_eir_base

  return(model)
}

#' @title Compute eir-pr scaling relationships
#' @description This function computes average annual values for the eir, the pr, and other
#' interesting terms and returns a table. It is computed for a model of class "human"
#' @param model a [list]
#' @param N an [integer], the size of the mesh on aEIR
#' @export
xde_scaling_Z = function(model, N=25){
  stopifnot(model$nPatches == 1)
  stopifnot(class(model$xde) == "human")

  model <- exDE::make_indices(model)
  y0 = exDE::get_inits(model)

  # Get the
  F_Z_base = model$MYZpar$Zf
  scl_a = with(model$Hpar,sum(wts_f*H)/sum(H))
  beta = with(model$Hpar, compute_beta(H, wts_f, TaR))
  scl_b = stats::integrate(model$F_EIR, 0, 365, pars=model, y=y0, beta=beta)$value
  model$MYZpar$scl = 1/scl_a/scl_b
  model$MYZpar$Zf = function(t, pars){with(pars$MYZpar,aeir*scl*F_Z_base(t, pars$MYZpar))}

  eir = 10^seq(-1,3, length.out=N)
  pr = 0*eir
  ni = 0*eir
  scaling = list()

  for(i in 1:N){
    model$MYZpar$aeir = eir[i]
    xde_tmp <- exDE::xde_stable_orbit(model)
    tmp <- xde_tmp$outputs$stable_orbit
    H = tmp$XH$H
    tot_pr <- rowSums(tmp$terms$pr*H)/rowSums(H)
    wts_f = model$Hpar$wts_f
    mean_ni <- with(tmp$terms, rowSums(ni*wts_f*H)/rowSums(wts_f*H))
    scaling[[i]] = with(tmp$terms, list(aeir=365*eir, eir=eir, pr=tot_pr, ni=mean_ni, pr_t = pr, ni_t = ni))
    pr[i] = mean(tot_pr)
    ni[i] = mean(mean_ni)
  }

  model$outputs$eirpr <- list(aeir=eir, eir=eir/365, pr=pr, ni=ni, scaling=scaling)
  model$MYZpar$Zf = F_Z_base

  return(model)
}

#' @title Compute lambda from an eirpr object using the Ross-Macdonald model
#' @description This function computes `m` and `lambda` for the output of one
#' of `xde_scaling_eir` or `xde_scaling_Z`.  The outputs are attached to eirpr
#' @param model a [list]
#' @export
xde_scaling_lambda = function(model){with(model,with(outputs,{
  stopifnot(exists("MYZss"))
  N = length(eirpr$eir)
  Z=Y=M=lambda=rep(0,N)
  for(i in 1:N){
    Z[i] = with(MYZss, eirpr$eir[i]/f/q)
    Y[i] = with(MYZss, exp(g*eip)*Z[i])
    M[i] = with(MYZss, f*q*eirpr$ni[i]/(f*q*eirpr$ni[i] + g)*Y[i])
    lambda[i] = M[i]/MYZss$g
  }
  model$outputs$eirpr$m=M
  model$outputs$eirpr$lambda=lambda
  return(model)
}))}


#' @title Set up the MYZss object for `xde_scaling_lambda`
#' @description This function computes several quantities that are require
#' @param model a [list]
#' @export
ssMYZ = function(model){with(model$MYZpar,{
  MYZss = list()
  Omega = make_Omega(g, sigma, calK, model$nPatches)
  MYZss$Omega = Omega
  MYZss$OmegaInv = solve(Omega)
  MYZss$Upsilon = expm(-Omega*eip)
  MYZss$UpsilonInv = expm(Omega*eip)
  beta = with(model$Hpar, compute_beta(H, wts_f, TaR))
  MYZss$beta = beta
  MYZss$betaInv = solve(beta)
  MYZss$f = f
  MYZss$q = q
  MYZss$g = g
  MYZss$eip = eip
  model$MYZss = MYZss
  return(model)
})}

#' @title Construct an eirpr object for an arbitary model
#' @description This takes a model and uses the XH component to define
#' the eirpr relationship using `xde_scaling_eir` then calls `xde_scaling_lambda`
#' @param model a [list]
#' @param N an [integer], the size of the mesh on aEIR
#' @param F_eir a [function], the size of the mesh on aEIR
#' @export
xde_scaling = function(model, N=25, F_eir=NULL){
  if(is.null(F_eir)) F_eir <- function(t, pars){0*t + 1}
  mod0 = exDE::xde_setup_cohort("mod0", F_eir)
  mod0$Xpar <- model$Xpar
  mod0$Hpar <- model$Hpar
  model$outputs$eirpr = xde_scaling_eir(mod0, N)$outputs$eirpr
  model <- ssMYZ(model)
  model <- pr2Lambda(model)
  return(model)
}

#' @title Using the eirpr matrix and a RM model, convert pr to Lambda
#' @description This takes a model and uses the XH component to define
#' the eirpr relationship using `xde_scaling_eir` then calls `xde_scaling_lambda`
#' @param pr a [numeric] vector
#' @param model a [list]
#' @param constrain a [logical]
#' @export
pr2Lambda = function(pr, model, constrain=TRUE){
  with(model, stopifnot(exists("MYZss")))
  with(model$MYZss,{
    eir = xde_pr2eir(pr, model, TRUE)$eir
    kappa = xde_pr2ni(pr, model, TRUE)$ni
    Z = (betaInv %*% eir)/f/q
    Y = OmegaInv %*% (UpsilonInv %*% (Omega %*% Z))
    M = diag(1/f/q/kappa, model$nPatches)%*%(diag(f*q*kappa, model$nPatches) + Omega) %*% Y
    Lambda = Omega %*% M
    if(constrain == TRUE) Lambda = pmax(Lambda,0)
    return(Lambda)
  })
}
