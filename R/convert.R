
#' @title Convert pr to eir
#' @description Use the outputs$eirpr table to convert a set of pr values into eir values
#' @param pr a [vector]
#' @param model a [list]
#' @param extend a [logical] option to determine whether to extend outside the range
#' @export
xde_pr2eir = function(pr, model, extend=FALSE){
  with(model$outputs, stopifnot(exists("eirpr")))
  PR = model$outputs$eirpr$pr
  EIR = model$outputs$eirpr$eir
  if(extend==TRUE){
    PR = c(0, PR, 1)
    EIR = c(0, EIR, 5*10^3/365)
  }
  ix = which(pr<=max(PR) & pr>=min(PR))
  get_eir = function(pr){
    if(pr == min(PR)) return(min(EIR))
    if(pr == max(PR)) return(max(EIR))
    ix = max(which(PR<pr))
    ff = (pr-PR[ix])/(PR[ix+1]-PR[ix])
    eir=EIR[ix] + ff*(EIR[ix+1]-EIR[ix])
    return(eir)
  }
  eir = 0*pr-1
  eir[ix] = sapply(pr[ix], get_eir)
  pr2eir = list(pr=pr[ix], eir=eir[ix])
  if(length(ix)>0) pr2eir$errors = c(pr=pr[-ix])
  return(pr2eir)
}

#' @title Convert eir to pr
#' @description Use the outputs$eirpr table to interpolate
#' @param eir a [vector]
#' @param model a [list]
#' @param extend a [logical] option to determine whether to extend outside the range
#' @export
xde_eir2pr = function(eir, model, extend=FALSE){
  with(model$outputs, stopifnot(exists("eirpr")))
  PR = model$outputs$eirpr$pr
  EIR = model$outputs$eirpr$eir
  if(extend==TRUE){
    PR = c(0, PR, 1)
    EIR = c(0, EIR, 5*10^3/365)
  }
  ix = which(eir<=max(EIR) & eir>=min(EIR))
  get_pr = function(eir){
    if(eir == min(EIR)) return(min(PR))
    if(eir == max(EIR)) return(max(PR))
    ix = max(which(EIR<eir))
    ff = (eir-EIR[ix])/(EIR[ix+1]-EIR[ix])
    pr=PR[ix] + ff*(PR[ix+1]-PR[ix])
    return(pr)
  }
  pr = 0*eir-1
  pr[ix] = sapply(eir[ix], get_pr)
  eir2pr = list(eir=eir[ix], pr=pr[ix])
  if(length(ix)>0) eir2pr$errors = c(pr=pr[-ix])
  return(eir2pr)
}

#' @title Convert eir to ni
#' @description Use the outputs$eirpr table to interpolate
#' @param eir a [vector]
#' @param model a [list]
#' @param extend a [logical] option to determine whether to extend outside the range
#' @export
xde_eir2ni = function(eir, model, extend=FALSE){
  with(model$outputs, stopifnot(exists("eirpr")))
  NI = model$outputs$eirpr$ni
  EIR = model$outputs$eirpr$eir
  if(extend==TRUE){
    NI = c(0, NI, 1)
    EIR = c(0, EIR, 5*10^3/365)
  }
  ix = which(eir<=max(EIR) & eir>=min(EIR))
  get_ni = function(eir){
    if(eir == min(EIR)) return(min(NI))
    if(eir == max(EIR)) return(max(NI))
    ix = max(which(EIR<eir))
    ff = (eir-EIR[ix])/(EIR[ix+1]-EIR[ix])
    ni=NI[ix] + ff*(NI[ix+1]-NI[ix])
    return(ni)
  }
  ni = 0*eir-1
  ni[ix] = sapply(eir[ix], get_ni)
  eir2ni = list(eir=eir[ix], ni=ni[ix])
  if(length(ix)>0) eir2ni$errors = c(ni=ni[-ix])
  return(eir2ni)
}

#' @title Convert pr to ni
#' @description Use the outputs$eirpr table to interpolate
#' @param pr a [vector]
#' @param model a [list]
#' @param extend a [logical] option to determine whether to extend outside the range
#' @export
xde_pr2ni = function(pr, model, extend=FALSE){
  with(model$outputs, stopifnot(exists("eirpr")))
  PR = model$outputs$eirpr$pr
  NI = model$outputs$eirpr$ni
  if(extend==TRUE){
    PR = c(0, PR, 1)
    NI = c(0, NI, max(NI))
  }
  ix = which(pr<=max(PR) & pr>=min(PR))
  get_ni = function(pr){
    if(pr == min(PR)) return(min(NI))
    if(pr == max(PR)) return(max(NI))
    ix = max(which(PR<pr))
    ff = (pr-PR[ix])/(PR[ix+1]-PR[ix])
    ni=NI[ix] + ff*(NI[ix+1]-NI[ix])
    return(ni)
  }
  ni = 0*pr-1
  ni[ix] = sapply(pr[ix], get_ni)
  pr2ni = list(pr=pr[ix], ni=ni[ix])
  if(length(ix)>0) pr2ni$errors = c(pr=pr[-ix])
  return(pr2ni)
}


#' @title Convert pr to mosquito density
#' @description Use the outputs$eirpr table to convert a set of pr values into scaled mosquito density values
#' @param pr a [vector]
#' @param model a [list]
#' @param extend a [logical] option to determine whether to extend outside the range
#' @export
xde_pr2m = function(pr, model, extend=FALSE){
  with(model$outputs, stopifnot(exists("eirpr")))
  PR = model$outputs$eirpr$pr
  m = model$outputs$eirpr$m
  if(extend==TRUE){
    PR = c(0, PR, 1)
    m = c(0, m, 10*max(m))
  }
  ix = which(pr<=max(PR) & pr>=min(PR))
  get_m = function(pr){
    if(pr == min(PR)) return(min(m))
    if(pr == max(PR)) return(max(m))
    ix = max(which(PR<pr))
    ff = (pr-PR[ix])/(PR[ix+1]-PR[ix])
    m=m[ix] + ff*(m[ix+1]-m[ix])
    return(m)
  }
  mm = 0*pr-1
  mm[ix] = sapply(pr[ix], get_m)
  pr2m = list(pr=pr[ix], m=mm[ix])
  if(length(ix)>0) pr2m$errors = c(pr=pr[-ix])
  return(pr2m)
}

#' @title Convert pr to lambda
#' @description Use the outputs$lambdapr table to convert a set of pr values into lambda values
#' @param pr a [vector]
#' @param model a [list]
#' @param extend a [logical] option to determine whether to extend outside the range
#' @export
xde_pr2lambda = function(pr, model, extend=FALSE){
  with(model$outputs, stopifnot(exists("eirpr")))
  PR = model$outputs$eirpr$pr
  lambda = model$outputs$eirpr$lambda
  if(extend==TRUE){
    PR = c(0, PR, 1)
    lambda = c(0, lambda, 5*10^3/365)
  }
  ix = which(pr<=max(PR) & pr>=min(PR))
  get_lambda = function(pr){
    if(pr == min(PR)) return(min(lambda))
    if(pr == max(PR)) return(max(lambda))
    ix = max(which(PR<pr))
    ff = (pr-PR[ix])/(PR[ix+1]-PR[ix])
    lambda=lambda[ix] + ff*(lambda[ix+1]-lambda[ix])
    return(lambda)
  }
  llambda = 0*pr-1
  llambda[ix] = sapply(pr[ix], get_lambda)
  pr2lambda = list(pr=pr[ix], lambda=llambda[ix])
  if(length(ix)>0) pr2lambda$errors = c(pr=pr[-ix])
  return(pr2lambda)
}
