
#' @title Convert pr to eir
#' @description Use the eirpr table to convert a set of pr values into eir values
#' @param pr a [vector]
#' @param model a [list]
#' @param extend a [logical] option to determine whether to extend outside the range
#' @export
xde_pr2eir = function(pr, model, extend=FALSE){
  model = with(model, if(exists("eirpr")){model}else{xde_scaling(model)})
  PR = model$eirpr$pr
  EIR = model$eirpr$eir
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
  eirpr = list(pr=pr[ix], eir=eir[ix])
  if(length(ix)>0) eirpr$errors = c(pr=pr[-ix])
  return(eirpr)
}

#' @title Convert eir to pr
#' @description Use the eirpr table to interpolate
#' @param eir a [vector]
#' @param model a [list]
#' @param extend a [logical] option to determine whether to extend outside the range
#' @export
xde_eir2pr = function(eir, model, extend=FALSE){
  model = with(model, if(exists("eirpr")){model}else{xde_scaling(model)})
  PR = model$eirpr$pr
  EIR = model$eirpr$eir
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
  preir = list(eir=eir[ix], pr=pr[ix])
  if(length(ix)>0) preir$errors = c(pr=pr[-ix])
  return(preir)
}

#' @title Convert eir to ni
#' @description Use the eirpr table to interpolate
#' @param eir a [vector]
#' @param model a [list]
#' @param extend a [logical] option to determine whether to extend outside the range
#' @export
xde_eir2ni = function(eir, model, extend=FALSE){
  model = with(model, if(exists("eirpr")){model}else{xde_scaling(model)})
  NI = model$eirpr$ni
  EIR = model$eirpr$eir
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
  nieir = list(eir=eir[ix], ni=ni[ix])
  if(length(ix)>0) nieir$errors = c(ni=ni[-ix])
  return(nieir)
}

#' @title Convert pr to ni
#' @description Use the eirpr table to interpolate
#' @param pr a [vector]
#' @param model a [list]
#' @param extend a [logical] option to determine whether to extend outside the range
#' @export
xde_pr2ni = function(pr, model, extend=FALSE){
  model = with(model, if(exists("eirpr")){model}else{xde_scaling(model)})
  PR = model$eirpr$pr
  NI = model$eirpr$ni
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
  nipr = list(pr=pr[ix], ni=ni[ix])
  if(length(ix)>0) nipr$errors = c(pr=pr[-ix])
  return(nipr)
}

