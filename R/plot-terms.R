
#' Basic plotting: plot the true EIR
#'
#' @param model a [list] specifying the model
#' @param clrs a [list] of colors
#' @param stable [logical]
#'
#' @export
plot_EIR <- function(model, clrs="black", stable=FALSE){
  vars = if(stable==TRUE){
    model$outputs$stable_orbits
  }else{
    model$outputs$orbits
  }
  with(vars$terms,{
    plot(time, eir, type = "l", xlab = "Time", ylab = "dEIR", ylim = range(0, eir), col = clrs)
    tk4 = pretty(eir)
    clr = "darkblue"
    axis(4, tk4, tk4, col=clr, col.axis =clr)
    mtext("aEIR", 4, -1.5, col = clr)
})}

#' Add lines for the true EIR
#'
#' @param terms a [list] with the outputs of exDE::parse_deout
#' @param model an [list]
#' @param clrs a [list]
#'
#' @export
lines_EIR <- function(terms, model, clrs="black"){
  with(terms,{
    if(model$nStrata==1) lines(time, eir, col=clrs)
    if(model$nStrata>1){
      if (length(clrs)==1) clrs=rep(clrs, nStrata)
      for(i in 1:model$nStrata)
        lines(time, eir[,i], col=clrs[i])
    }
  })
}

#' Basic plotting: plot the true aEIR
#'
#' @param model a [list] specifying the model
#' @param clrs a [list] of colors
#' @param stable [logical]
#'
#' @export
plot_aEIR <- function(model, clrs="black", stable=FALSE){
  vars = if(stable==TRUE){
    model$outputs$stable_orbits
  }else{
    model$outputs$orbits
  }
  with(vars$terms,{
    plot(time, eir, type = "l", xlab = "Time", ylab = "aEIR", ylim = range(0, eir), col = clrs)
    tk4 = pretty(365*eir)
    clr = "darkblue"
    axis(4, tk4, tk4, col=clr, col.axis =clr)
    mtext("aEIR", 4, -1.5, col = clr)
  })}

#' Add the orbits for the SIS model to a plot for models of human infection and immunity
#'
#' @param terms a [list] with the outputs of exDE::parse_deout
#' @param model an [list]
#' @param clrs a [list]
#'
#' @export
lines_aEIR <- function(terms, model, clrs="black"){
  with(terms,{
    if(model$nStrata==1) lines(time, eir, col=clrs)
    if(model$nStrata>1){
      if (length(clrs)==1) clrs=rep(clrs, nStrata)
      for(i in 1:model$nStrata)
        lines(time, eir[,i], col=clrs[i])
    }
  })
}

#' Plot the EIR vs. the PR
#'
#' @param model an [list]
#' @param clrs a [list]
#'
#' @export
plot_eirpr <- function(model, clrs = "black"){
  with(model$outputs$eirpr,
      plot(365*eir, pr, type = "l",
           xlab = "aEIR", ylab = "PR",
           xlim = range(0, 365*eir), ylim = c(0,1),
           col = clrs))
 }

#' Add lines for the EIR vs. the PR
#'
#' @param model an [list]
#' @param clrs a [list]
#'
#' @export
lines_eirpr <- function(model, clrs = "black"){
  with(model$outputs$eirpr, lines(365*eir, pr, col = clrs))
 }

#' Basic plotting: plot the true PR.
#'
#' @param model a [list] specifying the model
#' @param clrs a [list] of colors
#' @param stable [logical]
#'
#' @export
plot_PR = function(model, clrs, stable){
  vars = if(stable==TRUE){
    model$outputs$stable_orbits
  }else{
    model$outputs$orbits
  }
  with(vars$terms,
       plot(time, 0*time + 1, type = "n", ylim = c(0,1),
            ylab = "# Infected", xlab = "Time"))

  lines_PR(vars$terms, model$nStrata, clrs)
}

#' Add the orbits for the SIS model to a plot for models of human infection and immunity
#'
#' @param terms a [list] with the outputs of exDE::parse_deout
#' @param model an [list]
#' @param clrs a [list]
#'
#' @export
lines_PR = function(terms, model, clrs="black"){
  with(terms,{
    if(model$nStrata==1) lines(time, pr, col=clrs)
    if(model$nStrata>1){
      if (length(clrs)==1) clrs=rep(clrs, nStrata)
      for(i in 1:model$nStrata)
        lines(time, pr[,i], col=clrs[i])
    }
  })
}
