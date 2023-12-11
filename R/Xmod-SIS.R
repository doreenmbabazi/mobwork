#' Basic plotting: plot the density of infected humans for an SIS model
#'
#' @inheritParams plot_X
#' @export
plot_X.SIS = function(model, clrs="black", stable=FALSE){
  vars = if(stable==TRUE){
    model$outputs$stable_orbits
  }else{
    model$outputs$orbits
  }
  with(vars$XH,
      plot(time, 0*time + max(H), type = "n",
           ylab = "# Infected", xlab = "Time"))

  lines_X(vars$XH, model, clrs)
}

#' Add the orbits for the SIS model to a plot for models of human infection and immunity
#'
#' @inheritParams lines_X
#'
#' @export
lines_X.SIS = function(XH, model, clrs="black"){
  with(XH,{
    if(model$nStrata==1) lines(time, X, col=clrs)
    if(model$nStrata>1){
      if (length(clrs)==1) clrs=rep(clrs, nStrata)
      for(i in 1:model$nStrata){
        lines(time, X[,i], col=clrs[i])
      }
    }
})}
