
#' Basic plotting: plot mosquito population density for the "RM" model
#'
#' @inheritParams plot_M
#'
#' @export
plot_M.RM = function(model, clrs = "darkblue", stable=FALSE){
  vars = if(stable==TRUE){
    model$outputs$stable_orbits
  }else{
    model$outputs$orbits
  }

  with(vars$MYZ,
    plot(time, 0*time, type = "n", ylim = range(0,max(M)),
         ylab = "Mosquito Density", xlab = "Time"))

  lines_M(vars$MYZ, model, clrs)
}

#' Add the orbits for the SIS model to a plot for models of human infection and immunity
#'
#' @inheritParams lines_M
#'
#' @export
lines_M.RM = function(MYZ, model, clrs="darkblue"){
  with(MYZ,{
    if(model$nPatches==1) lines(time, M, col=clrs)
    if(model$nPatches>1){
      if (length(clrs)==1) clrs=rep(clrs, model$nPatches)
      for(i in 1:model$nPatches){
        lines(time, M[,i], col=clrs[i])
      }
    }
})}

#' Basic plotting: plot mosquito population density for the "RM" model
#'
#' @inheritParams plot_YZ
#'
#' @export
plot_YZ.RM = function(model, Yclrs = "purple", Zclrs = "darkred", stable=FALSE){
  vars = if(stable==TRUE){
    model$outputs$stable_orbits
  }else{
    model$outputs$orbits
  }

  with(vars$MYZ,
    plot(time, 0*time, type = "n", ylim = range(0,max(Y)),
         ylab = "Mosquito Density", xlab = "Time"))

  lines_YZ(vars$MYZ, model, Yclrs, Zclrs)
}

#' Add the orbits for the SIS model to a plot for models of human infection and immunity
#'
#' @inheritParams lines_YZ
#'
#' @export
lines_YZ.RM = function(MYZ, model, Yclrs="purple", Zclrs = "darkred"){
  with(MYZ,{
    if(model$nPatches==1){
      lines(time, Y, col=Yclrs)
      lines(time, Z, col=Zclrs)
    }
    if(model$nPatches>1){
      if (length(Yclrs)==1)
        Yclrs=rep(Yclrs, model$nPatches)
      if (length(Zclrs)==1)
        Zclrs=rep(Zclrs, model$nPatches)

      for(i in 1:model$nPatches){
        lines(time, Y[,i], col=Yclrs[i])
        lines(time, Z[,i], col=Zclrs[i])
      }
    }
})}


#' Basic plotting: plot the fraction infected and infective for the "RM" model
#'
#' @inheritParams plot_YZ_fracs
#'
#' @export
plot_YZ_fracs.RM = function(model, Yclrs = "purple", Zclrs = "darkred", stable=FALSE){
  vars = if(stable==TRUE){
    model$outputs$stable_orbits
  }else{
    model$outputs$orbits
  }

  with(vars$MYZ,
       plot(time, 0*time, type = "n", ylim = range(0,1),
            ylab = "Fraction Infected", xlab = "Time"))

  lines_YZ_fracs(vars$MYZ, model, Yclrs, Zclrs)
}

#' Add lines for the fraction infected and infective for the "RM" model
#'
#' @inheritParams lines_YZ_fracs
#'
#' @export
lines_YZ_fracs.RM = function(MYZ, model, Yclrs="purple", Zclrs="darkred"){
  with(MYZ,{
    if(model$nPatches==1) {
      lines(time, y, col=Yclrs)
      lines(time, z, col=Zclrs)
    }
    if(model$nPatches>1){
      if (length(Yclrs)==1)
        Yclrs=rep(Yclrs, model$nPatches)
      if (length(Zclrs)==1)
        Zclrs=rep(Zclrs, model$nPatches)

      for(i in 1:model$nPatches){
        lines(time, y[,i], col=Yclrs[i])
        lines(time, z[,i], col=Zclrs[i])
      }
    }
  })}
