## ----eval=F-------------------------------------------------------------------
#  dXdt.SIS <- function(t, y, pars, FoI) {
#    with(pars$Xpar, {
#      X <- y[X_ix]
#      H <- F_H(t, y, pars)
#  
#      dX <- FoI*(H - X) - r*X
#  
#      return(c(dX))
#    })
#  }

## ----eval=F-------------------------------------------------------------------
#  dXdt.SISdH <- function(t, y, pars, FoI) {
#    with(pars$Xpar, {
#  
#      H <- F_H(t, y, pars)
#      X <- y[X_ix]
#  
#      dX <- FoI*(H - X) - r*X + dHdt(t, X, pars)
#      dH <- Births(t, H, pars) + dHdt(t, H, pars)
#  
#      return(c(dX, dH))
#    })
#  }

## ----eval=F-------------------------------------------------------------------
#  setup_X.SIS = function(pars, Xname, Xopts=list()){
#    pars$Xname = "SIS"
#    pars = make_Xpar_SIS(pars, Xopts)
#    pars = make_Xinits_SIS(pars, Xopts)
#  
#    return(pars)
#  }

## ----eval=F-------------------------------------------------------------------
#  make_Xpar_SIS = function(pars, Xopts=list(),
#                           b=0.55, r=1/180, c=0.15){
#    with(Xopts,{
#      Xpar = list()
#      class(Xpar) <- c("SIS", "SISdX")
#  
#      Xpar$b = checkIt(b, pars$nStrata)
#      Xpar$c = checkIt(c, pars$nStrata)
#      Xpar$r = checkIt(r, pars$nStrata)
#  
#      pars$Xpar = Xpar
#      return(pars)
#  })}

## -----------------------------------------------------------------------------
make_Xinits_SIS = function(pars, Xopts = list(), X0=1){with(Xopts,{
  inits = list()
  inits$X0 = checkIt(X0, pars$nStrata)
  pars$Xinits = inits
  return(pars)
})}


