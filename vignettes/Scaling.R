## -----------------------------------------------------------------------------
suppressMessages(library(exDE))
suppressMessages(require(deSolve))
suppressMessages(require(rootSolve))
suppressMessages(require(mobwork))

## -----------------------------------------------------------------------------
F_eir = function(a, bday=0, scale=1){scale*(1.01 + sin(2*pi*(a+bday)/365))}

## -----------------------------------------------------------------------------
scl = integrate(F_eir, 0, 365)$value

## -----------------------------------------------------------------------------
make_F_eir = function(F_eir, bday, scale){return(function(a, scl=1){scl*F_eir(a, bday, scale)})}
F_eir1 = make_F_eir(F_eir, 0, 1/scl)

## -----------------------------------------------------------------------------
integrate(F_eir1, 0, 365, scl=2.3)$value

## -----------------------------------------------------------------------------
xde_setup_cohort(F_eir) -> sis

## -----------------------------------------------------------------------------
xde_solve(sis) -> sis

## -----------------------------------------------------------------------------
tt = seq(0:365)
p1 = list() 
p1$eirScale = 1 
plot(tt, F_eir(tt, scale=5), type = "l", xlab = "time (days)", ylab = expression(F[eir](t)))
lines(tt, F_eir(tt, scale=2)) 
lines(tt, F_eir(tt)) 
lines(tt, F_eir(tt, scale = 1/2)) 

## -----------------------------------------------------------------------------
xde_setup_cohort(F_eir) -> sis

## -----------------------------------------------------------------------------
xde_scaling_eir(sis, 25) -> sis

## -----------------------------------------------------------------------------
plot_eirpr(sis)

## -----------------------------------------------------------------------------
require(viridis)
clrs = turbo(25)

## ----fig.height=3.5, fig.width=5----------------------------------------------
plot_eirpr(sis)

with(sis$output$eirpr,{
  points(aeir, pr, col = clrs)
  lines(scaling[[5]]$aeir, scaling[[5]]$pr, col = clrs[5])
  lines(scaling[[10]]$aeir, scaling[[10]]$pr, col = clrs[10])
  lines(scaling[[15]]$aeir, scaling[[15]]$pr, col = clrs[15])
  lines(scaling[[20]]$aeir, scaling[[20]]$pr, col = clrs[20])
})

## -----------------------------------------------------------------------------
preir_i = xde_pr2eir(c(0.001, runif(25, 0, 1), 0.999), sis)

## -----------------------------------------------------------------------------
preir_i$errors

## -----------------------------------------------------------------------------
plot_eirpr(sis)
with(sis$outputs$eirpr, points(aeir, pr, pch = 15))
with(preir_i, points(365*eir, pr, pch = 19, col = "red"))

## ----eval=F-------------------------------------------------------------------
#  sis2 = split_stratum_by_biting(sis, 1, 1, .2, sqrt(10))
#  sis2 = split_stratum_by_biting(sis, 1, 1, .5, 1/sqrt(10))
#  sis2 <- xde_solve(sis2)
#  xde_scaling_eir(sis2, 25) -> sis2

## ----eval=F-------------------------------------------------------------------
#  plot(sis$outputs$eirpr$aeir, sis$outputs$eirpr$pr, type = "l", log = "x", xaxt= "n", xlab = "aEIR", ylab = "PR")
#  axis(1, 10^(-1:3), c(0.1, 1, 10, 100, 1000))
#  lines(sis2$outputs$eirpr$aeir, sis2$outputs$eirpr$pr, col = "darkred")

## -----------------------------------------------------------------------------
F_eir0 = function(a, bday=0, scale=1){scale*(0*a+1)}

## -----------------------------------------------------------------------------
sis0 = xde_setup_cohort(F_eir0)
xde_scaling_eir(sis0, 25) -> sis0

## -----------------------------------------------------------------------------
clrs = turbo(25)
with(sis$outputs$eirpr, plot(aeir, pr, type = "l", log = "x", xaxt= "n", xlab = "aEIR", ylab = "PR"))
axis(1, 10^(-1:3), c(0.1, 1, 10, 100, 1000))
lines(sis0$outputs$eirpr$aeir, sis0$outputs$eirpr$pr, col = "tomato", lwd=2) 

with(sis$outputs$eirpr, points(aeir, pr, col = clrs))
with(sis$outputs$eirpr, lines(scaling[[5]]$aeir, scaling[[5]]$pr, col = clrs[5]))
with(sis$outputs$eirpr, lines(scaling[[10]]$aeir, scaling[[10]]$pr, col = clrs[10]))
with(sis$outputs$eirpr, lines(scaling[[15]]$aeir, scaling[[15]]$pr, col = clrs[15]))
with(sis$outputs$eirpr, lines(scaling[[20]]$aeir, scaling[[20]]$pr, col = clrs[20]))

## -----------------------------------------------------------------------------
sip = xde_setup_cohort(F_eir0, Xname = "SIP")
sip$Xpar[[1]]$eta = 1/40
xde_scaling_eir(sip, 25) -> sip
sip1 = setup_exposure_nb(sip, 1/50)
xde_scaling_eir(sip1, 25) -> sip1

## -----------------------------------------------------------------------------
with(sis$outputs$eirpr, plot(aeir, pr, type = "l", log = "x", xaxt= "n", xlab = "aEIR", ylab = "PR"))
axis(1, 10^(-1:3), c(0.1, 1, 10, 100, 1000))
with(sip$outputs$eirpr, lines(aeir, pr, col = "darkorange"))
with(sip1$outputs$eirpr, lines(aeir, pr, col = "brown"))

## -----------------------------------------------------------------------------
#sis3 <- setup_exposure_nb(sis2, 1/50)
sis4 <- setup_exposure_nb(sis, 1/50)
#xde_scaling_eir(sis3, 25) -> sis3
xde_scaling_eir(sis4, 25) -> sis4

## -----------------------------------------------------------------------------
with(sis$outputs$eirpr, plot(aeir, pr, type = "l", log = "x", xaxt= "n", xlab = "aEIR", ylab = "PR"))
axis(1, 10^(-1:3), c(0.1, 1, 10, 100, 1000))
#with(sis2$outputs$eir, lines(aeir, pr, col = "blue"))
#with(sis3$outputs$eir, lines(aeir, pr, col = "purple"))
with(sis4$outputs$eir, lines(aeir, pr, col = "darkblue"))

## -----------------------------------------------------------------------------
sis5 <- setup_travel_foi(sis, delta_scale = 1/5/365)
xde_scaling_eir(sis5, 25) -> sis5

## -----------------------------------------------------------------------------
with(sis$outputs$eirpr, plot(aeir, pr, type = "l", log = "x", xaxt= "n", xlab = "aEIR", ylab = "PR"))
axis(1, 10^(-1:3), c(0.1, 1, 10, 100, 1000))
with(sis5$outputs$eir, lines(aeir, pr, col = "darkgreen")) 

