## ----echo=F, eval=F-----------------------------------------------------------
#  devtools::load_all()
#  devtools::load_all("../../exDE")

## -----------------------------------------------------------------------------
model <- xde_setup(MYZname="RMG", Lname="trace", Xname = "trace", nPatches=3)

## -----------------------------------------------------------------------------
model <- xde_solve(model)

## ----fig.height=7.5, fig.width=5.5--------------------------------------------
par(mfrow = c(2,1))
xde_plot_M(model)
xde_plot_YZ(model, add_axes = F)
xde_plot_YZ_fracs(model)

