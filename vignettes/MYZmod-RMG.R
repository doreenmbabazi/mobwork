## -----------------------------------------------------------------------------
model <- xde_setup(MYZname="RMG", Lname="trace", Xname = "trace", HPop = rep(1, 4), residence = 1:4, nPatches=4)

## -----------------------------------------------------------------------------
model <- xde_solve(model)

## ----fig.height=7.5, fig.width=5.5--------------------------------------------
par(mfrow = c(2,1))
xde_plot_M(model)
xde_plot_YZ(model, add_axes = F)
xde_plot_YZ_fracs(model)

