# hydro_bsc_project
Project to do hydrodynamics on cartesian and voronoi meshes. Includes Finite Volume MUSCL solver for shallow water and euler equations. Next steps are towards Discontinous Galerkin methods. Work in progress...

## Changes in last update
---
### Rayleigh Taylor instability
fixed initial conditions (pressure is now in hydrostatic equilibrium -> no more bouncing)

<p align="center">
  <img src="figures/RT_animation_425x425.gif" alt="1" width="45%">
  <img src="figures/image2D_1.png" alt="1" width="42%">
</p>

same triple plot tried for the quad shock simulation
<p align="center">
  <img src="figures/plot_idea.png" alt="1" width="75%">
</p>


---
## Discontinous Galerkin 1D scalar upwind advection

implemented 1d scalar upwind advection for arbitrary Np and monomial basisfunctions in C++ using Eigen library
even for very low resolution N = 20, the higher order schemes advect the initial condition quite well.

<p align="center">
  <img src="figures/12_DG_scalar_upwind_advection_20.png" alt="1" width="60%">
  <img src="figures/DG_scalar_upwind_advection.gif" alt="1" width="45%">
  <img src="figures/DG_scalar_upwind_advection2.gif" alt="1" width="45%">
</p>

Seems like we need slope limiting next

## Slope limiting
<p align="center">
  <img src="figures/img_slopelim1.png" alt="1" width="60%">
</p>


works quite well for discontinuities
<p align="center">
  <img src="figures/DG_bar_slopelimiting.gif" alt="1" width="60%">
</p>
however slope limiting does more bad than good to already smooth solutions
<p align="center">
  <img src="figures/DG_smooth_slopelimiting.gif" alt="1" width="60%">
</p>
therefore a more advanced slope limiting (e.g. better oscillation detection) would be nice but i guess thats for another project?

## L1-error comparison

<p align="center">
  <img src="figures/l1_formulas.png" alt="1" width="60%">
</p>
the integral here is approximated numerically by a sum where we can choose the resolution

### L1 error of advected step function
<p align="center">
  <img src="figures/DG_L1_over_time.png" alt="1" width="45%">
  <img src="figures/DG_L1_over_N.png" alt="1" width="45%">
</p>

- slope limiting again limits higher order methods to first order convergence at discontinuities and also limits the benefit of having higher order polynomials (e.g. compare Np = 1 vs Np = 2)
- interesting though that time evolution though is way less diffusive for the higher order methods? I guess we are not smoothing away structure after every step?


---
### next?
Question: Should we go 2D next? Or instead try to solve 1D SWE or Euler? I think 2D would be conceptually more interesting (e.g. nodal basis functions, more advanced mapping etc.). 1D on the other hand would allow for more comparison between Finite Volume and Discontinous Galerkin (L1 errors and so on)