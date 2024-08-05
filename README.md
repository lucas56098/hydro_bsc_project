# hydro_bsc_project
Project to do hydrodynamics on cartesian and voronoi meshes. Includes Finite Volume MUSCL solver for shallow water equations. Will eventually also include a second order FV solver for euler equations and maybe a DG solver

## Changes in last update
---
### Custom Boundaries
Basically one can now make arbitrarily shaped inner boundaries as long as the shape is fully connected. For that one needs the boundary points stored in order in a struct.csv file. File is generated in python (e.g. using geopandas for map data).

<p align="center">
  <img src="/figures/11_custom_boundary_plot.gif" alt="1" width="55%">
</p>

<p align="center">
  <img src="/figures/11_circular_boundary.gif" alt="1" width="49%">
  <img src="/figures/11_africa_boundary.gif" alt="1" width="49%">
</p>

---
# Euler-Equations

Euler equations out of the AM274 Lecture notes using $\rho$, u, v, E as variables:

<p align="center">
  <img src="/figures/eq1.png" alt="1" width="55%">
</p>

### 1st Order solver + Analytical test case (Shock tube)
Kept conceptually everything the same as in SWE solver. A typical test case with analytical solution is Sod's Shock Tube. For the analytical solution we use a [Package](https://github.com/ibackus/sod-shocktube) that basically does exactly what we need.


1D cartesian shock tube:

<p align="center">
  <img src="/figures/11_shock_tube_density.png" alt="1" width="32%">
  <img src="/figures/11_shock_tube_velocity.png" alt="1" width="32%">
  <img src="/figures/11_shock_tube_pressure.png" alt="1" width="32%">
</p>

2D cartesian:
- x vs y shock
- voronoi x-direction shock

<p align="center">
  <img src="/figures/11_shock_tube_xy_density.png" alt="1" width="45%">
  <img src="/figures/11_shock_tube_voronoi.png" alt="1" width="45%">
</p>

- diagonal cartesian shock
- diagonal voronoi shock

<p align="center">
  <img src="/figures/11_shock_tube_cartesian_diagonal.png" alt="1" width="45%">
  <img src="/figures/11_shock_tube_voronoi_diagonal.png" alt="1" width="45%">
</p>

-> first order seems to work? wasn't able to get KH-Instability to work with it though. Is it possible that a scheme is to diffusive or whatever that this stabilizes the situtation?

---
<p align="center">
  <img src="/figures/11_flow_against_boundary.gif" alt="1" width="45%">
  <img src="/figures/KH_with_wierd_shocks.gif" alt="1" width="45%">
</p>

---
### still to do
- further verify 1st order? what about kh instability, why does it not work, why those shocks or whatever on the right
- second order muscl + verification + convergence




