# hydro_bsc_project
project to eventually do hydrodynamics on different meshes, with maybe different physics and different solvers. still work in progress...

## Changes in last update
---
- visualization.ipynb -> vis_tk.py package

---
### Custom boundaries
Option to turn cells into boundary cells, int theory one could generate grids with arbitrary shapes here
<p align="center">
  <img src="/figures/2boundarystuff_diffustion.gif" alt="2boundarystuff_diffustion" width="45%">
  <img src="/figures/2boundarystuff_advection.gif" alt="2boundarystuff_advection" width="45%">
</p>

---
### Shallow Water Equations 2D

Using
<p align="center">
  <img src="/figures/0.png" alt="0" width="65%">
</p>

the equation reads

<p align="center">
  <img src="/figures/1.png" alt="1" width="25%">
</p>

Integrating over a Cell i leads to

<p align="center">
  <img src="/figures/2.png" alt="2" width="40%">
</p>

and applying Gauss's theorem we get

<p align="center">
  <img src="/figures/3.png" alt="3" width="40%">
</p>

For Voronoi we thus get

<p align="center">
  <img src="/figures/4.png" alt="4" width="50%">
</p>

which in the cartesian case simplifies to 

<p align="center">
  <img src="/figures/5.png" alt="5" width="60%">
</p>

### Flux approximations
But how do we approximate the fluxes?
Upwind? Can this even work since for example if $(\vec{F}_{i-\frac{1}{2}})_0 = - (\vec{F}_{i-\frac{1}{2}})_1$ how do i even define an upwind side. Componentwise?

<p align="center">
  <img src="/figures/swe_upwind_scheme_whatever.gif" alt="swe_upwind_scheme_whatever" width="45%">
</p>
Well that doesn't seem to work. What about Lax-Friedrichs as used before in Advection?

<p align="center">
  <img src="/figures/6.png" alt="6" width="45%">
</p>

<p align="center">
  <img src="/figures/swe_pure_lax_friedrich_aua.gif" alt="swe_pure_lax_friedrich_aua" width="45%">
</p>
Hmm. What if we just try FTCS (Finite Time Centered Space)

<p align="center">
  <img src="/figures/swe_first_try_ftcs.gif" alt="swe_first_try_ftcs" width="45%">
</p>
I mean it looks cool but out of one line there appear many and it needs additional vastly high resolution. Maybe some diffusion could solve this problem though less than Lax-Friedrichs. 
-> This is just a desperate try of me. Nothing foundated.
Try 

<p align="center">
  <img src="/figures/7.png" alt="7" width="50%">
</p>

Then we get
<p align="center">
  <img src="/figures/swe_first_try_lax_friedrich_less_diffusion.gif" alt="swe_first_try_lax_friedrich_less_diffusion" width="33%">
  <img src="/figures/swe_first_try_lax_friedrich_less_diffusion_u.gif" alt="swe_first_try_lax_friedrich_less_diffusion_u" width="33%">
  <img src="/figures/swe_first_try_lax_friedrich_less_diffusion_v.gif" alt="swe_first_try_lax_friedrich_less_diffusion_v" width="33%">
</p>

which at least looks somewhat like expected. I think we should update to some kind of Riemann solver as soon as possible or any other more suited Flux approximations for SWE e.g. like Roe/HLL. For now it seems to work by pure coincidence. I am quite sure the difficulties lie in the Flux approximations since i am quite confident about the initial update formulas for the FV scheme.

### Voronoi

Using

<p align="center">
  <img src="/figures/8.png" alt="8" width="50%">
</p>

and the adapted Lax-Friedrichs (but also tried all the other ones with similar results) we get

<p align="center">
  <img src="/figures/swe_voronoi.gif
  " alt="swe_voronoi_voronoi_alg" width="45%">
  <img src="/figures/swe_cartesian_voronoi_algorithm.gif" alt="swe_cartesian_voronoi_algorithm" width="45%">
</p>
while on the right we applied the general voronoi formula to the cartesian grid structure leading to the limit discussed before. To be honest i have no idea why the voronoi version doesn't work. I supect again the flux approximations here, since i am quite sure that the FV formula is right. But again a big mistery here.

---
### Boundary conditions for Shallow Water equations (cartesian)
1. Repeating (intrinsic in mesh structure)
2. Sink
<p align="center">
  <img src="/figures/9.png" alt="9" width="20%">
</p>
3. Reflecting
<p align="center">
  <img src="/figures/10.png" alt="10" width="20%">
</p>


<p align="center">
  <img src="/figures/swe_repeating_boundary_conditions.gif" alt="swe_repeating_boundary_conditions" width="45%">
  <img src="/figures/swe_sink_boundary.gif" alt="swe_sink_boundary" width="45%">
</p>
<p align="center">
  <img src="/figures/swe_boundary_hot.gif" alt="swe_boundary_hot" width="65%">
</p>

--- 
### Other examples
Box initial condition
<p align="center">
  <img src="/figures/hot_square_swe.gif" alt="hot_square_swe" width="65%">
</p>
Initial Flowing Box with no initial height
<p align="center">
  <img src="/figures/swe_flow_h.gif" alt="swe_flow_h" width="45%">
  <img src="/figures/swe_flow_u.gif" alt="swe_flow_u" width="45%">
</p>

---
still to do:

- Generalizing Boundary conditions. How?
- Find analytical solutions
- Get 2D Voronoi Mesh Version to work. How? Riemann Solver?
