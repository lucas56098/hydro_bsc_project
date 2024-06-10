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

$$ \vec{U} = \begin{pmatrix} h \\ hu \\ hv \end{pmatrix}, \;\;\; \vec{F} = \begin{pmatrix} hu \\ hu^2 + gh^2/2 \\ hvu \end{pmatrix}, \;\;\; \vec{G} = \begin{pmatrix} hv \\ hvu \\ hv^2 + gh^2/2 \end{pmatrix} $$

the equation reads

$$ \frac{\partial \vec{U}}{\partial t} + \frac{\partial \vec{F}}{\partial x} + \frac{\partial \vec{G}}{\partial y} = 0$$

Integrating over a Cell i leads to

$$ \frac{\partial}{\partial t} \int_{C_i}{dA \;\vec{U}} + \int_{C_i}{dA \; \vec{\nabla}(\vec{F}, \vec{G})} = 0$$

and applying Gauss's theorem we get

$$ A \frac{\partial <\vec{U_i}>}{\partial t} + \int_{\partial C_i}{ds \; (\vec{F}, \vec{G})\cdot \hat{n}} $$

For Voronoi we thus get

$$ \vec{U}_i^{n+1} = \vec{U}_i^{n} - \frac{\Delta t}{A} \sum_{j \in \partial C_i}{l_{i, j} \cdot (\vec{F}_{i,j} n_x + \vec{G}_{i,j}n_y)} $$

which in the cartesian case simplifies to 

$$ \vec{U}_i^{n+1} = \vec{U}_i^{n} + \frac{\Delta t}{\Delta x}(\vec{F}_{i-\frac{1}{2}} - \vec{F}_{i+\frac{1}{2}}) + \frac{\Delta t}{\Delta y}(\vec{G}_{i-\frac{1}{2}} - \vec{G}_{i+\frac{1}{2}}) $$

### Flux approximations
But how do we approximate the fluxes?
Upwind? Can this even work since for example if $(\vec{F}_{i-\frac{1}{2}})_0 = - (\vec{F}_{i-\frac{1}{2}})_1$ how do i even define an upwind side. Componentwise?

<p align="center">
  <img src="/figures/swe_upwind_scheme_whatever.gif" alt="swe_upwind_scheme_whatever" width="45%">
</p>
Well that doesn't seem to work. What about Lax-Friedrichs as used before in Advection?

$$ \vec{F}_{i-\frac{1}{2}} = \frac{1}{2} [\vec{F}_{i-1} + \vec{F}_{i}] - \frac{\Delta x}{2\Delta t} (\vec{U}_i - \vec{U}_{i-1}) $$

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

$$ \vec{F}_{i-\frac{1}{2}} = \frac{1}{2} [\vec{F}_{i-1} + \vec{F}_{i}] - \frac{1}{100} \cdot \frac{\Delta x}{2\Delta t} (\vec{U}_i - \vec{U}_{i-1}) $$

Then we get
<p align="center">
  <img src="/figures/swe_first_try_lax_friedrich_less_diffusion.gif" alt="swe_first_try_lax_friedrich_less_diffusion" width="33%">
  <img src="/figures/swe_first_try_lax_friedrich_less_diffusion_u.gif" alt="swe_first_try_lax_friedrich_less_diffusion_u" width="33%">
  <img src="/figures/swe_first_try_lax_friedrich_less_diffusion_v.gif" alt="swe_first_try_lax_friedrich_less_diffusion_v" width="33%">
</p>

which at least looks somewhat like expected. I think we should update to some kind of Riemann solver as soon as possible or any other more suited Flux approximations for SWE e.g. like Roe/HLL. For now it seems to work by pure coincidence. I am quite sure the difficulties lie in the Flux approximations since i am quite confident about the initial update formulas for the FV scheme.

### Voronoi

Using

$$ \vec{U}_i^{n+1} = \vec{U}_i^{n} - \frac{\Delta t}{A} \sum_{j \in \partial C_i}{l_{i, j} \cdot (\vec{F}_{i,j} n_x + \vec{G}_{i,j}n_y)} $$

and the adapted Lax-Friedrichs (but also tried all the other ones with similar results) we get

<p align="center">
  <img src="/figures/swe_voronoi_voronoi_alg.gif
  " alt="swe_voronoi_voronoi_alg" width="45%">
  <img src="/figures/swe_cartesian_voronoi_algorithm.gif" alt="swe_cartesian_voronoi_algorithm" width="45%">
</p>
while on the right we applied the general voronoi formula to the cartesian grid structure leading to the limit discussed before. To be honest i have no idea why the voronoi version doesn't work. I supect again the flux approximations here, since i am quite sure that the FV formula is right. But again a big mistery here.

---
### Boundary conditions for Shallow Water equations (cartesian)
1. Repeating (intrinsic in mesh structure)
2. Sink
$$ \begin{pmatrix} h^* \\ h^*u^* \\ h^*v^* \end{pmatrix} = \begin{pmatrix} h \\ hu \\ hv \end{pmatrix} $$
3. Reflecting
$$ \begin{pmatrix} h^* \\ h^*u^* \\ h^*v^* \end{pmatrix} = \begin{pmatrix} h \\ -hu \\ -hv \end{pmatrix} $$


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