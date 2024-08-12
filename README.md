# hydro_bsc_project
Project to do hydrodynamics on cartesian and voronoi meshes. Includes Finite Volume MUSCL solver for shallow water equations. Will eventually also include a second order FV solver for euler equations and maybe a DG solver

## Changes in last update
---
### Better circular boundary
added point in the middle of circle which makes mesh generation way more stable
<p align="center">
  <img src="/figures/output.png" alt="1" width="45%">
</p>

---
### MUSCL + another slope limiter
MUSCL scheme, slope limited for diagonal shock on voronoi:

<p align="center">
  <img src="/figures/_12_shock_tube_vornoi_diagonal.png" alt="1" width="45%">
</p>

First used slope limiter is not total variation diminishing, which leads to small oscillations at shocks. This leads to problems for e.g. KH instability

<p align="center">
  <img src="/figures/_12_KH_anim_with_shocks.gif" alt="1" width="45%">
</p>

Therefore implemented another slope limiter which holds TVD property:

first slope limiter used (e.g. arepo)
<p align="center">
  <img src="/figures/img1.png" alt="1" width="45%">
</p>

tvd slope limiter (by Duffell, P. C., and A. I. MacFadyen 2011 for TESS code)
<p align="center">
  <img src="/figures/img2.png" alt="1" width="40%">
</p>

It is however way more diffusive for shock tube. Therefore if possible: use the first one and if needed go to more robust slope limiter.

<p align="center">
  <img src="/figures/_12_shock_tube_1D_1st_MUSCL_TVD.png" alt="1" width="45%">
</p>

---
### Kelvin Helmholtz instability
Since the second slope limiter is too diffusive to get satisfying KH instability on my computer in reasonable time we need another solution. 
-> Adapted mesh resolution (e.g. more points around instability and less far away from it) + still use the non tvd slope limiter
Idea: lower resolution away from the initial instability is on larger scale than those shocks seen earlier. Therefore maybe they get smoothed away.

Point density given by uniform distribution in x direction and double gaussian in y direction + 5 lloyd iterations

<p align="center">
  <img src="/figures/_12_KH_adapted_mesh.png" alt="1" width="49%">
  <img src="/figures/_12_KH_adapted_mesh_final.gif" alt="1" width="49%">
</p>

The initial condition for this video includes no additional pertubation. So the pertubation from the voronoi mesh itself is enough.

---
### Rayleigh Taylor instability
Added gravity source term:

<p align="center">
  <img src="/figures/img3.png" alt="1" width="45%">
</p>

that in the FV update scheme can be just added here:
<p align="center">
  <img src="/figures/img4.png" alt="1" width="35%">
</p>

for a sin(x) initial condition with $\rho_1 > \rho_2, \; P_1 = P_2$ and no initial velocity we get:

<p align="center">
  <img src="/figures/_12_RT_voronoi.gif" alt="1" width="50%">
</p>

where the first flux limiter was used. One can clearly see a lot of oscillations that shouldn't happen. Also there is this bouncing because of the pressure gradient first building up first but im not sure how to get rid of that?

### "Quad-shock" (2D Riemann like problem)
Question: what happens in [your video](https://www.ita.uni-heidelberg.de/~dnelson/#group) in the bottom left corner? Is there some capturing of instabilies ongoing?
Given initial conditions:
<p align="center">
  <img src="/figures/img5.png" alt="1" width="45%">
</p>

For this simulation because of computation time there is also higher resolution in the bottom left quarter. Additionally to that the TVD slope limiter was used for the regions near the boundary because of instability reasons while the less diffusive non TVD slope limiter was used in the rest of the plot.
We get:
<p align="center">
  <img src="/figures/_12_quad_shock_mesh.png" alt="1" width="49%">
  <img src="/figures/_12_quad_shock_final.gif" alt="1" width="49%">
</p>

---
### other plots
low resolution KH, vortex shedding
<p align="center">
  <img src="/figures/_12_KH_lres2.gif" alt="1" width="49%">
  <img src="/figures/vortex_shredding_try_2-125 (verschoben).png" alt="1" width="49%">
</p>

---
### Discontinous Galerkin methods (intro)
-> intro into DG ([see pdf](./figures/DG_intro.pdf))

first rudimentary attempts (Np = 1, 1D, scalar upwind advection in python with manually calculated matrices)

<p align="center">
  <img src="/figures/_12_1D_DG_advection.gif" alt="1" width="49%">
  <img src="/figures/_12_1D_DG_advection_2.gif" alt="1" width="49%">
</p>

---
### next?
- how do we compare those euler FV animations to anything quantively? do we even want to do that? or instead go to discontinous galerkin fully?
- maybe use Eigen library in C++ since from now on we need many matrices and operations on them
- at least for 1D unclear how we can reuse any of the written program structure
- i think we should add a new DG_Solver class such that FV and DG are splitted cleanly
