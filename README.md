# hydro_bsc_project
Project to do hydrodynamics on cartesian and voronoi meshes. Includes Finite Volume MUSCL solver for shallow water and euler equations. Next steps are towards Discontinous Galerkin methods. Work in progress...

### Summary of last weeks meeting
---
#### Rayleigh Taylor instability
fixed initial conditions: P in hydrostatic equilibrium

<p align="center">
  <img src="figures/RT_animation_425x425.gif" alt="1" width="70%">
</p>


---
#### 1D DG scalar upwind advection
for arbitrary Np and monomial basisfunctions in C++ using Eigen library

<p align="center">
  <img src="figures/12_DG_scalar_upwind_advection_20.png" alt="1" width="60%">
</p>


#### Slope limiting
<p align="center">
  <img src="figures/img_slopelim1.png" alt="1" width="60%">
</p>

<p align="center">
  <img src="figures/DG_bar_slopelimiting.gif" alt="1" width="49%">
  <img src="figures/DG_smooth_slopelimiting.gif" alt="1" width="49%">
</p>

### L1-error of advected step function
<p align="center">
  <img src="figures/l1_formulas.png" alt="1" width="60%">
</p>
<p align="center">
  <img src="figures/DG_L1_over_N.png" alt="1" width="45%">
</p>

- slope limiting again limits higher order methods to first order convergence at discontinuities and also limits the benefit of having higher order polynomials (e.g. compare Np = 1 vs Np = 2)
---
### This week: 2D DG Scalar upwind advection
In general the integrals stay roughly the same
<p align="center">
  <img src="figures/int.png" alt="1" width="100%">
</p>
flux terms replaced with flux line integral. Partial derivative replaced with Gradient.
Also the projection here now is a (in our case diagonal) matrix instead of a coefficient. Need to respect that when changing into reference element.
<p align="center">
  <img src="figures/grad_jacob.png" alt="1" width="45%">
</p>

using basis functions as in [[Schaal, K. 2016]](http://www.kmschaal.de/PhDThesis_KevinSchaal.pdf). First few legendre polynomials are

<p align="center">
  <img src="figures/legendre.png" alt="1" width="70%">
</p>
use product of legendre polynomials as basis functions. This leads to 
<p align="center">
  <img src="figures/basis_P1_P2.png" alt="1" width="70%">
</p>

three different options implemented:
N_p = 0 : phi 0
N_p = 1: phi 0 to 3
N_p = 2: phi 0 to 5

#### Implementation
N_p = 0 and N_p = 1 of square without slope limiting on a 20x20 grid

<p align="center">
  <img src="figures/dg2d_animation_np0.gif" alt="1" width="45%">
  <img src="figures/dg2d_animation_np1.gif" alt="1" width="45%">
</p>

slope limited versions for N_p = 0, N_p = 1, N_p = 2 on a 100x100 grid
<p align="center">
  <img src="figures/0_P0_basis_100x100.png" alt="1" width="32%">
  <img src="figures/0_P1_basis_100x100_beta2.png" alt="1" width="32%">
  <img src="figures/0_P2_basis_100x100_beta2.png" alt="1" width="32%">
</p>
for the slope limiting it seems again that the higher order advantages get lost here (i.e almost no difference between second and third order) since effectively everything is limited to linear elements.

---
For smooth solutions without an active limiter the improvements are way more obvious since the limiting is not reducing the order here.

not slope limited gaussian for N_p = 0, N_p = 1, N_p = 2 on a 20x20 grid
<p align="center">
  <img src="figures/test6.png" alt="1" width="49%">
  <img src="figures/test5.png" alt="1" width="49%">
  <img src="figures/test4.png" alt="1" width="49%">
</p>

--- 
### Side project for FV euler: Airfoil
Thought that it just might be a cool idea to add an airfoil as another test scenario for the FV code (along KH, RT and Shock tube).

Plotted for different angles of attack. Colormap is pressure. Arrows give flow direction. Clearly for to high angles of attack the plane stalls.

<p align="center">
  <img src="figures/airfoil.png" alt="1" width="100%">
</p>

---
### other
- Who should we ask as a second examiner for the bachelor kolloquium?
- Decide on Plot style, x and y direction for the 2D plots or leave it completely?
<p align="center">
  <img src="figures/image2D_1.png" alt="1" width="50%">
  <img src="figures/KH_instab.png" alt="1" width="70%">
</p>



---
### next?
- add plot routines for DG from notebook to vis_tk
- change element wise calculation of M_ij and S_ij to one single global calculation
- Do the L1 error also for 2D and especially for the gaussian shape as well -> hopefully we see that for the smooth solutions the convergence is way better? would show that one needs more advanced slope limiters and better oscillation detection to benefit from higher order method
- For that we need to implement proper initial conditions using L2 projection and solving the integral there using gauss quadrature

further ideas?