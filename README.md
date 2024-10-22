# hydro_bsc_project
Project to do 2D hydrodynamics on Cartesian and Voronoi meshes using finite volume and discontinuous Galerkin methods. The project includes grid generation methods from the [vmp](https://github.com/lucas56098/voronoi_mesh_project) mesh generation code as well as solvers for the advection equation, the shallow water equations and Euler equations. The advection equation is solved in one and two dimensions using a simple FV upwind scheme, while the shallow water and Euler equations are solved using a first-order Godunov-type scheme and a second-order MUSCL-Hancock scheme with slope-limiting. To solve the Riemann problems a HLL approximate Riemann solver is used. For discontinous Galerkin methods a one and two dimensional Cartesian scalar upwind advection scheme is implemented using monomial basis functions and Legendre polynomials. Finally there exist various plotting and analysis routines in a python visualization toolkit. Note, that the program requires the Eigen library to run. This work was done as part of a bachelor thesis in physics at the Institute for theoretical Astrophysics in the Group of [Dr. Dylan Nelson](https://github.com/dnelson86), supervised by him and [Dr. Chris Byrohl](https://github.com/cbyrohl). The entire bachelor thesis can be accessed [here](https://github.com/lucas56098/hydro_bsc_project/blob/main/thesis.pdf). Otherwise feel free just to scroll through some of the plots produced using the program.

---
### 1D Advection
<p align="center">
  <img src="/figures/advected_step_function_evol_res100.png" alt="1" width="40%">
  <img src="/figures/advected_gaussian_evol_res100.png" alt="1" width="40%">
</p>

---
### Mesh generation with custom internal boundaries

<p align="center">
  <img src="/figures/custom_circle.png" alt="1" width="35%">
  <img src="/figures/custom_airfoil.png" alt="1" width="35%">
</p>

---
### 2D shallow water equations

<p align="center">
  <img src="/figures/box5.png" alt="1" width="80%">
  <img src="/figures/circle4.png" alt="1" width="80%">
  <img src="/figures/africa_europe4.png" alt="1" width="80%">
</p>

---
## Euler equations
- Sod's shock tube
<p align="center">
  <img src="/figures/euler_shock_tube_2nd_density.png" alt="1" width="40%">
  <img src="/figures/euler_shock_tube_2nd_velocity.png" alt="1" width="40%">
  <img src="/figures/euler_shock_tube_2nd_energy.png" alt="1" width="40%">
  <img src="/figures/euler_shock_tube_2nd_pressure.png" alt="1" width="40%">
</p>

- Kelvin-Helmholtz instability
<p align="center">
  <img src="/figures/img1_v2.png" alt="1" width="70%">
  <img src="/figures/img2_v2.png" alt="1" width="70%">
  <img src="/figures/img3_v2.png" alt="1" width="70%">
  <img src="/figures/img4_v2.png" alt="1" width="70%">
</p>

- Rayleigh-Taylor instability
<p align="center">
  <img src="/figures/RT_instab_quad2.png" alt="1" width="70%">
</p>

- 2D Riemann problem
<p align="center">
  <img src="/figures/quadshock_time_evolv3.png" alt="1" width="90%">
</p>

- vortex shedding
<p align="center">
  <img src="/figures/vshed_img2.png" alt="1" width="80%">
</p>

- flow around an airfoil
<p align="center">
  <img src="/figures/airfoil2.png" alt="1" width="80%">
</p>

---
### Discontinuous Galerkin scalar upwind advection
1D:
<p align="center">
  <img src="/figures/DG_slope_limited_50.png" alt="1" width="48%">
  <img src="/figures/DG_scalar_smooth_30.png" alt="1" width="48%">
</p>

2D:
<p align="center">
  <img src="/figures/DG2D_slbeta2_step_Nrow30.png" alt="1" width="90%">
  <img src="/figures/DG2D_gaussian.png" alt="1" width="90%">
</p>

spatial and temporal convergence of the 2D schemes:

<p align="center">
  <img src="/figures/advected_gaussian_spatial.png" alt="1" width="48%">
  <img src="/figures/advected_gaussian_time.png" alt="1" width="48%">
</p>
These show the potential higher-order convergence using DG, however the time convergence is only first-order because of the simple forward Euler time integration used. In future work one should implement a higher-order time integration...x