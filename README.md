# hydro_bsc_project
Project to do 2D hydrodynamics on Cartesian and static Voronoi meshes using finite volume and discontinuous Galerkin methods.

The project includes grid generation methods from the [vmp](https://github.com/lucas56098/voronoi_mesh_project) mesh generation code as well as solvers for the advection equation, the shallow water equations and Euler equations. The advection equation is solved in one and two dimensions using a simple FV upwind scheme, while the shallow water and Euler equations are solved using a first-order Godunov-type scheme and a second-order MUSCL-Hancock scheme with slope-limiting as in ([Springel 2010](https://arxiv.org/abs/0901.4107)). To solve the Riemann problems a HLL approximate Riemann solver is used. For discontinous Galerkin methods a one and two dimensional Cartesian scalar upwind advection scheme is implemented using monomial basis functions and Legendre polynomials following ([Schaal 2015](https://arxiv.org/abs/1506.06140)). Finally there exist various plotting and analysis routines in a python visualization toolkit. Note, that the program requires the Eigen library to run. 

This work was done as part of a bachelor thesis in physics at the Institute for theoretical Astrophysics in the Group of [Dr. Dylan Nelson](https://github.com/dnelson86), supervised by him and [Dr. Chris Byrohl](https://github.com/cbyrohl). The thesis can be accessed [here](https://github.com/lucas56098/hydro_bsc_project/blob/main/thesis.pdf). Otherwise feel free to just scroll through some of the plots produced using code.

---
### Advection
A one and two dimensional upwind scheme for the advection equation. At low resolution there is a lot of numerical diffusion.
<p align="center">
  <img src="/figures/advected_step_function_evol_res100.png" alt="1" width="40%">
  <img src="/figures/advected_gaussian_evol_res100.png" alt="1" width="40%">
</p>

---
### Mesh generation with custom internal boundaries
When generating a Voronoi mesh, we can place the seedpoints in specific ways to generate internal structures.
<p align="center">
  <img src="/figures/custom_circle.png" alt="1" width="35%">
  <img src="/figures/custom_airfoil.png" alt="1" width="35%">
</p>

---
### 2D shallow water equations
The shallow water equations can be solved using similar methods as for the euler equations and describe wave propagation, if the wavelength is way longer than the water depth.

<p align="center">
  <img src="/figures/box5.png" alt="1" width="80%">
  <img src="/figures/circle4.png" alt="1" width="80%">
  <img src="/figures/africa_europe4.png" alt="1" width="80%">
</p>

---
## Euler equations
Sod's shock tube is a classical test where the code can be compared to an analytical solution. Shown here is the second order MUSCL-scheme with slope limiting.
<p align="center">
  <img src="/figures/euler_shock_tube_2nd_density.png" alt="1" width="40%">
  <img src="/figures/euler_shock_tube_2nd_velocity.png" alt="1" width="40%">
  <img src="/figures/euler_shock_tube_2nd_energy.png" alt="1" width="40%">
  <img src="/figures/euler_shock_tube_2nd_pressure.png" alt="1" width="40%">
</p>

The Kelvin-Helmholtz instability is one of the most important fluid instabilities. The initial perturbation here is given by Voronoi mesh noise only. The size of the smallest possible eddies strongly depends on mesh resolution and the accuracy of the scheme.
<p align="center">
  <img src="/figures/img1_v2.png" alt="1" width="70%">
  <img src="/figures/img2_v2.png" alt="1" width="70%">
  <img src="/figures/img3_v2.png" alt="1" width="70%">
  <img src="/figures/img4_v2.png" alt="1" width="70%">
</p>

Another important fluid instability is the Rayleigh-Taylor instability shown here for two fluids under the influence of constant external gravity and a sinuodial perturbation.
<p align="center">
  <img src="/figures/RT_instab_quad2.png" alt="1" width="70%">
</p>

Additionally we looked at a 2D Riemann problem proposed by ([Kurganov and Tadmor, 2002](https://www.semanticscholar.org/paper/Solution-of-two%E2%80%90dimensional-Riemann-problems-for-Kurganov-Tadmor/a44da75f9a36ab879fb9073f2571801eb7bc74a3)).
<p align="center">
  <img src="/figures/quadshock_time_evolv3.png" alt="1" width="90%">
</p>

and observed vortex shedding
<p align="center">
  <img src="/figures/vshed_img2.png" alt="1" width="80%">
</p>

and the flow around an airfoil
<p align="center">
  <img src="/figures/airfoil2.png" alt="1" width="80%">
</p>

---
### Discontinuous Galerkin scalar upwind advection
In DG we propagate functions on the cells instead of cell averages as in FV. Slope limiting near discontinuities is really important here.

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

A convergence analysis shows higher order spatial convergence as expected for DG, the time convergence however is first order because of the forward Euler time integration.

<p align="center">
  <img src="/figures/advected_gaussian_spatial.png" alt="1" width="48%">
  <img src="/figures/advected_gaussian_time.png" alt="1" width="48%">
</p>
In future work one should implement a higher-order time integration... additionally extending DG to Euler equations and maybe Voronoi would be really exiting.