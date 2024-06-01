# hydro_bsc_project
project to eventually do hydrodynamics on different meshes, with maybe different physics and different solvers. still work in progress...

## Changes in last update
---
- General Code Cleanup
  - new folder structure
  - summed up advection solvers into one function
  - added support for negative $v_x$, $v_y$ in the cartesian advection
  - made visualization.ipynb more function oriented

- Initial conditions and L1-error estimates now defined using functions (analytic solutions)

---
### logQ plot to verify Q_diff plot

<p align="center">
  <img src="/figures/0hot_animation2D.gif" alt="0hot_animation2D" width="45%">
  <img src="/figures/0hot_loganimation2D.gif" alt="0hot_loganimation2D" width="45%">
</p>
<p align="center">
  <img src="/figures/0Q_conservation_with_analytic_line.png" alt="0hot_animation2D" width="60%">
</p>



---
### L1 error for 1D stepfunc advection on Vmesh

<p align="center">
  <img src="/figures/0L1_over_time_1D_adv.png" alt="0L1_over_time_1D_adv" width="45%">
  <img src="/figures/0L1_error_over_N_1D_adv.png" alt="delta_Q_total_float_precision" width="45%">
</p>

---
### L1 error for 2D circle advection on Vmesh

<p align="center">
  <img src="/figures/0L1_over_time2D_Vmesh_Circle.png" alt="0L1_over_time2D_Vmesh_Circle" width="45%">
  <img src="/figures/0L1_error_over_N_2D_Vmesh_Circle.png" alt="delta_Q_total_float_precision" width="45%">
</p>

---
### Angle dependence of L1 error on Voronoi vs. on Cartesian
<p align="center">
  <img src="/figures/0animation_vmesh_circle_upright.gif" alt="0animation_vmesh_circle_upright" width="45%">
  <img src="/figures/0animation_vmesh_circle_right.gif" alt="0animation_vmesh_circle_right" width="45%">
</p>

<p align="center">
  <img src="/figures/0animation_cmesh_circle_upright.gif" alt="0animation_cmesh_circle_upright" width="45%">
  <img src="/figures/0animation_cmesh_circle_right.gif" alt="0animation_cmesh_circle_right" width="45%">
</p>
<p align="center">
  <img src="/figures/0L1_error_over_time_up_right_comp_vmesh.png" alt="0L1_error_over_time_up_right_comp_vmesh" width="45%">
  <img src="/figures/0L1_error_over_time_up_right_comp_cmesh.png" alt="0L1_error_over_time_up_right_comp_cmesh" width="45%">
</p>
Higher angle dependence for Cartesian Mesh, error 45deg approx equal, for 0deg cmesh has lower error than vmesh

---
### Lloyd's Algorithm

<p align="center">
  <img src="/figures/0lloyds_algorithm.gif" alt="0animation_vmesh_circle_upright" width="60%">
</p>
<p align="center">
  <img src="/figures/0L1_error_over_time_up_right_comp_vmesh.png" alt="0L1_error_over_time_up_right_comp_vmesh" width="45%">
  <img src="/figures/0L1_error_over_time_lloyd.png" alt="0L1_error_over_time_up_right_comp_cmesh" width="45%">
</p>

---
### Repeating Boundary conditions (Cartesian)

<p align="center">
  <img src="/figures/0cartesian_repeating_boundary.gif" alt="0cartesian_repeating_boundary" width="45%">
  <img src="/figures/0delta_Q_total_repeating_boundary.png" alt="0delta_Q_total_repeating_boundary" width="45%">
</p>


Further ideas still to follow:
- Repeating Boundary Conditions on Vmesh/Cmesh
- Eventually Shallow Water Equation