# hydro_bsc_project
project to eventually do hydrodynamics on different meshes, with maybe different physics and different solvers. still work in progress...

## Changes in last update
---

- Another look at boundary conditions
- verification of HLL for dam break (diagonal/normal, voronoi/cartesian)
- verification of HLL for other initial conditions

---
### Another look at boundary conditions

Reflective boundary conditions:

For the flow of neighbour cell just set 

\[ 
\begin{bmatrix} \hat{u} \\ \hat{v} \end{bmatrix} = \begin{bmatrix} u \\ v \end{bmatrix} - 2 \cdot (\begin{bmatrix} u \\ v \end{bmatrix} \cdot \vec{n}) \cdot \vec{n}
\]

with n the normal of the face. This cancels out any flux component through the face, while keeping the flux component parallel to the face.


<p align="center">
  <img src="/figures/8_dam_break_45_deg_cartesian_reflective.gif" alt="6_2D_voronoi_low_res_50" width="45%">
  <img src="/figures/8_dam_break_45_deg_voronoi_reflective.gif" alt="6_2D_voronoi_low_res_50" width="45%">
</p>

Sink boundary. Here i thought just keeping u and v should be sufficient but apparently this "zero gradient" method also leads to some (even though way smaller) reflections. This is because for example at the top layer when just using the cell below we are systematically misjudging the cell that should be there for the ideal sink case. But i guess that is really hard to fix? I at least have no idea how to do this.
<p align='center'>
  <img src="/figures/8_dam_break_45_deg_cartesian_sink.gif" alt="6_2D_voronoi_low_res_50" width="45%">
  <img src="/figures/8_dam_break_45_deg_voronoi_sink.gif" alt="6_2D_voronoi_low_res_50" width="45%">
</p>

---
### Testing HLL on Dam Break 
normal dam break into x direction with voronoi and cartesian leads to expected results:
<p align='center'>
  <img src="/figures/8_dam_break_0_deg_v_c_analytic.gif" alt="6_2D_voronoi_low_res_50" width="60%">
</p>

when we do the diagonal dam break however we see effects of the "sink-reflections" (left). When we just look at a the points very close to the diagonal (dist < 0.05) we get the a filtered version (right) that shows the expected behaviour until at the end the boundary effects also reach the diagonal.
<p align='center'>
  <img src="/figures/8_dam_break_45_deg_v_c_analytic.gif" alt="6_2D_voronoi_low_res_50" width="45%">
  <img src="/figures/8_dam_break_45_deg_v_c_analytic_filtered.gif" alt="6_2D_voronoi_low_res_50" width="45%">
</p>

So i guess it works apart from that sink boundary.

--- 
### L1 error of 1D dam break 

As a function of time and over N

<p align='center'>
  <img src="/figures/8_L1_over_time.png" alt="6_2D_voronoi_low_res_50" width="45%">
  <img src="/figures/8_L1_error_over_N.png" alt="6_2D_voronoi_low_res_50" width="45%">
</p>

N_x here is the number of seedpoints in x direction. Therefore N_total = N_x * N_x. Not quite the 1/x dependence here?? Whatever that means

---
### Another analytical solution?
Have not found another analytical solution of a completely different type (e.g. no dam break) without changing dry/wet boundaries. Therefore i used a 1D Python [SWE_Solver](https://github.com/DanielCortild/SWE-Solver) to further counter check. Also slightly modified his plot_SWE function to return data in the form we need it here.

```pip install SWE_Solver```

With this solver we could counter check basically any inital condition. Lets take a look at a inital gaussian peak.

<p align='center'>
  <img src="/figures/8_compare_to_SWE_python_low_res.gif" alt="6_2D_voronoi_low_res_50" width="45%">
  <img src="/figures/8_compare_to_SWE_python_high_res.gif" alt="6_2D_voronoi_low_res_50" width="45%">
</p>

For low resolution and same N the solutions differ with regard to numerical diffusion (left) which is more or less expected since we use different solvers. For higher resultion (e.g. N = 1000, right) the **solutions converge to each other**. 

---
### Second Order (MUSCL)
next step i guess

---
### todo
next weeks:
- switch to euler equations and implement a muscl scheme there with hllc solver