# hydro_bsc_project
project to eventually do hydrodynamics on different meshes, with maybe different physics and different solvers. still work in progress...

## Changes in last update
---

- Another look at boundary conditions
- verification of HLL for dam break (diagonal/normal, voronoi/cartesian)
- verification of HLL for other initial conditions

- MUSCL scheme (work in progress)

---
### Another look at boundary conditions

Reflective boundary conditions:

For the flow of neighbour cell just set 

<p align="center">
  <img src="/figures/image_10.png" alt="6_2D_voronoi_low_res_50" width="30%">
</p>

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
### Second Order? (MUSCL)

Tried to implement second order muscl scheme following [arepo paper](https://wwwmpa.mpa-garching.mpg.de/~volker/arepo/arepo_paper.pdf) and Fundamentals of Simulation Methods script.

<p align='center'>
  <img src="/figures/image0.png" alt="6_2D_voronoi_low_res_50" width="45%">
</p>


We start by calculating the gradients for all components $\phi$ of $\vec{U}$.

<p align='center'>
  <img src="/figures/image2.png" alt="6_2D_voronoi_low_res_50" width="45%">
</p>

to those we apply slope limiters according to

<p align='center'>
  <img src="/figures/image3.png" alt="6_2D_voronoi_low_res_50" width="60%">
</p>

and then use the slope limited gradients to linearly extrapolate to the cell boundary (accounting also for half a time step).

<p align='center'>
  <img src="/figures/image1.png" alt="6_2D_voronoi_low_res_50" width="45%">
</p>

finally we use these states and our HLL Riemann solver to obtain an updated $U_i^{(n+1)}$

<p align='center'>
  <img src="/figures/image4.png" alt="6_2D_voronoi_low_res_50" width="33%">
</p>

---

### 1D cartesian 1st order against 2nd order scheme

First of all one can see that the 2nd-order approach is way more accurate (by roughly a magnitude) then the 1st-order approach. The L1-error also scales better but both algorithms apparently fail to provide real 1st or 2nd order scaling?? Reasons for that? 

<p align='center'>
  <img src="/figures/9first_second_order_dam_break.gif" alt="6_2D_voronoi_low_res_50" width="45%">
  <img src="/figures/9L1_error_over_N_problematic.png" alt="6_2D_voronoi_low_res_50" width="45%">
</p>

---

Here just to compare a version without the flux limiter.

<p align='center'>
  <img src="/figures/9second_order_flux_limiter.gif" alt="6_2D_voronoi_low_res_50" width="45%">
</p>

---
### still problems in 2D voronoi?
don't know the error yet

<p align='center'>
  <img src="/figures/9second_order_vmesh_cmesh.gif" alt="6_2D_voronoi_low_res_50" width="45%">
  <img src="/figures/9wierd_stuff_because_flux_limiter_not_tvd.gif" alt="6_2D_voronoi_low_res_50" width="45%">
</p>

---

### cartesian 2D example

Here one can see that at same grid resolution the 2nd order algorithm is way more accurate. But i wouldn't trust this result here as well until the problem for 2D voronoi is found

<p align='center'>
  <img src="/figures/91st_order_2D_cartesian.gif" alt="6_2D_voronoi_low_res_50" width="45%">
  <img src="/figures/92nd_order_2D_cartesian.gif" alt="6_2D_voronoi_low_res_50" width="45%">
</p>

---
### todo
- get 2D to work (probably a bug fix somewhere)
- solve scaling issues (How?)
- code cleanup
- then move on to euler

Because of exam preperation i won't be able to do much work in the next two weeks. I suggest we skip next week since i wont have enough progress by then. The week after that i have an exam somewhere tuesday morning. Think its best to just continue working on project after that.
