# hydro_bsc_project
project to eventually do hydrodynamics on different meshes, with maybe different physics and different solvers. still work in progress...

## Changes in last update
---
Last Week
- implemented roes approximate riemann solver plus additional diffusion
- dam break problem 1D with analytical solution

This Week
- comparison animation for 1D dam break problem
- turns out roe solver does not converge to analytical solution. Don't know where the mistake here is
- Implemented HLL solver does the job though. (Reflective Boundary condition needs to be adapted). HLL quite diffusife as expected

---
### 2D Shallow water equations with Roe's approximate Riemann Solver

From last time we had:

$$ \vec{U}_i^{n+1} = \vec{U}_i^n + \frac{\Delta t}{A} \sum_{j \in \partial C_i}{\vec{F}_{i, j} \cdot l_{i,j}} $$

However it was unclear how $\vec{F}_{i,j}$ should look. Roe's approximate Riemann solver estimates the flux as

$$ \vec{F}_{roe} = \frac{1}{2} \Bigl( (f_in_x + g_in_y) + (f_jn_x + g_jn_y) - |A|(U_i - U_j) \Bigr)$$

with $A = \frac{\partial F}{\partial U} |_{roe}$ the flux jacobian evaluated at the roe average state. The roe average state is

$$ \tilde{h} = \sqrt{h_i h_j}, \;\;\;\tilde{c} = \sqrt{\frac{1}{2} \bigl(gh_i + gh_j  \bigr)} $$
$$ \tilde{u} = \frac{u_i \sqrt{h_i} + u_j\sqrt{h_j}}{\sqrt{h_i} + \sqrt{h_j}}, \;\;\; \tilde{v} = \frac{v_i \sqrt{h_i} + v_j\sqrt{h_j}}{\sqrt{h_i} + \sqrt{h_j}} $$

Calculating |A| leads to 

$$|A| = \lambda_1 \cdot \lambda_2 \cdot \lambda_3$$

with $\lambda_1 = \tilde{u}n_x + \tilde{v}n_y,\;\;\;  \lambda_2 = \tilde{u}n_x + \tilde{v}n_y - \tilde{c} ,\;\;\; \lambda_3 = \tilde{u}n_x + \tilde{v}n_y + \tilde{c}$

with that flux estimate we now can use the FV-Scheme.

Only problem though is that it does not work apparently. Idk. Probably have to rewite it completely. Rather focus on HLLC instead or similar?

### 2D SWE using roes solver
looks quite nice, not physically correct though
<p align="center">
  <img src="/figures/6_2D_cartesian_low_res_50.gif" alt="6_2D_cartesian_low_res_50" width="45%">
  <img src="/figures/6_2D_voronoi_low_res_50.gif" alt="6_2D_voronoi_low_res_50" width="45%">
</p>

<p align="center">
  <img src="/figures/7_roe_looks_bad.gif" alt="6_2D_cartesian_low_res_50" width="65%">
</p>

Does not converge to solution with higher N...

--- 
### HLL Solver

\[
S_L = \min(u_L - \sqrt{gh_L}, u_R - \sqrt{gh_R})
\]

\[
S_R = \max(u_L + \sqrt{gh_L}, u_R + \sqrt{gh_R})
\]

\[
F_{i+\frac{1}{2}} =
\begin{cases} 
F_L & \text{if } S_L \geq 0 \\
F_R & \text{if } S_R \leq 0 \\
\frac{S_R F_L - S_L F_R + S_L S_R (U_R - U_L)}{S_R - S_L} & \text{if } S_L < 0 < S_R 
\end{cases}
\]

Beforehand spent some time thinking which side is the left/right side, because Cell-i/Cell-j is not always left or right. Otherwise flux wouldnt be conserved. (Basically you cant switch L and R states while expecting just a sign flip in the flux because velocities do not flip)

at the moment just sink or repeating boundary conditions... have to further check why this is the case
<p align="center">
  <img src="/figures/7_gauss_cartesian_sink_hll.gif" alt="6_2D_cartesian_low_res_50" width="45%">
  <img src="/figures/7_gauss_voronoi_sink.gif" alt="6_2D_voronoi_low_res_50" width="45%">
</p>

but generally the dam break problem looks promising

<p align="center">
  <img src="/figures/7_hll_1D_dam_break_sink.gif" alt="6_2D_voronoi_low_res_50" width="65%">
</p>

---
still to do:

- get reflective boundaries back to work
- rewrite Roe solver, alternatively just work on HLLC?

further ideas:
- What about second order? (MUSCL-Hancock?), using slope limiters here?
- Then decision on whether to progress with non constant ocean floor (e.g. source terms), start with Euler Equations or start DG already? I think Euler FV would make most sense but not sure