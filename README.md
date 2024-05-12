# hydro_bsc_project
project to eventually do hydrodynamics on different meshes, with maybe different physics and different solvers. still work in progress...

### changes in last update

-  delta_Q in faces replaced with Q in cell
-  solver class that does the updates of Q etc.
-  derived classes + virtual function for cell type handling (e.g. Q_Cell, Conway_Cell)


#### Conway's Game of Life

<p align="center">
  <img src="./figures/conway_cartesian.gif" alt="conway_cartesian" width="49%">
  <img src="./figures/conway_voronoi.gif" alt="conway_voronoi" width="49%">
</p>

Clearly does work better on cartesian mesh. System is really sensitive to rule changes

#### 1D - Solver
done by just generating "1D"-grids (eg. 45x1) and then using the same solver as before

#### 1D - Fake Diffiusion
time evolution on 1D-mesh:
<p align="center">
  <img src="./figures/1D_stepfunc.gif" alt="1Dstepfunc" width="80%">
</p>

time evolution on 2D-mesh binned and averaged to 1D:
<p align="center">
  <img src="./figures/2D_stepfunc_1Dbinned.gif" alt="2D_stepfunc_1Dbinned" width="80%">
</p>

#### Change in total Q
Plot shows difference between initial Q summed up over the mesh and the summed up Q at each simulation step. For Cartesian (left) and Voronoi (right). As one can see both are $\Delta Q < 10^{-11}$ and seem more or less random

<p align="center">
  <img src="./figures/delta_Q_total_c.png" alt="delta_Q_total_c" width="49%">
  <img src="./figures/delta_Q_total_v.png" alt="delta_Q_total_v" width="49%">
</p>

#### Advection

1st order finite difference upwind scheme for advection using

$$ u_i^{(n+1)} = u_i^{(n)} - v \frac{\Delta t}{\Delta x}\biggl[u_i^{(n)} - u_{i-1}^{(n)}\biggr] $$

with CFL timestep condition $\Delta t \leq \frac{\Delta x}{v}$

<p align="center">
  <img src="./figures/advection_anim.gif" alt="advection_anim" width="80%">
</p>

on the 1D-plots one can clearly see the numerical diffusion being more important for lower resolutions

<p align="center">
  <img src="./figures/advection1D_c_100cells_200step.gif" alt="advection1D_c_100cells_200step" width="32%">
  <img src="./figures/advection1D_c_1000cells_2000step.gif" alt="advection1D_c_1000cells_2000steps" width="32%">
  <img src="./figures/advection1D_c_10000cells_20000step.gif" alt="advection1D_c_10000cells_20000steps" width="32%">
</p>