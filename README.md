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

```cpp
        // calculate original mesh
        VoronoiMesh initial_vmesh(pts);
        initial_vmesh.do_point_insertion();
        
        // do multiple iterations of lloyds algorithm
        for (int i = 0; i<lloyd_iterations; i++) {

            // calculate centroids
            vector<Point> centroids;
            centroids.reserve(initial_vmesh.vcells.size());
            for (int i = 0; i<initial_vmesh.vcells.size(); i++) {
                centroids.push_back(initial_vmesh.vcells[i].get_centroid());
            }

            // replace original mesh seeds with centroids and calculate mesh again
            initial_vmesh = VoronoiMesh(centroids);
            initial_vmesh.do_point_insertion();

        }
```

```cpp
Point VoronoiCell::get_centroid() {

    double A = get_area();

    double sum_x = 0;
    double sum_y = 0;

    for (int i = 0; i<verticies.size(); i++) {
        sum_x += (verticies[i].x + verticies[(i+1)%verticies.size()].x)*(verticies[i].x * verticies[(i+1)%verticies.size()].y - verticies[(i+1)%verticies.size()].x * verticies[i].y);
        sum_y += (verticies[i].y + verticies[(i+1)%verticies.size()].y)*(verticies[i].x * verticies[(i+1)%verticies.size()].y - verticies[(i+1)%verticies.size()].x * verticies[i].y);
    }

    double C_x = -sum_x/(6*A);
    double C_y = -sum_y/(6*A);

    return Point(C_x, C_y);
}
```

<p align="center">
  <img src="/figures/0L1_error_over_time_up_right_comp_vmesh.png" alt="0L1_error_over_time_up_right_comp_vmesh" width="45%">
  <img src="/figures/0L1_error_over_time_lloyd.png" alt="0L1_error_over_time_up_right_comp_cmesh" width="45%">
</p>

---
### Repeating Boundary conditions

#### Cartesian

<p align="center">
  <img src="/figures/0cartesian_repeating_boundary.gif" alt="0cartesian_repeating_boundary" width="45%">
  <img src="/figures/0delta_Q_total_cartesian_repeating_boundary.png" alt="0delta_Q_total_repeating_boundary" width="45%">
</p>

```cpp
for (int j = 0; j<cells[i].edges.size(); j++) {
            if (cells[i].edges[j].is_boundary == false) {
                if (j == 0) {
                    cells[i].edges[j].neighbour = &cells[i - (i%n_hor) + ((i+n_hor - 1)%n_hor)];
                } else if (j == 1) {
                    cells[i].edges[j].neighbour = &cells[((((i - (i%n_hor))/n_hor)+1)%n_vert)*n_hor + (i%n_hor)];
                } else if (j == 2) {
                    cells[i].edges[j].neighbour = &cells[i - (i%n_hor) + ((i + 1)%n_hor)];
                } else if (j == 3) {
                    cells[i].edges[j].neighbour = &cells[((((i - (i%n_hor))/n_hor)+n_vert - 1)%n_vert)*n_hor + (i%n_hor)];
                }
                
            }
        }
```

#### Voronoi
<p align="center">
  <img src="/figures/0voronoi_repeating_boundary_low_res.gif" alt="0voronoi_repeating_boundary_low_res" width="45%">
  <img src="/figures/0voronoi_repeating_boundary_high_res.gif" alt="0voronoi_repeating_boundary_high_res" width="45%">
</p>
<p align="center">
  <img src="/figures/0delta_Q_total_vmesh_repeating_boundary.png" alt="0delta_Q_total_vmesh_repeating_boundary" width="45%">
</p>

```cpp
// preprocessing for repeating boundary conditions
    vector<Point> points_plus_ghost;
    int initial_pts_size = pts.size();
    if (repeating) {
        points_plus_ghost.reserve(pts.size() * 9);

        // put points into (middle/middle) block by shrinking them by a factor of 3
        for (int i = 0; i<pts.size(); i++) {
            points_plus_ghost.emplace_back((pts[i].x/3.0) + 1.0/3.0, (pts[i].y/3.0) + 1.0/3.0);
        }

        // add the same shrinked points again but shifted in all other 8 third blocks (up/middle/down, left/middle/right)
        vector<double> pos_X = {0., 1., 2., 0., 2., 0., 1., 2.};
        vector<double> pos_Y = {0., 0., 0., 1., 1., 2., 2., 2.};
        for (int i = 0; i<8; i++) {
            for (int j = 0; j<pts.size(); j++) {
                points_plus_ghost.emplace_back((pts[j].x/3.0) + pos_X[i] * 1.0/3.0, (pts[j].y/3.0) + pos_Y[i] * 1.0/3.0);
            }
        }
        
        // replace pts with pts + additional ghost cells (eg 8 times the pts all around)
        pts = points_plus_ghost;
    }

    // generate vmesh
    VoronoiMesh vmesh(pts);
    vmesh.do_point_insertion();

    // loop through all cells (of initial pts vector) to set everything but neighbour relations
    for (int i = 0; i<initial_pts_size; i++) {
      //...
    }
```
```cpp
neigbour_index = (vmesh.vcells[i].edges[j].index2)%initial_pts_size;
```
plus additional rescaling of quantities after all that

-  Question: What about L1 error calculation with repeated boundary conditions. Is this something we really need? Idk of a really simple way to do this yet



---
### other
just another example of using loyd's algorithm + repeating boundary conditions
<p align="center">
  <img src="/figures/0adv_vmesh_rep_bound_grid_no_loyd.gif" alt="0adv_vmesh_rep_bound_grid_no_loyd" width="45%">
  <img src="/figures/0adv_vmesh_rep_bound_grid_loyd.gif" alt="0adv_vmesh_rep_bound_grid_loyd" width="45%">
</p>

---
### Shallow Water Equations

conservative form:
$${\frac {\partial (\rho h )}{\partial t}}+{\frac {\partial (\rho h u)}{\partial x}}+{\frac {\partial (\rho h v)}{\partial y}}=0 $$
$${\frac {\partial (\rho h u)}{\partial t}}+{\frac {\partial }{\partial x}}\left(\rho h u^{2}+{\frac {1}{2}}\rho gh ^{2}\right)+{\frac {\partial (\rho h uv)}{\partial y}}=0$$ 
$${\frac {\partial (\rho h v)}{\partial t}}+{\frac {\partial }{\partial y}}\left(\rho h v^{2}+{\frac {1}{2}}\rho gh ^{2}\right)+{\frac {\partial (\rho h uv)}{\partial x}}=0 $$

with $\rho =$ density, h = fluid column height, (u, v) =  velocity averaged over column, g = gravitational acceleation. Assumptions: horizontal bed, neglible coriolis, friction, viscosity + wavelength >> water depth.

### 1D cartesian

$ v := 0,  \partial_y \to 0 $

$${\frac {\partial (\rho h )}{\partial t}}+{\frac {\partial (\rho h u)}{\partial x}}=0 $$
$${\frac {\partial (\rho h u)}{\partial t}}+{\frac {\partial }{\partial x}}\left(\rho h u^{2}+{\frac {1}{2}}\rho gh ^{2}\right)=0$$ 

Or written in an alternative way using
$$
U = \begin{bmatrix} h \\ hu\end{bmatrix},\;\; F = \begin{bmatrix} hu \\ hu^2 + \frac{1}{2}gh^2\end{bmatrix}
$$
$$
\frac{\partial U}{\partial t} + \frac{\partial F}{\partial x} = 0
$$
FV Update Scheme using Lax-Friedrichs Flux. Idk how upwind should work here?
$$
U_i^{n+1} = U_i^n - \frac{\Delta t}{A} \biggl[F_{i-\frac{1}{2}}^n l_y - F_{i+\frac{1}{2}}^n l_y \biggr]$$
$$
F_{i-\frac{1}{2}}^n = \frac{1}{2} \biggl[ F_{i-1}^n + F_i^n\biggr] - \frac{l_x}{2\Delta t} \biggl[U_i^n - U_{i-1}^n\biggr]
$$

#### Example 1

<p align="center">
  <img src="/figures/1_left_swe_2D.gif" alt="1left_swe_2D" width="45%">
  <img src="/figures/1_left_swe_2D_vel.gif" alt="1_left_swe_2D_vel" width="45%">
</p>
<p align="center">
  <img src="/figures/1_left_swe_1D.gif" alt="1_left_swe_1D" width="45%">
  <img src="/figures/1_left_swe_1D_vel.gif" alt="0adv_vmesh_rep_bound_grid_loyd" width="45%">
</p>

#### Example 2

<p align="center">
  <img src="/figures/1_mid_swe_2D.gif" alt="1_mid_swe_2D" width="45%">
  <img src="/figures/1_mid_swe_2D_vel.gif" alt="1_mid_swe_2D_vel" width="45%">
</p>
<p align="center">
  <img src="/figures/1_mid_swe_1D.gif" alt="0adv_vmesh_rep_bound_grid_no_loyd" width="45%">
  <img src="/figures/1_mid_swe_1D_vel.gif" alt="1_mid_swe_1D_vel" width="45%">
</p>
