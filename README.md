# hydro_bsc_project
project to eventually do hydrodynamics on different meshes, with maybe different physics and different solvers. still work in progress...

### Todos:
- further todos will be placed here

### Edges "none" or "face"
<p align="center">
  <img src="./figures/test_edgecolor_none.png" alt="none" width="45%">
  <img src="./figures/test_edgecolor_face.png" alt="face" width="45%">
</p>

- "none": with white borders leads to rendering issues for to many cells 
- "faces": no borders

### Cartesian Mesh

<p align="center">
  <img src="./figures/c_square_10k.gif" alt="square" width="30%">
  <img src="./figures/c_half_2k.gif" alt="half" width="30%">
  <img src="./figures/c_continuous_input.gif" alt="input" width="30%">
</p>

1. initial square with 10k cells
2. initial one half filled with 2k cells
3. continuous input in one cell on a 2k cell grid

### Vornoi Mesh
<p align="center">
  <img src="./figures/v_circle_2k.gif" alt="circle" width="30%">
  <img src="./figures/v_half_2k.gif" alt="half" width="30%">
  <img src="./figures/v_continuous_input.gif" alt="input" width="30%">
</p>

1. initial circle with 2k cells
2. initial one half filled with 2k cells
3. continuous input in one cell on a 2k cell grid