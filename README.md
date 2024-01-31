# CurvilinearGrids

[![Build Status](https://github.com/smillerc/CurvilinearGrids.jl/workflows/CI/badge.svg)](https://github.com/smillerc/CurvilinearGrids.jl/actions/workflows/CI.yml?query=branch%3Amaster) [![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)


`CurvilinearGrids.jl` is a Julia package that provides utilities for working with non-uniform curvilinear grids. The core function takes a grid and transforms it from $(x,y,z) \rightarrow (\xi,\eta,\zeta)$, where the transformed grid contains elements of unit length in each dimension. A common example of this is to use a body-fit mesh, e.g. a mesh around a wing, and transform it so that it becomes a uniform grid in $(\xi,\eta,\zeta)$. Then standard finite-difference stencils can be used on the uniform transformed grid. Below is an example of a cylindrical mesh in $(x,y)$ coordinates and the corresponding logical grid in $(\xi,\eta)$.

![Alt text](docs/image.png)

`CurvilinearGrids.jl` currently defines the `CurvilinearGrid2D` and `CurvilinearGrid3D` types. To construct these, you need to provide the functional form of the grid, e.g. `x(i,j), y(i,j)` for a 2D mesh. 

### Example
```julia
using CurvilinearGrids

"""Create a spherical grid as a function of (r,θ,ϕ)"""
function sphere_grid(nr, ntheta, nphi)
  r0, r1 = (1, 3) # min/max radius
  (θ0, θ1) = deg2rad.((35, 180 - 35)) # min/max polar angle
  (ϕ0, ϕ1) = deg2rad.((45, 360 - 45)) # min/max azimuthal angle

  # Linear spacing in each dimension
  # Sometimes (ξ, η, ζ) is used instead of (i, j, k), depending on preference
  r(ξ) = r0 + (r1 - r0) * ((ξ - 1) / (nr - 1))
  θ(η) = θ0 + (θ1 - θ0) * ((η - 1) / (ntheta - 1))
  ϕ(ζ) = ϕ0 + (ϕ1 - ϕ0) * ((ζ - 1) / (nphi - 1))

  # simple spherical to cartesian mapping
  x(i,j,k) = r(i) * sin(θ(j)) * cos(ϕ(k))
  y(i,j,k) = r(i) * sin(θ(j)) * sin(ϕ(k))
  z(i,j,k) = r(i) * cos(θ(j))

  return (x, y, z)
end

ni, nj, nk = (5, 9, 11) # number of nodes/vertices in each dimension
nhalo = 2 # halo cells needed for stencils (can be set to 0)

# Obtain the x, y, and z coordinate functions
x, y, z = sphere_grid(ni, nj, nk)

# Create the mesh
mesh = CurvilinearGrid3D(x, y, z, (ni, nj, nk), nhalo)
```
## Exported Functions

The API is still a work-in-progress, but for the moment, these functions are exported:

Here `idx` can be a `Tuple` or `CartesianIndex`, and mesh is an `AbstractCurvilinearGrid`. **Important:** The indices provided to these functions are aware of halo regions, so the functions do the offsets for you. This is by design, since fields attached to the mesh, like density or pressure for example, _will_ have halo regions, and loops through these fields typically have pre-defined limits that constrain the loop to only work in the non-halo cells. If you don't use halo cells, just set `nhalo=0` in the constructors.

The index can be a `Tuple`, scalar `Integer`, or `CartesianIndex`.
- `coord(mesh, idx)`: Get the $(x,y,z)$ coordinates at index `idx`. This can be 1, 2, or 3D.
- `centroid(mesh, idx)`:  Get the $(x,y,z)$ coordinates of the cell centroid at cell index `idx`. This can be 1, 2, or 3D.
- `metrics(mesh, idx)`: Get the cell metric information, e.g. $\xi_x, \xi_y$, etc.
- `conservative_metrics(mesh, idx)`: Get the metrics that are consistent with the Geometric Conservation Law (GCL)
- `jacobian(mesh, idx)`: Get the determinant of the Jacobian matrix of forward transformation, $J$
- `jacobian_matrix(mesh, idx)`: Get the Jacobian matrix of the forward transformation

These functions are primarily used to get the complete set of coordinates for plotting or post-processing. These do _not_ use halo regions, since there geometry is ill-defined here.
- `coords(mesh)` Get the: array of coordinates for the entire mesh (typically for writing to a .vtk for example)
- `centroids(mesh)` Get the: array of centroid coordinates for the entire mesh (typically for writing to a `.vtk` file)
## Grid Metrics

When solving equations such as the Navier-Stokes in transformed form (in $\xi,\eta,\zeta$), you need to include the grid metric terms. Providing these for the grid is the primary objective of `CurvilinearGrids.jl`. These conservative grid metrics satisfy the Geometric Conservation Law [(Thomas & Lombard 1979)](https://doi.org/10.2514/3.61273)

$$
\hat{\xi}_x, \hat{\xi}_y,...
$$

The subscript denotes a partial derivative, so $\xi_x = \partial \xi / \partial x$. 


## Jacobian matrices of transformation

Terminology can be somewhat confusing, but the "Jacobian matrix" is the matrix of partial derivatives that describe the forward or inverse transformation, and uses a bold-face $\textbf{J}$. The "Jacobian" then refers to the determinant of the Jacobian matrix, and is the non-bolded $J$. Some authors refer to the matrix as the "Jacobi matrix" as well. See [Wikipedia](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant) for more details.

Forward transformation, or $T: (\xi,\eta,\zeta) \rightarrow (x,y,z)$. These functions are what is provided to the `CurvilinearGrid` constructors. See the included examples above and in the unit tests.

$$
\textbf{J} = 
\begin{bmatrix}
x_\xi & x_\eta & x_\zeta \\
y_\xi & y_\eta & y_\zeta \\
z_\xi & z_\eta & z_\zeta
\end{bmatrix}
$$

$$
J = \det [\textbf{J}]
$$

Inverse transformation $T^{-1}$: $(x,y,z) \rightarrow (\xi,\eta,\zeta)$ : 

$$
\textbf{J}^{-1} = 
\begin{bmatrix}
\xi_x   & \xi_y   & \xi_z   \\
\eta_x  & \eta_y  & \eta_z  \\
\zeta_x & \zeta_y & \zeta_z
\end{bmatrix}
$$

$$
J^{-1} = \det [\textbf{J}^{-1}]
$$
