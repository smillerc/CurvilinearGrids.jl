# CurvilinearGrids

[![Build Status](https://github.com/smillerc/CurvilinearGrids.jl/workflows/CI/badge.svg)](https://github.com/smillerc/CurvilinearGrids.jl/actions/workflows/CI.yml?query=branch%3Amaster) [![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)


`CurvilinearGrids.jl` is a Julia package that provides utilities for working with non-uniform curvilinear grids. The core function takes a grid and transforms it from $(x,y,z) \rightarrow (\xi,\eta,\zeta)$, where the transformed grid contains elements of unit length in each dimension. A common example of this is to use a body-fit mesh, e.g. a mesh around a wing, and transform it so that it becomes a uniform grid in $(\xi,\eta,\zeta)$. Then standard finite-difference stencils can be used on the uniform transformed grid. Below is an example of a cylindrical mesh in $(x,y)$ coordinates and the corresponding logical grid in $(\xi,\eta)$.

![Alt text](docs/image.png)

`CurvilinearGrids.jl` currently defines the `CurvilinearMesh2D` and `CurvilinearMesh3D` types. To constructe these, you need to provide the functional form of the grid, e.g. `x(i,j), y(i,j)` for a 2D mesh. 

### Example
```julia
"""Create a sphereical grid as a function of (r,θ,ϕ)"""
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
  x(ξ, η, ζ) = r(ξ) * sin(θ(η)) * cos(ϕ(ζ))
  y(ξ, η, ζ) = r(ξ) * sin(θ(η)) * sin(ϕ(ζ))
  z(ξ, η, ζ) = r(ξ) * cos(θ(η))

  # Provide the 1-argument version as well, so ForwardDiff can
  # compute derivatives of f in the form of f(x)
  x((ξ, η, ζ)) = x(ξ, η, ζ)
  y((ξ, η, ζ)) = y(ξ, η, ζ)
  z((ξ, η, ζ)) = z(ξ, η, ζ)

  return (x, y, z)
end

ni, nj, nk = (5, 9, 11) # number of nodes/vertices in each dimension
nhalo = 2 # halo cells needed for stencils (can be set to 0)

# Obtain the x, y, and z coordinate functions
x, y, z = sphere_grid(ni, nj, nk)

# Create the mesh
mesh = CurvilinearMesh3D(x, y, z, (ni, nj, nk), nhalo)
```
## Exported Functions

The API is still a work-in-progress, but for the moment, these functions are exported:

Here `idx` can be a `Tuple` or `CartesianIndex`, and mesh is an `AbstractCurvilinearMesh`.

- `coord(mesh, idx)`: Get the $(x,y,z)$ coordinates at index `idx`. This can be 1, 2, or 3D.
- `centroid(mesh, idx)`:  Get the $(x,y,z)$ coordinates of the cell centroid at cell index `idx`. This can be 1, 2, or 3D.
- `coords(mesh)` Get the: array of coordinates for the entire mesh (typically for writing to a .vtk for example)
- `centroids(mesh)` Get the: array of centroid coordinates for the entire mesh (typically for writing to a `.vtk` file)
- `metrics(mesh, idx)`: 
- `jacobian(mesh, idx)`: 
- `jacobian_matrix(mesh, idx)`: 

## Grid Metrics

When solving equations such as the Navier-Stokes in transformed form (in $\xi,\eta,\zeta$), you need to include the grid metric terms. Providing these for the grid is the primary objective of `CurvilinearGrids.jl`. These conservative grid metrics satisfy the Geometric Conservation Law [(Thomas & Lombard 1979)](https://doi.org/10.2514/3.61273)

<!-- 
$$
\hat{\xi}_x = (y_\eta z)_\zeta − (y_\zeta z)_\eta\\
\hat{\xi}_y = (z_\eta x)_\zeta − (z_\zeta x)_\eta\\
\hat{\xi}_z = (x_\eta y)_\zeta − (x_\zeta y)_\eta\\
$$

The subscript denotes a partial derivative, so $\xi_x = \partial \xi / \partial x$. 
Jacobian matrices of transformation
Inverse transformation $T^{-1}$: $(\xi,\eta,\zeta) \rightarrow (x,y,z)$: $
\begin{bmatrix}
x_\xi & y_\xi & z_\xi \\
x_\eta & y_\eta & z_\eta \\
x_\zeta & y_\zeta & z_\zeta
\end{bmatrix}
$

Forward transformation $T$: $(x,y,z) \rightarrow (\xi,\eta,\zeta)$ : $\begin{bmatrix}
\xi_x & \eta_x & \zeta_x \\
\xi_y & \eta_y & \zeta_y \\
\xi_z & \eta_z & \zeta_z
\end{bmatrix}
$ -->
