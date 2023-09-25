# CurvilinearGrids

[![Build Status](https://github.com/smil/CurvilinearGrids.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/smil/CurvilinearGrids.jl/actions/workflows/CI.yml?query=branch%3Amain)


`CurvilinearGrids.jl` is a Julia package that provides utilities for working with non-uniform curvilinear grids. The core function takes a grid and transforms it from $(x,y,z) \rightarrow (\xi,\eta,\zeta)$, where the transformed grid contains elements of unit length in each dimension. A common example of this is to use a body-fit mesh, e.g. a mesh around a wing, and transform it so that it becomes a uniform grid in $(\xi,\eta,\zeta)$. Then standard finite-difference stencils can be used on the uniform transformed grid.

![Alt text](docs/image.png)


```julia
  function sphere_grid(nr, ntheta, nphi)
    r0, r1 = (1, 3)
    (θ0, θ1) = deg2rad.((35, 180 - 35))
    (ϕ0, ϕ1) = deg2rad.((45, 360 - 45))

    r(ξ) = r0 + (r1 - r0) * ((ξ - 1) / (nr - 1))
    θ(η) = θ0 + (θ1 - θ0) * ((η - 1) / (ntheta - 1))
    ϕ(ζ) = ϕ0 + (ϕ1 - ϕ0) * ((ζ - 1) / (nphi - 1))

    x(ξ, η, ζ) = r(ξ) * sin(θ(η)) * cos(ϕ(ζ))
    y(ξ, η, ζ) = r(ξ) * sin(θ(η)) * sin(ϕ(ζ))
    z(ξ, η, ζ) = r(ξ) * cos(θ(η))

    x((ξ, η, ζ)) = x(ξ, η, ζ)
    y((ξ, η, ζ)) = y(ξ, η, ζ)
    z((ξ, η, ζ)) = z(ξ, η, ζ)

    return (x, y, z)
  end

  ni, nj, nk = (5, 9, 11)
  nhalo = 0
  x, y, z = sphere_grid(ni, nj, nk)
```

## Exported Functions

```julia
mesh = CurvilinearMesh3D(...)

idx = (i,j,k) # or
idx = CartesianIndex((i,j,k))

coord(mesh, idx)
centroid(mesh, idx)

coords(mesh)
centroids(mesh)

metrics(mesh, idx)
jacobian(mesh, idx)
jacobian_matrix(mesh, idx)
```


## Grid Metrics Details

Jacobian matrix of transformation
$$
\begin{bmatrix}
x_\xi & y_\xi & z_\xi \\
x_\eta & y_\eta & z_\eta \\
x_\zeta & y_\zeta & z_\zeta
\end{bmatrix}
$$

$$
\begin{bmatrix}
\xi_x & \eta_x & \zeta_x \\
\xi_y & \eta_y & \zeta_y \\
\xi_z & \eta_z & \zeta_z
\end{bmatrix}
$$

Conservative grid metrics to satisfy the Geometric Conservation Law [(Thomas & Lombard 1979)](https://doi.org/10.2514/3.61273)

$$
\hat{\xi}_x = (y_\eta z)_\zeta − (y_\zeta z)_\eta\\
\hat{\xi}_y = (z_\eta x)_\zeta − (z_\zeta x)_\eta\\
\hat{\xi}_z = (x_\eta y)_\zeta − (x_\zeta y)_\eta\\
...
$$