# module CubedSphereMesh
using WriteVTK

# module CubedSphereMesh

# export CubeFace,
#        build_core_block_cube,
#        build_shell_block_conformal,
#        build_cubed_sphere_with_conformal_core

# ============================================================
# Face enumeration
# ============================================================

@enum CubeFace begin
  FaceXp   # +X
  FaceXm   # -X
  FaceYp   # +Y
  FaceYm   # -Y
  FaceZp   # +Z
  FaceZm   # -Z
end

# ============================================================
# Unit-sphere direction for each face (gnomonic form, but using X,Y directly)
# Here X,Y are the face-plane coordinates in [-1,1], not angles.
# ============================================================

@inline function cubed_sphere_direction(face::CubeFace, X, Y)
  s = inv(sqrt(1 + X^2 + Y^2))

  if face === FaceXp
    return s, X * s, Y * s
  elseif face === FaceXm
    return -s, -X * s, Y * s
  elseif face === FaceYp
    return -X * s, s, Y * s
  elseif face === FaceYm
    return X * s, -s, Y * s
  elseif face === FaceZp
    return -X * s, Y * s, s
  elseif face === FaceZm
    return X * s, Y * s, -s
  else
    error("Unknown CubeFace")
  end
end

# ============================================================
# Cartesian -> spherical helper (no allocations)
# ============================================================

@inline function cart2sph(x, y, z)
  r = sqrt(x * x + y * y + z * z)
  if r == 0.0
    return 0.0, 0.0, 0.0
  end
  θ = acos(clamp(z / r, -1.0, 1.0))     # polar angle
  ϕ = atan(y, x)                        # azimuth
  return r, θ, ϕ
end

# ============================================================
# Central core: a single cube block [-a,a]^3
# ============================================================

function build_core_block_cube(
  ξ::AbstractVector,   # [-1,1]
  η::AbstractVector,   # [-1,1]
  ζ::AbstractVector,   # [-1,1]
  a::Real,              # cube half-width
)
  Nx, Ny, Nz = length(ξ), length(η), length(ζ)

  x = Array{Float64}(undef, Nx, Ny, Nz)
  y = Array{Float64}(undef, Nx, Ny, Nz)
  z = Array{Float64}(undef, Nx, Ny, Nz)

  r = Array{Float64}(undef, Nx, Ny, Nz)
  θ = Array{Float64}(undef, Nx, Ny, Nz)
  ϕ = Array{Float64}(undef, Nx, Ny, Nz)

  @inbounds for k in 1:Nz, j in 1:Ny, i in 1:Nx
    xijk = a * ξ[i]
    yijk = a * η[j]
    zijk = a * ζ[k]

    x[i, j, k] = xijk
    y[i, j, k] = yijk
    z[i, j, k] = zijk

    rijk, θijk, ϕijk = cart2sph(xijk, yijk, zijk)
    r[i, j, k] = rijk
    θ[i, j, k] = θijk
    ϕ[i, j, k] = ϕijk
  end

  return (cartesian=(x=x, y=y, z=z), spherical=(r=r, θ=θ, ϕ=ϕ))
end

# ============================================================
# Shell face block: conformal inner boundary = cube face, outer boundary = sphere
#
# σ ∈ [0,1] is a blending coordinate:
#   σ=0 => point on cube surface (exact)
#   σ=1 => point on sphere of radius Rout (exact)
#
# Mapping is along the ray defined by the face unit direction n(X,Y):
#   x = λ(σ) * n
# where λ(0) is chosen so that the ray intersects the cube at the correct face.
# ============================================================

@inline function face_denom(face::CubeFace, nx, ny, nz)
  if face === FaceXp
    return nx
  elseif face === FaceXm
    return -nx
  elseif face === FaceYp
    return ny
  elseif face === FaceYm
    return -ny
  elseif face === FaceZp
    return nz
  elseif face === FaceZm
    return -nz
  else
    error("Unknown CubeFace")
  end
end

function build_shell_block_conformal(
  face::CubeFace,
  X::AbstractVector,     # face coordinate in [-1,1]
  Y::AbstractVector,     # face coordinate in [-1,1]
  σ::AbstractVector,     # blending coordinate in [0,1]
  a::Real,               # cube half-width (inner boundary)
  Rout::Real,             # outer sphere radius
)
  NX, NY, Nσ = length(X), length(Y), length(σ)

  x = Array{Float64}(undef, NX, NY, Nσ)
  y = Array{Float64}(undef, NX, NY, Nσ)
  z = Array{Float64}(undef, NX, NY, Nσ)

  r = Array{Float64}(undef, NX, NY, Nσ)
  θ = Array{Float64}(undef, NX, NY, Nσ)
  ϕ = Array{Float64}(undef, NX, NY, Nσ)

  @inbounds for k in 1:Nσ
    σk = σ[k]
    for j in 1:NY
      Yj = Y[j]
      for i in 1:NX
        Xi = X[i]

        # Unit direction from face parameter
        nx, ny, nz = cubed_sphere_direction(face, Xi, Yj)

        # Choose λ0 so that λ0*n lies exactly on the cube face for this block:
        # For example on +X face: want x=a => λ0 = a / nx (nx>0 on that face)
        d = face_denom(face, nx, ny, nz)      # positive on the chosen face
        λ0 = a / d

        # Blend from cube (λ0) to sphere (Rout)
        λ = (1 - σk) * λ0 + σk * Rout

        xijk = λ * nx
        yijk = λ * ny
        zijk = λ * nz

        x[i, j, k] = xijk
        y[i, j, k] = yijk
        z[i, j, k] = zijk

        rijk, θijk, ϕijk = cart2sph(xijk, yijk, zijk)
        r[i, j, k] = rijk
        θ[i, j, k] = θijk
        ϕ[i, j, k] = ϕijk
      end
    end
  end

  return (cartesian=(x=x, y=y, z=z), spherical=(r=r, θ=θ, ϕ=ϕ))
end

# ============================================================
# Full mesh: 1 core cube + 6 conformal shell faces
# ============================================================

function build_cubed_sphere_with_conformal_core(
  ξ::AbstractVector,
  η::AbstractVector,
  ζ::AbstractVector,
  X::AbstractVector,
  Y::AbstractVector,
  σ::AbstractVector,
  a::Real,
  Rout::Real,
)
  core = build_core_block_cube(ξ, η, ζ, a)

  faces = (FaceXp, FaceXm, FaceYp, FaceYm, FaceZp, FaceZm)
  shell = Dict{CubeFace,NamedTuple}()

  for face in faces
    shell[face] = build_shell_block_conformal(face, X, Y, σ, a, Rout)
  end

  return (core=core, shell=shell)
end

# Block 1 : Core cube
# Block 2 : +X shell face
# Block 3 : -X shell face
# Block 4 : +Y shell face
# Block 5 : -Y shell face
# Block 6 : +Z shell face
# Block 7 : -Z shell face

interface_connections = (

  # --------------------------------------------------
  # Core ↔ shell
  # --------------------------------------------------
  ((1, :ihi) => (2, :klo)),
  ((1, :ilo) => (3, :klo)),
  ((1, :jhi) => (4, :klo)),
  ((1, :jlo) => (5, :klo)),
  ((1, :khi) => (6, :klo)),
  ((1, :klo) => (7, :klo)),

  # --------------------------------------------------
  # Shell ↔ shell (edges)
  # --------------------------------------------------
  ((2, :ilo) => (4, :ihi)),
  ((2, :ihi) => (5, :ilo)),
  ((2, :jlo) => (6, :ihi)),
  ((2, :jhi) => (7, :ilo)),
  ((3, :ilo) => (5, :ihi)),
  ((3, :ihi) => (4, :ilo)),
  ((3, :jlo) => (6, :ilo)),
  ((3, :jhi) => (7, :ihi)),
  ((4, :ilo) => (3, :ihi)),
  ((4, :ihi) => (2, :ilo)),
  ((4, :jlo) => (6, :jhi)),
  ((4, :jhi) => (7, :jlo)),
  ((5, :ilo) => (2, :ihi)),
  ((5, :ihi) => (3, :ilo)),
  ((5, :jlo) => (6, :jlo)),
  ((5, :jhi) => (7, :jhi)),
  ((6, :ilo) => (3, :jlo)),
  ((6, :ihi) => (2, :jlo)),
  ((6, :jlo) => (5, :jlo)),
  ((6, :jhi) => (4, :jlo)),
  ((7, :ilo) => (2, :jhi)),
  ((7, :ihi) => (3, :jhi)),
  ((7, :jlo) => (4, :jhi)),
  ((7, :jhi) => (5, :jhi)),
)

# end # module

N = 33
Nr = 64

ξ = range(-1, 1; length=N)
η = range(-1, 1; length=N)
ζ = range(-1, 1; length=N)

# Use the SAME [-1,1] grids for shell face coordinates to be conformal
X = η
Y = ζ

σ = range(0, 1; length=Nr)

a = 0.5     # cube half-width
Rout = 2.0     # outer sphere radius

mesh = build_cubed_sphere_with_conformal_core(ξ, η, ζ, X, Y, σ, a, Rout)

# Access +X face
for (k, v) in mesh.shell
  #   coords = mesh[k]

  #   x = coords.cartesian.x
  #   y = coords.cartesian.y
  #   z = coords.cartesian.z

  #   θ = coords.spherical.θ
  #   ϕ = coords.spherical.ϕ

  @views vtk_grid("cubed_sphere_$(k)", (v.cartesian...)) do vtk
    dom = size(v.cartesian.x) .- 1
    indices = CartesianIndices(dom)

    vtk["index", VTKCellData(), component_names=["i", "j", "k"]] = (
      [idx.I[1] for idx in indices],
      [idx.I[2] for idx in indices],
      [idx.I[3] for idx in indices],
    )
  end
end

@views vtk_grid("cubed_sphere_core", (mesh.core.cartesian...)) do vtk
  dom = size(mesh.core.cartesian.x) .- 1
  indices = CartesianIndices(dom)

  vtk["index", VTKCellData(), component_names=["i", "j", "k"]] = (
    [idx.I[1] for idx in indices],
    [idx.I[2] for idx in indices],
    [idx.I[3] for idx in indices],
  )
end