module CurvilinearGridsMakieExt

using CurvilinearGrids
using Makie: Makie
import Makie: @recipe, mesh, mesh!, wireframe, wireframe!

const _GridLike = Union{
  CurvilinearGrids.AbstractUnifiedGrid,CurvilinearGrids.AbstractCurvilinearGrid
}
const _WireframeTarget = Union{_GridLike,CurvilinearGrids.MultiBlockMesh}
const _SurfaceTarget = _WireframeTarget
const GB = Makie.GeometryBasics

const _DEFAULT_MULTIBLOCK_COLORS = (
  :dodgerblue,
  :tomato,
  :seagreen,
  :orange,
  :mediumpurple,
  :goldenrod,
  :deeppink,
  :slategray,
  :teal,
  :brown,
)

@recipe(CurvilinearGridWireframe, mesh) do _scene
  Makie.Theme(; color=:black, linewidth=1.0, marker=nothing, markersize=8, colors=nothing)
end

function Makie.preferred_axis_type(plot::CurvilinearGridWireframe)
  mesh = Makie.to_value(plot[1])
  return _wireframe_dimension(mesh) == 3 ? Makie.Axis3 : Makie.Axis
end

function Makie.plot!(plot::CurvilinearGridWireframe)
  mesh = Makie.to_value(plot[1])
  color = _plot_attr(plot, :color, :black)
  linewidth = _plot_attr(plot, :linewidth, 1.0)
  marker = _plot_attr(plot, :marker, nothing)
  markersize = _plot_attr(plot, :markersize, 8)
  colors = _plot_attr(plot, :colors, nothing)

  _plot_wireframe!(
    plot,
    mesh;
    color=color,
    linewidth=linewidth,
    marker=marker,
    markersize=markersize,
    colors=colors,
  )
  return plot
end

wireframe(mesh::_WireframeTarget; kwargs...) = curvilineargridwireframe(mesh; kwargs...)
function wireframe!(ax, mesh::_WireframeTarget; kwargs...)
  curvilineargridwireframe!(ax, mesh; kwargs...)
end

@recipe(CurvilinearGridSurface, mesh) do _scene
  Makie.Theme(; color=:lightskyblue3, transparency=false, colors=nothing)
end

function Makie.preferred_axis_type(plot::CurvilinearGridSurface)
  mesh = Makie.to_value(plot[1])
  return _wireframe_dimension(mesh) == 3 ? Makie.Axis3 : Makie.Axis
end

function Makie.plot!(plot::CurvilinearGridSurface)
  mesh = Makie.to_value(plot[1])
  color = _plot_attr(plot, :color, :lightskyblue3)
  shading = _plot_attr(plot, :shading, true)
  transparency = _plot_attr(plot, :transparency, false)
  colors = _plot_attr(plot, :colors, nothing)

  _plot_surface!(
    plot, mesh; color=color, shading=shading, transparency=transparency, colors=colors
  )
  return plot
end

mesh(mesh::_SurfaceTarget; kwargs...) = curvilineargridsurface(mesh; kwargs...)
mesh!(ax, mesh::_SurfaceTarget; kwargs...) = curvilineargridsurface!(ax, mesh; kwargs...)

function _plot_wireframe!(
  scene,
  mb::CurvilinearGrids.MultiBlockMesh;
  color=:black,
  linewidth::Real=1.0,
  marker=nothing,
  markersize::Real=8,
  colors=nothing,
)
  block_colors = _resolve_block_colors(length(mb.blocks), colors)
  for i in eachindex(mb.blocks)
    _plot_wireframe!(
      scene,
      mb.blocks[i];
      color=block_colors[i],
      linewidth=linewidth,
      marker=marker,
      markersize=markersize,
      colors=nothing,
    )
  end
  return nothing
end

function _plot_wireframe!(
  scene,
  grid::_GridLike;
  color=:black,
  linewidth::Real=1.0,
  marker=nothing,
  markersize::Real=8,
  colors=nothing,
)
  cart = _cartesian_node_coordinates(grid)
  N = length(cart)

  if N == 1
    x = cart[1]
    y = zero.(x)
    Makie.lines!(scene, x, y; color=color, linewidth=linewidth)
    if marker !== nothing
      Makie.scatter!(scene, x, y; color=color, marker=marker, markersize=markersize)
    end
  elseif N == 2
    x, y = cart
    _plot_wireframe_2d!(scene, x, y; color=color, linewidth=linewidth)
    if marker !== nothing
      Makie.scatter!(
        scene, vec(x), vec(y); color=color, marker=marker, markersize=markersize
      )
    end
  elseif N == 3
    x, y, z = cart
    _plot_wireframe_3d_boundaries!(scene, x, y, z; color=color, linewidth=linewidth)
    if marker !== nothing
      xb, yb, zb = _boundary_point_cloud(x, y, z)
      Makie.scatter!(scene, xb, yb, zb; color=color, marker=marker, markersize=markersize)
    end
  else
    throw(ArgumentError("Unsupported mesh dimension $N for wireframe plotting."))
  end

  return nothing
end

function _plot_wireframe_2d!(scene, x::AbstractMatrix, y::AbstractMatrix; color, linewidth)
  for i in axes(x, 1)
    Makie.lines!(
      scene, vec(@view(x[i, :])), vec(@view(y[i, :])); color=color, linewidth=linewidth
    )
  end
  for j in axes(x, 2)
    Makie.lines!(
      scene, vec(@view(x[:, j])), vec(@view(y[:, j])); color=color, linewidth=linewidth
    )
  end
  return nothing
end

function _plot_wireframe_3d_boundaries!(
  scene,
  x::AbstractArray{<:Real,3},
  y::AbstractArray{<:Real,3},
  z::AbstractArray{<:Real,3};
  color,
  linewidth,
)
  i0, i1 = firstindex(x, 1), lastindex(x, 1)
  j0, j1 = firstindex(x, 2), lastindex(x, 2)
  k0, k1 = firstindex(x, 3), lastindex(x, 3)

  for i in (i0, i1)
    for j in axes(x, 2)
      Makie.lines!(
        scene,
        vec(@view(x[i, j, :])),
        vec(@view(y[i, j, :])),
        vec(@view(z[i, j, :]));
        color=color,
        linewidth=linewidth,
      )
    end
    for k in axes(x, 3)
      Makie.lines!(
        scene,
        vec(@view(x[i, :, k])),
        vec(@view(y[i, :, k])),
        vec(@view(z[i, :, k]));
        color=color,
        linewidth=linewidth,
      )
    end
  end

  for j in (j0, j1)
    for i in axes(x, 1)
      Makie.lines!(
        scene,
        vec(@view(x[i, j, :])),
        vec(@view(y[i, j, :])),
        vec(@view(z[i, j, :]));
        color=color,
        linewidth=linewidth,
      )
    end
    for k in axes(x, 3)
      Makie.lines!(
        scene,
        vec(@view(x[:, j, k])),
        vec(@view(y[:, j, k])),
        vec(@view(z[:, j, k]));
        color=color,
        linewidth=linewidth,
      )
    end
  end

  for k in (k0, k1)
    for i in axes(x, 1)
      Makie.lines!(
        scene,
        vec(@view(x[i, :, k])),
        vec(@view(y[i, :, k])),
        vec(@view(z[i, :, k]));
        color=color,
        linewidth=linewidth,
      )
    end
    for j in axes(x, 2)
      Makie.lines!(
        scene,
        vec(@view(x[:, j, k])),
        vec(@view(y[:, j, k])),
        vec(@view(z[:, j, k]));
        color=color,
        linewidth=linewidth,
      )
    end
  end

  return nothing
end

function _plot_surface!(
  scene,
  mesh::CurvilinearGrids.MultiBlockMesh;
  color=:lightskyblue3,
  shading::Bool=true,
  transparency::Bool=false,
  colors=nothing,
)
  block_colors = _resolve_block_colors(length(mesh.blocks), colors)
  for i in eachindex(mesh.blocks)
    _plot_surface!(
      scene,
      mesh.blocks[i];
      color=block_colors[i],
      shading=shading,
      transparency=transparency,
      colors=nothing,
    )
  end
  return nothing
end

function _plot_surface!(
  scene,
  mesh::_GridLike;
  color=:lightskyblue3,
  shading::Bool=true,
  transparency::Bool=false,
  colors=nothing,
)
  cart = _cartesian_node_coordinates(mesh)
  N = length(cart)

  if N == 1
    throw(ArgumentError("`mesh` surface plotting is not defined for 1D grids."))
  elseif N == 2
    x, y = cart
    smesh = _structured_surface_mesh_2d(x, y)
    Makie.mesh!(scene, smesh; color=color, shading=shading, transparency=transparency)
  elseif N == 3
    x, y, z = cart
    _plot_surface_3d_boundaries!(
      scene, x, y, z; color=color, shading=shading, transparency=transparency
    )
  else
    throw(ArgumentError("Unsupported mesh dimension $N for surface plotting."))
  end

  return nothing
end

@inline _structured_index(i::Int, j::Int, ni::Int) = i + (j - 1) * ni

function _structured_surface_mesh_2d(x::AbstractMatrix, y::AbstractMatrix)
  ni, nj = size(x)
  size(y) == (ni, nj) ||
    throw(ArgumentError("2D coordinate arrays must have matching sizes."))

  points = Vector{GB.Point3f}(undef, ni * nj)
  for j in 1:nj, i in 1:ni
    k = _structured_index(i, j, ni)
    points[k] = GB.Point3f(Float32(x[i, j]), Float32(y[i, j]), 0.0f0)
  end

  faces = Vector{GB.GLTriangleFace}(undef, 2 * (ni - 1) * (nj - 1))
  f = 1
  for j in 1:(nj - 1), i in 1:(ni - 1)
    v11 = _structured_index(i, j, ni)
    v21 = _structured_index(i + 1, j, ni)
    v12 = _structured_index(i, j + 1, ni)
    v22 = _structured_index(i + 1, j + 1, ni)
    faces[f] = GB.GLTriangleFace(v11, v21, v12)
    faces[f + 1] = GB.GLTriangleFace(v21, v22, v12)
    f += 2
  end
  return GB.Mesh(points, faces)
end

function _structured_surface_mesh_3d(
  x::AbstractMatrix, y::AbstractMatrix, z::AbstractMatrix
)
  ni, nj = size(x)
  size(y) == (ni, nj) == size(z) ||
    throw(ArgumentError("3D face coordinate arrays must have matching sizes."))

  points = Vector{GB.Point3f}(undef, ni * nj)
  for j in 1:nj, i in 1:ni
    k = _structured_index(i, j, ni)
    points[k] = GB.Point3f(Float32(x[i, j]), Float32(y[i, j]), Float32(z[i, j]))
  end

  faces = Vector{GB.GLTriangleFace}(undef, 2 * (ni - 1) * (nj - 1))
  f = 1
  for j in 1:(nj - 1), i in 1:(ni - 1)
    v11 = _structured_index(i, j, ni)
    v21 = _structured_index(i + 1, j, ni)
    v12 = _structured_index(i, j + 1, ni)
    v22 = _structured_index(i + 1, j + 1, ni)
    faces[f] = GB.GLTriangleFace(v11, v21, v12)
    faces[f + 1] = GB.GLTriangleFace(v21, v22, v12)
    f += 2
  end
  return GB.Mesh(points, faces)
end

function _plot_surface_3d_boundaries!(
  scene,
  x::AbstractArray{<:Real,3},
  y::AbstractArray{<:Real,3},
  z::AbstractArray{<:Real,3};
  color,
  shading,
  transparency,
)
  i0, i1 = firstindex(x, 1), lastindex(x, 1)
  j0, j1 = firstindex(x, 2), lastindex(x, 2)
  k0, k1 = firstindex(x, 3), lastindex(x, 3)

  for i in (i0, i1)
    smesh = _structured_surface_mesh_3d(
      @view(x[i, :, :]), @view(y[i, :, :]), @view(z[i, :, :])
    )
    Makie.mesh!(scene, smesh; color=color, shading=shading, transparency=transparency)
  end
  for j in (j0, j1)
    smesh = _structured_surface_mesh_3d(
      @view(x[:, j, :]), @view(y[:, j, :]), @view(z[:, j, :])
    )
    Makie.mesh!(scene, smesh; color=color, shading=shading, transparency=transparency)
  end
  for k in (k0, k1)
    smesh = _structured_surface_mesh_3d(
      @view(x[:, :, k]), @view(y[:, :, k]), @view(z[:, :, k])
    )
    Makie.mesh!(scene, smesh; color=color, shading=shading, transparency=transparency)
  end
  return nothing
end

function _boundary_point_cloud(
  x::AbstractArray{<:Real,3}, y::AbstractArray{<:Real,3}, z::AbstractArray{<:Real,3}
)
  i0, i1 = firstindex(x, 1), lastindex(x, 1)
  j0, j1 = firstindex(x, 2), lastindex(x, 2)
  k0, k1 = firstindex(x, 3), lastindex(x, 3)

  xb = vcat(
    vec(@view(x[i0, :, :])),
    vec(@view(x[i1, :, :])),
    vec(@view(x[:, j0, :])),
    vec(@view(x[:, j1, :])),
    vec(@view(x[:, :, k0])),
    vec(@view(x[:, :, k1])),
  )
  yb = vcat(
    vec(@view(y[i0, :, :])),
    vec(@view(y[i1, :, :])),
    vec(@view(y[:, j0, :])),
    vec(@view(y[:, j1, :])),
    vec(@view(y[:, :, k0])),
    vec(@view(y[:, :, k1])),
  )
  zb = vcat(
    vec(@view(z[i0, :, :])),
    vec(@view(z[i1, :, :])),
    vec(@view(z[:, j0, :])),
    vec(@view(z[:, j1, :])),
    vec(@view(z[:, :, k0])),
    vec(@view(z[:, :, k1])),
  )

  return xb, yb, zb
end

function _wireframe_dimension(mesh)
  if mesh isa CurvilinearGrids.MultiBlockMesh
    return _wireframe_dimension(first(mesh.blocks))
  end
  cart = _cartesian_node_coordinates(mesh)
  return length(cart) == 3 ? 3 : 2
end

function _cartesian_node_coordinates(mesh::_GridLike)
  raw = CurvilinearGrids.coords(mesh)
  q = _normalize_coordinate_arrays(raw)
  return _to_cartesian(q, _coordinate_system_trait(mesh))
end

@inline _normalize_coordinate_arrays(x::AbstractVector) = (x,)

function _normalize_coordinate_arrays(raw::Tuple)
  N = length(raw)
  if N == 2
    return _normalize_coordinate_arrays_2(raw[1], raw[2])
  elseif N == 3
    return _normalize_coordinate_arrays_3(raw[1], raw[2], raw[3])
  end
  throw(ArgumentError("Unsupported coordinate tuple length $N."))
end

function _normalize_coordinate_arrays_2(x::AbstractVector, y::AbstractVector)
  xx = [x[i] for i in eachindex(x), _ in eachindex(y)]
  yy = [y[j] for _ in eachindex(x), j in eachindex(y)]
  return (xx, yy)
end

function _normalize_coordinate_arrays_2(x::AbstractArray, y::AbstractArray)
  size(x) == size(y) ||
    throw(ArgumentError("2D coordinate arrays must have matching sizes."))
  return (x, y)
end

function _normalize_coordinate_arrays_3(
  x::AbstractVector, y::AbstractVector, z::AbstractVector
)
  xx = [x[i] for i in eachindex(x), _ in eachindex(y), _ in eachindex(z)]
  yy = [y[j] for _ in eachindex(x), j in eachindex(y), _ in eachindex(z)]
  zz = [z[k] for _ in eachindex(x), _ in eachindex(y), k in eachindex(z)]
  return (xx, yy, zz)
end

function _normalize_coordinate_arrays_3(
  x::AbstractArray, y::AbstractArray, z::AbstractArray
)
  size(x) == size(y) == size(z) ||
    throw(ArgumentError("3D coordinate arrays must have matching sizes."))
  return (x, y, z)
end

@inline _to_cartesian(q::NTuple{1,Any}, ::CurvilinearGrids.CoordinateSystemTrait) = q

@inline _to_cartesian(q::NTuple{2,Any}, ::CurvilinearGrids.CartesianCS) = q
@inline _to_cartesian(q::NTuple{2,Any}, ::CurvilinearGrids.CurvilinearCS) = q
@inline _to_cartesian(q::NTuple{2,Any}, ::CurvilinearGrids.AxisymmetricCS{:x}) = q
@inline _to_cartesian(q::NTuple{2,Any}, ::CurvilinearGrids.AxisymmetricCS{:y}) = q
@inline _to_cartesian(q::NTuple{2,Any}, ::CurvilinearGrids.CylindricalCS) = q

function _to_cartesian((r, theta)::NTuple{2,Any}, ::CurvilinearGrids.SphericalCS)
  x = @. r * sin(theta)
  y = @. r * cos(theta)
  return (x, y)
end

@inline _to_cartesian(q::NTuple{3,Any}, ::CurvilinearGrids.CartesianCS) = q
@inline _to_cartesian(q::NTuple{3,Any}, ::CurvilinearGrids.CurvilinearCS) = q

function _to_cartesian((r, theta, z)::NTuple{3,Any}, ::CurvilinearGrids.CylindricalCS)
  x = @. r * cos(theta)
  y = @. r * sin(theta)
  return (x, y, z)
end

function _to_cartesian((r, theta, phi)::NTuple{3,Any}, ::CurvilinearGrids.SphericalCS)
  st = @. sin(theta)
  x = @. r * st * cos(phi)
  y = @. r * st * sin(phi)
  z = @. r * cos(theta)
  return (x, y, z)
end

@inline _to_cartesian(q::Tuple, ::CurvilinearGrids.CoordinateSystemTrait) = q

@inline _coordinate_system_trait(grid::CurvilinearGrids.AbstractUnifiedGrid) = CurvilinearGrids.coordinate_system(
  grid
)

@inline _coordinate_system_trait(::CurvilinearGrids.CylindricalGrid1D) = CurvilinearGrids.CylindricalCS()
@inline _coordinate_system_trait(::CurvilinearGrids.CylindricalOrthogonalGrid1D) = CurvilinearGrids.CylindricalCS()
@inline _coordinate_system_trait(::CurvilinearGrids.SphericalGrid1D) = CurvilinearGrids.SphericalCS()
@inline _coordinate_system_trait(::CurvilinearGrids.SphericalGrid3D) = CurvilinearGrids.SphericalCS()
@inline _coordinate_system_trait(::CurvilinearGrids.SphericalOrthogonalGrid1D) = CurvilinearGrids.SphericalCS()
@inline _coordinate_system_trait(::CurvilinearGrids.SphericalBasisCurvilinearGrid3D) = CurvilinearGrids.CartesianCS()
@inline _coordinate_system_trait(::CurvilinearGrids.AxisymmetricOrthogonalGrid2D) = CurvilinearGrids.AxisymmetricCS{
  :y
}()

@inline function _coordinate_system_trait(grid::CurvilinearGrids.AxisymmetricGrid2D)
  return if grid.rotational_axis === :x
    CurvilinearGrids.AxisymmetricCS{:x}()
  else
    CurvilinearGrids.AxisymmetricCS{:y}()
  end
end

@inline _coordinate_system_trait(::CurvilinearGrids.AbstractCurvilinearGrid) = CurvilinearGrids.CartesianCS()

function _resolve_block_colors(n::Int, colors)
  if colors === nothing
    n <= length(_DEFAULT_MULTIBLOCK_COLORS) || throw(
      ArgumentError(
        "Default multiblock color palette supports at most $(length(_DEFAULT_MULTIBLOCK_COLORS)) blocks. Pass `colors` with one color per block.",
      ),
    )
    return collect(_DEFAULT_MULTIBLOCK_COLORS[1:n])
  end

  if colors isa Tuple || colors isa AbstractVector
    length(colors) >= n ||
      throw(ArgumentError("`colors` must provide at least one color per block (need $n)."))
    return collect(colors[1:n])
  end

  throw(
    ArgumentError("`colors` must be `nothing` or a tuple/vector with one color per block.")
  )
end

@inline function _plot_attr(plot, key::Symbol, default)
  haskey(plot.attributes, key) || return default
  return Makie.to_value(plot.attributes[key])
end

end
