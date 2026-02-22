module VTKOutput

using WriteVTK

using ..GridTypes

export save_vtk

@inline _matrix_component_names(::Val{1}) = ["x11"]
@inline _matrix_component_names(::Val{2}) = ["x11", "x12", "x21", "x22"]
@inline _matrix_component_names(::Val{3}) = [
  "x11", "x12", "x13", "x21", "x22", "x23", "x31", "x32", "x33"
]
@inline _gcl_component_names(::Val{1}) = ["I1"]
@inline _gcl_component_names(::Val{2}) = ["I1", "I2"]
@inline _gcl_component_names(::Val{3}) = ["I1", "I2", "I3"]
@inline _index_component_names(::Val{1}) = ["i"]
@inline _index_component_names(::Val{2}) = ["i", "j"]
@inline _index_component_names(::Val{3}) = ["i", "j", "k"]
@inline _axis_name(::Val{1}, axis::Int) =
  axis == 1 ? "i" : throw(ArgumentError("Invalid axis $axis for 1D."))
@inline function _axis_name(::Val{2}, axis::Int)
  return if axis == 1
    "i"
  elseif axis == 2
    "j"
  else
    throw(ArgumentError("Invalid axis $axis for 2D."))
  end
end
@inline function _axis_name(::Val{3}, axis::Int)
  return if axis == 1
    "i"
  elseif axis == 2
    "j"
  elseif axis == 3
    "k"
  else
    throw(ArgumentError("Invalid axis $axis for 3D."))
  end
end

function _matrix_components(data, domain::CartesianIndices{N}, ::Val{N}) where {N}
  components = ntuple(Val(N * N)) do c
    i = ((c - 1) ÷ N) + 1
    j = ((c - 1) % N) + 1
    getindex.(data[domain], i, j)
  end
  return components
end

function _write_metric_matrix_field!(
  vtk, field_name::AbstractString, data, domain::CartesianIndices{N}, ::Val{N}
) where {N}
  vtk[field_name, VTKCellData(), component_names = _matrix_component_names(Val(N))] = _matrix_components(
    data, domain, Val(N)
  )
  return nothing
end

function _write_unified_metrics!(
  vtk,
  mesh::Union{MappedGrid{N},DiscreteGrid{N}},
  domain::CartesianIndices{N};
  refresh_metrics::Bool=false,
) where {N}
  cm = cell_metrics(mesh; refresh=refresh_metrics)
  fm = face_metrics(mesh; refresh=refresh_metrics)

  vtk["J", VTKCellData()] = getproperty.(cm.forward[domain], :J)
  _write_metric_matrix_field!(vtk, "cell_forward", cm.forward, domain, Val(N))
  _write_metric_matrix_field!(vtk, "cell_inverse", cm.inverse, domain, Val(N))

  for axis in 1:N
    axis_name = _axis_name(Val(N), axis)
    vtk["face_$(axis_name)_forward_J", VTKCellData()] = getproperty.(
      fm[axis].forward[domain], :J
    )
    vtk["face_$(axis_name)_inverse_J", VTKCellData()] = getproperty.(
      fm[axis].inverse[domain], :J
    )
    _write_metric_matrix_field!(
      vtk, "face_$(axis_name)_forward", fm[axis].forward, domain, Val(N)
    )
    _write_metric_matrix_field!(
      vtk, "face_$(axis_name)_inverse", fm[axis].inverse, domain, Val(N)
    )
    _write_metric_matrix_field!(
      vtk, "face_$(axis_name)_conserved", fm[axis].conserved, domain, Val(N)
    )
  end

  gcl_residual = GridTypes.gcl(fm, domain)
  if N == 1
    vtk["GCL", VTKCellData()] = gcl_residual[domain]
  else
    vtk["GCL", VTKCellData(), component_names = _gcl_component_names(Val(N))] = ntuple(
      d -> gcl_residual[d][domain], N
    )
  end

  return nothing
end

function _write_cell_indices!(vtk, domain::CartesianIndices{N}) where {N}
  indices = collect(domain)

  if N == 1
    vtk["index", VTKCellData()] = [idx.I[1] for idx in indices]
  else
    vtk["index", VTKCellData(), component_names = _index_component_names(Val(N))] = ntuple(
      d -> [idx.I[d] for idx in indices], N
    )
  end

  return nothing
end

function _write_extra_cell_data!(vtk, domain::CartesianIndices, extra_cell_data)
  if isnothing(extra_cell_data)
    return nothing
  end

  for (key, value) in pairs(extra_cell_data)
    vtk[String(key), VTKCellData()] = value[domain]
  end

  return nothing
end

function _save_vtk_unified(
  mesh::Union{MappedGrid{N},DiscreteGrid{N}},
  fn::AbstractString;
  include_metrics::Bool=true,
  refresh_metrics::Bool=false,
  extra_cell_data=nothing,
) where {N}
  @info "Writing to $fn.vti"

  if include_metrics && mesh.metric_caches === nothing
    throw(
      ArgumentError(
        "`save_vtk(...; include_metrics=true)` requires metric storage. Recreate the grid without `compute_metrics=false, cache_mode=:off` or call with `include_metrics=false`.",
      ),
    )
  end

  xyz_n = coords(mesh)
  domain = mesh.iterators.cell.domain

  @views vtk_grid(fn, xyz_n) do vtk
    if include_metrics
      _write_unified_metrics!(vtk, mesh, domain; refresh_metrics=refresh_metrics)
    end

    _write_cell_indices!(vtk, domain)
    _write_extra_cell_data!(vtk, domain, extra_cell_data)
  end

  return nothing
end

"""
    save_vtk(mesh::AbstractCurvilinearGrid3D, fn=\"mesh\")

Write a 3D curvilinear grid and its metrics to a `.vti` VTK file named `fn`.
"""
function save_vtk(mesh::AbstractCurvilinearGrid3D, fn="mesh")
  @info "Writing to $fn.vti"

  xyz_n = coords(mesh)
  domain = mesh.iterators.cell.domain

  @views vtk_grid(fn, xyz_n) do vtk
    vtk["J", VTKCellData()] = mesh.cell_center_metrics.J[domain]

    # vtk["volume", VTKCellData()] = cellvolume.(Ref(mesh), domain)

    I1, I2, I3 = GridTypes.gcl(mesh.edge_metrics, domain)

    vtk["GCL", VTKCellData(), component_names = ["I1", "I2", "I3"]] = (
      I1[domain], I2[domain], I3[domain]
    )

    indices = collect(domain)

    vtk["index", VTKCellData(), component_names = ["i", "j", "k"]] = (
      [idx.I[1] for idx in indices],
      [idx.I[2] for idx in indices],
      [idx.I[3] for idx in indices],
    )

    for (name, edge) in pairs(mesh.edge_metrics)
      for dim in (:ξ̂, :η̂, :ζ̂)
        var = edge[dim]
        vtk["$(name)_$(dim)", VTKCellData(), component_names = ["x1", "x2", "x3"]] = (
          var.x₁[domain], var.x₂[domain], var.x₃[domain]
        )
      end
    end

    vtk["xi", VTKCellData(), component_names = ["x1", "x2", "x3", "t"]] = (
      mesh.cell_center_metrics.ξ.x₁[domain],
      mesh.cell_center_metrics.ξ.x₂[domain],
      mesh.cell_center_metrics.ξ.x₃[domain],
      mesh.cell_center_metrics.ξ.t[domain],
    )

    vtk["eta", VTKCellData(), component_names = ["x1", "x2", "x3", "t"]] = (
      mesh.cell_center_metrics.η.x₁[domain],
      mesh.cell_center_metrics.η.x₂[domain],
      mesh.cell_center_metrics.η.x₃[domain],
      mesh.cell_center_metrics.η.t[domain],
    )

    vtk["zeta", VTKCellData(), component_names = ["x1", "x2", "x3", "t"]] = (
      mesh.cell_center_metrics.ζ.x₁[domain],
      mesh.cell_center_metrics.ζ.x₂[domain],
      mesh.cell_center_metrics.ζ.x₃[domain],
      mesh.cell_center_metrics.ζ.t[domain],
    )

    vtk["dx_di", VTKCellData(), component_names = ["xi", "eta", "zeta"]] = (
      mesh.cell_center_metrics.x₁.ξ[domain],
      mesh.cell_center_metrics.x₁.η[domain],
      mesh.cell_center_metrics.x₁.ζ[domain],
    )

    vtk["dy_di", VTKCellData(), component_names = ["xi", "eta", "zeta"]] = (
      mesh.cell_center_metrics.x₂.ξ[domain],
      mesh.cell_center_metrics.x₂.η[domain],
      mesh.cell_center_metrics.x₂.ζ[domain],
    )

    vtk["dz_di", VTKCellData(), component_names = ["xi", "eta", "zeta"]] = (
      mesh.cell_center_metrics.x₃.ξ[domain],
      mesh.cell_center_metrics.x₃.η[domain],
      mesh.cell_center_metrics.x₃.ζ[domain],
    )
  end
end

"""
    save_vtk(mesh::Union{MappedGrid{2},DiscreteGrid{2}}, fn="mesh"; kwargs...)
    save_vtk(mesh::Union{MappedGrid{3},DiscreteGrid{3}}, fn="mesh"; kwargs...)

Write unified mapped/discrete grids to VTK. Metric fields are included by
default when metric storage is available.

# Keywords
  - `include_metrics`: Write cell/face metric and GCL fields. Default: `true`.
  - `refresh_metrics`: Force metric-cache refresh before writing. Default: `false`.
  - `extra_cell_data`: Optional named tuple or dictionary of additional cell-data
    arrays indexed by `mesh.iterators.cell.domain`.
"""
function save_vtk(
  mesh::Union{MappedGrid{2},DiscreteGrid{2}},
  fn="mesh";
  include_metrics::Bool=true,
  refresh_metrics::Bool=false,
  extra_cell_data=nothing,
)
  _save_vtk_unified(
    mesh,
    fn;
    include_metrics=include_metrics,
    refresh_metrics=refresh_metrics,
    extra_cell_data=extra_cell_data,
  )
end

function save_vtk(
  mesh::Union{MappedGrid{3},DiscreteGrid{3}},
  fn="mesh";
  include_metrics::Bool=true,
  refresh_metrics::Bool=false,
  extra_cell_data=nothing,
)
  _save_vtk_unified(
    mesh,
    fn;
    include_metrics=include_metrics,
    refresh_metrics=refresh_metrics,
    extra_cell_data=extra_cell_data,
  )
end

"""
    save_vtk(mesh::SphericalGrid3D, fn=\"mesh\"; extra_cell_data=nothing)

Write a spherical grid to VTK, optionally including additional `extra_cell_data` arrays.
"""
function save_vtk(mesh::SphericalGrid3D, fn="mesh"; extra_cell_data=nothing)
  @info "Writing to $fn.vti"

  x = @view mesh.cartesian_node_coordinates.x[mesh.iterators.node.domain]
  y = @view mesh.cartesian_node_coordinates.y[mesh.iterators.node.domain]
  z = @view mesh.cartesian_node_coordinates.z[mesh.iterators.node.domain]

  domain = mesh.iterators.cell.domain

  indices = collect(domain)

  @views vtk_grid(fn, (x, y, z)) do vtk
    vtk["volume", VTKCellData()] = mesh.cell_volumes[domain]

    vtk["index", VTKCellData(), component_names = ["i", "j", "k"]] = (
      [idx.I[1] for idx in indices],
      [idx.I[2] for idx in indices],
      [idx.I[3] for idx in indices],
    )

    vtk["face_area_i₊½", VTKCellData()] = mesh.face_areas.i₊½[domain]
    vtk["face_area_j₊½", VTKCellData()] = mesh.face_areas.j₊½[domain]
    vtk["face_area_k₊½", VTKCellData()] = mesh.face_areas.k₊½[domain]
    if !isnothing(extra_cell_data)
      for (key, value) in pairs(extra_cell_data)
        vtk[String(key), VTKCellData()] = value[domain]
      end
    end
  end
end

"""
    save_vtk((x, y, z)::NTuple{3,AbstractArray}, fn=\"mesh\")

Write raw `(x, y, z)` coordinates to a `.vti` VTK file named `fn`.
"""
function save_vtk((x, y, z)::NTuple{3,AbstractArray{T,3}}, fn="mesh") where {T}
  @info "Writing to $fn.vti"

  vtk_grid(fn, x, y, z) do vtk
  end
end

"""
    save_vtk((x, y)::NTuple{2,AbstractArray}, fn=\"mesh\")

Write raw `(x, y)` coordinates to a `.vti` VTK file named `fn`.
"""
function save_vtk((x, y)::NTuple{2,AbstractMatrix{N}}, fn="mesh") where {N}
  @info "Writing to $fn.vti"

  vtk_grid(fn, x, y) do vtk
  end
end

"""
    save_vtk(mesh::AbstractCurvilinearGrid2D, fn=\"mesh\")

Write a 2D curvilinear grid and its metrics to a `.vti` VTK file named `fn`.
"""
function save_vtk(mesh::AbstractCurvilinearGrid2D, fn="mesh")
  @info "Writing to $fn.vti"

  xyz_n = coords(mesh)
  domain = mesh.iterators.cell.domain

  @views vtk_grid(fn, xyz_n) do vtk
    vtk["J", VTKCellData()] = mesh.cell_center_metrics.J[domain]

    # vtk["volume", VTKCellData()] = cellvolume.(Ref(mesh), domain)

    I1, I2 = GridTypes.gcl(mesh.edge_metrics, domain)

    vtk["GCL", VTKCellData(), component_names = ["I1", "I2"]] = (I1[domain], I2[domain])

    indices = collect(domain)

    vtk["index", VTKCellData(), component_names = ["i", "j"]] = (
      [idx.I[1] for idx in indices], [idx.I[2] for idx in indices]
    )

    for (name, edge) in pairs(mesh.edge_metrics)
      for dim in (:ξ̂, :η̂)
        var = edge[dim]
        vtk["$(name)_$(dim)", VTKCellData(), component_names = ["x1", "x2"]] = (
          var.x₁[domain], var.x₂[domain]
        )
      end
    end

    vtk["xi", VTKCellData(), component_names = ["x1", "x2", "t"]] = (
      mesh.cell_center_metrics.ξ.x₁[domain],
      mesh.cell_center_metrics.ξ.x₂[domain],
      mesh.cell_center_metrics.ξ.t[domain],
    )

    vtk["eta", VTKCellData(), component_names = ["x1", "x2", "t"]] = (
      mesh.cell_center_metrics.η.x₁[domain],
      mesh.cell_center_metrics.η.x₂[domain],
      mesh.cell_center_metrics.η.t[domain],
    )

    vtk["dx_di", VTKCellData(), component_names = ["xi", "eta"]] = (
      mesh.cell_center_metrics.x₁.ξ[domain], mesh.cell_center_metrics.x₁.η[domain]
    )

    vtk["dy_di", VTKCellData(), component_names = ["xi", "eta"]] = (
      mesh.cell_center_metrics.x₂.ξ[domain], mesh.cell_center_metrics.x₂.η[domain]
    )
  end
end

end
