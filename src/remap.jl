using Interpolations

"""Get interpolation functions of the mesh coordinates"""
function mesh_coord_interpolators(mesh::CurvilinearGrid1D)
  x = CurvilinearGrids.coords(mesh)
  x_itp = interpolate(x, BSpline(Linear()))

  return (; x=x_itp)
end

"""Get interpolation functions of the mesh coordinates"""
function mesh_coord_interpolators(mesh::CurvilinearGrid2D)
  x, y = CurvilinearGrids.coords(mesh)
  x_itp = interpolate(x, BSpline(Linear()))
  y_itp = interpolate(y, BSpline(Linear()))

  return (; x=x_itp, y=y_itp)
end

"""Get interpolation functions of the mesh coordinates"""
function mesh_coord_interpolators(mesh::CurvilinearGrid3D)
  x, y, z = CurvilinearGrids.coords(mesh)
  x_itp = interpolate(x, BSpline(Linear()))
  y_itp = interpolate(y, BSpline(Linear()))
  z_itp = interpolate(z, BSpline(Linear()))

  return (; x=x_itp, y=y_itp, z=z_itp)
end

"""
    change_resolution(mesh::CurvilinearGrid1D, ni::Int)

Re-interpolate the resolution of the given `mesh` to have the new dimensions `ni`.
"""
function change_resolution(mesh::CurvilinearGrid1D, ni::Int)
  itps = mesh_coord_interpolators(mesh)

  nx, ny = size(itps.x.coefs)
  new_x = itps.x(range(1; stop=nx, length=ni), range(1; stop=ny, length=nj))

  new_mesh = CurvilinearGrid1D(new_x, mesh.discretization_scheme_name)

  return new_mesh
end

"""
    change_resolution(mesh::CurvilinearGrid2D, (ni, nj)::NTuple{2,Int})

Re-interpolate the resolution of the given `mesh` to have the new dimensions `(ni, nj)`.
"""
function change_resolution(mesh::CurvilinearGrid2D, (ni, nj)::NTuple{2,Int})
  itps = mesh_coord_interpolators(mesh)

  nx, ny = size(itps.x.coefs)
  new_x = itps.x(range(1; stop=nx, length=ni), range(1; stop=ny, length=nj))
  new_y = itps.y(range(1; stop=nx, length=ni), range(1; stop=ny, length=nj))

  new_mesh = CurvilinearGrid2D(new_x, new_y, mesh.discretization_scheme_name)

  return new_mesh
end

"""
    change_resolution(mesh::CurvilinearGrid3D, (ni, nj, nk)::NTuple{3,Int})

Re-interpolate the resolution of the given `mesh` to have the new dimensions `(ni, nj)`.
"""
function change_resolution(mesh::CurvilinearGrid3D, (ni, nj, nk)::NTuple{3,Int})
  itps = mesh_coord_interpolators(mesh)

  nx, ny, nz = size(itps.x.coefs)
  new_x = itps.x(
    range(1; stop=nx, length=ni), #
    range(1; stop=ny, length=nj), #
    range(1; stop=nz, length=nk), #
  )

  new_y = itps.y(
    range(1; stop=nx, length=ni), #
    range(1; stop=ny, length=nj), #
    range(1; stop=nz, length=nk), #
  )

  new_z = itps.z(
    range(1; stop=nx, length=ni), #
    range(1; stop=ny, length=nj), #
    range(1; stop=nz, length=nk), #
  )

  new_mesh = CurvilinearGrid3D(new_x, new_y, new_z, mesh.discretization_scheme_name)

  return new_mesh
end

"""
    scale_resolution(mesh::CurvilinearGrid2D, (α, β)::NTuple{2,Real})

Scale / re-interpolate the resolution of the given `mesh` by a factor of (α,β), where 
α is along the first dimension, β along the 2nd dimension
"""
function scale_resolution(mesh::CurvilinearGrid1D, α::Real)
  itps = mesh_coord_interpolators(mesh)

  ni, nj = size(itps.x.coefs)
  new_x = itps.x(range(1; stop=ni, step=inv(α)), range(1; stop=nj, step=inv(β)))
  new_y = itps.y(range(1; stop=ni, step=inv(α)), range(1; stop=nj, step=inv(β)))

  new_mesh = CurvilinearGrid2D(new_x, new_y, mesh.discretization_scheme_name)

  return new_mesh
end

"""
    scale_resolution(mesh::CurvilinearGrid2D, (α, β)::NTuple{2,Real})

Scale / re-interpolate the resolution of the given `mesh` by a factor of (α,β), where 
α is along the first dimension, β along the 2nd dimension
"""
function scale_resolution(mesh::CurvilinearGrid2D, (α, β)::NTuple{2,Real})
  itps = mesh_coord_interpolators(mesh)

  ni, nj = size(itps.x.coefs)
  new_x = itps.x(range(1; stop=ni, step=inv(α)), range(1; stop=nj, step=inv(β)))
  new_y = itps.y(range(1; stop=ni, step=inv(α)), range(1; stop=nj, step=inv(β)))

  new_mesh = CurvilinearGrid2D(new_x, new_y, mesh.discretization_scheme_name)

  return new_mesh
end

"""
    scale_resolution(mesh::CurvilinearGrid3D, (α, β, γ)::NTuple{3,Real})

Scale / re-interpolate the resolution of the given `mesh` by a factor of (α,β,γ), where 
α is along the first dimension, β along the 2nd dimension
"""
function scale_resolution(mesh::CurvilinearGrid3D, (α, β, γ)::NTuple{3,Real})
  itps = mesh_coord_interpolators(mesh)

  ni, nj, nk = size(itps.x.coefs)
  new_x = itps.x(
    range(1; stop=ni, step=inv(α)), # 
    range(1; stop=nj, step=inv(β)), # 
    range(1; stop=nk, step=inv(γ)), #
  )

  new_y = itps.y(
    range(1; stop=ni, step=inv(α)), # 
    range(1; stop=nj, step=inv(β)), # 
    range(1; stop=nk, step=inv(γ)), #
  )

  new_z = itps.z(
    range(1; stop=ni, step=inv(α)), # 
    range(1; stop=nj, step=inv(β)), # 
    range(1; stop=nk, step=inv(γ)), #
  )

  new_mesh = CurvilinearGrid3D(new_x, new_y, new_z, mesh.discretization_scheme_name)

  return new_mesh
end

"""
    remap_cell_data(mesh_1::CurvilinearGrid1D, mesh_2::CurvilinearGrid1D, data_on_1::AbstractArray{T,1}) where {T}

Remap (e.g. coarsen or refine) cell-based data sized for `mesh_1` onto `mesh_2`. This currently uses
BSpline(Linear()) interpolation from Interpolations.jl
"""
function remap_cell_data(
  mesh_1::CurvilinearGrid1D, mesh_2::CurvilinearGrid1D, data_on_1::AbstractArray{T,1}
) where {T}
  dims_1 = size(mesh_1.iterators.cell.domain)
  dims_2 = size(mesh_2.iterators.cell.domain)

  @assert size(data_on_1) == size(mesh_1.iterators.cell.domain) "The input data size ($(size(data_on_1))) is incorrect. The data extents must not include halo regions, only the inner region, which has size ($(size(mesh_1.iterators.cell.domain)))."

  data_itp = interpolate(data_on_1, BSpline(Linear()))

  data_on_2 = data_itp(range(1; stop=dims_1[1], length=dims_2[1]))

  return data_on_2
end
"""
    remap_cell_data(mesh_1::CurvilinearGrid2D, mesh_2::CurvilinearGrid2D, data_on_1::AbstractArray{T,2}) where {T}

Remap (e.g. coarsen or refine) cell-based data sized for `mesh_1` onto `mesh_2`. This currently uses
BSpline(Linear()) interpolation from Interpolations.jl
"""
function remap_cell_data(
  mesh_1::CurvilinearGrid2D, mesh_2::CurvilinearGrid2D, data_on_1::AbstractArray{T,2}
) where {T}
  dims_1 = size(mesh_1.iterators.cell.domain)
  dims_2 = size(mesh_2.iterators.cell.domain)

  @assert size(data_on_1) == size(mesh_1.iterators.cell.domain) "The input data size ($(size(data_on_1))) is incorrect. The data extents must not include halo regions, only the inner region, which has size ($(size(mesh_1.iterators.cell.domain)))."

  data_itp = interpolate(data_on_1, BSpline(Linear()))

  data_on_2 = data_itp(
    range(1; stop=dims_1[1], length=dims_2[1]), # 
    range(1; stop=dims_1[2], length=dims_2[2]), #
  )

  return data_on_2
end
"""
    remap_cell_data(mesh_1::CurvilinearGrid3D, mesh_2::CurvilinearGrid3D, data_on_1::AbstractArray{T,3}) where {T}

Remap (e.g. coarsen or refine) cell-based data sized for `mesh_1` onto `mesh_2`. This currently uses
BSpline(Linear()) interpolation from Interpolations.jl
"""
function remap_cell_data(
  mesh_1::CurvilinearGrid3D, mesh_2::CurvilinearGrid3D, data_on_1::AbstractArray{T,3}
) where {T}
  dims_1 = size(mesh_1.iterators.cell.domain)
  dims_2 = size(mesh_2.iterators.cell.domain)

  @assert size(data_on_1) == size(mesh_1.iterators.cell.domain) "The input data size ($(size(data_on_1))) is incorrect. The data extents must not include halo regions, only the inner region, which has size ($(size(mesh_1.iterators.cell.domain)))."

  data_itp = interpolate(data_on_1, BSpline(Linear()))

  data_on_2 = data_itp(
    range(1; stop=dims_1[1], length=dims_2[1]), # 
    range(1; stop=dims_1[2], length=dims_2[2]), #
    range(1; stop=dims_1[3], length=dims_2[3]), #
  )

  return data_on_2
end
