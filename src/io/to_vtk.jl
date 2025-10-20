module VTKOutput

using WriteVTK

using ..GridTypes

export save_vtk

"""Write the mesh to .VTK format"""
function save_vtk(mesh::AbstractCurvilinearGrid3D, fn="mesh")
  @info "Writing to $fn.vti"

  xyz_n = coords(mesh)
  domain = mesh.iterators.cell.domain

  @views vtk_grid(fn, xyz_n) do vtk
    vtk["J", VTKCellData()] = mesh.cell_center_metrics.J[domain]

    # vtk["volume", VTKCellData()] = cellvolume.(Ref(mesh), domain)

    for (name, edge) in pairs(mesh.edge_metrics)
      for dim in (:ξ̂, :η̂, :ζ̂)
        var = edge[dim]
        vtk["$(name)_$(dim)", VTKCellData(), component_names=["x1", "x2", "x3"]] = (
          var.x₁[domain], var.x₂[domain], var.x₃[domain]
        )
      end
    end

    vtk["xi", VTKCellData(), component_names=["x1", "x2", "x3", "t"]] = (
      mesh.cell_center_metrics.ξ.x₁[domain],
      mesh.cell_center_metrics.ξ.x₂[domain],
      mesh.cell_center_metrics.ξ.x₃[domain],
      mesh.cell_center_metrics.ξ.t[domain],
    )

    vtk["eta", VTKCellData(), component_names=["x1", "x2", "x3", "t"]] = (
      mesh.cell_center_metrics.η.x₁[domain],
      mesh.cell_center_metrics.η.x₂[domain],
      mesh.cell_center_metrics.η.x₃[domain],
      mesh.cell_center_metrics.η.t[domain],
    )

    vtk["zeta", VTKCellData(), component_names=["x1", "x2", "x3", "t"]] = (
      mesh.cell_center_metrics.ζ.x₁[domain],
      mesh.cell_center_metrics.ζ.x₂[domain],
      mesh.cell_center_metrics.ζ.x₃[domain],
      mesh.cell_center_metrics.ζ.t[domain],
    )

    vtk["dx_di", VTKCellData(), component_names=["xi", "eta", "zeta"]] = (
      mesh.cell_center_metrics.x₁.ξ[domain],
      mesh.cell_center_metrics.x₁.η[domain],
      mesh.cell_center_metrics.x₁.ζ[domain],
    )

    vtk["dy_di", VTKCellData(), component_names=["xi", "eta", "zeta"]] = (
      mesh.cell_center_metrics.x₂.ξ[domain],
      mesh.cell_center_metrics.x₂.η[domain],
      mesh.cell_center_metrics.x₂.ζ[domain],
    )

    vtk["dz_di", VTKCellData(), component_names=["xi", "eta", "zeta"]] = (
      mesh.cell_center_metrics.x₃.ξ[domain],
      mesh.cell_center_metrics.x₃.η[domain],
      mesh.cell_center_metrics.x₃.ζ[domain],
    )
  end
end

"""Write (x,y,z) coordinates to .vtk"""
function save_vtk((x, y, z)::NTuple{3,AbstractArray{T,3}}, fn="mesh") where {T}
  @info "Writing to $fn.vti"

  vtk_grid(fn, x, y, z) do vtk
  end
end

"""Write (x,y) coordinates to .vtk"""
function save_vtk((x, y)::NTuple{2,AbstractMatrix{N}}, fn="mesh") where {N}
  @info "Writing to $fn.vti"

  vtk_grid(fn, x, y) do vtk
  end
end

function save_vtk(mesh::AbstractCurvilinearGrid2D, fn="mesh")
  @info "Writing to $fn.vti"

  xyz_n = coords(mesh)
  domain = mesh.iterators.cell.domain

  @views vtk_grid(fn, xyz_n) do vtk
    vtk["J", VTKCellData()] = mesh.cell_center_metrics.J[domain]

    # vtk["volume", VTKCellData()] = cellvolume.(Ref(mesh), domain)

    vtk["xi", VTKCellData(), component_names=["x1", "x2", "t"]] = (
      mesh.cell_center_metrics.ξ.x₁[domain],
      mesh.cell_center_metrics.ξ.x₂[domain],
      mesh.cell_center_metrics.ξ.t[domain],
    )

    vtk["eta", VTKCellData(), component_names=["x1", "x2", "t"]] = (
      mesh.cell_center_metrics.η.x₁[domain],
      mesh.cell_center_metrics.η.x₂[domain],
      mesh.cell_center_metrics.η.t[domain],
    )

    vtk["dx_di", VTKCellData(), component_names=["xi", "eta"]] = (
      mesh.cell_center_metrics.x₁.ξ[domain], mesh.cell_center_metrics.x₁.η[domain]
    )

    vtk["dy_di", VTKCellData(), component_names=["xi", "eta"]] = (
      mesh.cell_center_metrics.x₂.ξ[domain], mesh.cell_center_metrics.x₂.η[domain]
    )
  end
end

end
