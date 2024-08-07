module VTKOutput

using WriteVTK

using ..GridTypes

export save_vtk

"""Write the mesh to .VTK format"""
function save_vtk(mesh::CurvilinearGrid3D, fn="mesh")
  @info "Writing to $fn.vti"

  xyz_n = coords(mesh)
  domain = mesh.iterators.cell.domain

  @views vtk_grid(fn, xyz_n) do vtk
    vtk["J", VTKCellData()] = mesh.cell_center_metrics.J[domain]
  end
end

function save_vtk(mesh::AbstractCurvilinearGrid2D, fn="mesh")
  @info "Writing to $fn.vti"

  xyz_n = coords(mesh)
  domain = mesh.iterators.cell.domain

  @views vtk_grid(fn, xyz_n) do vtk
    vtk["J", VTKCellData()] = mesh.cell_center_metrics.J[domain]

    vtk["volume", VTKCellData()] = cellvolume.(Ref(mesh), domain)

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
  end
end

end
