module VTKOutput

using WriteVTK

using ..GridTypes

export save_vtk

"""Write the mesh to .VTK format"""
function to_vtk(mesh::AbstractCurvilinearGrid, filename)
  xyz = coords(mesh)

  vtk_grid(filename, xyz) do vtk
    vtk["volume"] = mesh.J
  end
end

function save_vtk(mesh::CurvilinearGrid2D)
  fn = "mesh"
  @info "Writing to $fn.vti"

  xyz_n = coords(mesh)
  domain = mesh.iterators.cell.domain

  @views vtk_grid(fn, xyz_n) do vtk
    vtk["J", VTKCellData()] = mesh.cell_center_metrics.J[domain]

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
