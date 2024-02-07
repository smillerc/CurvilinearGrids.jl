module VTKOutput

using WriteVTK

using ..GridTypes

export to_vtk

"""Write the mesh to .VTK format"""
function to_vtk(mesh::AbstractCurvilinearGrid, filename)
  xyz = coords(mesh)

  vtk_grid(filename, xyz) do vtk
    vtk["volume"] = mesh.J
  end
end

function save_vtk(mesh)
  fn = "wavy"
  @info "Writing to $fn.vti"

  xyz_n = CurvilinearGrids.coords(mesh)
  domain = mesh.iterators.cell.domain

  @views vtk_grid(fn, xyz_n) do vtk
    for (k, v) in pairs(mesh.cell_center_metrics)
      vtk["$k"] = v[domain]
    end

    for (edge_name, data) in pairs(mesh.edge_metrics)
      for (k, v) in pairs(data)
        vtk["$(k)_$(edge_name)"] = v[domain]
      end
    end
  end
end

end
