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

end
