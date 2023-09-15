module VTKOutput

using WriteVTK

using ..GridTypes

export to_vtk

"""Write the mesh to .VTK format"""
function to_vtk(mesh::AbstractCurvilinearMesh, filename)
  xyz = coords(mesh)
  vtk_grid(filename, xyz) do vtk
    # add datasets...
  end
end

end
