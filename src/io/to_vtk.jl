module VTKOutput

using WriteVTK

using ..GridTypes

export to_vtk

"""Write the mesh to .VTK format"""
function to_vtk(mesh::AbstractCurvilinearMesh, filename)
  xyz = coords(mesh)
  cell_centered_data = zeros((mesh.nnodes .- 1)...)

  @inbounds for I in CartesianIndices(cell_centered_data)
    cell_centered_data[I] = abs(jacobian(mesh, I))
  end

  vtk_grid(filename, xyz) do vtk
    vtk["volume"] = cell_centered_data
  end
end

end
