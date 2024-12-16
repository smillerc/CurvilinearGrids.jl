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
    vtk["J", VTKCellData()] = mesh.cell_center_metrics.forward.J[domain]

    # vtk["volume", VTKCellData()] = cellvolume.(Ref(mesh), domain)

    vtk["xi", VTKCellData(), component_names=["x1", "x2", "x3", "t"]] = (
      mesh.cell_center_metrics.inverse.ξ.x₁[domain],
      mesh.cell_center_metrics.inverse.ξ.x₂[domain],
      mesh.cell_center_metrics.inverse.ξ.x₃[domain],
      mesh.cell_center_metrics.inverse.ξ.t[domain],
    )

    vtk["eta", VTKCellData(), component_names=["x1", "x2", "x3", "t"]] = (
      mesh.cell_center_metrics.inverse.η.x₁[domain],
      mesh.cell_center_metrics.inverse.η.x₂[domain],
      mesh.cell_center_metrics.inverse.η.x₃[domain],
      mesh.cell_center_metrics.inverse.η.t[domain],
    )

    vtk["zeta", VTKCellData(), component_names=["x1", "x2", "x3", "t"]] = (
      mesh.cell_center_metrics.inverse.ζ.x₁[domain],
      mesh.cell_center_metrics.inverse.ζ.x₂[domain],
      mesh.cell_center_metrics.inverse.ζ.x₃[domain],
      mesh.cell_center_metrics.inverse.ζ.t[domain],
    )

    vtk["dx_di", VTKCellData(), component_names=["xi", "eta", "zeta"]] = (
      mesh.cell_center_metrics.forward.x₁.ξ[domain],
      mesh.cell_center_metrics.forward.x₁.η[domain],
      mesh.cell_center_metrics.forward.x₁.ζ[domain],
    )

    vtk["dy_di", VTKCellData(), component_names=["xi", "eta", "zeta"]] = (
      mesh.cell_center_metrics.forward.x₂.ξ[domain],
      mesh.cell_center_metrics.forward.x₂.η[domain],
      mesh.cell_center_metrics.forward.x₂.ζ[domain],
    )

    vtk["dz_di", VTKCellData(), component_names=["xi", "eta", "zeta"]] = (
      mesh.cell_center_metrics.forward.x₃.ξ[domain],
      mesh.cell_center_metrics.forward.x₃.η[domain],
      mesh.cell_center_metrics.forward.x₃.ζ[domain],
    )
  end
end

function save_vtk(mesh::AbstractCurvilinearGrid2D, fn="mesh")
  @info "Writing to $fn.vti"

  xyz_n = coords(mesh)
  domain = mesh.iterators.cell.domain

  @views vtk_grid(fn, xyz_n) do vtk
    vtk["J", VTKCellData()] = mesh.cell_center_metrics.forward.J[domain]

    # vtk["volume", VTKCellData()] = cellvolume.(Ref(mesh), domain)

    vtk["xi", VTKCellData(), component_names=["x1", "x2", "t"]] = (
      mesh.cell_center_metrics.inverse.ξ.x₁[domain],
      mesh.cell_center_metrics.inverse.ξ.x₂[domain],
      mesh.cell_center_metrics.inverse.ξ.t[domain],
    )

    vtk["eta", VTKCellData(), component_names=["x1", "x2", "t"]] = (
      mesh.cell_center_metrics.inverse.η.x₁[domain],
      mesh.cell_center_metrics.inverse.η.x₂[domain],
      mesh.cell_center_metrics.inverse.η.t[domain],
    )
  end
end

end
