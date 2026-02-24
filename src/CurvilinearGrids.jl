module CurvilinearGrids

include("discretization_schemes/DiscretizationSchemes.jl")
using .DiscretizationSchemes
export MonotoneExplicitGradientScheme
export SecondOrder, FourthOrder, SixthOrder, EighthOrder
export compute_first_derivatives!, compute_second_derivatives!
export cell_center_derivatives!, interpolate_to_edge!

include("metric_schemes/MetricDiscretizationSchemes.jl")
using .MetricDiscretizationSchemes

include("RectilinearArrays.jl")
using .RectilinearArrays

include("grids/GridTypes.jl")
using .GridTypes

include("multiblock/multiblock.jl")
using .MultiBlockMeshes: BlockFace, BlockInterface
using .MultiBlockMeshes:
  validate_multiblock!, build_interface_caches!, invalidate_interface_caches!
using .MultiBlockMeshes: exchange_interface!, exchange_all_interfaces!
using .MultiBlockMeshes: MultiBlockMesh

include("remapping/remapping.jl")
using .RemappingSchemes: RemapCache, build_remap_cache, remap_scalar, remap_scalar!
using .RemappingSchemes: source_overlap_mass, validate_remap_cache

export AbstractCurvilinearGrid
export AbstractCurvilinearGrid1D
export AbstractCurvilinearGrid2D
export AbstractCurvilinearGrid3D
export CurvilinearGrid1D, CurvilinearGrid2D, CurvilinearGrid3D
export RectilinearGrid2D, RectilinearGrid3D
export UniformGrid1D, UniformGrid2D, UniformGrid3D
export CylindricalGrid1D, SphericalGrid1D
export AxisymmetricGrid2D
export SphericalGrid3D, SphericalBasisCurvilinearGrid3D
export CartesianOrthogonalGrid1D,
  CylindricalOrthogonalGrid1D, SphericalOrthogonalGrid1D, AxisymmetricOrthogonalGrid2D

export AbstractUnifiedGrid
export MappedGrid, DiscreteGrid, OrthogonalGrid
export Metric, ConservedMetric

export CoordinateSystemTrait
export CartesianCS, CylindricalCS, SphericalCS, AxisymmetricCS, CurvilinearCS

export BasisTrait
export CartesianBasis, SphericalBasis
export EdgeInterpolationSchemeTrait
export EdgeInterpolationOrder1, EdgeInterpolationOrder2, EdgeInterpolationOrder3
export coordinate_system, basis_trait

export rectilinear_grid,
  rtheta_grid, rthetaphi_grid, rectilinear_cylindrical_grid, rectilinear_spherical_grid
export axisymmetric_rectilinear_grid, axisymmetric_rtheta_grid

export update!
export cellsize, cellsize_withhalo
export coord, coords
export centroid, centroids, cartesian_centroid
export metrics, jacobian, jacobian_matrix
export conservative_metrics
export metrics_with_jacobian
export cell_metrics, face_metrics, cell_indices
export radius, centroid_radius, centroid_radii
export cellvolume, cellvolumes
export face_area, outward_face_normal

export forward_cell_metrics, inverse_cell_metrics
export InverseCoordinateResult, computational_coordinate
export BlockFace, BlockInterface
export MultiBlockMesh
export validate_multiblock!, build_interface_caches!, invalidate_interface_caches!
export exchange_interface!, exchange_all_interfaces!
export invalidate_cell_metrics!, invalidate_face_metrics!
export refresh_cell_metrics!, refresh_face_metrics!
export RemapCache, build_remap_cache, remap_scalar, remap_scalar!
export source_overlap_mass, validate_remap_cache

include("operators/Operators.jl")
using .Operators
export cell_center_derivative, edge_derivative
# export cell_center_curl, edge_curl # not working yet
export cell_center_divergence, edge_divergence
export cell_center_gradient, edge_gradient

include("adapt_to_gpu.jl")

include("io/to_vtk.jl")
using .VTKOutput
export save_vtk

include("io/to_h5.jl")
using .H5Output
export write_coordinates, read_coordinates

include("mesh_functions/stretching_functions.jl")
export one_sided_stretch, double_sided_stretch, one_sided_with_initial_spacing

include("grids/surface_mesh.jl")
export SurfaceGrid, extract_surface_mesh

include("mesh_functions/estimate_yplus.jl")
export estimate_wall_distance

include("mesh_functions/rigid_body_transformations.jl")
export translate!, rotate!, scale!

include("remap.jl")
export change_resolution, scale_resolution, remap_cell_data

end
