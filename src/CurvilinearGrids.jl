module CurvilinearGrids

include("indexing_utils.jl")
using .IndexingUtils

include("metrics.jl")
using .MetricTypes

include("metric_schemes/metric_schemes.jl")
using .MetricDiscretizationSchemes
export MEG6Scheme, update_metrics!
export J
export ξx, ηx, ζx, ξy, ηy, ζy, ξz, ηz, ζz
export xξ, xη, xζ, yξ, yη, yζ, zξ, zη, zζ
export ξ̂x, η̂x, ζ̂x, ξ̂y, η̂y, ζ̂y, ξ̂z, η̂z, ζ̂z
export ξ̂xᵢ₊½, η̂xᵢ₊½, ζ̂xᵢ₊½, ξ̂yᵢ₊½, η̂yᵢ₊½, ζ̂yᵢ₊½, ξ̂zᵢ₊½, η̂zᵢ₊½, ζ̂zᵢ₊½
export ξ̂xⱼ₊½, η̂xⱼ₊½, ζ̂xⱼ₊½, ξ̂yⱼ₊½, η̂yⱼ₊½, ζ̂yⱼ₊½, ξ̂zⱼ₊½, η̂zⱼ₊½, ζ̂zⱼ₊½
export ξ̂xₖ₊½, η̂xₖ₊½, ζ̂xₖ₊½, ξ̂yₖ₊½, η̂yₖ₊½, ζ̂yₖ₊½, ξ̂zₖ₊½, η̂zₖ₊½, ζ̂zₖ₊½

include("metric_schemes/finitediff_metrics.jl")
using .FiniteDifferenceMetrics

include("grids/GridTypes.jl")
using .GridTypes
export AbstractCurvilinearGrid
export AbstractCurvilinearGrid1D
export AbstractCurvilinearGrid2D
export AbstractCurvilinearGrid3D
export CurvilinearGrid1D, CurvilinearGrid2D, CurvilinearGrid3D
export CylindricalGrid1D, SphericalGrid1D
export AxisymmetricGrid2D
export RectlinearGrid, RThetaGrid, RThetaPhiGrid
export AxisymmetricRectlinearGrid, AxisymmetricRThetaGrid

export update!
export cellsize, cellsize_withhalo
export coord, coords
export centroid, centroids
export metrics, jacobian, jacobian_matrix
export conservative_metrics
export metrics_with_jacobian
export cell_metrics, cell_indices
export radius
export cellvolume

include("adapt_to_gpu.jl")

include("io/to_vtk.jl")
using .VTKOutput
export save_vtk

end
