"""
Multi-block connectivity support for unified grids.
"""
module MultiBlockMeshes

using LinearAlgebra
using StaticArrays

import ..GridTypes: AbstractMappedOrDiscreteGrid
import ..GridTypes: CoordinateSystemTrait, CartesianCS, CurvilinearCS, SphericalCS
import ..GridTypes: BasisTrait, CartesianBasis, SphericalBasis
import ..GridTypes: coordinate_system, basis_trait
import ..GridTypes: coord, centroid, forward_cell_metrics
import ..GridTypes: computational_coordinate

export BlockFace
export BlockInterface
export MultiBlockMesh
export validate_multiblock!, build_interface_caches!, invalidate_interface_caches!
export exchange_interface!, exchange_all_interfaces!

include("types.jl")
include("validation.jl")
include("cache.jl")
include("transfer.jl")

end
