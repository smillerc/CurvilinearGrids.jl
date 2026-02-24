"""
Conservative remapping support for unified mapped/discrete grids.
"""
module RemappingSchemes

using LinearAlgebra
using SparseArrays
using StaticArrays

import ..GridTypes: AbstractMappedOrDiscreteGrid
import ..GridTypes:
  CoordinateSystemTrait,
  CartesianCS,
  CurvilinearCS,
  CylindricalCS,
  AxisymmetricCS,
  SphericalCS
import ..GridTypes:
  coordinate_system, coord, computational_coordinate, cellvolumes, cellsize

export RemapCache
export build_remap_cache
export remap_scalar, remap_scalar!
export source_overlap_mass, validate_remap_cache

include("types.jl")
include("conservative_scalar.jl")

end
