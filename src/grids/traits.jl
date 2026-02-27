#
# Grid trait system shared by both orthogonal and unified grids.
#

abstract type CoordinateSystemTrait end
struct CartesianCS <: CoordinateSystemTrait end
struct CylindricalCS <: CoordinateSystemTrait end
struct SphericalCS <: CoordinateSystemTrait end
struct AxisymmetricCS{Axis} <: CoordinateSystemTrait end
struct CurvilinearCS <: CoordinateSystemTrait end

abstract type BasisTrait end
struct CartesianBasis <: BasisTrait end
struct SphericalBasis <: BasisTrait end
