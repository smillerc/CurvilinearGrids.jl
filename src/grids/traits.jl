#
# Grid trait system shared by both orthogonal and unified grids.
#

abstract type CoordinateSystemTrait end
struct CartesianCS <: CoordinateSystemTrait end
struct CylindricalCS <: CoordinateSystemTrait end
struct SphericalCS <: CoordinateSystemTrait end
struct AxisymmetricCS{Axis} <: CoordinateSystemTrait end

"""Grid geometry with no certified separation into orthogonal coordinate axes."""
struct GeneralGeometry end

"""Grid geometry with independent coordinates along native orthogonal axes."""
struct OrthogonalGeometry end

@inline _orthogonal_chart(::CoordinateSystemTrait, ::Val) = GeneralGeometry()
@inline _orthogonal_chart(::Union{CartesianCS, CylindricalCS, SphericalCS}, ::Union{Val{1}, Val{2}, Val{3}}) =
    OrthogonalGeometry()
@inline _orthogonal_chart(::Union{AxisymmetricCS{:x}, AxisymmetricCS{:y}}, ::Val{2}) =
    OrthogonalGeometry()

abstract type BasisTrait end
struct CartesianBasis <: BasisTrait end
struct SphericalBasis <: BasisTrait end
