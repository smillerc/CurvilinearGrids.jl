# Positive abscissas and weights of the eight-point Gauss-Legendre rule.
const _FACE_EDGE_GAUSS_LEGENDRE_RULE = (
    (0.1834346424956498, 0.362683783378362),
    (0.525532409916329, 0.3137066458778873),
    (0.7966664774136267, 0.2223810344533745),
    (0.9602898564975363, 0.1012285362903763),
)

# Convert on the host and map the rule from [-1, 1] to a unit-width face.
function _face_edge_gauss_legendre_rule(::Type{T}) where {T}
    return map(_FACE_EDGE_GAUSS_LEGENDRE_RULE) do (node, weight)
        (T(node / 2), T(weight / 2))
    end
end

_uses_weighted_face_metrics(cs, basis, dimension) = false
_uses_weighted_face_metrics(::SphericalCS, ::SphericalBasis, ::Union{Val{2}, Val{3}}) = true
_uses_weighted_face_metrics(::AxisymmetricCS, ::CartesianBasis, ::Val{2}) = true

@inline function _spherical_mapped_point(maps, t, ξ::SVector{N, T}, params) where {N, T}
    return SVector{N, T}(map(f -> f(t, ξ..., params), values(maps)))
end

@inline function _spherical_edge_sample(maps, jacobian, t, ξ::SVector{2, T}, params, axis) where {T}
    r, θ = _spherical_mapped_point(maps, t, ξ, params)
    F = _as_smatrix(Val(2), jacobian(t, ξ..., params))
    sθ = _spherical_face_polar_sine(θ)
    cθ = cos(θ)
    # In (R,Z) = (r sin θ, r cos θ), the radial edge integrand is -R dZ.
    return -r * sθ * (cθ * F[1, axis] - r * sθ * F[2, axis])
end

@inline function _spherical_edge_sample(maps, jacobian, t, ξ::SVector{3, T}, params, axis) where {T}
    r, θ, ϕ = _spherical_mapped_point(maps, t, ξ, params)
    F = _as_smatrix(Val(3), jacobian(t, ξ..., params))
    sθ = _spherical_face_polar_sine(θ)
    cθ = cos(θ)
    sϕ, cϕ = sincos(ϕ)
    radial = SVector(sθ * cϕ, sθ * sϕ, cθ)
    polar = SVector(cθ * cϕ, cθ * sϕ, -sθ)
    azimuthal = SVector(-sϕ, cϕ, zero(T))
    tangent = radial * F[1, axis] + polar * (r * F[2, axis]) + azimuthal * (r * sθ * F[3, axis])
    return cross(r * radial, tangent) / T(2)
end

@inline function _spherical_edge_integral(maps, jacobian, t, ξ::SVector{N, T}, params, axis, rule) where {N, T}
    integral = N == 2 ? zero(T) : zero(ξ)
    for (offset, weight) in rule
        lower = setindex(ξ, ξ[axis] - offset, axis)
        upper = setindex(ξ, ξ[axis] + offset, axis)
        integral += weight * (
            _spherical_edge_sample(maps, jacobian, t, lower, params, axis) +
                _spherical_edge_sample(maps, jacobian, t, upper, params, axis)
        )
    end
    return integral
end

@inline function _spherical_cylindrical_radius(q::SVector{2, T}) where {T}
    r, θ = q
    return r * _spherical_face_polar_sine(θ)
end

@inline function _spherical_face_vector(maps, jacobian, t, ξ::SVector{2, T}, params, axis, rule) where {T}
    tangent = 3 - axis
    lower = setindex(ξ, ξ[tangent] - T(0.5), tangent)
    upper = setindex(ξ, ξ[tangent] + T(0.5), tangent)
    r_lower = _spherical_cylindrical_radius(_spherical_mapped_point(maps, t, lower, params))
    r_upper = _spherical_cylindrical_radius(_spherical_mapped_point(maps, t, upper, params))
    radial = _spherical_edge_integral(maps, jacobian, t, ξ, params, tangent, rule)
    axial = ((r_upper - r_lower) * (r_upper + r_lower)) / T(2)
    sign = axis == 1 ? one(T) : -one(T)
    return (sign * T(2π)) * SVector(radial, axial)
end

@inline function _spherical_face_vector(maps, jacobian, t, ξ::SVector{3, T}, params, axis, rule) where {T}
    b = mod1(axis + 1, 3)
    c = mod1(axis + 2, 3)
    b_lower = setindex(ξ, ξ[b] - T(0.5), b)
    b_upper = setindex(ξ, ξ[b] + T(0.5), b)
    c_lower = setindex(ξ, ξ[c] - T(0.5), c)
    c_upper = setindex(ξ, ξ[c] + T(0.5), c)
    # Shared edge potentials cancel in the sum of Cartesian face vectors.
    return _spherical_edge_integral(maps, jacobian, t, b_upper, params, c, rule) -
        _spherical_edge_integral(maps, jacobian, t, b_lower, params, c, rule) -
        (
        _spherical_edge_integral(maps, jacobian, t, c_upper, params, b, rule) -
            _spherical_edge_integral(maps, jacobian, t, c_lower, params, b, rule)
    )
end

@inline function _spherical_conserved_row(maps, jacobian, t, ξ::SVector{N, T}, params, axis, rule) where {N, T}
    q = _spherical_mapped_point(maps, t, ξ, params)
    vector = _spherical_face_vector(maps, jacobian, t, ξ, params, axis, rule)
    basis = _basis_to_cartesian_matrix(SphericalCS(), Val(N), q, T)
    native = transpose(basis) * vector
    scale = conservation_face_metric_component_scale(SphericalCS(), SphericalBasis(), q)
    # Define the chart value at a true collapse by its zero physical extension.
    return map((v, w) -> iszero(v) && iszero(w) ? zero(T) : v / w, native, scale)
end

@kernel function _fill_spherical_face_rows_kernel!(
        conserved, maps, jacobian, local_faces, global_cells, nhalo, t, params,
        axis, rule, ::Val{N}, ::Type{T},
    ) where {N, T}
    i = @index(Global, Linear)
    I = local_faces[i]
    conserved[I] = _spherical_conserved_metric(
        maps, jacobian, t, global_cells[I], nhalo, params, axis, rule, Val(N), T
    )
end

@inline function _spherical_conserved_metric(
        maps, jacobian, t, global_cell, nhalo, params, axis, rule, ::Val{N}, ::Type{T},
    ) where {N, T}
    ξ = SVector{N, T}(
        ntuple(d -> global_cell[d] - nhalo + (d == axis ? one(T) : T(0.5)), Val(N))
    )
    row = _spherical_conserved_row(maps, jacobian, t, ξ, params, axis, rule)
    matrix = SMatrix{N, N, T, N * N}(
        ntuple(k -> mod1(k, N) == axis ? row[cld(k, N)] : zero(T), Val(N * N))
    )
    return ConservedMetric(matrix)
end

function _discrete_spherical_face_collapsed(nodes, I::CartesianIndex{N}, axis) where {N}
    origin = ntuple(d -> I[d] + (d == axis ? 1 : 0), Val(N))
    corners = CartesianIndices(ntuple(d -> d == axis ? (0:0) : (0:1), Val(N)))
    for component in 1:2
        coordinates = nodes[component]
        value = coordinates[origin...]
        if iszero(value) || (component == 2 && value == oftype(value, π))
            collapsed = all(corners) do corner
                index = ntuple(d -> origin[d] + corner[d], Val(N))
                coordinates[index...] == value
            end
            collapsed && return true
        end
    end
    return false
end

function _fill_spherical_conserved_faces!(storage, grid, t, params, dim::Val{N}, ::Type{T}) where {N, T}
    N in (2, 3) || return nothing
    coordinate_system(grid) isa SphericalCS || return nothing
    basis_trait(grid) isa SphericalBasis || return nothing

    owned = grid.iterators.cell.domain.indices
    global_cells = grid.iterators.global_domain.cell.full
    jacobian = grid.metric_functions_cache.forward.jacobian
    rule = _face_edge_gauss_legendre_rule(T)
    for axis in 1:N
        lower = max(1, first(owned[axis]) - 1)
        faces = CartesianIndices(setindex(owned, lower:last(owned[axis]), axis))
        conserved = storage[axis].physical_conserved
        copyto!(conserved, storage[axis].conserved)
        if grid.backend isa KernelAbstractions.CPU
            Threads.@threads for I in faces
                if grid isa DiscreteGrid && _discrete_spherical_face_collapsed(grid.node_coordinates, I, axis)
                    # A multilinear face lies in the exact nodal plane even when
                    # floating-point interpolation perturbs its constant coordinate.
                    conserved[I] = zero(eltype(conserved))
                else
                    conserved[I] = _spherical_conserved_metric(
                        grid.mapping_functions, jacobian, t, global_cells[I], grid.nhalo, params, axis, rule, dim, T
                    )
                end
            end
        else
            _fill_spherical_face_rows_kernel!(grid.backend)(
                conserved, grid.mapping_functions, jacobian, faces, global_cells,
                grid.nhalo, t, params, axis, rule, dim, T;
                ndrange = length(faces),
            )
        end
    end
    KernelAbstractions.synchronize(grid.backend)
    return nothing
end
