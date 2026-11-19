@inline function _axisymmetric_axial_integrand(radial_map, jacobian, t, ξ, params, axial, tangent)
    radius = radial_map(t, ξ..., params)
    F = _as_smatrix(Val(2), jacobian(t, ξ..., params))
    return abs(radius) * F[axial, tangent]
end

@inline function _axisymmetric_conserved_face(
        maps, jacobian, quadrature, t, params, global_cell, nhalo, axis, ::Val{R}, ::Type{T},
    ) where {R, T}
    half = T(0.5)
    ξ = SVector{2, T}(
        ntuple(d -> global_cell[d] - nhalo + (d == axis ? one(T) : half), Val(2))
    )
    tangent = 3 - axis
    axial = 3 - R
    radial_map = maps[R]

    lower = setindex(ξ, ξ[tangent] - half, tangent)
    upper = setindex(ξ, ξ[tangent] + half, tangent)
    radius_lower = radial_map(t, lower..., params)
    radius_upper = radial_map(t, upper..., params)
    integral_r = half * (radius_upper * abs(radius_upper) - radius_lower * abs(radius_lower))

    integral_z = zero(T)
    for (offset, weight) in quadrature
        lower = setindex(ξ, ξ[tangent] - offset, tangent)
        upper = setindex(ξ, ξ[tangent] + offset, tangent)
        integral_z += weight * (
            _axisymmetric_axial_integrand(radial_map, jacobian, t, lower, params, axial, tangent) +
                _axisymmetric_axial_integrand(radial_map, jacobian, t, upper, params, axial, tangent)
        )
    end

    # Positive computational face vectors are (|R|dy, -|R|dx) on i and its
    # negative on j. The factor 2π cancels against the face scale.
    row = R == 1 ? SVector(integral_z, -integral_r) : SVector(integral_r, -integral_z)
    radius = T(radial_map(t, ξ..., params))
    if iszero(radius)
        all(iszero, row) ||
            throw(DomainError(radius, "An axisymmetric face with zero radius has nonzero integrated measure"))
    else
        row = (axis == 1 ? one(T) : -one(T)) * row / abs(radius)
    end
    matrix = axis == 1 ? SMatrix{2, 2, T, 4}(row[1], zero(T), row[2], zero(T)) :
        SMatrix{2, 2, T, 4}(zero(T), row[1], zero(T), row[2])
    return ConservedMetric(matrix)
end

@kernel function _fill_axisymmetric_face_rows_kernel!(
        conserved, maps, jacobian, quadrature, faces, global_cells, nhalo, t, params, axis, radial, ::Type{T},
    ) where {T}
    i = @index(Global, Linear)
    I = faces[i]
    conserved[I] = _axisymmetric_conserved_face(
        maps, jacobian, quadrature, t, params, global_cells[I], nhalo, axis, radial, T,
    )
end

function _fill_axisymmetric_conserved_faces!(storage, grid, t, params, ::Val{N}, ::Type{T}) where {N, T}
    N == 2 || return nothing
    cs = coordinate_system(grid)
    cs isa Union{AxisymmetricCS{:x}, AxisymmetricCS{:y}} || return nothing
    basis_trait(grid) isa CartesianBasis || return nothing

    radial = cs isa AxisymmetricCS{:x} ? Val(2) : Val(1)
    owned = grid.iterators.cell.domain.indices
    global_cells = grid.iterators.global_domain.cell.full
    maps = grid.mapping_functions
    jacobian = grid.metric_functions_cache.forward.jacobian
    quadrature = _face_edge_gauss_legendre_rule(T)
    for axis in 1:2
        lower = max(1, first(owned[axis]) - 1)
        faces = CartesianIndices(setindex(owned, lower:last(owned[axis]), axis))
        conserved = storage[axis].physical_conserved
        copyto!(conserved, storage[axis].conserved)
        if grid.backend isa CPU
            Threads.@threads for i in eachindex(faces)
                I = faces[i]
                conserved[I] = _axisymmetric_conserved_face(
                    maps, jacobian, quadrature, t, params, global_cells[I], grid.nhalo, axis, radial, T,
                )
            end
        else
            _fill_axisymmetric_face_rows_kernel!(grid.backend)(
                conserved, maps, jacobian, quadrature, faces, global_cells, grid.nhalo, t, params, axis, radial, T;
                ndrange = length(faces),
            )
        end
    end
    KernelAbstractions.synchronize(grid.backend)
    return nothing
end
