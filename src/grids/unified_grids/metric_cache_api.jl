#
# Metric cache API
#

"""
    invalidate_cell_metrics!(grid::AbstractMappedOrDiscreteGrid)

Mark cell-metric cache entries as stale.

# Arguments
  - `grid`: Target mapped/discrete grid.

# Returns
`nothing`.
"""
function invalidate_cell_metrics!(grid::AbstractMappedOrDiscreteGrid)
    if !_has_metric_storage(grid)
        return nothing
    end
    grid.metric_caches.cell.valid = false
    return nothing
end

"""
    invalidate_face_metrics!(grid::AbstractMappedOrDiscreteGrid)

Mark face-metric cache entries as stale.

# Arguments
  - `grid`: Target mapped/discrete grid.

# Returns
`nothing`.
"""
function invalidate_face_metrics!(grid::AbstractMappedOrDiscreteGrid)
    if !_has_metric_storage(grid)
        return nothing
    end
    grid.metric_caches.face.valid = false
    return nothing
end

"""
    refresh_cell_metrics!(grid::AbstractMappedOrDiscreteGrid; include_halo_region=false)

Recompute and return cell metrics for a mapped/discrete unified grid.

# Arguments
  - `grid`: Target mapped/discrete grid.

# Keywords
  - `include_halo_region`: Reserved compatibility flag. Default: `false`.

# Returns
Cell metric cache payload.
"""
function refresh_cell_metrics!(
        grid::AbstractMappedOrDiscreteGrid; include_halo_region::Bool = false
    )
    _require_metric_storage(grid, "refresh_cell_metrics!")
    return _refresh_cell_metrics!(grid; include_halo_region = include_halo_region)
end

"""
    refresh_face_metrics!(grid::AbstractMappedOrDiscreteGrid; include_halo_region=false)

Recompute and return face metrics for a mapped/discrete unified grid.

# Arguments
  - `grid`: Target mapped/discrete grid.

# Keywords
  - `include_halo_region`: Reserved compatibility flag. Default: `false`.

# Returns
Face metric cache payload.
"""
function refresh_face_metrics!(
        grid::AbstractMappedOrDiscreteGrid; include_halo_region::Bool = false
    )
    _require_metric_storage(grid, "refresh_face_metrics!")
    data = _refresh_face_metrics!(grid; include_halo_region = include_halo_region)
    return map(_chart_face_metrics, data)
end

"""
    cell_metrics(grid::AbstractMappedOrDiscreteGrid; refresh=false)

Access cell metric cache data.

# Arguments
  - `grid`: Target mapped/discrete grid.

# Keywords
  - `refresh`: Recompute before returning data. Default: `false`.

# Returns
Cell metric cache payload.
"""
function cell_metrics(grid::AbstractMappedOrDiscreteGrid; refresh::Bool = false)
    _require_metric_storage(grid, "cell_metrics")
    if refresh || !grid.metric_caches.cell.valid
        return refresh_cell_metrics!(grid)
    end
    return grid.metric_caches.cell.data
end

"""
    face_metrics(grid::AbstractMappedOrDiscreteGrid; refresh=false)

Access face metric cache data.
The conserved rows in this payload use the computational chart. Use
`physical_conserved_face_metrics` for the rows that define physical vector areas.

# Arguments
  - `grid`: Target mapped/discrete grid.

# Keywords
  - `refresh`: Recompute before returning data. Default: `false`.

# Returns
Face metric cache payload.
"""
function face_metrics(grid::AbstractMappedOrDiscreteGrid; refresh::Bool = false)
    _require_metric_storage(grid, "face_metrics")
    if refresh || !grid.metric_caches.face.valid
        return refresh_face_metrics!(grid)
    end
    return map(_chart_face_metrics, grid.metric_caches.face.data)
end

@inline _chart_face_metrics(face) =
    (; face.forward, face.inverse, face.conserved)

"""
    physical_conserved_face_metrics(grid)

Borrow the active conserved face rows that define physical vector areas.
Multiply each row by `conservation_face_metric_component_scale` at that
face to get its vector area in the grid's native basis. This accessor
refreshes a stale face cache.
"""
function physical_conserved_face_metrics(grid::AbstractMappedOrDiscreteGrid)
    _require_metric_storage(grid, "physical_conserved_face_metrics")
    if !grid.metric_caches.face.valid
        _refresh_face_metrics!(grid)
    end
    return map(grid.metric_caches.face.data) do face
        isnothing(face.physical_conserved) ? face.conserved : face.physical_conserved
    end
end

function cell_metrics(::AbstractOrthogonalGrid; refresh::Bool = false)
    throw(ArgumentError("`cell_metrics` is undefined for `AbstractOrthogonalGrid`."))
end

function face_metrics(::AbstractOrthogonalGrid; refresh::Bool = false)
    throw(ArgumentError("`face_metrics` is undefined for `AbstractOrthogonalGrid`."))
end
