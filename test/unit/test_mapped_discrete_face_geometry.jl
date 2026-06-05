function _mapped_discrete_face_geometry_cases()
  return (
    (; name="1D CurvilinearCS", dim=Val(1), coordinate_system=CurvilinearCS(), basis=CartesianBasis()),
    (; name="1D CartesianCS", dim=Val(1), coordinate_system=CartesianCS(), basis=CartesianBasis()),
    (; name="1D CylindricalCS", dim=Val(1), coordinate_system=CylindricalCS(), basis=CartesianBasis()),
    (; name="1D SphericalCS", dim=Val(1), coordinate_system=SphericalCS(), basis=SphericalBasis()),
    (; name="2D CurvilinearCS", dim=Val(2), coordinate_system=CurvilinearCS(), basis=CartesianBasis()),
    (; name="2D CartesianCS", dim=Val(2), coordinate_system=CartesianCS(), basis=CartesianBasis()),
    (; name="2D CylindricalCS", dim=Val(2), coordinate_system=CylindricalCS(), basis=CartesianBasis()),
    (;
      name="2D AxisymmetricCS{:x}",
      dim=Val(2),
      coordinate_system=AxisymmetricCS{:x}(),
      basis=CartesianBasis(),
    ),
    (;
      name="2D AxisymmetricCS{:y}",
      dim=Val(2),
      coordinate_system=AxisymmetricCS{:y}(),
      basis=CartesianBasis(),
    ),
    (; name="2D SphericalCS", dim=Val(2), coordinate_system=SphericalCS(), basis=SphericalBasis()),
    (; name="3D CurvilinearCS", dim=Val(3), coordinate_system=CurvilinearCS(), basis=CartesianBasis()),
    (; name="3D CartesianCS", dim=Val(3), coordinate_system=CartesianCS(), basis=CartesianBasis()),
    (; name="3D SphericalCS", dim=Val(3), coordinate_system=SphericalCS(), basis=SphericalBasis()),
  )
end

_mapped_discrete_face_geometry_celldims(::Val{1}) = ((8,), (16,), (32,))
_mapped_discrete_face_geometry_celldims(::Val{2}) = ((8, 8), (16, 16), (32, 32))
_mapped_discrete_face_geometry_celldims(::Val{3}) = ((4, 4, 4), (8, 8, 8), (12, 12, 12))

_mapped_discrete_face_locations(::Val{1}) = (:ihi,)
_mapped_discrete_face_locations(::Val{2}) = (:ihi, :jhi)
_mapped_discrete_face_locations(::Val{3}) = (:ihi, :jhi, :khi)

function _mapped_discrete_grid_pair(
  ::Val{1}, celldims, coordinate_system::CoordinateSystemTrait, basis::BasisTrait
)
  ni, = celldims
  nhalo = 2
  params = (; ni, q0=1.1, dq=0.8 / ni, amp=0.02)
  qmap(t, i, p) = p.q0 + (i - 1) * p.dq + p.amp * sinpi(2 * (i - 1) / p.ni)

  mapped = MappedGrid(
    qmap,
    params,
    celldims,
    nhalo;
    backend=CPU(),
    coordinate_system=coordinate_system,
    basis=basis,
    cache_mode=:eager,
  )
  qnodes = mapped.node_coordinates[1][mapped.iterators.node.full]
  discrete = DiscreteGrid(
    qnodes,
    nhalo;
    coordinate_system=coordinate_system,
    basis=basis,
    halo_coords_included=true,
    cache_mode=:eager,
  )
  return mapped, discrete
end

function _mapped_discrete_grid_pair(
  ::Val{2}, celldims, coordinate_system::CoordinateSystemTrait, basis::BasisTrait
)
  ni, nj = celldims
  nhalo = 2
  params = (; ni, nj, q10=1.0, q20=0.7, dq1=0.8 / ni, dq2=0.6 / nj, a1=0.025, a2=0.018)
  q1map(t, i, j, p) = p.q10 + (i - 1) * p.dq1 + p.a1 * sinpi(2 * (j - 1) / p.nj)
  q2map(t, i, j, p) = p.q20 + (j - 1) * p.dq2 + p.a2 * sinpi(2 * (i - 1) / p.ni)

  mapped = MappedGrid(
    q1map,
    q2map,
    params,
    celldims,
    nhalo;
    backend=CPU(),
    coordinate_system=coordinate_system,
    basis=basis,
    cache_mode=:eager,
  )
  q1nodes = mapped.node_coordinates[1][mapped.iterators.node.full]
  q2nodes = mapped.node_coordinates[2][mapped.iterators.node.full]
  discrete = DiscreteGrid(
    q1nodes,
    q2nodes,
    nhalo;
    coordinate_system=coordinate_system,
    basis=basis,
    halo_coords_included=true,
    cache_mode=:eager,
  )
  return mapped, discrete
end

function _mapped_discrete_grid_pair_3d(
  q1map, q2map, q3map, params, celldims, coordinate_system::CoordinateSystemTrait, basis::BasisTrait
)
  nhalo = 2
  mapped = MappedGrid(
    q1map,
    q2map,
    q3map,
    params,
    celldims,
    nhalo;
    backend=CPU(),
    coordinate_system=coordinate_system,
    basis=basis,
    diff_backend=AutoForwardDiff(),
    conserved_metric_scheme=ADThomasLombardMetric(),
    cache_mode=:eager,
  )
  q1nodes = mapped.node_coordinates[1][mapped.iterators.node.full]
  q2nodes = mapped.node_coordinates[2][mapped.iterators.node.full]
  q3nodes = mapped.node_coordinates[3][mapped.iterators.node.full]
  discrete = DiscreteGrid(
    q1nodes,
    q2nodes,
    q3nodes,
    nhalo;
    coordinate_system=coordinate_system,
    basis=basis,
    halo_coords_included=true,
    cache_mode=:eager,
  )
  return mapped, discrete
end

function _mapped_discrete_grid_pair(
  ::Val{3}, celldims, coordinate_system::SphericalCS, basis::SphericalBasis
)
  ni, nj, nk = celldims
  params = (;
    ni,
    nj,
    nk,
    r0=1.2,
    theta0=0.45,
    phi0=0.25,
    dr=0.4 / ni,
    dtheta=0.35 / nj,
    dphi=0.45 / nk,
    ar=0.012,
    atheta=0.008,
    aphi=0.006,
  )
  q1map(t, i, j, k, p) =
    p.r0 + (i - 1) * p.dr + p.ar * sinpi(2 * (j - 1) / p.nj) * sinpi((k - 1) / p.nk)
  q2map(t, i, j, k, p) =
    p.theta0 + (j - 1) * p.dtheta + p.atheta * sinpi(2 * (i - 1) / p.ni)
  q3map(t, i, j, k, p) =
    p.phi0 + (k - 1) * p.dphi + p.aphi * sinpi(2 * (i - 1) / p.ni) * sinpi(2 * (j - 1) / p.nj)

  return _mapped_discrete_grid_pair_3d(
    q1map, q2map, q3map, params, celldims, coordinate_system, basis
  )
end

function _mapped_discrete_grid_pair(
  ::Val{3}, celldims, coordinate_system::CoordinateSystemTrait, basis::BasisTrait
)
  ni, nj, nk = celldims
  params = (;
    ni,
    nj,
    nk,
    q10=0.6,
    q20=0.8,
    q30=1.0,
    dq1=0.7 / ni,
    dq2=0.6 / nj,
    dq3=0.5 / nk,
    a1=0.014,
    a2=0.011,
    a3=0.009,
  )
  q1map(t, i, j, k, p) =
    p.q10 + (i - 1) * p.dq1 + p.a1 * sinpi(2 * (j - 1) / p.nj) * sinpi(2 * (k - 1) / p.nk)
  q2map(t, i, j, k, p) =
    p.q20 + (j - 1) * p.dq2 + p.a2 * sinpi(2 * (i - 1) / p.ni) * sinpi(2 * (k - 1) / p.nk)
  q3map(t, i, j, k, p) =
    p.q30 + (k - 1) * p.dq3 + p.a3 * sinpi(2 * (i - 1) / p.ni) * sinpi(2 * (j - 1) / p.nj)

  return _mapped_discrete_grid_pair_3d(
    q1map, q2map, q3map, params, celldims, coordinate_system, basis
  )
end

function _mapped_discrete_max_gcl(grid)
  domain = grid.iterators.cell.domain
  residuals = CurvilinearGrids.GridTypes.gcl(face_metrics(grid), domain)
  residual_tuple = residuals isa Tuple ? residuals : (residuals,)
  max_residual = 0.0
  for residual in residual_tuple
    max_residual = max(max_residual, maximum(abs, residual[domain]))
  end
  return max_residual
end

function _mapped_discrete_face_geometry_errors(case, celldims)
  mapped, discrete = _mapped_discrete_grid_pair(
    case.dim, celldims, case.coordinate_system, case.basis
  )

  @test discrete.iterators.cell.domain == mapped.iterators.cell.domain
  @test _mapped_discrete_max_gcl(mapped) < 1e-12
  @test _mapped_discrete_max_gcl(discrete) < 1e-12

  max_metric = 0.0
  max_area = 0.0
  max_normal = 0.0
  max_coordinate = 0.0
  for I in mapped.iterators.cell.domain
    for loc in _mapped_discrete_face_locations(case.dim)
      mapped_geom = face_flux_geometry(mapped, I.I, loc)
      discrete_geom = face_flux_geometry(discrete, I.I, loc)
      max_metric = max(max_metric, norm(discrete_geom.metric_vector - mapped_geom.metric_vector))
      max_area = max(max_area, abs(discrete_geom.area - mapped_geom.area))
      max_normal = max(max_normal, norm(discrete_geom.normal - mapped_geom.normal))
      max_coordinate = max(max_coordinate, norm(discrete_geom.coordinate - mapped_geom.coordinate))
    end
  end

  return (; celldims, max_metric, max_area, max_normal, max_coordinate)
end

function _mapped_discrete_refinement_orders(errors, field)
  return map(zip(errors[1:(end - 1)], errors[2:end])) do (coarse, fine)
    resolution_ratio = first(fine.celldims) / first(coarse.celldims)
    log(getproperty(coarse, field) / getproperty(fine, field)) / log(resolution_ratio)
  end
end

function _test_mapped_discrete_refines(errors, field)
  values = map(error -> getproperty(error, field), errors)
  @test all(isfinite, values)
  if maximum(values) < 1e-12
    @test true
    return nothing
  end

  @test values[2] < values[1]
  @test values[3] < values[2]
  @test values[end] < 0.35 * values[1]
  observed_order = minimum(_mapped_discrete_refinement_orders(errors, field))
  @test observed_order > 1.0
  return nothing
end

@testset "MappedGrid vs DiscreteGrid face geometry across coordinate systems" begin
  for case in _mapped_discrete_face_geometry_cases()
    @testset "$(case.name)" begin
      errors = map(
        celldims -> _mapped_discrete_face_geometry_errors(case, celldims),
        _mapped_discrete_face_geometry_celldims(case.dim),
      )

      _test_mapped_discrete_refines(errors, :max_metric)
      _test_mapped_discrete_refines(errors, :max_area)
      _test_mapped_discrete_refines(errors, :max_normal)
      _test_mapped_discrete_refines(errors, :max_coordinate)
    end
  end
end
