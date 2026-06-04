using CurvilinearGrids
using KernelAbstractions
using Test

const _METAL = @static if Sys.isapple()
  try
    pkgid = Base.PkgId(Base.UUID("dde4c033-4e86-420c-a63e-0dd931031962"), "Metal")
    Base.require(pkgid)
    Base.loaded_modules[pkgid]
  catch
    nothing
  end
else
  nothing
end

function _unified_test_backends()
  backends = Tuple{Symbol,Any}[(:cpu, CPU())]

  if !isnothing(_METAL)
    try
      metal_backend = _METAL.MetalBackend()
      probe = KernelAbstractions.zeros(metal_backend, Float32, 1)
      fill!(probe, 1.0f0)
      KernelAbstractions.synchronize(metal_backend)
      push!(backends, (:metal, metal_backend))
    catch err
      @warn "Skipping Metal backend in unified backend harness." error = sprint(
        showerror, err
      )
    end
  end

  return backends
end

const _UNIFIED_TEST_BACKENDS = _unified_test_backends()

function _host_face_metrics(face_storage::Tuple)
  ntuple(
    axis -> (;
      forward=Array(face_storage[axis].forward),
      inverse=Array(face_storage[axis].inverse),
      conserved=Array(face_storage[axis].conserved),
    ),
    length(face_storage),
  )
end

function _max_gcl(face_storage::Tuple, domain::CartesianIndices{2})
  I1, I2 = CurvilinearGrids.GridTypes.gcl(face_storage, domain)
  return max(maximum(abs, I1[domain]), maximum(abs, I2[domain]))
end

@testset "Unified mapped/discrete backend harness" begin
  for (backend_name, backend) in _UNIFIED_TEST_BACKENDS
    @testset "backend=$(backend_name)" begin
      number_type = backend_name === :metal ? Float32 : Float64
      ωx = number_type(0.2)
      ωy = number_type(0.3)
      xmap(t, ξ, η, p) = ξ + p.ax * sin(ωx * η) + p.ct * t
      ymap(t, ξ, η, p) = η + p.ay * cos(ωy * ξ) - p.ct * t

      params = (; ax=number_type(0.08), ay=number_type(0.05), ct=number_type(0.01))
      mapped = MappedGrid(
        xmap,
        ymap,
        params,
        (8, 9),
        2;
        backend=backend,
        cache_mode=:eager,
        T=number_type,
        t=zero(number_type),
      )

      mapped_domain = mapped.iterators.cell.domain
      mapped_cm = cell_metrics(mapped)
      mapped_fm = face_metrics(mapped)
      mapped_cm_host = (;
        forward=Array(mapped_cm.forward), inverse=Array(mapped_cm.inverse)
      )
      mapped_fm_host = _host_face_metrics(mapped_fm)

      mapped_J = getproperty.(mapped_cm_host.forward[mapped_domain], :J)
      @test all(isfinite, mapped_J)
      @test minimum(mapped_J) > 0

      mapped_tol = backend_name === :cpu ? 1e-12 : 1e-7
      @test _max_gcl(mapped_fm_host, mapped_domain) < mapped_tol

      mapped_prev = Array(mapped.node_coordinates[1])
      update!(mapped, number_type(1), params)
      mapped_next = Array(mapped.node_coordinates[1])
      @test !all(mapped_prev .== mapped_next)

      x_nodes = Array(mapped.node_coordinates[1])
      y_nodes = Array(mapped.node_coordinates[2])
      if backend isa CPU
        discrete = DiscreteGrid(
          x_nodes,
          y_nodes,
          mapped.nhalo;
          backend=backend,
          T=number_type,
          halo_coords_included=true,
          interpolation=:linear,
          cache_mode=:eager,
        )

        discrete_domain = discrete.iterators.cell.domain
        discrete_cm = cell_metrics(discrete)
        discrete_fm = face_metrics(discrete)
        discrete_cm_host = (;
          forward=Array(discrete_cm.forward), inverse=Array(discrete_cm.inverse)
        )
        discrete_fm_host = _host_face_metrics(discrete_fm)

        discrete_J = getproperty.(discrete_cm_host.forward[discrete_domain], :J)
        @test all(isfinite, discrete_J)
        @test minimum(discrete_J) > 0
        @test _max_gcl(discrete_fm_host, discrete_domain) < 1e-12
      else
        @test_throws ArgumentError DiscreteGrid(
          x_nodes,
          y_nodes,
          mapped.nhalo;
          backend=backend,
          T=number_type,
          halo_coords_included=true,
          interpolation=:linear,
          cache_mode=:eager,
        )
      end
    end
  end
end
