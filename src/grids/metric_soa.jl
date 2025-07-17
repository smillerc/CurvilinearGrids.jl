
function get_metric_soa(celldims::NTuple{3,Int}, backend, T)
  cell_center_metrics = (
    J=KernelAbstractions.zeros(backend, T, celldims),
    ξ=StructArray((
      x₁=KernelAbstractions.zeros(backend, T, celldims),
      x₂=KernelAbstractions.zeros(backend, T, celldims),
      x₃=KernelAbstractions.zeros(backend, T, celldims),
      t=KernelAbstractions.zeros(backend, T, celldims),
    )),
    η=StructArray((
      x₁=KernelAbstractions.zeros(backend, T, celldims),
      x₂=KernelAbstractions.zeros(backend, T, celldims),
      x₃=KernelAbstractions.zeros(backend, T, celldims),
      t=KernelAbstractions.zeros(backend, T, celldims),
    )),
    ζ=StructArray((
      x₁=KernelAbstractions.zeros(backend, T, celldims),
      x₂=KernelAbstractions.zeros(backend, T, celldims),
      x₃=KernelAbstractions.zeros(backend, T, celldims),
      t=KernelAbstractions.zeros(backend, T, celldims),
    )),
    ξ̂=StructArray((
      x₁=KernelAbstractions.zeros(backend, T, celldims),
      x₂=KernelAbstractions.zeros(backend, T, celldims),
      x₃=KernelAbstractions.zeros(backend, T, celldims),
      t=KernelAbstractions.zeros(backend, T, celldims),
    )),
    η̂=StructArray((
      x₁=KernelAbstractions.zeros(backend, T, celldims),
      x₂=KernelAbstractions.zeros(backend, T, celldims),
      x₃=KernelAbstractions.zeros(backend, T, celldims),
      t=KernelAbstractions.zeros(backend, T, celldims),
    )),
    ζ̂=StructArray((
      x₁=KernelAbstractions.zeros(backend, T, celldims),
      x₂=KernelAbstractions.zeros(backend, T, celldims),
      x₃=KernelAbstractions.zeros(backend, T, celldims),
      t=KernelAbstractions.zeros(backend, T, celldims),
    )),
    x₁=StructArray((
      ξ=KernelAbstractions.zeros(backend, T, celldims),
      η=KernelAbstractions.zeros(backend, T, celldims),
      ζ=KernelAbstractions.zeros(backend, T, celldims),
    )),
    x₂=StructArray((
      ξ=KernelAbstractions.zeros(backend, T, celldims),
      η=KernelAbstractions.zeros(backend, T, celldims),
      ζ=KernelAbstractions.zeros(backend, T, celldims),
    )),
    x₃=StructArray((
      ξ=KernelAbstractions.zeros(backend, T, celldims),
      η=KernelAbstractions.zeros(backend, T, celldims),
      ζ=KernelAbstractions.zeros(backend, T, celldims),
    )),
  )

  edge_metrics = (
    i₊½=(
      J=KernelAbstractions.zeros(backend, T, celldims),
      ξ=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
        x₃=KernelAbstractions.zeros(backend, T, celldims),
      )),
      η=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
        x₃=KernelAbstractions.zeros(backend, T, celldims),
      )),
      ζ=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
        x₃=KernelAbstractions.zeros(backend, T, celldims),
      )),
      ξ̂=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
        x₃=KernelAbstractions.zeros(backend, T, celldims),
        t=KernelAbstractions.zeros(backend, T, celldims),
      )),
      η̂=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
        x₃=KernelAbstractions.zeros(backend, T, celldims),
        t=KernelAbstractions.zeros(backend, T, celldims),
      )),
      ζ̂=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
        x₃=KernelAbstractions.zeros(backend, T, celldims),
        t=KernelAbstractions.zeros(backend, T, celldims),
      )),
    ),
    j₊½=(
      J=KernelAbstractions.zeros(backend, T, celldims),
      ξ=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
        x₃=KernelAbstractions.zeros(backend, T, celldims),
      )),
      η=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
        x₃=KernelAbstractions.zeros(backend, T, celldims),
      )),
      ζ=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
        x₃=KernelAbstractions.zeros(backend, T, celldims),
      )),
      ξ̂=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
        x₃=KernelAbstractions.zeros(backend, T, celldims),
        t=KernelAbstractions.zeros(backend, T, celldims),
      )),
      η̂=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
        x₃=KernelAbstractions.zeros(backend, T, celldims),
        t=KernelAbstractions.zeros(backend, T, celldims),
      )),
      ζ̂=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
        x₃=KernelAbstractions.zeros(backend, T, celldims),
        t=KernelAbstractions.zeros(backend, T, celldims),
      )),
    ),
    k₊½=(
      J=KernelAbstractions.zeros(backend, T, celldims),
      ξ=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
        x₃=KernelAbstractions.zeros(backend, T, celldims),
      )),
      η=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
        x₃=KernelAbstractions.zeros(backend, T, celldims),
      )),
      ζ=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
        x₃=KernelAbstractions.zeros(backend, T, celldims),
      )),
      ξ̂=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
        x₃=KernelAbstractions.zeros(backend, T, celldims),
        t=KernelAbstractions.zeros(backend, T, celldims),
      )),
      η̂=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
        x₃=KernelAbstractions.zeros(backend, T, celldims),
        t=KernelAbstractions.zeros(backend, T, celldims),
      )),
      ζ̂=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
        x₃=KernelAbstractions.zeros(backend, T, celldims),
        t=KernelAbstractions.zeros(backend, T, celldims),
      )),
    ),
  )

  return cell_center_metrics, edge_metrics
end

function get_metric_soa_rectilinear3d(celldims::NTuple{3,Int}, backend, T)
  cell_center_metrics = (
    J=KernelAbstractions.zeros(backend, T, celldims),
    ξ=StructArray((
      x₁=RectilinearArrays.zeros(T, backend, (2,3), celldims),
      x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
    )),
    η=StructArray((
      x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      x₂=RectilinearArrays.zeros(T, backend, (1,3), celldims),
      x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
    )),
    ζ=StructArray((
      x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      x₃=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
    )),
    ξ̂=StructArray((
      x₁=RectilinearArrays.zeros(T, backend, (1,), celldims),
      x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
    )),
    η̂=StructArray((
      x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      x₂=RectilinearArrays.zeros(T, backend, (2,), celldims),
      x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
    )),
    ζ̂=StructArray((
      x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      x₃=RectilinearArrays.zeros(T, backend, (3,), celldims),
      t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
    )),
    x₁=StructArray((
      ξ=RectilinearArrays.zeros(T, backend, (2,3), celldims),
      η=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      ζ=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
    )),
    x₂=StructArray((
      ξ=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      η=RectilinearArrays.zeros(T, backend, (1,3), celldims),
      ζ=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
    )),
    x₃=StructArray((
      ξ=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      η=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      ζ=RectilinearArrays.zeros(T, backend, (1,2), celldims),
    )),
  )

  edge_metrics = (
    i₊½=(
      J=KernelAbstractions.zeros(backend, T, celldims),
      ξ=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      η=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      ζ=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      )),
      ξ̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      η̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (2,), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      ζ̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (3,), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
    ),
    j₊½=(
      J=KernelAbstractions.zeros(backend, T, celldims),
      ξ=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      η=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      ζ=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      )),
      ξ̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      η̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (2,), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      ζ̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (3,), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
    ),
    k₊½=(
      J=KernelAbstractions.zeros(backend, T, celldims),
      ξ=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      η=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      ζ=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      )),
      ξ̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      η̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (2,), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      ζ̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (3,), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
    ),
  )

  return cell_center_metrics, edge_metrics
end

function get_metric_soa_uniform3d(celldims::NTuple{3,Int}, backend, T)
  cell_center_metrics = (
    J=KernelAbstractions.zeros(backend, T, celldims),
    ξ=StructArray((
      x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
    )),
    η=StructArray((
      x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
    )),
    ζ=StructArray((
      x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
    )),
    ξ̂=StructArray((
      x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
    )),
    η̂=StructArray((
      x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
    )),
    ζ̂=StructArray((
      x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
    )),
    x₁=StructArray((
      ξ=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      η=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      ζ=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
    )),
    x₂=StructArray((
      ξ=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      η=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      ζ=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
    )),
    x₃=StructArray((
      ξ=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      η=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      ζ=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
    )),
  )

  edge_metrics = (
    i₊½=(
      J=KernelAbstractions.zeros(backend, T, celldims),
      ξ=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      η=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      ζ=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      ξ̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      η̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      ζ̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
    ),
    j₊½=(
      J=KernelAbstractions.zeros(backend, T, celldims),
      ξ=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      η=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      ζ=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      ξ̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      η̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      ζ̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
    ),
    k₊½=(
      J=KernelAbstractions.zeros(backend, T, celldims),
      ξ=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      η=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      ζ=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      ξ̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      η̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
      ζ̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        x₃=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2,3), celldims),
      )),
    ),
  )

  return cell_center_metrics, edge_metrics
end

function get_metric_soa(celldims::NTuple{2,Int}, backend, T)
  cell_center_metrics = (
    J=KernelAbstractions.zeros(backend, T, celldims),
    ξ=StructArray((
      x₁=KernelAbstractions.zeros(backend, T, celldims),
      x₂=KernelAbstractions.zeros(backend, T, celldims),
      t=KernelAbstractions.zeros(backend, T, celldims),
    )),
    η=StructArray((
      x₁=KernelAbstractions.zeros(backend, T, celldims),
      x₂=KernelAbstractions.zeros(backend, T, celldims),
      t=KernelAbstractions.zeros(backend, T, celldims),
    )),
    ξ̂=StructArray((
      x₁=KernelAbstractions.zeros(backend, T, celldims),
      x₂=KernelAbstractions.zeros(backend, T, celldims),
      t=KernelAbstractions.zeros(backend, T, celldims),
    )),
    η̂=StructArray((
      x₁=KernelAbstractions.zeros(backend, T, celldims),
      x₂=KernelAbstractions.zeros(backend, T, celldims),
      t=KernelAbstractions.zeros(backend, T, celldims),
    )),
    x₁=StructArray((
      ξ=KernelAbstractions.zeros(backend, T, celldims),
      η=KernelAbstractions.zeros(backend, T, celldims),
    )),
    x₂=StructArray((
      ξ=KernelAbstractions.zeros(backend, T, celldims),
      η=KernelAbstractions.zeros(backend, T, celldims),
    )),
  )

  edge_metrics = (
    i₊½=(
      J=KernelAbstractions.zeros(backend, T, celldims),
      ξ=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
      )),
      η=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
      )),
      ξ̂=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
        t=KernelAbstractions.zeros(backend, T, celldims),
      )),
      η̂=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
        t=KernelAbstractions.zeros(backend, T, celldims),
      )),
    ),
    j₊½=(
      J=KernelAbstractions.zeros(backend, T, celldims),
      ξ=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
      )),
      η=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
      )),
      ξ̂=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
        t=KernelAbstractions.zeros(backend, T, celldims),
      )),
      η̂=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        x₂=KernelAbstractions.zeros(backend, T, celldims),
        t=KernelAbstractions.zeros(backend, T, celldims),
      )),
    ),
  )

  return cell_center_metrics, edge_metrics
end

function get_metric_soa_rectilinear2d(celldims::NTuple{2,Int}, backend, T)
  cell_center_metrics = (
    J=KernelAbstractions.zeros(backend, T, celldims),
    ξ=StructArray((
      x₁=RectilinearArrays.zeros(T, backend, (2,), celldims),
      x₂=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      t=RectilinearArrays.zeros(T, backend, (1,2), celldims),
    )),
    η=StructArray((
      x₁=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      x₂=RectilinearArrays.zeros(T, backend, (1,), celldims),
      t=RectilinearArrays.zeros(T, backend, (1,2), celldims),
    )),
    ξ̂=StructArray((
      x₁=RectilinearArrays.zeros(T, backend, (1,), celldims),
      x₂=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      t=RectilinearArrays.zeros(T, backend, (1,2), celldims),
    )),
    η̂=StructArray((
      x₁=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      x₂=RectilinearArrays.zeros(T, backend, (2,), celldims),
      t=RectilinearArrays.zeros(T, backend, (1,2), celldims),
    )),
    x₁=StructArray((
      ξ=RectilinearArrays.zeros(T, backend, (2,), celldims),
      η=RectilinearArrays.zeros(T, backend, (1,2), celldims),
    )),
    x₂=StructArray((
      ξ=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      η=RectilinearArrays.zeros(T, backend, (1,), celldims),
    )),
  )

  edge_metrics = (
    i₊½=(
      J=KernelAbstractions.zeros(backend, T, celldims),
      ξ=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (2,), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      )),
      η=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,), celldims),
      )),
      ξ̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      )),
      η̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (2,), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      )),
    ),
    j₊½=(
      J=KernelAbstractions.zeros(backend, T, celldims),
      ξ=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (2,), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      )),
      η=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,), celldims),
      )),
      ξ̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      )),
      η̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (2,), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      )),
    ),
  )

  return cell_center_metrics, edge_metrics
end

function get_metric_soa_uniform2d(celldims::NTuple{2,Int}, backend, T)
  cell_center_metrics = (
    # As long as halo cells don't matter
    # J=KernelAbstractions.zeros(backend, T, celldims),
    J=RectilinearArrays.zeros(T, backend, (1,2), celldims),
    ξ=StructArray((
      x₁=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      x₂=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      t=RectilinearArrays.zeros(T, backend, (1,2), celldims),
    )),
    η=StructArray((
      x₁=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      x₂=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      t=RectilinearArrays.zeros(T, backend, (1,2), celldims),
    )),
    ξ̂=StructArray((
      x₁=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      x₂=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      t=RectilinearArrays.zeros(T, backend, (1,2), celldims),
    )),
    η̂=StructArray((
      x₁=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      x₂=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      t=RectilinearArrays.zeros(T, backend, (1,2), celldims),
    )),
    x₁=StructArray((
      ξ=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      η=RectilinearArrays.zeros(T, backend, (1,2), celldims),
    )),
    x₂=StructArray((
      ξ=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      η=RectilinearArrays.zeros(T, backend, (1,2), celldims),
    )),
  )

  edge_metrics = (
    i₊½=(
      J=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      ξ=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      )),
      η=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      )),
      ξ̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      )),
      η̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      )),
    ),
    j₊½=(
      J=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      ξ=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      )),
      η=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      )),
      ξ̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      )),
      η̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,2), celldims),
        x₂=RectilinearArrays.zeros(T, backend, (1,2), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,2), celldims),
      )),
    ),
  )

  return cell_center_metrics, edge_metrics
end

function get_metric_soa(celldims::NTuple{1,Int}, backend, T)
  cell_center_metrics = (
    J=KernelAbstractions.zeros(backend, T, celldims),
    ξ=StructArray((
      x₁=KernelAbstractions.zeros(backend, T, celldims),
      t=KernelAbstractions.zeros(backend, T, celldims),
    )),
    ξ̂=StructArray((
      x₁=KernelAbstractions.zeros(backend, T, celldims),
      t=KernelAbstractions.zeros(backend, T, celldims),
    )),
    x₁=StructArray((ξ=KernelAbstractions.zeros(backend, T, celldims),)),
  )

  edge_metrics = (
    i₊½=(
      J=KernelAbstractions.zeros(backend, T, celldims),
      ξ=StructArray((x₁=KernelAbstractions.zeros(backend, T, celldims),)),
      ξ̂=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        t=KernelAbstractions.zeros(backend, T, celldims),
      )),
    ),
  )

  return cell_center_metrics, edge_metrics
end

function get_metric_soa_uniform1d(celldims::NTuple{1,Int}, backend::Backend, ::Type{T}) where {T}
  cell_center_metrics = (
    J=RectilinearArrays.zeros(T, backend, (1,), celldims),
    ξ=StructArray((
      x₁=RectilinearArrays.zeros(T, backend, (1,), celldims),
      t=RectilinearArrays.zeros(T, backend, (1,), celldims),
    )),
    ξ̂=StructArray((
      x₁=RectilinearArrays.zeros(T, backend, (1,), celldims),
      t=RectilinearArrays.zeros(T, backend, (1,), celldims),
    )),
    x₁=StructArray((ξ=RectilinearArrays.zeros(T, backend, (1,), celldims),)),
  )

  edge_metrics = (
    i₊½=(
      J=RectilinearArrays.zeros(T, backend, (1,), celldims),
      ξ=StructArray((x₁=RectilinearArrays.zeros(T, backend, (1,), celldims),)),
      ξ̂=StructArray((
        x₁=RectilinearArrays.zeros(T, backend, (1,), celldims),
        t=RectilinearArrays.zeros(T, backend, (1,), celldims),
      )),
    ),
  )

  return cell_center_metrics, edge_metrics
end
