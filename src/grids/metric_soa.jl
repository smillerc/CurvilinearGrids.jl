
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

function get_metric_soa_rectlinear2d(celldims::NTuple{2,Int}, backend, T)
  cell_center_metrics = (
    J=KernelAbstractions.zeros(backend, T, celldims),
    ξ=StructArray((
      x₁=RectlinearArrays.zeros(T, backend, (2,), celldims),
      x₂=RectlinearArrays.zeros(T, backend, (1,2), celldims),
      t=RectlinearArrays.zeros(T, backend, (1,2), celldims),
    )),
    η=StructArray((
      x₁=RectlinearArrays.zeros(T, backend, (1,2), celldims),
      x₂=RectlinearArrays.zeros(T, backend, (1,), celldims),
      t=RectlinearArrays.zeros(T, backend, (1,2), celldims),
    )),
    ξ̂=StructArray((
      x₁=RectlinearArrays.zeros(T, backend, (1,), celldims),
      x₂=RectlinearArrays.zeros(T, backend, (1,2), celldims),
      t=RectlinearArrays.zeros(T, backend, (1,2), celldims),
    )),
    η̂=StructArray((
      x₁=RectlinearArrays.zeros(T, backend, (1,2), celldims),
      x₂=RectlinearArrays.zeros(T, backend, (2,), celldims),
      t=RectlinearArrays.zeros(T, backend, (1,2), celldims),
    )),
    x₁=StructArray((
      ξ=RectlinearArrays.zeros(T, backend, (2,), celldims),
      η=RectlinearArrays.zeros(T, backend, (1,2), celldims),
    )),
    x₂=StructArray((
      ξ=RectlinearArrays.zeros(T, backend, (1,2), celldims),
      η=RectlinearArrays.zeros(T, backend, (1,), celldims),
    )),
  )

  edge_metrics = (
    i₊½=(
      ξ=StructArray((
        x₁=RectlinearArrays.zeros(T, backend, (2,), celldims),
        x₂=RectlinearArrays.zeros(T, backend, (1,2), celldims),
      )),
      η=StructArray((
        x₁=RectlinearArrays.zeros(T, backend, (1,2), celldims),
        x₂=RectlinearArrays.zeros(T, backend, (1,), celldims),
      )),
      ξ̂=StructArray((
        x₁=RectlinearArrays.zeros(T, backend, (1,), celldims),
        x₂=RectlinearArrays.zeros(T, backend, (1,2), celldims),
        t=RectlinearArrays.zeros(T, backend, (1,2), celldims),
      )),
      η̂=StructArray((
        x₁=RectlinearArrays.zeros(T, backend, (1,2), celldims),
        x₂=RectlinearArrays.zeros(T, backend, (2,), celldims),
        t=RectlinearArrays.zeros(T, backend, (1,2), celldims),
      )),
    ),
    j₊½=(
      ξ=StructArray((
        x₁=RectlinearArrays.zeros(T, backend, (2,), celldims),
        x₂=RectlinearArrays.zeros(T, backend, (1,2), celldims),
      )),
      η=StructArray((
        x₁=RectlinearArrays.zeros(T, backend, (1,2), celldims),
        x₂=RectlinearArrays.zeros(T, backend, (1,), celldims),
      )),
      ξ̂=StructArray((
        x₁=RectlinearArrays.zeros(T, backend, (1,), celldims),
        x₂=RectlinearArrays.zeros(T, backend, (1,2), celldims),
        t=RectlinearArrays.zeros(T, backend, (1,2), celldims),
      )),
      η̂=StructArray((
        x₁=RectlinearArrays.zeros(T, backend, (1,2), celldims),
        x₂=RectlinearArrays.zeros(T, backend, (2,), celldims),
        t=RectlinearArrays.zeros(T, backend, (1,2), celldims),
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
      ξ=StructArray((x₁=KernelAbstractions.zeros(backend, T, celldims),)),
      ξ̂=StructArray((
        x₁=KernelAbstractions.zeros(backend, T, celldims),
        t=KernelAbstractions.zeros(backend, T, celldims),
      )),
    ),
  )

  return cell_center_metrics, edge_metrics
end
