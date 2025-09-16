function get_node_cell_iterators(
  x::AbstractArray{T,3}, y::AbstractArray{T,3}, z::AbstractArray{T,3}, nhalo
) where {T}

  #
  @assert size(x) == size(y) == size(z)
  nnodes = size(x)
  ni, nj, nk = nnodes

  ncells = nnodes .- 1
  ni_cells = ni - 1
  nj_cells = nj - 1
  nk_cells = nk - 1
  lo = nhalo + 1

  limits = (
    node=(ilo=lo, ihi=ni + nhalo, jlo=lo, jhi=nj + nhalo, klo=lo, khi=nk + nhalo),
    cell=(
      ilo=lo,
      ihi=ni_cells + nhalo,
      jlo=lo,
      jhi=nj_cells + nhalo,
      klo=lo,
      khi=nk_cells + nhalo,
    ),
  )

  nodeCI = CartesianIndices(nnodes .+ 2nhalo)
  cellCI = CartesianIndices(ncells .+ 2nhalo)

  node = (
    full=nodeCI,
    domain=nodeCI[
      (begin + nhalo):(end - nhalo),
      (begin + nhalo):(end - nhalo),
      (begin + nhalo):(end - nhalo),
    ],
  )
  cell = (
    full=cellCI,
    domain=cellCI[
      (begin + nhalo):(end - nhalo),
      (begin + nhalo):(end - nhalo),
      (begin + nhalo):(end - nhalo),
    ],
  )
  return limits, (; node, cell)
end

function get_node_cell_iterators(
  x::AbstractArray{T,2}, y::AbstractArray{T,2}, nhalo
) where {T}

  #
  @assert size(x) == size(y)

  nnodes = size(x)
  ni, nj = nnodes

  ncells = nnodes .- 1
  ni_cells, nj_cells = ncells
  lo = nhalo + 1
  limits = (
    node=(ilo=lo, ihi=ni + nhalo, jlo=lo, jhi=nj + nhalo),
    cell=(ilo=lo, ihi=ni_cells + nhalo, jlo=lo, jhi=nj_cells + nhalo),
  )

  nodeCI = CartesianIndices(nnodes .+ 2nhalo)
  cellCI = CartesianIndices(ncells .+ 2nhalo)

  node = (
    full=nodeCI, domain=nodeCI[(begin + nhalo):(end - nhalo), (begin + nhalo):(end - nhalo)]
  )
  cell = (
    full=cellCI, domain=cellCI[(begin + nhalo):(end - nhalo), (begin + nhalo):(end - nhalo)]
  )
  return limits, (; node, cell)
end

function get_node_cell_iterators(x::AbstractVector{T}, nhalo) where {T}
  nnodes = size(x)
  ni, = nnodes

  ncells = nnodes .- 1
  ni_cells, = ncells
  lo = nhalo + 1
  limits = (node=(ilo=lo, ihi=ni + nhalo), cell=(ilo=lo, ihi=ni_cells + nhalo))

  nodeCI = CartesianIndices(nnodes .+ 2nhalo)
  cellCI = CartesianIndices(ncells .+ 2nhalo)

  node = (full=nodeCI, domain=nodeCI[(begin + nhalo):(end - nhalo)])
  cell = (full=cellCI, domain=cellCI[(begin + nhalo):(end - nhalo)])

  return limits, (; node, cell)
end
