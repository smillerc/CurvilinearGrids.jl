function get_iterators(
  nodedims::NTuple{N,Int}, halo_coords_included::Bool, nhalo::Int
) where {N}

  #
  if halo_coords_included
    cellCI = CartesianIndices(nodedims .- 1 .+ 2nhalo)
    nodeCI = CartesianIndices(nodedims .+ 2nhalo)
  else
    cellCI = CartesianIndices(nodedims .- 1)
    nodeCI = CartesianIndices(nodedims)
  end

  node = (full=nodeCI, domain=expand(nodeCI, -nhalo))
  cell = (full=cellCI, domain=expand(cellCI, -nhalo))
  limits = domain_limits(cell.domain)
  return limits, (; node, cell)
end

function domain_limits(cell_domain::CartesianIndices{1})
  ilo = first(cell_domain.indices[1])
  ihi = last(cell_domain.indices[1])

  limits = (
    node=(; ilo=ilo, ihi=ihi + 1), #
    cell=(; ilo=ilo, ihi=ihi),     #
  )

  return limits
end

function domain_limits(cell_domain::CartesianIndices{2})
  ilo = first(cell_domain.indices[1])
  ihi = last(cell_domain.indices[1])
  jlo = first(cell_domain.indices[2])
  jhi = last(cell_domain.indices[2])

  limits = (
    node=(; ilo=ilo, ihi=ihi + 1, jlo=jlo, jhi=jhi + 1), #
    cell=(; ilo=ilo, ihi=ihi, jlo=jlo, jhi=jhi),         #
  )

  return limits
end

function domain_limits(cell_domain::CartesianIndices{3})
  ilo = first(cell_domain.indices[1])
  ihi = last(cell_domain.indices[1])
  jlo = first(cell_domain.indices[2])
  jhi = last(cell_domain.indices[2])
  klo = first(cell_domain.indices[3])
  khi = last(cell_domain.indices[3])

  limits = (
    node=(; ilo=ilo, ihi=ihi + 1, jlo=jlo, jhi=jhi + 1, klo=klo, khi=khi + 1), #
    cell=(; ilo=ilo, ihi=ihi, jlo=jlo, jhi=jhi, klo=klo, khi=khi),             #
  )

  return limits
end