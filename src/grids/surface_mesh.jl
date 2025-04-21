using CartesianDomains

function extract_surface_mesh(mesh::CurvilinearGrid2D, loc::Symbol)
  i, j, = (1, 2)
  full = mesh.iterators.node.domain

  if loc === :ilo
    domain = extract_from_lower(full, i, 1)
  elseif loc === :jlo
    domain = extract_from_lower(full, j, 1)
  elseif loc === :ihi
    domain = extract_from_upper(full, i, 1)
  elseif loc === :jhi
    domain = extract_from_upper(full, j, 1)
  end

  x = Array(mesh.node_coordinates.x[domain])
  y = Array(mesh.node_coordinates.y[domain])

  return (x, y)
end

function extract_surface_mesh(mesh::CurvilinearGrid3D, loc::Symbol)
  i, j, k = (1, 2, 3)
  full = mesh.iterators.node.domain

  if loc === :ilo
    domain = extract_from_lower(full, i, 1)
  elseif loc === :jlo
    domain = extract_from_lower(full, j, 1)
  elseif loc === :klo
    domain = extract_from_lower(full, k, 1)
  elseif loc === :ihi
    domain = extract_from_upper(full, i, 1)
  elseif loc === :jhi
    domain = extract_from_upper(full, j, 1)
  elseif loc === :khi
    domain = extract_from_upper(full, k, 1)
  end

  x = Array(mesh.node_coordinates.x[domain])
  y = Array(mesh.node_coordinates.y[domain])
  z = Array(mesh.node_coordinates.z[domain])

  return (x, y, z)
end