function get_node_cell_iterators(
  nodeCI::CartesianIndices{3}, cellCI::CartesianIndices{3}, nhalo
)
  _iterators = (
    node=(
      full=nodeCI,
      domain=nodeCI[
        (begin + nhalo):(end - nhalo),
        (begin + nhalo):(end - nhalo),
        (begin + nhalo):(end - nhalo),
      ],
      ilo_halo=nodeCI[
        begin:(begin + nhalo), (begin + nhalo):(end - nhalo), (begin + nhalo):(end - nhalo)
      ],
      jlo_halo=nodeCI[
        (begin + nhalo):(end - nhalo), begin:(begin + nhalo), (begin + nhalo):(end - nhalo)
      ],
      klo_halo=nodeCI[
        (begin + nhalo):(end - nhalo), (begin + nhalo):(end - nhalo), begin:(begin + nhalo)
      ],
      ihi_halo=nodeCI[
        (end - nhalo):end, (begin + nhalo):(end - nhalo), (begin + nhalo):(end - nhalo)
      ],
      jhi_halo=nodeCI[
        (begin + nhalo):(end - nhalo), (end - nhalo):end, (begin + nhalo):(end - nhalo)
      ],
      khi_halo=nodeCI[
        (begin + nhalo):(end - nhalo), (begin + nhalo):(end - nhalo), (end - nhalo):end
      ],
    ),
    cell=(
      full=cellCI,
      domain=cellCI[
        (begin + nhalo):(end - nhalo),
        (begin + nhalo):(end - nhalo),
        (begin + nhalo):(end - nhalo),
      ],
      ilo_halo=cellCI[
        begin:(begin + nhalo - 1),
        (begin + nhalo):(end - nhalo),
        (begin + nhalo):(end - nhalo),
      ],
      jlo_halo=cellCI[
        (begin + nhalo):(end - nhalo),
        begin:(begin + nhalo - 1),
        (begin + nhalo):(end - nhalo),
      ],
      klo_halo=cellCI[
        (begin + nhalo):(end - nhalo),
        (begin + nhalo):(end - nhalo),
        begin:(begin + nhalo - 1),
      ],
      ihi_halo=cellCI[
        (end - nhalo + 1):end, (begin + nhalo):(end - nhalo), (begin + nhalo):(end - nhalo)
      ],
      jhi_halo=cellCI[
        (begin + nhalo):(end - nhalo), (end - nhalo + 1):end, (begin + nhalo):(end - nhalo)
      ],
      khi_halo=cellCI[
        (begin + nhalo):(end - nhalo), (begin + nhalo):(end - nhalo), (end - nhalo + 1):end
      ],
    ),
  )

  return _iterators
end

function get_node_cell_iterators(
  nodeCI::CartesianIndices{2}, cellCI::CartesianIndices{2}, nhalo
)
  _iterators = (
    node=(
      full=nodeCI,
      domain=nodeCI[(begin + nhalo):(end - nhalo), (begin + nhalo):(end - nhalo)],
      ilo_halo=nodeCI[begin:(begin + nhalo), (begin + nhalo):(end - nhalo)],
      ihi_halo=nodeCI[(end - nhalo):end, (begin + nhalo):(end - nhalo)],
      jlo_halo=nodeCI[(begin + nhalo):(end - nhalo), begin:(begin + nhalo)],
      jhi_halo=nodeCI[(begin + nhalo):(end - nhalo), (end - nhalo):end],
    ),
    cell=(
      full=cellCI,
      domain=cellCI[(begin + nhalo):(end - nhalo), (begin + nhalo):(end - nhalo)],
      ilo_halo=cellCI[begin:(begin + nhalo - 1), (begin + nhalo):(end - nhalo)],
      jlo_halo=cellCI[(begin + nhalo):(end - nhalo), begin:(begin + nhalo - 1)],
      ihi_halo=cellCI[(end - nhalo + 1):end, (begin + nhalo):(end - nhalo)],
      jhi_halo=cellCI[(begin + nhalo):(end - nhalo), (end - nhalo + 1):end],
    ),
  )
  return _iterators
end

function get_node_cell_iterators(
  nodeCI::CartesianIndices{1}, cellCI::CartesianIndices{1}, nhalo
)
  _iterators = (
    node=(
      full=nodeCI,
      domain=nodeCI[(begin + nhalo):(end - nhalo)],
      ilo_halo=nodeCI[begin:(begin + nhalo)],
      ihi_halo=nodeCI[(end - nhalo):end],
    ),
    cell=(
      full=cellCI,
      domain=cellCI[(begin + nhalo):(end - nhalo)],
      ilo_halo=cellCI[begin:(begin + nhalo - 1)],
      ihi_halo=cellCI[(end - nhalo + 1):end],
    ),
  )
  return _iterators
end
