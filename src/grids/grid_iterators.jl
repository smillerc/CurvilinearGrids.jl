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
    ),
    cell=(
      full=cellCI,
      domain=cellCI[
        (begin + nhalo):(end - nhalo),
        (begin + nhalo):(end - nhalo),
        (begin + nhalo):(end - nhalo),
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
    ),
    cell=(
      full=cellCI,
      domain=cellCI[(begin + nhalo):(end - nhalo), (begin + nhalo):(end - nhalo)],
    ),
  )
  return _iterators
end

function get_node_cell_iterators(
  nodeCI::CartesianIndices{1}, cellCI::CartesianIndices{1}, nhalo
)
  _iterators = (
    node=(full=nodeCI, domain=nodeCI[(begin + nhalo):(end - nhalo)]),
    cell=(full=cellCI, domain=cellCI[(begin + nhalo):(end - nhalo)]),
  )
  return _iterators
end
