"""
    RemapCache{N,T,...}

Reusable sparse overlap map for conservative scalar remapping between two unified
grids with dimension `N`.

# Fields
  - `source_shape`: Source interior cell shape.
  - `destination_shape`: Destination interior cell shape.
  - `overlap_volume`: Sparse matrix with entries
    `overlap_volume[dst_linear, src_linear] = |Omega_src ∩ Omega_dst|`.
  - `source_cell_volumes`: Source cell volumes in source linear ordering.
  - `destination_cell_volumes`: Destination cell volumes in destination linear ordering.
  - `source_overlap_volumes`: Per-source overlap volumes (`sum(overlap_volume; dims=1)`).
  - `destination_overlap_volumes`: Per-destination covered volumes (`sum(overlap_volume; dims=2)`).
  - `source_state_token`: Hash token of source grid state when cache was built.
  - `destination_state_token`: Hash token of destination grid state when cache was built.
  - `quadrature_order`: Tensor-product quadrature order used for overlap assembly.
  - `sample_count`: Number of samples per source cell used during assembly.
"""
struct RemapCache{
  N,
  T,
  Ssrc<:NTuple{N,Int},
  Sdst<:NTuple{N,Int},
  M<:SparseMatrixCSC{T,Int},
  Vsrc<:AbstractVector{T},
  Vdst<:AbstractVector{T},
  Vosrc<:AbstractVector{T},
  Vodst<:AbstractVector{T},
}
  source_shape::Ssrc
  destination_shape::Sdst
  overlap_volume::M
  source_cell_volumes::Vsrc
  destination_cell_volumes::Vdst
  source_overlap_volumes::Vosrc
  destination_overlap_volumes::Vodst
  source_state_token::UInt
  destination_state_token::UInt
  quadrature_order::Int
  sample_count::Int
end
