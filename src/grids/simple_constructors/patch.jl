function coons_patch(c0, c1, d0, d1)
  Lc(s, t) = (1 - t) * c0(s) + t * c1(s)
  Ld(s, t) = (1 - s) * d0(t) + s * d1(t)

  function B(s, t)
    return c0(0) * (1 - s) * (1 - t) +
           c0(1) * s * (1 - t) +
           c1(0) * (1 - s) * t +
           c1(1) * s * t
  end

  C(s, t) = Lc(s, t) + Ld(s, t) - B(s, t)

  return C
end

#             x⃗2
#        ───────────►     
#      3              2  
#      ┌──────────────┐ ▲ 
#    ▲ │              │ │ 
#    │ │              │ │ 
#    │ │              │ │ 
# x⃗3 │ │              │ │ x⃗1
#    │ │              │ │ 
#    │ │              │ │ 
#    │ │              │ │ 
#      └──────────────┘   
#      0              1   
#        ──────────►    
#             x⃗0
function PatchGrid(x⃗0, x⃗1, x⃗2, x⃗3)
  @assert length(x⃗0) == length(x⃗2)
  @assert length(x⃗1) == length(x⃗3)

  @assert first(x⃗0) ≈ first(x⃗3)
  @assert first(x⃗1) ≈ last(x⃗1)
  @assert last(x⃗1) ≈ last(x⃗2)
  @assert first(x⃗3) ≈ last(x⃗3)
end