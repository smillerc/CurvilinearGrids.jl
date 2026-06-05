using Documenter
using CurvilinearGrids

DocMeta.setdocmeta!(
  CurvilinearGrids, :DocTestSetup, :(using CurvilinearGrids); recursive=true
)

makedocs(;
  modules=[CurvilinearGrids],
  sitename="CurvilinearGrids.jl",
  authors="CurvilinearGrids contributors",
  format=Documenter.HTML(; prettyurls=get(ENV, "CI", "false") == "true"),
  checkdocs=:none,
  doctest=false,
  pages=[
    "Home" => "index.md",
    "Getting Started" => "getting-started.md",
    "Metric Schemes and GCL" => "metric-schemes.md",
    "Unified Grid API" => "unified-grid-api.md",
    "Solver Interface" => "solver-interface.md",
    "Mesh Quality and Testing" => "mesh-quality.md",
    "API Reference" => "api-reference.md",
  ],
)
