
# using MPI
using CurvilinearGrids

# MPI.Init()
# const comm = MPI.COMM_WORLD
# const rank = MPI.Comm_rank(comm)
# const nprocs = MPI.Comm_size(comm)

# begin
#   pg1 = PartitionedRectlinearGrid(0.0, 1.1, 11, 2, 2, 1)
#   @show coords(pg1)
#   pg2 = PartitionedRectlinearGrid(0.0, 1.1, 11, 2, 2, 2)
#   @show coords(pg2)
#   nothing
# end

begin
  for rank in 1:6
    #     @show rank
    println("\nrank: $rank")
    # rank = 2
    nhalo = 2
    split_frac = (3, 2)
    pg = PartitionedRectlinearGrid(
      (0.0, 0.0), (2.1, 1.0), (21, 10), nhalo, split_frac, rank
    )

    @show pg.onbc
    # save_vtk(pg, "mesh_rank_$(rank)")
    #     x, y = coords(pg)
    #     display(pg.node_coordinates.x)
    #     display(pg.node_coordinates.y)
  end
  #   nothing
end

# MPI.Finalize()
