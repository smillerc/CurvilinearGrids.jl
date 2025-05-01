using BenchmarkTools
import CurvilinearGrids

const SUITE = BenchmarkGroup()

SUITE["rectlinear"] = BenchmarkGroup()

# 2D grid benchmark
ni, nj = (40, 80)
x0, x1 = (0, 2)
y0, y1 = (1, 3)

SUITE["rectlinear"]["2D"] = @benchmarkable(CurvilinearGrids.rectlinear_grid(($x0, $y0), ($x1, $y1), ($ni, $nj), :MEG6))

# 3D grid benchmark
x2, x3 = (0.0, 2.0)
y2, y3 = (1, 3)
z2, z3 = (-1, 2)
ni1, nj1, nk1 = (40, 80, 120)

SUITE["rectlinear"]["3D"] = @benchmarkable(CurvilinearGrids.rectlinear_grid(($x2, $y2, $z2), ($x3, $y3, $z3), ($ni1, $nj1, $nk1), :meg6))
