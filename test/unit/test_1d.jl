
using CurvilinearGrids

ni = 10
nhalo = 2
x(i) = i^2
m = CurvilinearGrid1D(x, (ni,), nhalo)