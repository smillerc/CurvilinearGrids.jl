
using Plots, CurvilinearGrids

function polar_mesh(ni, nj, nk)
  r0, r1 = (0, 2)
  # θ0, θ1 = (0, 2pi)
  θ0, θ1 = (0, 2)# times π
  z0, z1 = (0, 3)

  r(i) = r0 + (r1 - r0) * ((i - 1) / (ni - 1))
  θ(j) = θ0 + (θ1 - θ0) * ((j - 1) / (nj - 1))
  # θ(j) = 0
  _z(k) = z0 + (z1 - z0) * ((k - 1) / (nk - 1))

  x(i, j, k) = r(i) * cospi(θ(j))
  y(i, j, k) = r(i) * sinpi(θ(j))
  # x(i, j, k) = r(i) * cos(θ(j))
  # y(i, j, k) = r(i) * sin(θ(j))
  z(i, j, k) = _z(k)

  return (x, y, z)
end

function polar_mesh2d(ni, nj)
  r0, r1 = (0, 2)
  # θ0, θ1 = (0, 1e-10)
  z0, z1 = (0, 3)

  R(i) = r0 + (r1 - r0) * ((i - 1) / (ni - 1))
  Z(j) = z0 + (z1 - z0) * ((k - 1) / (nk - 1))

  x(i, j) = R(i) * cos(θ(j))
  # y(i, j, k) = r(i) * sin(θ(j))
  y(i, j) = _z(k)

  return (x, y)
end

ni, nj, nk = (5, 2, 11)
nhalo = 0
x, y, z = polar_mesh(ni, nj, nk)
mesh = CurvilinearGrid3D(x, y, z, (ni, nj, nk), nhalo)

xyz = CurvilinearGrids.coords(mesh)
xyz_c = CurvilinearGrids.centroids(mesh)

xn = @view xyz[1, :, :, :]
yn = @view xyz[2, :, :, :]
zn = @view xyz[3, :, :, :]

xc = @view xyz_c[1, :, :, :]
yc = @view xyz_c[2, :, :, :]
zc = @view xyz_c[3, :, :, :]

jacobian_matrix(mesh, (2, 1, 20))

p = plot(; xlabel="R", ylabel="Z", aspect_ratio=:equal)

# J = zeros(cellsize_withhalo(mesh))
# for idx in CartesianIndices(J)
#   i, j, k = idx.I
#   J[idx] = jacobian(mesh, (i + 0.5, j + 0.5, k + 0.5))
# end

# xc1d = xc[:, 1, 1]
# zc1d = zc[1, 1, :]
# size(J[:, 1, :])
# heatmap!(xc1d, zc1d, J[:, 1, :]'; colormap=:Spectral)
# heatmap(J[:, 1, :]'; colormap=:Spectral)

for k in axes(xn, 3)
  plot!(xn[:, 1, k], zn[:, 1, k]; color=:black, marker=:dot, label=nothing)
end

# for k in axes(xc, 3)
#   scatter!(xc[:, 1, k], zc[:, 1, k]; color=:black, marker=:cross, label=nothing)
# end

for i in axes(xn, 1)
  plot!(xn[i, 1, :], zn[i, 1, :]; color=:black, marker=:dot, label=nothing)
end
p
# begin
#   ijk = (5.5, 0, 2.5)
#   # given (i,j,k) what is the xyz position?
#   x⃗ = coord(mesh, ijk) # → (1.25, 0, 0.75)

#   J_ijk = jacobian_matrix(mesh, ijk)
#   J⁻¹_ijk = inv(J_ijk)

#   # given (x,y,z), what is (i,j,k)?
#   ijk_computed = J⁻¹_ijk * x⃗ #.+ (mesh.nhalo + 1)
# end
# J_ijk
# # J
# p

begin
  ξx = 2.0
  ξy = 0.0
  ξz = 0.0
  ηx = 0.0
  ηy = 1.4335402848056648e15
  ηz = 0.0
  ζx = 0.0
  ζy = 0.0
  ζz = 3.333333333333333
  ξt = 0.0
  ηt = 0.0
  ζt = 0.0
  J = 1.0463605494025897e-16

  ((J * ξx)^2 + (J * ξz)^2) / J
end

ξₓ¹
ξₓ²
ξₓ³