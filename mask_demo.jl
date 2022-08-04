## Libraries:
using GaussQuadrature
using LinearAlgebra
using Plots
using Revise
using StaticArrays

using PathIntersections

## Curve parameters ------------------------------------------------------------------
pacman = PresetGeometries.Pacman(R=0.5, x0=0.13, y0=0.2)
pizza_slice = PresetGeometries.Pacman(R=0.5, x0=0, y0=0, first_jaw=-pi/8, second_jaw=pi/6)
rect = PresetGeometries.Rectangle(Lx=0.85, Ly=0.6, x0=0, y0=0, theta0=pi/4)
ellipse = PresetGeometries.Ellipse(Rx=0.5, Ry=0.25, x0=0., y0=0., theta0=pi/6)
circle = PresetGeometries.Circle(R=0.95)

curve = pacman
# Plot random points according to whether they are inside the curve
nx = 50
ny = 50
x_coords = LinRange(-1,1, nx)
y_coords = LinRange(-1,1, ny)

test_pts = [(x_coords[i], y_coords[j]) for i=1:nx, j=1:ny]
# test_pts = [(0.25, -0.05), (0.25, 0.05)]
ds = 1e-3

# Plot the interior points
# interior_pts = zeros(Bool, size(test_pts))
# exterior_pts = zeros(Bool, size(test_pts))
# @. interior_pts = is_contained(test_pts, curve, ds=ds)
# @. exterior_pts = !interior_pts
# scatter(getindex.(test_pts[interior_pts], 1), getindex.(test_pts[interior_pts],2))

interior_pts = zeros(Bool, size(test_pts))
exterior_pts = zeros(Bool, size(test_pts))
@. interior_pts = is_contained(curve, test_pts)
@. exterior_pts = !interior_pts
scatter(getindex.(test_pts[interior_pts], 1), getindex.(test_pts[interior_pts],2))

# Plot the exterior/boundary points
scatter!(getindex.(test_pts[exterior_pts], 1), getindex.(test_pts[exterior_pts],2))

# Plot the curve's boundary
s = LinRange(0,1, 100)
boundary_pts = curve.(s)
scatter!(getindex.(boundary_pts, 1), getindex.(boundary_pts,2))

plot!(leg=false)
