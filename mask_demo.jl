## Libraries:
using GaussQuadrature
using LinearAlgebra
using Plots
using Revise
using StaticArrays

using PathIntersections

## Curve parameters ------------------------------------------------------------------
pacman = PresetGeometries.Pacman(R=0.5, x0=0.13, y0=0.2)
pizza_slice = PresetGeometries.Pacman(R=0.5, x0=0, y0=0, first_jaw=-pi/8, second_jaw=pi/8)
rect = PresetGeometries.Rectangle(Lx=0.85, Ly=0.6, x0=-0.5, y0=0.5)
ellipse = PresetGeometries.Ellipse(Rx=0.5, Ry=0.25, x0=0., y0=0., theta0=pi/6)
ellipse_dg = PresetGeometries.Ellipse(Rx=0.5, Ry=0.25, x0=1e-3, y0=1e-3, theta0=pi/6)

curve = pacman
# Plot the boundary of the curve
s = 0:1e-3:1
boundary_pts = curve.(s)
scatter(getindex.(boundary_pts,1), getindex.(boundary_pts,2))


# Plot random points according to whether they are inside the curve
nx = 100
ny = 100
x_coords = LinRange(-1,1, nx)
y_coords = LinRange(-1,1, ny)

test_pts = [(x_coords[i], y_coords[j]) for i=1:nx, j=1:ny]
# test_pts = [(0.131313, 0.333)]
ds = 1/1000

# Plot the interior points
interior_pts = @. is_contained(test_pts, curve, ds=ds)
scatter!(getindex.(test_pts[interior_pts], 1), getindex.(test_pts[interior_pts],2))

# Plot the exterior/boundary points
scatter!(getindex.(test_pts[.!interior_pts], 1), getindex.(test_pts[.!interior_pts],2))

plot!(leg=false)
