## Libraries:
using LinearAlgebra
using Plots
using StaticArrays
using Revise

using PathIntersections

## Curve parameters ------------------------------------------------------------------
# curve = PresetGeometries.Rectangle(Lx=0.5, Ly=0.3, x0=0, y0=0)
curve = PresetGeometries.Ellipse(Rx=0.5, Ry=0.3, x0=0, y0=0)
## Intersection search parameters ----------------------------------------------------
# Mesh parameters
# Domain: [-1,1]^2
x_coords = LinRange(-1, 1, 5)
y_coords = LinRange(-1, 1, 5)

coords = [x_coords, y_coords]

arc_tol = 1e-8
corner_tol = 1e-8
ds = 1/100

intersections = find_mesh_intersections(coords, curve, ds, arc_tol, corner_tol)


# # Plot the mesh
# pts = intersections[1].pt
plot(leg=false)
pts = [ intersections[i].pt for i = 1:length(intersections) ]
scatter!(getindex.(pts, 1), getindex.(pts, 2))
# plot x lines
for i in eachindex(x_coords)
    plot!(SVector{2}(x_coords[i], x_coords[i]), SVector{2}(extrema(y_coords)...))
    plot!(SVector{2}(extrema(x_coords)...), SVector{2}(y_coords[i], y_coords[i]))
end
plot!(leg=false)