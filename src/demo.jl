## Libraries:
using LinearAlgebra
using Plots
using StaticArrays
using Revise

using PathIntersections

## Curve parameters ------------------------------------------------------------------
circle = PresetGeometries.Circle(R=0.2, x0=-0.5, y0=0.5)
ellipse = PresetGeometries.Ellipse(Rx=0.4, Ry=0.2, x0=0.5, y0=0.5, theta0=pi/4, orientation=-1)
pacman = PresetGeometries.Pacman(R=0.5, y0=-0.3, first_jaw=pi/4, second_jaw=7pi/4, orientation=1)

## Intersection search parameters ----------------------------------------------------
# Mesh parameters
# Domain: [-1,1]^2
x_coords = LinRange(-1, 1, 50)
y_coords = LinRange(-1, 1, 50)

coords = [x_coords, y_coords]

arc_tol = 1e-8
corner_tol = 1e-8
ds = [1/100, 1/100]#, 1/300]
# curves = [circle, ellipse, pacman];
curves = [circle, ellipse]#, pacman];

intersections_by_curve = find_mesh_intersections(coords, curves, ds, arc_tol, corner_tol)


# # Plot the mesh
# pts = intersections[1].pt
plot(leg=false)
for intersections in intersections_by_curve
    # TODO: fix for change from StructArrays
    pts = [ intersections[i].pt for i = 1:length(intersections) ]
    scatter!(getindex.(pts, 1), getindex.(pts, 2))
    # plot x lines
    for i in eachindex(x_coords)
        plot!(SVector{2}(x_coords[i], x_coords[i]), SVector{2}(extrema(y_coords)...))
        plot!(SVector{2}(extrema(x_coords)...), SVector{2}(y_coords[i], y_coords[i]))
    end
end
plot!(leg=false)