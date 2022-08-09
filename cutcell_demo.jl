## Libraries:
using GaussQuadrature
using LinearAlgebra
using Plots
using StaticArrays

using PathIntersections

## Curve parameters ------------------------------------------------------------------
pacman = PresetGeometries.Pacman(R=0.5, x0=0.13, y0=0.2)
pizza_slice = PresetGeometries.Pacman(R=0.5, x0=0, y0=0, first_jaw=-pi/8, second_jaw=pi/8)
rect = PresetGeometries.Rectangle(Lx=0.85, Ly=0.6, x0=-0.5, y0=0.5)
ellipse = PresetGeometries.Ellipse(Rx=0.5, Ry=0.25, x0=0., y0=0., theta0=pi/6)

ellipse_dg = PresetGeometries.Ellipse(Rx=0.5, Ry=0.25, x0=1e-3, y0=1e-3, theta0=pi/6)
all_curves = [pacman]
curve = pacman
test_pts = [ (0.2, 0.0), (0.5, 0.2), (-0.9, -0.5), (0.13, 0.5)]

## Intersection search parameters ----------------------------------------------------
# Mesh parameters
# Domain: [-1,1]^2
x_coords = LinRange(-1, 1, 21)
y_coords = LinRange(-1, 1, 21)

coords = [x_coords, y_coords]

arc_tol = 1e-12
corner_tol = 1e-12
ds = 1/100

ref_pts, ref_wts = legendre(5) # Legendre-Gauss quadrature
ref_quad = (ref_pts, ref_wts)
intersections = find_mesh_intersections(coords, all_curves, ds, arc_tol, corner_tol)
# regions, cutcell_i, cutcell_quad, cutcells = get_cutcell_nodes(coords, all_curves, ref_quad)


# # Plot the mesh
plot(leg=false)
for i in eachindex(x_coords)
    plot!(SVector{2}(x_coords[i], x_coords[i]), SVector{2}(extrema(y_coords)...))
end

for i in eachindex(y_coords)
    plot!(SVector{2}(extrema(x_coords)...), SVector{2}(y_coords[i], y_coords[i]))
end


# # Plot points along the curves
# s = LinRange(0, 1, 50)
# pts = @. rect(s)
# scatter!(getindex.(pts, 1), getindex.(pts, 2))

# pts = @. ellipse(s)
# scatter!(getindex.(pts, 1), getindex.(pts, 2))

# # Plot points along a specific cutcell
# s = LinRange(0, 1, 50)
# i1, i2 =  31, 31
# pts = @. cutcells[i1](s)
# scatter!(getindex.(pts, 1), getindex.(pts, 2))
# # pts = @. cutcells[i2](s)
# scatter!(getindex.(pts, 1), getindex.(pts, 2))


# Plot the quadrature nodes on a single cell
# i_cell = 3
# for i_cell = 9:16
# for i_face in 1:length(cutcell_quad[:pts][i_cell])
#     scatter!(getindex.(cutcell_quad[:pts][i_cell][i_face], 1), getindex.(cutcell_quad[:pts][i_cell][i_face], 2))
# end
# end

# # Plot the quadrature nodes on all cutcells
# for cell_pts in cutcell_quad[:pts]
#     for face in 1:length(cell_pts)
#         scatter!(getindex.(cell_pts[face], 1), getindex.(cell_pts[face], 2))
#     end
# end

# Plot stop points
c=1
pts = [ intersections[c][i].pt for i = 1:length(intersections[c]) ]
scatter!(getindex.(pts, 1), getindex.(pts, 2))

# Plot the extra point
is_inside = @. is_contained(test_pts, curve, ds = 1/25)
scatter!(getindex.(test_pts, 1), getindex.(test_pts,2))


plot!(leg=false)
