## Libraries:
using GaussQuadrature
using LinearAlgebra
using Plots
using Revise
using StaticArrays

using PathIntersections

## Curve parameters ------------------------------------------------------------------
rect = PresetGeometries.Rectangle(Lx=0.85, Ly=0.6, x0=0, y0=0)
ellipse = PresetGeometries.Ellipse(Rx=0.25, Ry=0.15, x0=0.5, y0=-0.5, theta0=pi/6)

## Intersection search parameters ----------------------------------------------------
# Mesh parameters
# Domain: [-1,1]^2
x_coords = LinRange(-1, 1, 11)
y_coords = LinRange(-1, 1, 11)

coords = [x_coords, y_coords]

arc_tol = 1e-8
corner_tol = 1e-8
ds = 1/100

ref_pts, ref_wts = legendre(5) # Legendre-Gauss quadrature
ref_quad = (ref_pts, ref_wts)
# intersections = find_mesh_intersections(coords, curve, ds, arc_tol, corner_tol)
regions, cutcell_i, cutcell_quad, cutcells = get_cutcell_nodes(coords, [rect], ref_quad)


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
# i = 3
# pts = @. cutcells[i](s)
# scatter!(getindex.(pts, 1), getindex.(pts, 2))
# pts = @. cutcells[i+4](s)
# scatter!(getindex.(pts, 1), getindex.(pts, 2))


# # Plot the quadrature nodes on a single cell
# i_cell
# for i_face in 1:length(cutcell_quad[i_cell][1])
#     scatter!(getindex.(cutcell_quad[i_cell][1][i_face], 1), getindex.(cutcell_quad[i_cell][1][i_face], 2))
# end

# Plot the quadrature nodes on all cutcells
for cell_quad in cutcell_quad
    for face in 1:length(cell_quad[1])
        scatter!(getindex.(cell_quad[1][face], 1), getindex.(cell_quad[1][face], 2))
    end
end

# # Plot stop points
# pts = [ intersections[i].pt for i = 1:length(intersections) ]
# scatter!(getindex.(pts, 1), getindex.(pts, 2))
# println("Finished plotting intersections and stop pts")
plot!(leg=false)