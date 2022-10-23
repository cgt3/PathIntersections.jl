using SparseArrays
using ForwardDiff
using .PresetGeometries

const region_schemes = Set([:binary, :unique])

# Note: not using neighbors for now, just assuming curves define hard boundaries
# struct Neighbor{T_region, T_indices, Bool, T_curves}
#     region::T_region
#     global_index::T_indices
#     is_cartesian::Bool
#     shared_curves::T_curves
# end

# struct BoundaryCell{T_neighbors}
#     boundary::PiecewiseCurve
#     neighbors::T_neighbors
# end

# Rules:
# region_scheme = :binary && return_exterior = true, 
#  - only the exterior elements are returned
# region_scheme = :binary && return_exterior = false,
#  - only the interior elements are returned
# region_scheme = :unique,
#  - each curve gets its own region whose index corresponds to the curve's index
#  - the regions returns are the curves' interiors
#  - the region not occupied by any curve is the "exterior", which has index 0? (N+1)? 
#  - if regions/curves overlap curves that come later will occupy their intersection
#    

# TODO: supports multiple regions but does not compute their interior element information
# Assumes curves do not overlap and also do not share elements
function find_regions(mesh_coords, curves, stop_pts_by_curve, region_scheme::Symbol=:binary; return_exterior=true)
    if length(curves) != length(stop_pts_by_curve)
        error("find_regions: number of curves does not match stop pts list")
    end
    
    # For sizing the sparse matrices
    n,m = length(mesh_coords[1]), length(mesh_coords[2])
    num_curves = length(curves)
    if region_scheme == :binary
        num_region = 2
    else
        num_regions = num_curves
    end

    elements_by_region = [ spzeros(n,m) for i = 1:num_regions ]
    boundary_cells = SparseMatrixCSC{PiecewiseCurve, Int64}[] # Will work?
    all_elements = spzeros(n,m)
    for c = 1:num_curves
        dC_ds(s) = ForwardDiff.derivative(curves[c], s)

        if region_scheme == :binary
            region = 1
        else
            region = c
        end
        region_elements = spzeros(n,m)
        region_boundary_cells = spzeros(n,m)

        # 1) Find boundary cells for curve c
        for i = 1:length(stop_pts_by_curve[c])
            dC_ds(s) = ForwardDiff.derivative(curve, s)
            start_pt = stop_pts_by_curve[c][i]
            ij_start = stop_pt.indices

            # If no other curve has marked this cell yet, mark it
            if all_elements[mesh_indices] === nothing
                all_elements[mesh_indices] = region

                # TODO: if the intersection was a corner the cell indicated by the indices
                # may not be a cut cell; check on either side of the s-value
                
                # Compute the cut cell's boundary's stop points (starts at the first point in 
                # stop_pts_by_curve in this cell)
                cutcell_stop_pts = [0]
                cutcell_sub_bounds = [] # TODO: goodbye type stability
                cutcell_subcurves = []

                # Mark the portions of the curve in the cell
                sub_bound = [start_pt.s, 0]
                i_start = i
                start_node = [ mesh_coords[dim][start_pt.indices[dim]] for dim=1:length(start_pt.indices) ]
                n_start = outward_normal(dC_ds, start_pt.s)
                i_end = i_start + 1
                num_stop_pts = 1
                while start_pt.indices == stop_pts_by_curve[c][i_end].indices
                    sub_bound[2] = stop_pts_by_curve[c][i_end].s
                    push!(cutcell_sub_bounds, sub_bound)
                    push!(cutcell_subcurves, curves[c])
                    sub_bound[1] = sub_bound[2]
                    num_stop_pts += 1
                    i_end += 1
                end
                end_pt = stop_pts_by_curve[c][i_end]
                sub_bound[2] = end_pt.s
                push!(cutcell_sub_bounds, sub_bound)
                push!(cutcell_subcurves, curves[c])

                # Add the portions of the boundary coming from the Cartesian mesh
                n_end = outward_normal(dC_ds, end_pt.s)
                ij_end = end_pt.indices
                # n_mid = outward_normal(dC_ds, 0.5*(start_pt.s + end_pt.s))
                
                # Find the first mesh node on the boundary
                is_corner = (dim[1] == dim[2])
                if !is_corner
                    step = (!end_pt.dim).*sign.(n_end)
                    target_dim = dim[1] ? 1 : 2
                    mesh_node = [ mesh_coords[dim][end_pt.indices[dim]] for dim=1:length(end_pt.indices) ]
                    
                    num_stop_pts += 1
                    push!(cutcell_sub_bounds, [0,1]) # Line segment
                    if step[target_dim] < 0
                        push!(cutcell_subcurves, Line(end_pt.pt, mesh_node))
                    else
                        mesh_node[target_dim] = mesh_coords[target_dim][end_pt.indices[target_dim]+1]
                        push!(cutcell_subcurves, Line(end_pt.pt, mesh_node))
                    end

                    # Process the other mesh nodes (if they exist)
                    if step[1] != 0 # Always follows the CW orientation due to the use of right hand rule in def of normal?
                        step = -step
                    end
                    step = !step

                    while mesh_node != 
                        
                        push!(cutcell_subcurves, Line(end_pt.pt, mesh_node))
                    end
                else # corner
                    target_dim = 1
                end

                if !end_is_corner
                else
                end

                # Add the Cartesian nodes to the cell's boundaries
                # Find the face/corner we're on
                # if we're on a corner, consider the x direction first
                index_step = end_pt.dim.*n_end
                # + index_step[1]
                # + index_step[2]
                # - index_step[1]
                # - index_step[2]
                if stop_pts_by_curve[c][i_end].dim[]

                    # Figure out the direction to go in from the outward normal/tangents

                push!(cutcell_stop_pts, stop_pts_by_curve[c][i_end].s)



                # add the rest of the cell
                subcurves = 
                sub_bounds = 
                boundary_curve = PiecewiseCurve(cutcell_stop_pts, cutcell_subcurves, cutcell_sub_bounds)
                region_boundary_cells[mesh_indices] = BoundaryCell(boundary_curve, neighbors)
            end
        end
        # 2) Mark the interior cells for curve c

        push!(elements_by_region, region_elements)
        push!(boundary_cells, region_boundary_cells)
    end

    # TODO: do a sweep on the ultimate boundary and set those cells?    
    return all_elements, elements_by_region, boundary_cells
end