using ForwardDiff
using SparseArrays

# Note: this function only works properly for points on the boundaries
function get_element_index(entry_pt, exit_pt, tan_entry, mesh_coords; tol=1e-12)
    diff = entry_pt.indices .- exit_pt.indices
    abs_diff = abs.(diff)

    # If the indices are the same, return their common index
    if Tuple(entry_pt.indices) == Tuple(exit_pt.indices)
        return Tuple(entry_pt.indices)
    end

    # If the indices are not on the same element return a invalid index
    for i = 1:length(abs_diff)
        if abs_diff[i] > 1
            @warn "Point indices are not on the same element: entry=$(entry_pt.indices), exit=$(exit_pt.indices)"
            return (-1,-1)
        end
    end
    i_element, j_element = 0, 0
    if sum(abs_diff) == 2
        i_element = min(entry_pt.indices[1], exit_pt.indices[1])
        j_element = min(entry_pt.indices[2], exit_pt.indices[2])
    elseif abs_diff[1] == 1
        i_element = min(entry_pt.indices[1], exit_pt.indices[1])
        if entry_pt.pt[1] != mesh_coords[1][entry_pt.indices[1]] ||
           exit_pt.pt[1]  != mesh_coords[1][entry_pt.indices[1]] ||
           tan_entry[1] > 0
            j_element = entry_pt.indices[2]
        else
            j_element = entry_pt.indices[2] - 1
        end
    elseif abs_diff[2] == 1
        if entry_pt.pt[2] != mesh_coords[2][entry_pt.indices[2]] ||
           exit_pt.pt[2]  != mesh_coords[2][entry_pt.indices[2]] ||
           tan_entry[2] > 0
            i_element = entry_pt.indices[1]
        else
            i_element = entry_pt.indices[1] - 1
        end
        j_element = min(entry_pt.indices[2], exit_pt.indices[2])
    end
    return (i_element, j_element)
end

function get_face_index(point, I_element)
    # Note: the order of the if/elseif matters for correctly classifying corners
    # Bottom face
    if point.indices[2] == I_element[2] && point.indices[1] != I_element[1] + 1 && point.dim[2] == true
        i_face = 1
    # Left face
    elseif point.indices[1] == I_element[1] && point.dim[1] == true
        i_face = 2
    # Top face
    elseif point.indices[2] == I_element[2] + 1
        i_face = 3
    # Right face
    else
        i_face = 4
    end
    # Check if the point is on the corner of the face
    is_corner = false
    if point.dim[1] && point.dim[2]
        is_corner = true
    end
    return i_face, is_corner
end

function define_regions(mesh_coords, curves, stop_pts; binary_regions=false, edge_tol=1e-12)
    Line = PresetGeometries.Line

    nx = length(mesh_coords[1])
    ny = length(mesh_coords[2])

    regions_by_element = zeros(Int, nx-1, ny-1)
    cutcell_indices = spzeros(Int, nx-1, ny-1)
    cutcells = PiecewiseCurve[]

    num_curves = length(curves)
    
    # For incrementing the element index to corner/face indices
    # i.e. corner_index = element_index + pt_incr[face index]
    pt_incr = ((0,0),(0,1),(1,1),(1,0))
    cutcell_index = 1
    # For each curve
    for c in 1:num_curves
        if binary_regions == true
            region = 1
        else
            region = c
        end
        num_stop_pts = length(stop_pts[c])
        tangent(s) = ForwardDiff.derivative(curves[c], s)


        # For filling in the curve's region later
        i_min, i_max = nx, 1
        j_min, j_max = ny, 1
        region_mask = zeros(Bool, nx-1, ny-1)

        # For each stop point on that curve
        i_prev = 0
        i = 1
        I_last = (0,0) # for the indices of the last element to be marked
        i0_was_skipped = false
        while (i < num_stop_pts && !i0_was_skipped) || (i > i_prev && i0_was_skipped)
            i_prev = i
            i_entry = i
            entry_pt = stop_pts[c][i_entry]
            tan_entry = tangent(entry_pt.s)

            cutcell_curves = Function[]
            cutcell_sub_bounds = Tuple[] # TODO: can make type stable?
            sub_bounds = [1.0*entry_pt.s, entry_pt.s]

            # 0) Make sure we have not left the domain
            if  (abs(entry_pt.pt[1] - mesh_coords[1][1] ) < edge_tol && tan_entry[1] < 0) || # x0
                (abs(entry_pt.pt[1] - mesh_coords[1][nx]) < edge_tol && tan_entry[1] > 0) || # xN
                (abs(entry_pt.pt[2] - mesh_coords[2][1] ) < edge_tol && tan_entry[2] < 0) || # y0
                (abs(entry_pt.pt[2] - mesh_coords[2][ny]) < edge_tol && tan_entry[2] > 0)  # yN
                i += 1 # move to the next point if we've left the domain
            # Make sure we don't start inside an element
            elseif sum(entry_pt.dim) == 0
                if i == 1
                    i0_was_skipped = true
                end
                i += 1
            else # We're on the boundary of an element
                # 1) Find the exit point for this element and add any interior stop points
                num_subcurves = 0
                j = i
                j_next = j < num_stop_pts ? j+1 : 2
                found_exit = false
                while found_exit == false
                    # If we wrap all the way around the curve make sure to use the proper sub-bounds
                    if stop_pts[c][j_next].s < stop_pts[c][j].s #j_next == num_stop_pts 
                        sub_bounds[1] = 0 - (1 - stop_pts[c][j].s)
                    else
                        sub_bounds[1] = sub_bounds[2]
                    end
                    sub_bounds[2] = stop_pts[c][j_next].s
                    push!(cutcell_curves, curves[c])
                    push!(cutcell_sub_bounds, Tuple(sub_bounds))
                    num_subcurves += 1

                    # WARNING: If tunneling is to be allowed this condition won't
                    # suffice
                    if sum(stop_pts[c][j_next].dim) != 0
                        found_exit = true
                    end
                    j = j_next
                    j_next = j < num_stop_pts ? j+1 : 2
                end
                i_exit = j
                exit_pt = stop_pts[c][i_exit]
                
                # 2) Identify the element we're on based on the entry and exit point
                #    and mark said element as cut.
                I_element = get_element_index(entry_pt, exit_pt, tan_entry, mesh_coords)
                if regions_by_element[I_element...] != 0
                    @warn "define_regions: tunneling has occured; element $I_element (cutcell $(length(cutcells)+1)) has been assigned more than one cutcell"
                end
                regions_by_element[I_element...] = region
                region_mask[I_element...] = true
                i_min, i_max = minimum([I_element[1], i_min]), maximum([I_element[1], i_max])
                j_min, j_max = minimum([I_element[2], j_min]), maximum([I_element[2], j_max])


                # 3) Determine the faces the entry and exit points are on
                f_entry, _ = get_face_index(entry_pt, I_element)
                f_exit, exit_is_corner = get_face_index(exit_pt, I_element)
                
                f_entry -= 1 # don't add the Cartesian node on the entry pt's face
                if exit_is_corner
                    f_exit += 1
                end

                if f_entry < f_exit
                    f_entry += 4
                end

                # 4) Add the Cartesian nodes to the boundary curves
                curr_pt = exit_pt.pt
                next_pt = curr_pt # Just so this var persists past the while loop
                for f = f_exit:f_entry
                    next_pt = ( mesh_coords[1][I_element[1] + pt_incr[(f+3) % 4 + 1][1]], 
                                mesh_coords[2][I_element[2] + pt_incr[(f+3) % 4 + 1][2]] )
                    push!(cutcell_curves, Line(Tuple(curr_pt), Tuple(next_pt)))
                    push!(cutcell_sub_bounds, (0,1))
                    num_subcurves += 1
                    curr_pt = next_pt
                end
                # 5) Close the curve
                push!(cutcell_curves, Line(Tuple(next_pt), Tuple(entry_pt.pt)))
                push!(cutcell_sub_bounds, (0,1))
                num_subcurves += 1

                cutcell_stop_pts = [k/num_subcurves for k=0:num_subcurves]
                push!(cutcells, PiecewiseCurve(cutcell_stop_pts, cutcell_curves, cutcell_sub_bounds))
                
                # 6) Add the cut cell and its index to their respective arrays
                cutcell_indices[I_element...] = cutcell_index
                cutcell_index += 1

                # 7) Increment i for the next element
                i = i_exit 
                I_last = I_element
            end # if: whether we're in the domain or not
        end # while num_stop_pts[c]

        # 8) Mark elements inside the curve as excluded
        for j = j_min:j_max
            # Find where the interior starts
            i_start = i_min
            while i_start+1 < i_max && !(region_mask[i_start,j] == true && region_mask[i_start+1,j] == false) 
                i_start += 1
            end
            
            i_end = i_max
            while i_end-1 > i_min && !(region_mask[i_end,j] == true && region_mask[i_end-1,j] == false) 
                i_end -= 1
            end

            for i = i_start:i_end # Note is i_end < i_start this should have no effect
                region_mask[i,j] = true
            end
        end

        # # WARNING: This if-statement may not catch every case
        # # Check the direction of the normal to see if the mask needs to be flipped
        # n_out = outward_normal(tangent, 1)
        # if  0 < (sign(n_out[1]) + i_last) <= nx && 0 < (sign(n_out[2]) + j_last) <= ny &&
        #     (!region_mask[sign(n_out[1]) + i_last, j_last] || !region_mask[i_last, sign(n_out[2]) + j_last] )
        #     region_mask = !region_mask
        # end

        # Add the interior to the region map, but only overwrite zeros
        interior = @. (regions_by_element == 0) && region_mask
        regions_by_element[interior] .= -region
        regions_by_element

    end # for: c in 1:num_curves

    return regions_by_element, cutcell_indices, cutcells
end # function

function get_cutcell_nodes(mesh_coords, curves, ref_quad; 
    ds=DEFAULT_DS,
    arc_tol=DEFAULT_ARC_TOL,
    corner_tol=DEFAULT_CORNER_TOL,
    binary_regions=DEFAULT_BINARY_REGIONS,
    ref_domain=DEFAULT_REF_DOMAIN,
    normalization=DEFAULT_NORMALIZATION,
    closure_tol=1e-12,
    normals_wrt_cell=true )

    # 1) Get mesh intersections and curve stop points
    stop_pts = find_mesh_intersections(mesh_coords, curves, ds, arc_tol, corner_tol,
        closed_list=true, closure_tol=closure_tol)

    # 2) Calculate cutcells
    regions_by_element, cutcell_indices, cutcells = define_regions(mesh_coords, curves,
        stop_pts, binary_regions=binary_regions)

    # 3) Compute quadrature on the cutcell boundaries
    all_pts = Vector{Vector{Vector{Float64}}}[]
    all_wts = Vector{Vector{Float64}}[]
    all_n =   Vector{Vector{Vector{Float64}}}[]
    for c = 1:length(cutcells)
        cell_pts, cell_wts, cell_n = map_line_quadrature(ref_quad, cutcells[c], cutcells[c].stop_pts,
            ref_domain=ref_domain, normalization=normalization )

        if normals_wrt_cell == true
            cell_n = @. -cell_n
        end
        push!(all_pts, cell_pts)
        push!(all_wts, cell_wts)
        push!(all_n, cell_n)
    end
    all_quadratures = Dict(:pts => all_pts, :wts => all_wts, :n => all_n)
    
    return regions_by_element, cutcell_indices, all_quadratures, cutcells
end