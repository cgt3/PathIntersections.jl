using ForwardDiff
using SparseArrays

function get_element_index(I_entry, I_next, tan_entry; tol=1e-12)
    diff = I_entry .- I_next
    abs_diff = abs.(diff)

    # If the indices are the same, return their common index
    if Tuple(I_entry) == Tuple(I_next)
        return Tuple(I_entry)
    end

    # If the indices are not on the same element return a invalid index
    for i = 1:length(abs_diff)
        if abs_diff[i] > 1
            return (-1,-1)
        end
    end

    i_element, j_element = 0, 0
    if sum(abs_diff) == 2
        i_element = minimum([I_entry[1], I_next[1]])
        j_element = minimum([I_entry[2], I_next[2]])
    elseif abs_diff[1] == 1
        i_element = minimum([I_entry[1], I_next[1]])
        j_element = I_entry[2] + (tan_entry[2] >= -tol ? 0 : -1)
    elseif abs_diff[2] == 1
        i_element = I_entry[1] + (tan_entry[1] >= -tol ? 0 : -1)
        j_element = minimum([I_entry[2], I_next[2]])
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

function define_regions(mesh_coords, curves, stop_pts; binary_regions=false)
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
        region_mask = zeros(Bool, nx, ny)

        # For each stop point on that curve
        i_prev = 0
        i = 1
        I_last = (0,0) # for the indices of the last element to be marked
        while i > i_prev && i < num_stop_pts
            i_prev = i
            i_entry = i
            entry_pt = stop_pts[c][i_entry]
            tan_entry = tangent(entry_pt.s)

            cutcell_curves = Function[]
            cutcell_sub_bounds = Tuple[] # TODO: can make type stable?
            sub_bounds = [entry_pt.s, 0.0]

            # 0) Make sure we have not left the domain
            if (entry_pt.indices[1] == 1 && tan_entry[1] < 0) || # x0
                (entry_pt.indices[1] == nx && tan_entry[1] > 0) || # xN
                (entry_pt.indices[2] == 1  && tan_entry[2] < 0) || # y0
                (entry_pt.indices[2] == ny && tan_entry[2] > 0)  # yN
                i += 1 # move to the next point if we've left the domain
            # Make sure we don't start inside an element
            elseif sum(entry_pt.dim) == 0
                i += 1
            else # We're inside an element
                # 1) Identify the element we're on based on the entry point and the next stop point
                #    and mark said element as cut.
                I_element = get_element_index(entry_pt.indices, stop_pts[c][i+1].indices, tan_entry)
                regions_by_element[I_element...] = region
                region_mask[I_element...] = true
                i_min, i_max = minimum([I_element[1], i_min]), maximum([I_element[1], i_max])
                j_min, j_max = minimum([I_element[2], j_min]), maximum([I_element[2], j_max])

                # 2) Find the exit point for this element and add any interior stop points
                num_subcurves = 0
                j = i
                j_next = j+1 <= num_stop_pts ? j+1 : 2
                while ( j == i || sum(stop_pts[c][j].dim) == 0 ) && I_element == get_element_index(entry_pt.indices, stop_pts[c][j_next].indices, tan_entry)
                    j = j_next
                    sub_bounds[2] = stop_pts[c][j].s
                    push!(cutcell_curves, curves[c])
                    push!(cutcell_sub_bounds, Tuple(sub_bounds))
                    if j == num_stop_pts
                        sub_bounds[1] = stop_pts[c][1].s
                    else
                        sub_bounds[1] = sub_bounds[2]
                    end
                    num_subcurves += 1
                    j_next = j+1 <= num_stop_pts ? j+1 : 2
                end
                i_exit = j
                exit_pt = stop_pts[c][i_exit]

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
            while i+1 < i_max && !(region_mask[i,j] == true && region_mask[i+1,j] == false) 
                i_start += 1
            end
            
            # Find where the interior ends
            i_end = i_max
            while i-1 < i_min && !(region_mask[i,j] == true && region_mask[i-1,j] == false) 
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

        # # Add the interior to the region map, but only overwrite zeros
        # @. interior = (regions_by_element == 0) && region_mask
        # @. regions_by_element[interior] =  -region

    end # for: c in 1:num_curves

    return regions_by_element, cutcell_indices, cutcells
end # function

function get_cutcell_nodes(mesh_coords, curves, ref_quad; 
    ds=DEFAULT_DS,
    arc_tol=DEFAULT_ARC_TOL,
    corner_tol=DEFAULT_CORNER_TOL,
    binary_regions=DEFAULT_BINARY_REGIONS,
    ref_domain=DEFAULT_REF_DOMAIN,
    normalization=DEFAULT_NORMALIZATION )

    # 1) Get mesh intersections and curve stop points
    stop_pts = find_mesh_intersections(mesh_coords, curves, ds, arc_tol, corner_tol)

    # 2) Calculate cutcells
    regions_by_element, cutcell_indices, cutcells = define_regions(mesh_coords, curves,
        stop_pts, binary_regions=binary_regions)

    # 3) Compute quadrature on the cutcell boundaries
    cutcell_quadratures = Tuple[]
    for c = 1:length(cutcells)
        cutcell_quad = map_line_quadrature(ref_quad, cutcells[c], cutcells[c].stop_pts,
            ref_domain=ref_domain, normalization=normalization )
        push!(cutcell_quadratures, cutcell_quad)
    end
    
    return regions_by_element, cutcell_indices, cutcell_quadratures, cutcells
end