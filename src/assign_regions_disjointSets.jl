function cutCartesianShareSide(I_cart, I_cut, mesh_coords, cutcell_indices, cutcells)
    # I_cart and I_cut should be off by one in one direction
    d_diff = I_cart[1] == I_cut[1] ? 2 : 1

    if I_cart[d_diff] < I_cut[d_diff]
        target_coord_d = mesh_coords[d_diff][I_cut[d_diff]]
    else
        target_coord_d = mesh_coords[d_diff][I_cart[d_diff]]
    end

    # Check the faces of the cut cell that are lines
    # To be neighbors the cut cell has to have an uncut Cartesian face
    cut_element_faces = cutcells[cutcell_indices[I_cut...]].subcurves
    for f in eachindex(cut_element_faces)
        if cut_element_faces[f] isa PresetGeometries.Line
            if abs(cut_element_faces[f].start_pt[d_diff] - cut_element_faces[f].end_pt[d_diff]) < 1e-10 &&
               abs(cut_element_faces[f].start_pt[d_diff] - target_coord_d) < 1e-10
               return true
            end
        end 
    end

    return false
end

function cutCutShareSide(I1, I2, mesh_coords, cutcell_indices, cutcells)
        # I1 and I2 should be off by one in one direction
        d_diff = I1[1] == I2[1] ? 2 : 1

        if I1[d_diff] < I2[d_diff]
            target_coord_d = mesh_coords[d_diff][I2[d_diff]]
        else
            target_coord_d = mesh_coords[d_diff][I1[d_diff]]
        end
    
        # Check the faces of the cut cell that are lines
        cut_element_faces1 = cutcells[cutcell_indices[I1...]].subcurves
        for f1 in eachindex(cut_element_faces1)
            if cut_element_faces1[f1] isa PresetGeometries.Line
                if abs(cut_element_faces1[f1].start_pt[d_diff] - cut_element_faces1[f1].end_pt[d_diff]) < 1e-10 &&
                    abs(cut_element_faces1[f1].start_pt[d_diff] - target_coord_d) < 1e-10
                    
                    cut_element_faces2 = cutcells[cutcell_indices[I2...]].subcurves
                    for f2 in eachindex(cut_element_faces2)
                        if cut_element_faces2[f2] isa PresetGeometries.Line
                            if abs(cut_element_faces2[f2].start_pt[d_diff] - cut_element_faces2[f2].end_pt[d_diff]) < 1e-10 &&
                               abs(cut_element_faces2[f2].start_pt[d_diff] - target_coord_d) < 1e-10
                               return true
                            end
                        end
                    end
                end
            end 
        end
    
        return false
end


function is_adjacent(I, I_nbr, mesh_coords, cutcell_indices, cutcells)
    # Both elements are Cartesian
    if cutcell_indices[I...] == 0 && cutcell_indices[I_nbr...] == 0
        return true
    # I is Cartesian, I_nbr is cut
    elseif cutcell_indices[I...] == 0 && cutcell_indices[I_nbr...] != 0
        return cutCartesianShareSide(I, I_nbr, mesh_coords, cutcell_indices, cutcells)
    # I is cut, I_nbr is Cartesian
    elseif cutcell_indices[I...] != 0 && cutcell_indices[I_nbr...] == 0
        return cutCartesianShareSide(I_nbr, I, mesh_coords, cutcell_indices, cutcells)
    # Both elements are cut
    else
        return cutCutShareSide(I, I_nbr, mesh_coords, cutcell_indices, cutcells)
    end
end


function assign_regions!(regions_by_element, mesh_coords, cutcell_indices, cutcells; binary_regions=true)
    size_mesh = size(cutcell_indices)

    set = Matrix{Vector{Int64}}(undef, size_mesh)

    rep_value = Vector{Int64}[]
    count = 0
    for j in axes(set,2)
        for i in axes(set,1)
            # Check preceeding neighbors
            join_top_neighbor = false
            if (i-1) > 0 && is_adjacent((;i,j), (;i=i-1, j), mesh_coords, cutcell_indices, cutcells)
                join_top_neighbor = true
            end
            
            join_left_neighbor = false
            if (j-1) > 0 && is_adjacent((;i,j), (;i, j=j-1), mesh_coords, cutcell_indices, cutcells)
                join_left_neighbor = true
            end

            if join_top_neighbor
                if join_left_neighbor
                    # Merge the top and left neighbor sets
                    # When peforming a union take the top set's rep value
                    set[i,j] = set[i-1,j]
                    if set[i-1,j][1] != set[i,j-1][1]
                        overwritten_val = set[i,j-1][1]
                        for k in eachindex(rep_value)
                            if rep_value[k][1] == overwritten_val
                                rep_value[k][1] = set[i-1,j][1]
                            end
                        end
                    end
                else
                    # Join the top neigbor's set only
                    set[i,j] = set[i-1, j]
                end
            else
                if join_left_neighbor
                    # Join the left neighbor's set only
                    set[i,j] = set[i,j-1]
                else
                    # Create a new set
                    count += 1
                    push!(rep_value, [count])
                    set[i,j] = rep_value[count]
                end
            end #join_top_neighbor
        end # for j
    end  # for i

    # Use the disjoint sets to determine the regions
    # Note: this only works for binary regions (i.e. embedded boundaries, not interfaces)
    # Determine the active region by looking at a cut cell
    i_cut, j_cut, _ = findnz(cutcell_indices)
    active_region = set[i_cut[1], j_cut[1]][1]
    for i in axes(regions_by_element, 1)
        for j in axes(regions_by_element, 2)
            if set[i,j][1] != active_region && regions_by_element[i,j] == 0
                regions_by_element[i,j] = -1
            end
        end
    end
end