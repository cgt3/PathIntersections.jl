## Helper Functions -------------------------------------------------------------
function get_ds(ds::Real, s)
    return ds
end

function get_ds(ds::Function, s)
    return ds(s)
end
function update_intersection_bounds(pt, dims, coords, nearby_indices)
    indices = copy(nearby_indices)
    for d = 1:length(dims)
        # Goal: ensure the indices not involved in the intersection point to 
        #       mesh values <= the pt's coordinates (i.e. coords[d] <= pt[d])
        if dims[d] == false # this index was not involved in the intersection
            # Bring the index up if it is too low
            while indices[d] < length(coords[d]) && coords[d][indices[d]+1] <= pt[d]
                indices[d] += 1 
            end
            # Bring the index down if it is too high
            while indices[d] > 1 && coords[d][indices[d]] > pt[d]
                indices[d] -= 1
            end
        end
    end
    return indices
end

function tighten_bounds!(pt_new, dim, coords, indices_lb, indices_ub)
    # Make sure the indices are as low as possible
    while indices_lb[dim] > 1 && pt_new[dim] <= coords[dim][indices_lb[dim]]
        indices_lb[dim] -= 1
    end
    while indices_ub[dim] > 1 && pt_new[dim] < coords[dim][indices_ub[dim] - 1]
        indices_ub[dim] -= 1
    end
    # Make sure the indices are as high as possible
    while indices_lb[dim] < length(coords[dim]) && pt_new[dim] > coords[dim][indices_lb[dim] + 1]
        indices_lb[dim] += 1
    end
    while indices_ub[dim] < length(coords[dim]) && pt_new[dim] >= coords[dim][indices_ub[dim]]
        indices_ub[dim] += 1
    end
end


function secant_single_dim(intercept, targetDim, curve, s_lb, s_ub, arc_tol)
    pt_lb = curve(s_lb)
    pt_ub = curve(s_ub)

    # So these variables persist past the while loop
    s_new = s_lb
    pt_new = pt_lb

    pt_dist = arc_tol + 1;
    while pt_dist > arc_tol
        # Find the secant intercept from s_lb and s_ub:
        s_new = (intercept - pt_lb[targetDim]) / (pt_ub[targetDim] - pt_lb[targetDim]) * (s_ub - s_lb) + s_lb;
        pt_new = curve(s_new)

        # Update one of the bounds with the new point
        if pt_new[targetDim] < intercept
            pt_lb = pt_new
            s_lb = s_new
        elseif pt_new[targetDim] > intercept
            pt_ub = pt_new
            s_ub = s_new
        else    
            return s_new, pt_new
        end

        # Update the other bound to be closer to avoid stalling under
        # the arc length condition
        s_avg = 0.5*(s_ub + s_lb)
        pt_avg = curve(0.5*(s_ub + s_lb))
        if s_ub != s_new && pt_avg[targetDim] > intercept
            s_ub = s_avg
            pt_ub = pt_avg
        elseif s_lb != s_new && pt_avg[targetDim] < intercept
            s_lb = s_avg
            pt_lb = pt_avg
        elseif pt_avg[targetDim] == intercept
            return s_avg, pt_avg
        end

        # Calculate the Euclidean distance between the two bounds to check
        # against tol
        pt_dist = norm(pt_ub .- pt_lb,2)
    end

    return s_new, pt_new
end

function outward_normal(dC_ds, s; normalization=false)
    tangent = dC_ds(s)
    # Note: need array if normalization is selected because tuples are immutable
    normal = SVector{2}(tangent[2], -tangent[1])
    if normalization == true
        return normal ./ norm(normal)
    end

    return normal
end

function insert_sorted!(intersections, new_pt; include_duplicates=false, s_tol=1e-8)
    num_pts = length(intersections)
    i = num_pts
    if num_pts == 0 || intersections[end].s < new_pt.s
        push!(intersections, new_pt)
        return
    end

    if include_duplicates == true
        if intersections[i].s <= new_pt.s
            push!(intersections, new_pt)
        else
            push!(intersections, intersections[i])
            while i > 1 && new_pt.s < intersections[i].s
                intersections[i] = intersections[i-1]
                i -= 1
            end
            intersections[i] = new_pt
        end
    else # include_duplicate == true
        # check if this point already exists
        while i > 1 && new_pt.s < intersections[i].s 
            i -= 1
        end
        if abs(intersections[i].s - new_pt.s) < s_tol
            return
        else
            push!(intersections, intersections[num_pts])
            for j = num_pts:-1:i+1
                intersections[j] = intersections[j-1]
            end
            intersections[i] = new_pt
        end
    end
end