## Helper Functions -------------------------------------------------------------
function get_ds(ds::Real, s::Real)
    return ds
end

function get_ds(ds::Function, s::Real)
    return ds(s)
end

function tighten_bounds(pt_new, dim, coords, indices_lb, indices_ub)
    while indices_lb[dim] > 1 && pt_new[dim] <= coords[dim][indices_lb[dim]]
        indices_lb[dim] -= 1
    end
    while indices_ub[dim] > 1 && pt_new[dim] < coords[dim][indices_ub[dim] - 1]
        indices_ub[dim] -= 1
    end

    while indices_lb[dim] < length(coords[dim]) && pt_new[dim] > coords[dim][indices_lb[dim] + 1]
        indices_lb[dim] += 1
    end
    while indices_ub[dim] < length(coords[dim]) && pt_new[dim] >= coords[dim][indices_ub[dim]]
        indices_ub[dim] += 1
    end
end


function secant_single_dim(intercept::Real, targetDim::Integer, curve::Function, s_lb::Real, s_ub::Real, arc_tol::Real)
    pt_lb = curve(s_lb)
    pt_ub = curve(s_ub)

    # So these variables persist past the while loop
    s_new = s_lb
    pt_new = pt_lb

    pt_dist = arc_tol + 1;
    while pt_dist > arc_tol
        # Find the secant intercept from s_lb and s_ub:
        s_new = (intercept - pt_lb[targetDim]) / (pt_ub[targetDim] - pt_lb[targetDim]) * (s_ub-s_lb) + s_lb;
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
        pt_dist = norm(pt_ub - pt_lb,2)
    end

    return s_new, pt_new
end
