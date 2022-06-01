# Function for finding multiple dim intersections against a single curve
function find_mesh_intersections(coords, curve::Function,
     ds=DEFAULT_DS,
    arc_tol=100*eps(eltype(coords[1])),
    corner_tol=100*eps(eltype(coords[1])),
    closed_list=true,
    closure_tol=1e-12)

    s_tol = 1e-12
    numDim = length(coords)
    # Default step size to something (approximately) smaller than the smallest mesh size
    if ds == DEFAULT_DS
        # Find the minimum mesh size
        min_by_dim = zeros(numDim, 1)
        for d = 1:numDim
            min_by_dim[d] = minimum( coords[d][2:end] - coords[d][1:end-1] )
        end
        dx_min = minimum(min_by_dim)

        # Sample the derivative of the curve to approximate its max deriv value
        s_sample = LinRange(0,1,1000)
        dcurve_ds(s) = ForwardDiff.derivative(curve, s)
        deriv_sample = @. dcurve_ds(s_sample)
        deriv_mag = @. norm(deriv_sample)
        max_deriv_mag = maximum(deriv_mag)

        ds = minimum([1/1000, dx_min / max_deriv_mag])
    end

    # Start walking along the curve at s = 0, watching for the curve's stop points
    if :stop_pts in fieldnames(typeof(curve))
        stop_pts = curve.stop_pts
        i_stop_pts = 1
    else
        i_stop_pts = -1
    end

    s = 0
    pt_curr = curve(s)
    
    # Find the indices to the closest mesh values in each dimension
    indices_ub = ones(Int, numDim)
    indices_lb = ones(Int, numDim)
    for d = 1:numDim
        # TODO: can change to binary search since the coord vectors are sorted
        while indices_ub[d] < length(coords[d]) && coords[d][indices_ub[d]] < pt_curr[d]
            indices_ub[d] += 1 
        end
        
        indices_lb[d] = indices_ub[d]
        while indices_lb[d] > 1 && coords[d][indices_lb[d]] >= pt_curr[d]
            indices_lb[d] -= 1
        end
    end

    # Walk along the curve, sensing for intersections in each dim
    intersections = MeshIntersection[]
    intersectionIndex = 0
    dim = Array{Bool,1}(undef, numDim)

    s_new = s
    pt_new = pt_curr
    while s <= 1
        # Check if we are currently in the domain or crossing into or out of it
        curr_isInsideDomain = true
        new_isInsideDomain = true
        for d = 1:numDim
            if pt_curr[d] < coords[d][1] || pt_curr[d] > coords[d][end]
                curr_isInsideDomain = false
            end
            if pt_new[d] < coords[d][1] || pt_new[d] > coords[d][end]
                new_isInsideDomain = false
            end
        end
        isInsideDomain = curr_isInsideDomain | new_isInsideDomain

        compoundIntersection = false
        wasOnBoundary = false
        for d = 1:numDim
            # WARNING: lines of intersection will be counted as a series of intersections
            #          that are ds distance apart

            # If the last intersection did not involve this dimension, check it for an
            # intersection
            if compoundIntersection == false || intersections[intersectionIndex].dim[d] == false
                dim[d] = false # Since using undef doesn't guarantee dim is initialized to false
                intersectionOccurred = false
                wasOnBoundary = false
                indices = copy(indices_lb)

                # If we're still inside the domain, check for an intersection; if not
                # just update the mesh bounds
                if isInsideDomain
                    # The new point is on the lower bound
                    if pt_new[d] == coords[d][indices_lb[d]]
                        s_intercept, pt_intercept = s_new, pt_new
                        intersectionOccurred = true
                        wasOnBoundary = true
                    # The new point is on the upper bound
                    elseif pt_new[d] == coords[d][indices_ub[d]]
                        s_intercept, pt_intercept = s_new, pt_new
                        indices[d] = indices_ub[d]
                        intersectionOccurred = true
                        wasOnBoundary = true
                    # The new point crossed the lower bound
                    elseif pt_new[d] < coords[d][indices_lb[d]]
                        s_intercept, pt_intercept = secant_single_dim(coords[d][indices_lb[d]], d, curve, s_new, s, arc_tol)
                        intersectionOccurred = true
                    # The new point crossed the upper bound
                    elseif pt_new[d] > coords[d][indices_ub[d]]
                        s_intercept, pt_intercept = secant_single_dim(coords[d][indices_ub[d]], d, curve, s, s_new, arc_tol)
                        indices[d] = indices_ub[d]
                        intersectionOccurred = true
                    end # Otherwise, if no boundaries are crossed keep walking
                
                    if intersectionOccurred
                        # If an intersection was encountered, check if it crossed boundaries in
                        # the other dimensions before adding it
                        for i = 1:numDim
                            if abs(pt_intercept[i] - coords[i][indices_lb[i]]) < corner_tol 
                                dim[i] = true
                                indices[i] = indices_lb[i]
                                if i != d
                                    compoundIntersection = true
                                end
                            elseif abs(pt_intercept[i] - coords[i][indices_ub[i]]) < corner_tol 
                                dim[i] = true
                                indices[i] = indices_ub[i]
                                if i != d
                                    compoundIntersection = true
                                end
                            else
                                dim[i] = false
                            end
                        end
                        
                        push!(intersections, MeshIntersection(s_intercept, pt_intercept, copy(dim), copy(indices)) )
                        intersectionIndex += 1

                        # TODO: if the new point went past multiple bounds, find missing intersections?

                        # For compound intersections: update the bounds in each dimension involved
                        for i = 1:numDim
                            # If the intersection was ON a boundary (rather than crossing it),
                            # the bounds need to be updated with the next point and not the
                            # current one
                            if wasOnBoundary == false && dim[i] == true # this dim participated in the compound intersection
                                tighten_bounds(pt_new, i, coords, indices_lb, indices_ub)
                            end
                        end # Update bounds for-loop
                    end # if intersection occurred
                else # is outside domain; keep bounds up to date
                    tighten_bounds(pt_new, d, coords, indices_lb, indices_ub)
                end
            end # if intersection hasn't already been encountered
        end # d for-loop
        
        # If the most recently discovered intersection was a stop point,
        # increment the stop point index
        if i_stop_pts > 0 && i_stop_pts <= length(stop_pts) && intersectionIndex > 0 && abs(intersections[intersectionIndex].s - stop_pts[i_stop_pts]) <= s_tol 
            i_stop_pts += 1
        # If the next stop point is between the current point and the next point or is the end point, add it now
        elseif i_stop_pts > 0 && i_stop_pts <= length(stop_pts) && 
            ( s + ds > stop_pts[i_stop_pts] || (s == 1 && stop_pts[i_stop_pts] == 1) )
            push!(intersections, MeshIntersection(stop_pts[i_stop_pts], curve(stop_pts[i_stop_pts]), false*dim, copy(indices_lb))) 
            intersectionIndex += 1
            i_stop_pts += 1
        end
        s = s_new
        s_new = s + get_ds(ds,s)

        # Make sure to test the endpoint
        if s_new > 1 && s != 1
            s_new = 1
        end
        
        pt_curr = pt_new
        pt_new = curve(s_new)
        
        if wasOnBoundary == true
            for d = 1:numDim
                if dim[d] == true
                    tighten_bounds(pt_new, d, coords, indices_lb, indices_ub)
                end
            end
        end
    end # s while-loop

    if closed_list == true && sum(intersections[1].pt - intersections[end].pt .> closure_tol) != 0
        push!(intersections, intersections[1])
    end
    return intersections
end


function find_mesh_intersections(coords, curves::Union{Array, Tuple}, ds=DEFAULT_DS, 
    arc_tol=100*eps(eltype(coords[1])), corner_tol=100*eps(eltype(coords[1])))
    numCurves = length(curves)
    intersectionsByCurve = []

    # Allow scalar arguments without changing functionality
    if typeof(ds) <:Real
        ds = ds*ones(numCurves)
    end

    if typeof(arc_tol) <:Real
        arc_tol = arc_tol*ones(numCurves)
    end

    if typeof(corner_tol) <:Real
        corner_tol = corner_tol*ones(numCurves)
    end

    # Process each curve
    if typeof(ds) <: Function
        for c = 1:numCurves
            intersections = find_mesh_intersections(coords, curves[c], ds, arc_tol[c], corner_tol[c])
            push!(intersectionsByCurve, intersections)
        end
    else
        for c = 1:numCurves
            intersections = find_mesh_intersections(coords, curves[c], ds[c], arc_tol[c], corner_tol[c])
            push!(intersectionsByCurve, intersections)
        end
    end
    return intersectionsByCurve
end
