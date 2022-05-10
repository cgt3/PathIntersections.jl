# Author: Christina G Taylor
# Date:   April 2022
# Description:
#   Code for identifying intersections between a multi-dimensional Cartesian 
#   grid and parameterized paths. Intersection points are solved using a 
#   bracketted secant method. 
#
# Missing Functionality:
#  1) Non-simple intersections are not yet fully supported; some may be
#     captured, others may not depending on the step size (ds) used.
#
#  2) Lines of intersection will return multiple points of intersection
#     all ds distance from one another. The endpoints may not be returned.
#
#  3) Curves with sufficiently varying arc length per step may result in
#     missing intersections
#
# Inputs/User-defined objects:
#   coords                : An array of arrays of coordiantes. Each subarray
#                           contains the Cartesian mesh coordinates in that
#                           dimension. Each subarray must be in increasing
#                           order and contain the endpoints of the domain.
#                            Indexing: coords[dim][coord index]
#   curves                : An array of functions. Each function parameterizes a
#                           single object. Functions are called using:
#                            pt = curve[curve index](s, params[curve index])
#   ds                    : Array of step sizes to use when searching each curve.
#                           Must have same dimension as 'curves'.
#   arc_tol                : Array of tolerances used to find intersection points.
#                           The Euclidean distance from the approximate intersection
#                           point and the true intersection point is guaranteed to
#                           be less than this value. Must have same dimension as 
#                           'curves'.
# corner_tol               : The tolerance used to determine if a given dimension
#                           participated in a previously found intersection.
#
# Outputs/Results:
#   intersections         : A vector of (s,pt,dim,I) tuples describing where the
#                           curves intersect the mesh where
#                               s = the approximate parameter value where the curve
#                                   intersects the mesh
#                              pt = the point returned by the curve at s
#                             dim = tuple of booleans indicating which dimension(s)
#                                   the intersection occured in where
#                                      true = that dim participated in the intersection
#                                     false = that dim did not participate
#                               I = a tuple of indices into the coordinate
#                                   arrays corresponding to the intersected value
#                                   (if that dimension participated in the 
#                                   intersection) or the closest coord value less
#                                   than the intersection point
#
# #########################################################################################



module PathIntersections
using LinearAlgebra
using StructArrays

export Intersection

export find_mesh_intersections

mutable struct Intersection
    s::Real
    pt::AbstractArray{Real}
    dim::Vector{Bool}
    I::AbstractArray{Real}
end

## Helper Functions -------------------------------------------------------------
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


# Function for finding multiple dim intersections against a single curve
function find_mesh_intersections(coords, curve::Function, ds::Real, arc_tol::Real, corner_tol::Real)
    numDim = length(coords)
    # Start walking along the curve at s = 0
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
    intersections = StructArray{Intersection}(undef, 0)
    intersectionIndex = 0
    dim = Array{Bool,1}(undef, numDim)

    # TODO: make sure to check the end point (should be the same as the initial point for
    #       closed curves)
    s_new = s
    pt_new = pt_curr
    while s <= 1
        # Check is we are currently in the domain or crossing into or out of it
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
                        
                        push!(intersections, Intersection(s_intercept, pt_intercept, copy(dim), copy(indices)) )
                        intersectionIndex += 1

                        # TODO: if the new point went past multiple bounds, find missing intersections?

                        # For compound intersections: update the bounds in each dimension involved
                        
                        # tighten_bounds(pt_new, d, coords, indices_lb, indices_ub)
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
        
        # TODO: add adaptive stepping in s for curves with non-uniform arc-length per step?
        s = s_new
        s_new = s + ds

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

    return intersections
end


function find_mesh_intersections(coords, curves::Array, ds, arc_tol, corner_tol)
    numCurves = length(curves)
    intersectionsByCurve = []

    for c = 1:numCurves
        intersections = find_mesh_intersections(coords, curves[c], ds[c], arc_tol[c], corner_tol[c])
        push!(intersectionsByCurve, intersections)
    end
    return intersectionsByCurve
end

end # module