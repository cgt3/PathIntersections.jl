
struct ConstFunction
    value::Real
end

function (f::ConstFunction)(x)
    return f.value
end

struct PiecewiseFunction # Callable object
    stop_pts::Array # Warning: these must be in increasing order and include s = 0,1
    subfunctions::Array
    sub_x_bounds::Array
    
    is_continuous::Bool

    function PiecewiseFunction(input_stop_pts, input_subfunctions, input_x_bounds; continuity_tol=DEFAULT_CONTINUITY_TOL)
        if length(input_stop_pts)-1 != length(input_subfunctions)
            error("PiecewiseFunction: inconsistent number of stop points and subfunctions")
        end
        if length(input_stop_pts)-1 != length(input_x_bounds)
            error("PiecewiseFunction: inconsistent number of stop points and bounds")
        end

        # Calculate arc_length of the segments and total arc length
        num_functions = length(input_subfunctions)
        is_continuous = true
        for i = 1:num_functions
            # Update continuity if needed
            if i < num_functions && norm(input_subfunctions[i](input_x_bounds[i][2]) - input_subfunctions[i+1](input_x_bounds[i+1][1]) ) >= continuity_tol
                is_continuous = false
            end

            # TODO: Add an option for enforcing left-continuity at stop points instead of right?
            # Check if the subfunctions are defined at their start points (the end point is assumed to belong to the next curve)
            if isnan( input_subfunctions[i](input_x_bounds[i][1]) )
                is_continuous = false
            end
        end
        return new(input_stop_pts, input_subfunctions, input_x_bounds, is_continuous)
    end

    function PiecewiseFunction(input_stop_pts, input_subfunctions; continuity_tol=DEFAULT_CONTINUITY_TOL)
        num_curves = length(input_stop_pts) - 1
        input_x_bounds = @. [ [input_stop_pts[i], input_stop_pts[i+1]] for i=1:num_curves ]
        return PiecewiseFunction(input_stop_pts, input_subfunctions, input_x_bounds, continuity_tol = continuity_tol)
    end
end # PiecewiseFunction struct def

function (f::PiecewiseFunction)(x)
    # Find the subfunction to call
    i = 1 # Index to the closest stop point <= x
    # Note: left endpoints below to the current subfunction, right endpoints to the next function
    while i < length(f.stop_pts) && x >= f.stop_pts[i+1]
        i += 1
    end 

    # Map x to the subfunctions's proper x value
    if i == length(f.stop_pts)
        i -= 1
        sub_x = f.sub_x_bounds[i][2]
    else
        sub_x = (f.sub_x_bounds[i][2] - f.sub_x_bounds[i][1]) / (f.stop_pts[i+1] - f.stop_pts[i])*(x - f.stop_pts[i]) + f.sub_x_bounds[i][1]
    end
    return f.subfunctions[i](sub_x)
end 


struct PiecewiseCurve # Callable object for piecewise defined curves
    stop_pts::Array # Warning: these must be in increasing order and include s = 0,1
    subcurves::Array
    sub_s_bounds::Array
    
    is_continuous::Bool
    is_closed::Bool
    arc_lengths::Array
    total_arc_length::Real

    function PiecewiseCurve(input_stop_pts, input_subcurves, input_s_bounds; continuity_tol=DEFAULT_CONTINUITY_TOL, steps_per_segment=1000)
        if input_stop_pts[1] != 0
            error("PiecewiseCurve: stop_pts must start at 0")
        end
        if input_stop_pts[end] != 1
            error("PiecewiseCurve: stop_pts must end at 1")
        end

        if length(input_stop_pts)-1 != length(input_subcurves)
            error("PiecewiseCurve: inconsistent number of stop points and subcurves")
        end
        if length(input_stop_pts)-1 != length(input_s_bounds)
            error("PiecewiseCurve: inconsistent number of stop points and bounds")
        end

        # Calculate arc_length of the segments and total arc length
        arc_lengths = zeros(length(input_subcurves))
        num_curves = length(input_subcurves)
        is_continuous = true
        for i = 1:num_curves
            # Update continuity if needed
            if i < num_curves && norm(input_subcurves[i](input_s_bounds[i][2]) - input_subcurves[i+1](input_s_bounds[i+1][1]) ) >= continuity_tol
                is_continuous = false
            end

            # Compute the arc length of this segment
            sub_ds = (input_s_bounds[i][2] - input_s_bounds[i][1]) / steps_per_segment
            sub_s_prev = input_s_bounds[i][1] 
            pt_prev = input_subcurves[i](sub_s_prev)
            for j = 1:steps_per_segment
                sub_s = sub_s_prev + sub_ds
                pt = input_subcurves[i](sub_s)
                arc_lengths[i] += norm(pt - pt_prev)

                sub_s_prev = sub_s
                pt_prev = pt
            end
        end
        total_arc_length = sum(arc_lengths)

        # If the curve is continuous and its end points match, mark it as closed
        is_closed = false
        if is_continuous == true && norm(input_subcurves[end](input_s_bounds[end][2]) - input_subcurves[1](input_s_bounds[1][1]) ) < continuity_tol
            is_closed = true
        end

        return new(copy(input_stop_pts), copy(input_subcurves), copy(input_s_bounds), is_continuous, is_closed, arc_lengths, total_arc_length)
    end
    
    function PiecewiseCurve(input_stop_pts, input_subcurves; continuity_tol=DEFAULT_CONTINUITY_TOL)
        input_s_bounds = @. [ (input_stop_pts[i], input_stop_pts[i+1]) for i=1:length(input_stop_pts)-1 ]
        return PiecewiseFunction(input_stop_pts, input_subcurves, input_s_bounds, continuity_tol = continuity_tol)
    end

    # TODO: Should these be outside the struct definition?
    function ds_by_num_steps(curve::PiecewiseCurve, num_steps::Integer; continuity_tol=DEFAULT_CONTINUITY_TOL)
        # Determine how much of the curve's length comes from each of its segments
        # and use this fraction to allot steps.
        length_fractions = curve.arc_lengths / curve.total_arc_length
        steps_by_segment = trunc.(num_steps * length_fractions + 0.5)

        # Modify the number of steps on the last segment to make sure the correct
        # number of steps is returned.
        sum_steps = sum(steps_by_segment)
        extra_steps = sum_steps - num_steps
        steps_by_segment[end] -= extra_steps

        # calculate ds for each segment
        subfunctions = ConstFunction[]
        for i = 1:length(curve.subcurves)
            push!(subfunctions, ConstFunction(curve.arc_length[i]) / steps_by_segment[i])
        end

        return PiecewiseFunction(curve.stop_pts, subfunctions, continuity_tol)
    end

    function ds_by_arc_length(C::PiecewiseCurve, step_arc_length::Real; continuity_tol::Real=DEFAULT_CONTINUITY_TOL)
        subfunctions = ConstFunction[]
        for i = 1:length(curve.subcurves)
            subvalue = (curve.sub_s_bounds[i][2] - curve.sub_s_bounds[i][1])*(curve.arc_lengths[i] / step_arc_length)
            push!(subfunctions, ConstFunction(subvalue))
        end
        return PiecewiseFunction(stop_pts, subfunctions, continuity_tol)
    end
end # struct PiecwiseCurve def

function (curve::PiecewiseCurve)(s::Real)
    # Wrap s back to [0,1]
    if s > 1 || s < -1 # => bring s to [-1, 1]
        s = s - trunc(Int,s)
    end
    if s < 0 # => brings s to [0,1]
        s = 1 + s
    end

    # Find the curve to operate on
    i = 1 # Index to the closest stop point <= s
    while i < length(curve.stop_pts) && s >= curve.stop_pts[i+1]
        i += 1
    end 

    # Map s to the subcurve's proper s value
    if i == length(curve.stop_pts)
        i = i-1
        sub_s = curve.sub_s_bounds[i][2]
    else
        sub_s = (curve.sub_s_bounds[i][2] - curve.sub_s_bounds[i][1]) / (curve.stop_pts[i+1] - curve.stop_pts[i])*(s - curve.stop_pts[i]) + curve.sub_s_bounds[i][1]
    end
    return curve.subcurves[i](sub_s)
end
