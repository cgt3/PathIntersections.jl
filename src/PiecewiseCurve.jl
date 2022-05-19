
struct ConstFunction <:Function
    value::Real
end
function (f::ConstFunction)(x)
    return f.value
end

struct PiecewiseCurve{T_points, T_curves, T_bounds, 
    T_continuity, T_closed, T_subArcLengths, T_arcLength} <: Function # Callable object for piecewise defined curves

    stop_pts::T_points # Warning: these must be in increasing order and include s = 0,1
    subcurves::T_curves
    sub_bounds::T_bounds
    
    is_continuous::T_continuity
    is_closed::T_closed
    arc_lengths::T_subArcLengths
    total_arc_length::T_arcLength
end # struct PiecwiseCurve def

function PiecewiseCurve(input_stop_pts, input_subcurves, input_s_bounds; continuity_tol=DEFAULT_CONTINUITY_TOL, steps_per_segment=1000)
    # TODO: Remove this for PiecewiseFunctions
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
        if i < num_curves && norm(input_subcurves[i](input_s_bounds[i][2]) .- input_subcurves[i+1](input_s_bounds[i+1][1]) ) >= continuity_tol
            is_continuous = false
        end

        # TODO: arc lengths are not needed for PiecewiseFunction
        # Compute the arc length of this segment
        sub_ds = (input_s_bounds[i][2] - input_s_bounds[i][1]) / steps_per_segment
        sub_s_prev = input_s_bounds[i][1] 
        pt_prev = input_subcurves[i](sub_s_prev)
        for j = 1:steps_per_segment
            sub_s = sub_s_prev + sub_ds
            pt = input_subcurves[i](sub_s)
            arc_lengths[i] += norm(pt .- pt_prev)

            sub_s_prev = sub_s
            pt_prev = pt
        end
    end
    total_arc_length = sum(arc_lengths)

    # TODO: is_closed is not needed for PiecewiseFunction
    # If the curve is continuous and its end points match, mark it as closed
    is_closed = false
    if is_continuous == true && norm(input_subcurves[end](input_s_bounds[end][2]) .- input_subcurves[1](input_s_bounds[1][1]) ) < continuity_tol
        is_closed = true
    end

    return PiecewiseCurve(input_stop_pts, input_subcurves, input_s_bounds, is_continuous, is_closed, arc_lengths, total_arc_length)
end

function PiecewiseCurve(input_stop_pts, input_subcurves; continuity_tol=DEFAULT_CONTINUITY_TOL)
    input_s_bounds = @. [ (input_stop_pts[i], input_stop_pts[i+1]) for i=1:length(input_stop_pts)-1 ]
    return PiecewiseFunction(input_stop_pts, input_subcurves, input_s_bounds, continuity_tol = continuity_tol)
end

# TODO: PiecewiseFunction is currently broken
const PiecwiseFunction{T_pts, T_curves, T_bounds, T_continuity} = PiecewiseCurve{T_pts, T_curves, T_bounds, T_continuity, nothing, nothing, nothing}
function PiecewiseFunction(input_stop_pts, input_subcurves, input_s_bounds; continuity_tol=DEFAULT_CONTINUITY_TOL)
    return PiecewiseFunction(input_stop_pts, input_subcurves, input_s_bounds; continuity_tol=DEFAULT_CONTINUITY_TOL, num_steps=1000)
end

function (curve::PiecewiseCurve)(s)
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
        sub_s = curve.sub_bounds[i][2]
    else
        sub_s = (curve.sub_bounds[i][2] - curve.sub_bounds[i][1]) / (curve.stop_pts[i+1] - curve.stop_pts[i])*(s - curve.stop_pts[i]) + curve.sub_bounds[i][1]
    end
    return curve.subcurves[i](sub_s)
end

function ds_by_num_steps(curve::PiecewiseCurve, num_steps; continuity_tol=DEFAULT_CONTINUITY_TOL)
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

function ds_by_arc_length(curve::PiecewiseCurve, step_arc_length::Real; continuity_tol::Real=DEFAULT_CONTINUITY_TOL)
    subfunctions = ConstFunction[]
    for i = 1:length(curve.subcurves)
        subvalue = (curve.sub_bounds[i][2] - curve.sub_bounds[i][1]) / (curve.arc_lengths[i] / step_arc_length)
        push!(subfunctions, ConstFunction(subvalue))
    end
    return PiecewiseFunction(curve.stop_pts, subfunctions, continuity_tol)
end