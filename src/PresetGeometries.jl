module PresetGeometries

import ..PiecewiseCurve

# Note: these functions return functions or PiecewiseCurves
export customCircle
export customEllipse
export customPacman

# Note: these functions are curves
export circle, ellipse, line
export unitCircle
export defaultPacman

# Non-exported functions
function circle(s::Real, radius::Real, x0::Real, y0::Real)
    return [radius*cos(2*pi*s) + x0, radius*sin(2*pi*s) + y0]
end

function ellipse(s::Real, rx::Real, ry::Real, x0::Real, y0::Real, offset_angle::Real)
    return [rx*cos(2*pi*s+offset_angle) + x0, ry*sin(2*pi*s+offset_angle) + y0]
end

function line(s::Real, start_pt::Union{Array, Tuple}, end_pt::Union{Array, Tuple})
    return [(end_pt[1] - start_pt[1])*s + start_pt[1], (end_pt[2] - start_pt[2])*s + start_pt[2] ] 
end

function customCircle(; radius::Real=1, x0::Real=0, y0::Real=0)
    return s->circle(s, radius, x0, y0)
end

function customEllipse(; rx::Real=1, ry::Real=1, x0::Real=0, y0::Real=0, offset_angle::Real=0)
    return s->ellipse(s, rx, ry, x0, y0, offset_angle)
end


function customPacman(; radius::Real=1, x0::Real=0, y0::Real=0, upper_jaw::Real=pi/4, lower_jaw::Real=-pi/4)
    # Enforce upper_jaw < lower_jaw
    while lower_jaw < 0
        lower_jaw += 2*pi
    end
    while upper_jaw < 0
        upper_jaw += 2*pi
    end

    if upper_jaw == lower_jaw
        lower_jaw = upper_jaw + 2*pi
    elseif upper_jaw > lower_jaw
        swap = upper_jaw
        upper_jaw = lower_jaw
        lower_jaw = swap
    end

    # Compute the s values at the corners by balancing arc-length
    total_arc_length = 2*pi*radius / (2*pi - (upper_jaw-lower_jaw)) + 2*radius
    s1 = radius / total_arc_length
    s2 = (total_arc_length - radius) / total_arc_length

    # Compute the jaw corner point locations
    upper_corner = radius*cos(upper_jaw)+x0, radius*sin(upper_jaw)+y0
    lower_corner = radius*cos(lower_jaw)+x0, radius*sin(lower_jaw)+y0

    # Create arrays for the PiecewiseCurve
    stop_pts = [0, s1, s2, 1]
    subcurves = [(s->line(s, [x0,y0], upper_corner)), 
                 (s->circle(s, radius, x0, y0)), 
                 (s->line(s, lower_corner, [x0,y0]))
                ];
    s_bounds = [ [0,1], [upper_jaw/(2*pi), lower_jaw/(2*pi)], [0,1] ]
    # return PiecewiseCurve(copy(stop_pts), copy(subcurves), copy(s_bounds))
    
    pacman = PiecewiseCurve(stop_pts, subcurves, s_bounds)
    return pacman
end

unitCircle = customCircle()
defaultPacman = customPacman()
    
end # module