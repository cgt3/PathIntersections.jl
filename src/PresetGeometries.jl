module PresetGeometries

import ..PiecewiseCurve

# Note: these are callable structs
export Line, Ellipse, Circle, Pacman

export CurveOrientation
@enum CurveOrientation pos=1 neg=-1

abstract type PresetGeometry end

struct Line{T_start, T_end} <: Function
    start_pt::T_start
    end_pt::T_end
end

struct Ellipse{T_Rx, T_Ry, T_x0, T_y0, T_theta0, T_orientation} <: Function
    Rx::T_Rx
    Ry::T_Ry
    x0::T_x0
    y0::T_y0
    theta0::T_theta0
    orientation::T_orientation #CurveOrientation
end

# Needed to default arguments
function Ellipse(; Rx=1, Ry=1, x0=0, y0=0, theta0=0, orientation=1)
    return Ellipse(Rx, Ry, x0, y0, theta0, orientation)
end

const Circle{T_R, T_x0, T_y0, T_theta0, T_orientation} = Ellipse{T_R, T_R, T_x0, T_y0, T_theta0, T_orientation}
function Circle(; R=1, x0=0, y0=0, theta0=0, orientation=1)
    return Ellipse(R, R, x0, y0, theta0, orientation)
end

struct Pacman{T_radius, T_first, T_second, T_x0, T_y0, T_orientation} <: Function
    R::T_radius
    first_jaw::T_first
    second_jaw::T_second
    x0::T_x0
    y0::T_y0
    orientation::T_orientation #CurveOrientation
    func # should not be set by users; calculated by the constructor
end

# Constructor to enforce first/second_jaw in [0,2pi] for constructing func,
# the PiecewiseCurve function
function Pacman(; R=1, first_jaw=pi/4, second_jaw=7*pi/4, x0=0, y0=0, orientation=1)
    # Enforce first_jaw in [0, 2pi]
    first_jaw_bounded = first_jaw
    while first_jaw_bounded < 0
        first_jaw_bounded += 2*pi
    end
    while first_jaw_bounded > 2*pi
        first_jaw_bounded -= 2*pi
    end

    # Enforce second_jaw in [0, 2pi]
    second_jaw_bounded = second_jaw
    while second_jaw_bounded < 0
        second_jaw_bounded += 2*pi
    end
    while second_jaw_bounded > 2*pi
        second_jaw_bounded -= 2*pi
    end

    # Compute the s values at the corners by balancing arc-length
    jaw_diff = second_jaw_bounded - first_jaw_bounded
    if orientation*sign(jaw_diff) == -1
        if sign(jaw_diff) == 1
            jaw_diff = 2*pi - (second_jaw_bounded - first_jaw_bounded)
        else
            jaw_diff = 2*pi + (second_jaw_bounded - first_jaw_bounded)
        end
    end

    total_arc_length = R*(2*pi - abs(jaw_diff)) + 2*R
    s1 = R / total_arc_length
    s2 = (total_arc_length - R) / total_arc_length

    # Compute the jaw corner point locations for the lines
    circle_ref = Circle(R=R, x0=x0, y0=y0)
    first_corner = circle_ref(first_jaw_bounded/(2*pi))
    second_corner = circle_ref(second_jaw_bounded/(2*pi))
    
    # Create arrays for the PiecewiseCurve
    stop_pts = (0, s1, s2, 1)
    subcurves = ( Line((x0,y0), first_corner), 
                  Circle(R=R, x0=x0, y0=y0, orientation=orientation), 
                  Line(second_corner, (x0,y0))
                );

    s_bounds = ( (0,1), (orientation*first_jaw_bounded/(2*pi), (orientation*first_jaw_bounded + abs(jaw_diff))/(2*pi)), (0,1) )
    
    func = PiecewiseCurve(stop_pts, subcurves, s_bounds)
    return Pacman(R, first_jaw_bounded, second_jaw_bounded, x0, y0, orientation, func)
end

# Make the structs callable
function (E::Ellipse)(s) # Also does circles
    return ( E.Rx*cos(E.orientation*2*pi*s + E.theta0) + E.x0, E.Ry*sin(E.orientation*2*pi*s + E.theta0) + E.y0)
end

function (L::Line)(s)
    return ( (L.end_pt[1] - L.start_pt[1])*s + L.start_pt[1], (L.end_pt[2] - L.start_pt[2])*s + L.start_pt[2] ) 
end

function (p::Pacman)(s)
    return p.func(s)
end
    
end # module