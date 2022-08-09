module PresetGeometries

using StaticArrays
import LinearAlgebra.norm

import ..PiecewiseCurve

# Note: these are callable structs
export Line, Ellipse, Circle, Pacman

struct Line{T_pt1, T_pt2} <: Function
    start_pt::T_pt1
    end_pt::T_pt2
end

struct Rectangle{T_Lx, T_Ly, T_x0, T_y0, T_theta0, T_orientation, T_pts, T_function} <: Function
    Lx::T_Lx
    Ly::T_Ly
    x0::T_x0
    y0::T_y0
    theta0::T_theta0
    orientation::T_orientation #CurveOrientation
    # These should not be set by users; calculated by the constructor
    stop_pts::T_pts
    func::T_function
end

struct Ellipse{T_Rx, T_Ry, T_x0, T_y0, T_theta0, T_orientation} <: Function
    Rx::T_Rx
    Ry::T_Ry
    x0::T_x0 # For the center of the rectangle
    y0::T_y0 # For the center of the rectangle
    theta0::T_theta0
    orientation::T_orientation
end

struct Pacman{T_radius, T1_first, T1_second, T2_first, T2_second, T_x0, T_y0, T_orientation, T_pts, T_function} <: Function
    R::T_radius
    first_jaw::T1_first
    second_jaw::T1_second
    theta_lb::T2_first
    theta_ub::T2_second
    x0::T_x0
    y0::T_y0
    orientation::T_orientation
    # These should not be set by users; calculated by the constructor
    stop_pts::T_pts
    func::T_function 
end

# Needed to default arguments
function Rectangle(; Lx=1, Ly=1, x0=0, y0=0, theta0=0, orientation=1)

    if orientation == 1
        right  = Line(( Lx/2, -Ly/2), ( Lx/2,  Ly/2))
        top    = Line(( Lx/2,  Ly/2), (-Lx/2,  Ly/2))
        left   = Line((-Lx/2,  Ly/2), (-Lx/2, -Ly/2))
        bottom = Line((-Lx/2, -Ly/2), ( Lx/2, -Ly/2))
        edges = SVector(right, top, left, bottom)

        first_corner = Ly / (2*(Lx + Ly))
    else
        bottom = Line(( Lx/2, -Ly/2), (-Lx/2, -Ly/2))
        left   = Line((-Lx/2, -Ly/2), (-Lx/2,  Ly/2))
        top    = Line((-Lx/2,  Ly/2), ( Lx/2,  Ly/2))
        right  = Line(( Lx/2,  Ly/2), ( Lx/2, -Ly/2))
        edges= SVector(bottom, left, top, right)

        first_corner = Lx / (2*(Lx + Ly))
    end
    stop_pts = (0, first_corner, 0.5, 0.5+first_corner, 1)
    sub_bounds = ((0,1), (0,1), (0,1), (0,1))
    func = PiecewiseCurve(stop_pts, edges, sub_bounds)

    return Rectangle(Lx, Ly, x0, y0, theta0, orientation, stop_pts, func)
end

function Ellipse(; Rx=1, Ry=1, x0=0, y0=0, theta0=0, orientation=1)
    return Ellipse(Rx, Ry, x0, y0, theta0, orientation)
end

# Simplified geometries
const Square{T_L, T_x0, T_y0, T_theta0, T_orientation} = Rectangle{T_L, T_L, T_x0, T_y0, T_theta0, T_orientation}
function Square(; L=1, x0=0, y0=0, theta0=0, orientation=1)
    return Rectangle(Lx=L, Ly=L, x0=x0, y0=y0, theta0=theta0, orientation=orientation)
end

const Circle{T_R, T_x0, T_y0, T_theta0, T_orientation} = Ellipse{T_R, T_R, T_x0, T_y0, T_theta0, T_orientation}
function Circle(; R=1, x0=0, y0=0, theta0=0, orientation=1)
    return Ellipse(R, R, x0, y0, theta0, orientation)
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
            jaw_diff = 2*pi - jaw_diff
        else
            jaw_diff = 2*pi + jaw_diff
        end
    end
    jaw_diff_mag = abs(jaw_diff)

    total_arc_length = R*(2*pi - jaw_diff_mag) + 2*R
    s1 = R / total_arc_length
    s2 = (total_arc_length - R) / total_arc_length

    # Compute the jaw corner point locations for the lines
    circle_ref = Circle(R=R, x0=x0, y0=y0)
    first_corner = circle_ref(first_jaw_bounded/(2*pi))
    second_corner = circle_ref(second_jaw_bounded/(2*pi))
    
    # Define the jaws so that they are the upper and lower bounds of an interval
    if orientation > 0 # positive orientation
        theta_lb = first_jaw_bounded
        theta_ub = first_jaw_bounded + jaw_diff_mag
    else
        theta_lb = second_jaw_bounded
        theta_ub = second_jaw_bounded + jaw_diff_mag
    end

    if theta_ub > 2*pi || theta_lb > 2*pi
        theta_lb -= 2*pi
        theta_ub -= 2*pi
    elseif  theta_ub < -2*pi || theta_lb < -2*pi
        theta_lb += 2*pi
        theta_ub += 2*pi
    end

    # Create arrays for the PiecewiseCurve
    stop_pts = (0, s1, s2, 1)
    subcurves = ( Line((x0,y0), first_corner), 
                  Circle(R=R, x0=x0, y0=y0, orientation=orientation), 
                  Line(second_corner, (x0,y0))
                );
    s_bounds = ( (0,1), (orientation*first_jaw_bounded/(2*pi), (orientation*first_jaw_bounded + jaw_diff_mag)/(2*pi)), (0,1) )
    
    func = PiecewiseCurve(stop_pts, subcurves, s_bounds)
    return Pacman(R, first_jaw_bounded, second_jaw_bounded, theta_lb, theta_ub, x0, y0, orientation, stop_pts, func)
end

# Make the structs callable
function (L::Line)(s)
    return SVector( (L.end_pt[1] - L.start_pt[1])*s + L.start_pt[1], (L.end_pt[2] - L.start_pt[2])*s + L.start_pt[2] )
end

function (R::Rectangle)(s) # Also does squares
    (x,y) = R.func(s)
    r = norm((x,y))
    theta = angle(x + y*im)
    return SVector(r*cos(theta + R.theta0) + R.x0, r*sin(theta + R.theta0) + R.y0)
end

function (E::Ellipse)(s) # Also does circles
    (x,y) = (E.Rx*cos(E.orientation*2*pi*s), E.Ry*sin(E.orientation*2*pi*s))
    r = norm((x,y))
    theta = angle(x + y*im)
    return SVector(r*cos(theta + E.theta0) + E.x0, r*sin(theta + E.theta0) + E.y0)
end

function (p::Pacman)(s)
    return p.func(s)
end
    
end # module