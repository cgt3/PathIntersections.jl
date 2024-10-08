module PresetGeometries

using StaticArrays
import LinearAlgebra.norm

import ..PiecewiseCurve

# Note: these are callable structs
export Line, Ellipse, Circle, Pacman, Fish, BiconvexAirfoil

struct Line{T_pt1, T_pt2} <: Function
    start_pt::T_pt1
    end_pt::T_pt2
end

struct Rectangle{T_Lx, T_Ly, T_x0, T_y0, T_theta0, T_orientation, T_pts, T_function} <: Function
    Lx::T_Lx
    Ly::T_Ly
    x0::T_x0 # for the center of the rectangle
    y0::T_y0 # for the center of the rectangle
    theta0::T_theta0
    orientation::T_orientation #CurveOrientation
    # These should not be set by users; calculated by the constructor
    stop_pts::T_pts
    func::T_function
end

struct Ellipse{T_Rx, T_Ry, T_x0, T_y0, T_theta0, T_orientation} <: Function
    Rx::T_Rx
    Ry::T_Ry
    x0::T_x0 # For the center of the ellipse
    y0::T_y0 # For the center of the ellipse
    theta0::T_theta0
    orientation::T_orientation
end

struct Pacman{T_radius, T1_first, T1_second, T2_first, T2_second, T_x0, T_y0, T_orientation, T_pts, T_function} <: Function
    R::T_radius
    first_jaw::T1_first   # The user-defined angle of the first jaw/radii
    second_jaw::T1_second # The user-defiend angle of the second jaw/radii
    theta1::T2_first      # The angle of the first jaw/radii; guaranteed to be in [0, 2pi]
    theta2::T2_second     # The angle of the second jaw/radii; > theta1 for pos orientation, < theta1 for neg
    x0::T_x0
    y0::T_y0
    orientation::T_orientation
    # These should not be set by users; calculated by the constructor
    stop_pts::T_pts
    func::T_function 
end

const FISH_T_START = 0.08;

struct Fish{T_x0, T_y0, float, T_pts, T_func} <: Function
    x0::T_x0
    y0::T_y0
    scale::float
    stop_pts::T_pts
    func::T_func
end

struct BiconvexAirfoil{Float64, T_x0, T_y0, T_pts, T_function} <: Function
    R_top::Float64
    R_btm::Float64
    L::Float64 # chord length
    t::Float64 # maximum thickness
    AoA::Float64 # angle of attack, measured CLOCKWISE from the -x axis
    x0::T_x0 # for the center of the chord line
    y0::T_y0 # for the center of the chord line
    scale::Float64
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
    stop_pts = (0.0, first_corner, 0.5, 0.5+first_corner, 1.0)
    sub_bounds = ((0.0, 1.0), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0))
    func = PiecewiseCurve(stop_pts, edges, sub_bounds)

    return Rectangle(Lx, Ly, x0, y0, theta0, orientation, stop_pts, func)
end

function Ellipse(; Rx=1.0, Ry=1.0, x0=0.0, y0=0.0, theta0=0.0, orientation=1)
    return Ellipse(Rx, Ry, x0, y0, theta0, orientation)
end

# Simplified geometries
const Square{T_L, T_x0, T_y0, T_theta0, T_orientation} = Rectangle{T_L, T_L, T_x0, T_y0, T_theta0, T_orientation}
function Square(; L=1.0, x0=0.0, y0=0.0, theta0=0.0, orientation=1)
    return Rectangle(Lx=L, Ly=L, x0=x0, y0=y0, theta0=theta0, orientation=orientation)
end

const Circle{T_R, T_x0, T_y0, T_theta0, T_orientation} = Ellipse{T_R, T_R, T_x0, T_y0, T_theta0, T_orientation}
function Circle(; R=1.0, x0=0.0, y0=0.0, theta0=0.0, orientation=1)
    return Ellipse(R, R, x0, y0, theta0, orientation)
end

# Constructor to enforce first/second_jaw in [0,2pi] for constructing func,
# the PiecewiseCurve function
function Pacman(; R=1, first_jaw=pi/4, second_jaw=7*pi/4, x0=0.0, y0=0.0, orientation=1)
    # Enforce theta1 (the first jaw) in [0, 2pi]
    theta1 = first_jaw
    while theta1 < 0
        theta1 += 2*pi
    end
    while theta1 > 2*pi
        theta1 -= 2*pi
    end

    # Start theta2 (the second jaw) in [0, 2pi]
    theta2 = second_jaw
    while theta2 < 0
        theta2 += 2*pi
    end
    while theta2 > 2*pi
        theta2 -= 2*pi
    end

    # Enforce theta2 (the second jaw) > theta1 for pos orientation
    # or theta2 < theta1 for neg orientation
    while orientation > 0 && theta2 < theta1
        theta2 += 2*pi
    end
    while orientation < 0 && theta2 > theta1
        theta2 -= 2*pi
    end

    # Compute the s values at the corners by balancing arc-length
    jaw_diff = theta2 - theta1

    # Compute the jaw corner point locations for the lines
    circle_ref = Circle(R=R, x0=x0, y0=y0)
    first_corner = circle_ref(theta1/(2*pi))
    second_corner = circle_ref(theta2/(2*pi))
    
    # Define the stop points for balanced arc length
    total_arc_length = R*abs(jaw_diff) + 2*R
    s1 = R / total_arc_length
    s2 = (total_arc_length - R) / total_arc_length

    # Create arrays for the PiecewiseCurve
    stop_pts = (0, s1, s2, 1)
    subcurves = ( Line((x0,y0), first_corner), 
                  Circle(R=R, x0=x0, y0=y0), # Note: the orientation is reversed by the subbounds when necessary 
                  Line(second_corner, (x0,y0))
                );
    sub_bounds = ( (0.0, 1.0), (theta1/(2*pi), theta2/(2*pi)), (0.0, 1.0) )
    
    func = PiecewiseCurve(stop_pts, subcurves, sub_bounds)
    return Pacman(R, first_jaw, second_jaw, theta1, theta2, x0, y0, orientation, stop_pts, func)
end

function Fish(; scale=1, x0=0.0, y0=0.0)
    Rx = 2.5/5
    Ry = 1.25/5

    p1 = (x0 + scale*Rx*cos(2*pi*FISH_T_START), y0 + scale*Ry*sin(2*pi*FISH_T_START))
    p2 = (x0 + scale*(Rx + 0.5*Ry), y0 + scale*Ry)
    p3 = (p2[1], y0 - scale*Ry)
    p4 = (p1[1], y0 - scale*Ry*sin(2*pi*FISH_T_START))

    # Based on arc length for even step sizes
    t1 = 0.77
    t2 = 0.81
    t3 = 0.95

    stop_pts = (0.0, t1, t2, t3, 1.0)
    subcurves = ( Ellipse(Rx=scale*Rx, Ry=scale*Ry, x0=x0, y0=y0), 
                  Line(p4, p3),
                  Line(p3, p2),
                  Line(p2, p1) )
    sub_bounds = ( (FISH_T_START, 1 - FISH_T_START), (0.0,1.0), (0.0,1.0), (0.0,1.0))
    func = PiecewiseCurve(stop_pts, subcurves, sub_bounds)
    return Fish(x0, y0, scale, stop_pts, func)
end

# Default is the biconvex airfoil from Spring 2021 MECH 591 Gas Dynamics Exam I
# normalized to a chord length of 1.
function BiconvexAirfoil(; R1=40.0, R2=55.0, L=12.0, AoA=0.0, scale=1.0, x0=0.0, y0=0.0, normalize=true)
    if L > R1 || L > R2
        throw(error("PresentGeometries.jl: BiconvexAirfoil constructor: chord length must be smaller than both radii."))
    end

    # Compute the angles swept within each of the circles
    theta1_trailing = acos(0.5*L/R1)
    theta1_leading  = pi - theta1_trailing
    
    theta2_leading  =   pi + acos(0.5*L/R2)
    theta2_trailing = 2*pi - acos(0.5*L/R2)

    if normalize == true
        R1 = R1 / L
        R2 = R2 / L
        L = 1.0
    end
    
    # Compute the thickness
    t_above = R1*(1-sin(theta1_trailing))
    t_below = R2*(1-sin(acos(0.5*L/R2)))
    t = t_above + t_below

    # Calculate the circle centers relative to (x0, y0)
    y1_0 = - (R1 - t_above)
    y2_0 =   (R2 - t_below)

    # Calculate arc length of each surface
    arcL1 = R1*(theta1_leading - theta1_trailing)
    arcL2 = R2*(theta2_trailing - theta2_leading)
    stop_pts = (0.0, arcL1 / (arcL1 + arcL2), 1.0)

    subcurves = (
        Circle(R=scale*R1, x0=0.0, y0=scale*y1_0),
        Circle(R=scale*R2, x0=0.0, y0=scale*y2_0)
    )

    sub_bounds = ( (theta1_trailing, theta1_leading)./(2*pi), (theta2_leading, theta2_trailing)./(2*pi) )

    func = PiecewiseCurve(stop_pts, subcurves, sub_bounds)
    return BiconvexAirfoil(R1, R2, L, t, AoA, x0, y0, scale, stop_pts, func)
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

function (P::Pacman)(s)
    return P.func(s)
end

function(F::Fish)(s)
    return F.func(s)
end

function(B::BiconvexAirfoil)(s)
    (x,y) = B.func(s)
    r = norm((x,y))
    theta = angle(x + y*im)
    return SVector(r*cos(theta - B.AoA) + B.x0, r*sin(theta - B.AoA) + B.y0)
end
    
end # module