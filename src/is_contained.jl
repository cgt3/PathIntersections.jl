# WARNING: if ds is too small you may see lines of mislabelled points near 
#          aligned with concave corners and/or the start-stop point (s=0,1)
function is_contained(curve, test_pt; ds=1e-4)
    numDim = length(test_pt)
    # 1) Count the number of times rays drawn from the pt cross the curve
    ray_intersections = zeros(Int, 2, length(test_pt))
    above, below = 1,2

    # pre-allocate
    is_less_than = zeros(Bool, numDim)
    is_less_than_prev = zeros(Bool, numDim)
    is_not_equal = zeros(Bool, numDim)
    ray_was_crossed = zeros(Bool, numDim)
    above_indices = zeros(Bool, numDim)
    below_indices = zeros(Bool, numDim)

    pt_curve = curve(1-ds)
    is_less_than_prev = pt_curve .< test_pt
    s = 0
    is_at_end = false
    while is_at_end == false
        pt_curve = curve(s)
        # Use both < and > to ensure equality is not counted
        @. is_less_than = pt_curve < test_pt
        @. is_not_equal = pt_curve != test_pt
        @. ray_was_crossed = (is_less_than != is_less_than_prev) && is_not_equal

        # Note: reversing only works for 2D
        reverse!(ray_was_crossed)
        @. above_indices = pt_curve > test_pt && ray_was_crossed
        @. below_indices = pt_curve < test_pt && ray_was_crossed

        @. ray_intersections[above, above_indices] += ray_was_crossed[above_indices]
        @. ray_intersections[below, below_indices] += ray_was_crossed[below_indices]

        # Get ready for the next iteration
        @. is_less_than_prev[is_not_equal] = is_less_than[is_not_equal]
        s += ds
        if s >= 1
            s = 1
            is_at_end = true
        end
    end

    # 2) Determine if an intersection occurred or not
    is_enclosed = true
    for i = 1:numDim
        # For a point to be inside the curve ALL its rays must have an
        # odd number of intersections (not including non-simple intersections)
        if ray_intersections[above,i] % 2 == 0 || ray_intersections[below,i] % 2 == 0
            is_enclosed = false
            break
        end
    end

    # TODO: take the curve's orientation into account
    return is_enclosed
end


function is_contained(E::PresetGeometries.Ellipse, pt)
    # Remove the rotation
    r = sqrt((pt[1]-E.x0)^2 + (pt[2]-E.y0)^2)
    theta = angle((pt[2]-E.x0) + (pt[2]-E.y0)*im)
    x,y = r*cos(theta - E.theta0), r*sin(theta - E.theta0)

    if E.orientation == 1
        return (x/E.Rx)^2 + (y/E.Ry)^2 <= 1
    else
        return (x/E.Rx)^2 + (y/E.Ry)^2 >= 1
    end
end


function is_contained(R::PresetGeometries.Rectangle, pt)
    # Remove the rotation
    r = sqrt((pt[1]-R.x0)^2 + (pt[2]-R.y0)^2)
    theta = angle((pt[1]-R.x0) + (pt[2]-R.y0)*im)
    x,y = r*cos(theta - R.theta0), r*sin(theta - R.theta0)

    return abs(x) < R.Lx/2 && abs(y) <= R.Ly/2
end

function is_contained(P::PresetGeometries.Pacman, pt)
    # Calculate the points radius and angle from the center of the Pacman
    x = pt[1] - P.x0
    y = pt[2] - P.y0

    r = sqrt(x^2 + y^2)
    theta = angle(x + y*im)
    theta_lb = (P.orientation == 1) ? P.theta1 : P.theta2
    theta_ub = (P.orientation == 1) ? P.theta2 : P.theta1
    if theta > theta_lb + 2*pi
        theta -= 2*pi
    elseif theta < theta_ub - 2*pi
        theta += 2*pi
    end

    return r <= P.R && theta_lb <= theta <= theta_ub
end


function is_contained(F::PresetGeometries.Fish, pt)
    # Check the ellipse
    if is_contained(F.func.subcurves[1], pt)
        return true
    end

    # Check if we are past the tail in x
    if pt[1] < F.func.subcurves[2].start_pt[1] || 
       pt[1] > F.func.subcurves[2].end_pt[1]
        return false
    end
    
    # Coarse check of whether we are past the tail in y0
    if pt[2] < F.func.subcurves[2].end_pt[2] || 
       pt[2] > F.func.subcurves[3].end_pt[2]
         return false
     end

    # Fine tuned check if we are past the tail in y
    t = (pt[1] - F.func.subcurves[2].start_pt[1]) / (F.func.subcurves[2].end_pt[1] - F.func.subcurves[2].start_pt[1] )
    y_fish = F.func.subcurves[2](t)[2]

    if pt[2] < y_fish || pt[2] > 2*F.y0 - y_fish
        return false
    end

    return true
end


function is_contained(B::PresetGeometries.BiconvexAirfoil, pt)
    # Remove the rotation from the AoA
    r = sqrt((pt[1]-B.x0)^2 + (pt[2]-B.y0)^2)
    theta = angle((pt[1]-B.x0) + (pt[2]-B.y0)*im)
    x_new, y_new = r*cos(theta + B.AoA), r*sin(theta + B.AoA)

    # Check the two circles defining the airfoil
    if is_contained(B.func.subcurves[1], (x_new,y_new)) && 
       is_contained(B.func.subcurves[2], (x_new,y_new))
        return true
    else
        return false
    end
end

# generalize `is_contained` to multiple curves
is_contained(curves::Tuple, pt; kwargs...) = all(map(c -> is_contained(c, pt; kwargs...), curves))
