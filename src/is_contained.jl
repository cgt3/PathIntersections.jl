# WARNING: if ds is too small you may see lines of mislabelled points near 
#          aligned with concave corners and/or the start-stop point (s=0,1)
function is_contained(test_pt, curve; ds=1e-4)
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

        # Note: reverse only works for 2D
        ray_was_crossed = reverse(ray_was_crossed)
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