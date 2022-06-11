# WARNING: if ds is too small you may see lines of mislabelled points near 
#          aligned with concave corners and/or the start-stop point (s=0,1)
function is_contained(test_pt, curve; ds=1e-4)
    # 1) Count the number of times rays drawn from the pt cross the curve
    ray_intersections = zeros(Int, 2, length(test_pt))
    above, below = 1,2

    pt_curve = curve(0)
    is_less_than_prev = pt_curve .< test_pt
    s = ds
    is_at_end = false
    while is_at_end == false
        pt_curve = curve(s)
        # Use both < and > to ensure equality is not counted
        is_less_than = pt_curve .< test_pt
        is_not_equal = pt_curve .!= test_pt
        ray_was_crossed = @. (is_less_than != is_less_than_prev) && is_not_equal

        # Note: reverse only works for 2D
        ray_was_crossed = reverse(ray_was_crossed)
        above_indices = pt_curve .> test_pt .&& ray_was_crossed
        below_indices = pt_curve .< test_pt .&& ray_was_crossed

        ray_intersections[above, above_indices] .+= ray_was_crossed[above_indices]
        ray_intersections[below, below_indices] .+= ray_was_crossed[below_indices]

        # Get ready for the next iteration
        is_less_than_prev[is_not_equal] = is_less_than[is_not_equal]
        s += ds
        if s >= 1
            s = 1
            is_at_end = true
        end
    end

    # 2) Determine if an intersection occurred or not
    for i = 1:length(test_pt)
        # For a point to be inside the curve ALL its rays must have an
        # odd number of intersections (not including non-simple intersections)
        if ray_intersections[above,i] % 2 == 0 || ray_intersections[below,i] % 2 == 0
            # 3) Consider the curve's orientation too
            return false
        end
    end

    # 3) 
    return true
end