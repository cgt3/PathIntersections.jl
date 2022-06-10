using ForwardDiff
using LinearAlgebra
using Revise

function map_line_quadrature(ref_quadr, curve::Function, s_domain::Tuple;
    ref_domain=DEFAULT_REF_DOMAIN,
    normalization=DEFAULT_NORMALIZATION)

    # Map the reference quadr points to s values
    ref_pts, ref_weights = ref_quadr
    x_lb, x_ub = ref_domain
    s_lb, s_ub = s_domain
    scaling = (s_ub - s_lb)/(x_ub - x_lb)
    s_pts = @. scaling*(ref_pts - x_lb) + s_lb


    # Calculate the derivative and normal values at each of the s values
    dc_ds(s) = ForwardDiff.derivative(curve, s)
    normals = outward_normal.(dc_ds, s_pts, normalization=normalization)

    # Calculate the evaluation points and adjusted weights
    line_pts = @. curve(s_pts)
    line_weights = @. ref_weights.*norm(dc_ds(s_pts))*scaling

    return line_pts, line_weights, normals
end



function map_line_quadrature(ref_quadr, curve, consectutive_stop_points::AbstractArray;
    ref_domain=DEFAULT_REF_DOMAIN,
    normalization=DEFAULT_NORMALIZATION )

    pts_by_segment = []
    wts_by_segment = []
    normals_by_segment = []
    for i = 2:length(consectutive_stop_points)
        segment_pts, segment_wts, segment_normals = map_line_quadrature(ref_quadr, curve,
            (consectutive_stop_points[i-1], consectutive_stop_points[i]), ref_domain=ref_domain,
             normalization=normalization)
        push!(pts_by_segment, segment_pts)
        push!(wts_by_segment, segment_wts)
        push!(normals_by_segment, segment_normals)
    end

    return pts_by_segment, wts_by_segment, normals_by_segment
end



function map_line_quadrature_nonconsecutive(ref_quadr, curve, sub_bounds;
    ref_domain=DEFAULT_REF_DOMAIN,
    normalization=DEFAULT_NORMALIZATION )

    pts_by_segment = []
    wts_by_segment = []
    normals_by_segment = []
    for i = 1:length(sub_bounds)
        segment_pts, segment_wts, segment_normals = map_line_quadrature(ref_quadr, curve, sub_bounds[i], 
            ref_domain=ref_domain, normalization=normalization)
        push!(pts_by_segment, segment_pts)
        push!(wts_by_segment, segment_wts)
        push!(normals_by_segment, segment_normals)
    end

    return pts_by_segment, wts_by_segment, normals_by_segment
end
