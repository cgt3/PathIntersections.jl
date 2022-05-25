using ForwardDiff
using LinearAlgebra
using Revise

function map_line_quadrature(ref_quadr, curve::Function, s_domain; ref_domain=(-1,1), normalization=false)

    # Map the reference quadr points to s values
    ref_pts, ref_weights = ref_quadr
    x_lb, x_ub = ref_domain
    s_lb, s_ub = s_domain
    s_pts = @. (s_ub - s_lb)/(x_ub - x_lb)*(ref_pts - x_lb) + s_lb

    # Calculate the derivative and normal values at each of the s values
    dc_ds(s) = ForwardDiff.derivative(curve, s)
    tangents = @. dc_ds(s_pts)


    normals = [ [ tangents[i][2], -tangents[i][1] ] for i in 1:length(tangents) ]
    if normalization == true
        for i = 1:length(normals)
            normals[i] = normals[i] / norm(normals[i])
        end
    end

    # Calculate the evaluation points and adjusted weights
    line_pts = @. curve(s_pts)
    line_weights = @. ref_weights.*norm(dc_ds(s_pts))

    return line_pts, line_weights, normals
end


function map_line_quadrature_consecutive(ref_quadr,
    curve,
    stop_points;
    ref_domain_all=(-1,1),
    normalization_all=false )

    pts_by_segment = []
    wts_by_segment = []
    normals_by_segment = []
    for i = 2:length(stop_points)
        segment_pts, segment_wts, segment_normals = map_line_quadrature(ref_quadr, curve, (stop_points[i-1], stop_points[i]), 
            ref_domain=ref_domain_all, normalization=normalization_all)
        push!(pts_by_segment, segment_pts)
        push!(wts_by_segment, segment_wts)
        push!(normals_by_segment, segment_normals)
    end

    return pts_by_segment, wts_by_segment, normals_by_segment
end



function map_line_quadrature_multiple(ref_quadr,
    curve,
    sub_bounds;
    ref_domain_all=(-1,1),
    normalization_all=false )

    pts_by_segment = []
    wts_by_segment = []
    normals_by_segment = []
    for i = 1:length(sub_bounds)
        segment_pts, segment_wts, segment_normals = map_line_quadrature(ref_quadr, curve, sub_bounds[i], 
            ref_domain=ref_domain_all, normalization=normalization_all)
        push!(pts_by_segment, segment_pts)
        push!(wts_by_segment, segment_wts)
        push!(normals_by_segment, segment_normals)
    end

    return pts_by_segment, wts_by_segment, normals_by_segment
end
